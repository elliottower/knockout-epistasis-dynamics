"""Run one ODE configuration on one network, and write what it produced.

    uv run python -m scripts.run_ode_arm --network lambda_phage \
        --construction hillcube_normalized --hill-n 10 \
        --solver experiments/ode_v2/solver_primary.json \
        --classifier experiments/ode_v2/classifier_primary.json \
        --n-init 32 --seed 42 --out results/ode_v2

Solver and classifier settings are read from JSON files rather than chosen here, so the
files the registration freezes are the files the run used. The legacy audit passes no
classifier.

A run has three stages, so that a large network can be spread over many machines
(`scripts/modal_ode_arm.py`); run from this command line, they execute in sequence:

  1. shards      each shard sweeps a contiguous range of coalitions into its own file,
                 checkpointed as it goes and resumable
  2. replicate   a deterministic sample of coalitions is recomputed from scratch, on a
                 different machine when run on Modal, in chunks with their own files
  3. finalize    the shards are checked to tile every coalition exactly once under one
                 configuration, then merged; the replicates are compared; the record is written

Files under --out, for <stem> = <network>__<construction>__n<hill_n>__k<hill_k>__c<clamp_value>:

    <stem>.shard_<start>-<end>.npz   per-coalition and per-trajectory arrays for one range
    <stem>.replicate_<i>.npz         one chunk of the replicate sample, with full tails
    <stem>.json                      status and class totals, replicate agreement, scores, provenance

A sweep with any coalition lacking an accepted value is written with its totals and
`"scored": false`, and no spectrum. It is never scored on the coalitions that happened to
succeed.
"""

import argparse
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np

from grn_coalition_sweep import extract_rule_fourier
from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.ode_engine import (
    CONSTRUCTIONS,
    STATUSES,
    TRAJECTORY_CLASSES,
    ClassifierConfig,
    CompiledNetwork,
    DynamicsConfig,
    ODEEngineError,
    SolverAttempt,
    SolverConfig,
    compile_model,
    initial_states,
    prepare,
    provenance,
    simulate_coalition,
    sweep,
)
from scripts.attractors import network_basins
from scripts.paths import RESULTS
from scripts.walsh_estimators import (
    ENERGY,
    ESTIMATORS,
    EstimatorError,
    fractions,
    higher_order_ratio,
    order3_fraction,
)

BOOLEAN_REFERENCE = RESULTS / "paper_number_reconciliation.json"
NULL_BAND = 0.005  # 0.5 percentage points of spectral energy, the paper's null band
REPLICATE_FRACTION = 0.05
# Per replicate chunk: full tails for this many coalitions in sample order, and for up to this
# many more that contain an oscillating, extended or failed trajectory.
TAILS_IN_SAMPLE_ORDER = 4
TAILS_FLAGGED = 12
# Two runs of the same trajectory agree when |x - y| <= REPLICATE_ATOL + REPLICATE_RTOL * max(|x|, |y|).
# The legacy engine differs from itself by up to 1.14e-12 relative on one platform
# (results/audit/solver_determinism.json).
REPLICATE_ATOL = 1e-12
REPLICATE_RTOL = 1e-9
# Fingerprint fields that describe the machine rather than the configuration. Shards run on
# different machines of one image, so these may differ between shards only if the image does.
PLATFORM_FIELDS = ("python", "machine", "system", "numpy", "scipy", "blas", "container_image")
# Fields only a replicate file's fingerprint carries.
REPLICATE_FIELDS = ("role", "chunk", "coalitions_sha256")


class ShardError(ODEEngineError):
    """The shard files on disk do not add up to one complete run of one configuration."""


def load_solver(path: Path) -> SolverConfig:
    spec = json.loads(path.read_text())
    try:
        attempts = tuple(SolverAttempt(a["method"], int(a["max_rhs_evals"]), a["when"]) for a in spec["attempts"])
        return SolverConfig(attempts=attempts, rtol=float(spec["rtol"]), atol=float(spec["atol"]),
                            max_step=float(spec["max_step"]), t_max=float(spec["t_max"]),
                            t_tail=float(spec["t_tail"]), n_tail_samples=int(spec["n_tail_samples"]))
    except KeyError as err:
        raise ODEEngineError(f"{path}: solver file is missing {err}") from err


def load_classifier(path: Path) -> ClassifierConfig:
    spec = json.loads(path.read_text())
    try:
        return ClassifierConfig(
            fixed_tail_range_max=float(spec["fixed_tail_range_max"]),
            fixed_derivative_max=float(spec["fixed_derivative_max"]),
            oscillatory_tail_range_min=float(spec["oscillatory_tail_range_min"]),
            halves_envelope_tol=float(spec["halves_envelope_tol"]),
            bounded_margin=float(spec["bounded_margin"]),
            extension_horizons=tuple(float(h) for h in spec["extension_horizons"]))
    except KeyError as err:
        raise ODEEngineError(f"{path}: classifier file is missing {err}") from err


@dataclass(frozen=True)
class ArmSpec:
    """Everything that defines one run. Every stage rebuilds its inputs from this alone."""

    network: str
    construction: str
    hill_n: float
    solver_path: Path
    classifier_path: Path | None
    n_init: int
    seed: int
    clamp_value: float = 0.0
    keep_states: bool = False
    hill_k: float = 0.5

    def __post_init__(self) -> None:
        if self.network not in ALL_MODELS:
            raise ODEEngineError(f"unknown network {self.network!r}")
        if self.construction not in CONSTRUCTIONS:
            raise ODEEngineError(f"unknown construction {self.construction!r}")
        if self.classifier_path is None and self.construction != "operatorwise_legacy":
            raise ODEEngineError("only the legacy audit runs without a classifier")
        if self.n_init < 2 or self.n_init % 2:
            raise ODEEngineError(f"n_init must be even and at least 2 for the split-half estimator, got {self.n_init}")

    @property
    def stem(self) -> str:
        return f"{self.network}__{self.construction}__n{self.hill_n:g}__k{self.hill_k:g}__c{self.clamp_value:g}"

    def net(self) -> CompiledNetwork:
        info = ALL_MODELS[self.network]
        return compile_model(info["rules"], info["output_nodes"])

    def dynamics(self) -> DynamicsConfig:
        return DynamicsConfig(self.construction, hill_n=self.hill_n, hill_k=self.hill_k)

    def solver(self) -> SolverConfig:
        return load_solver(self.solver_path)

    def classifier(self) -> ClassifierConfig | None:
        return load_classifier(self.classifier_path) if self.classifier_path else None


def plan_shards(n_nodes: int, coalitions_per_shard: int) -> list[tuple[int, int]]:
    total = 2**n_nodes
    return [(s, min(s + coalitions_per_shard, total)) for s in range(0, total, coalitions_per_shard)]


def shard_file(out: Path, spec: ArmSpec, start: int, end: int) -> Path:
    return out / f"{spec.stem}.shard_{start:07d}-{end:07d}.npz"


def replicate_file(out: Path, spec: ArmSpec, chunk: int) -> Path:
    return out / f"{spec.stem}.replicate_{chunk:03d}.npz"


def run_shard(spec: ArmSpec, start: int, end: int, out: Path, checkpoint_every: int = 64,
              on_checkpoint: Callable[[], None] | None = None) -> str:
    """Sweep coalitions [start, end) into their shard file. Resumes from the file if present."""
    out.mkdir(parents=True, exist_ok=True)
    run = sweep(spec.net(), spec.dynamics(), spec.solver(), network=spec.network, n_init=spec.n_init,
                seed=spec.seed, clamp_value=spec.clamp_value, classifier=spec.classifier(),
                coalition_range=(start, end), checkpoint=shard_file(out, spec, start, end),
                checkpoint_every=checkpoint_every, keep_states=spec.keep_states, on_checkpoint=on_checkpoint)
    return run["fingerprint"]


def configuration_of(fingerprint: dict) -> dict:
    """A fingerprint without its coalition range, machine fields and replicate fields: what every
    shard and replicate of one run must share."""
    return {k: v for k, v in fingerprint.items()
            if k != "coalition_range" and k not in PLATFORM_FIELDS and k not in REPLICATE_FIELDS}


def platform_of(fingerprint: dict) -> dict:
    return {k: fingerprint[k] for k in PLATFORM_FIELDS}


def expected_fingerprint(spec: ArmSpec) -> dict:
    """The fingerprint, less its coalition range, that this spec produces with the network and
    code as they are on disk now."""
    prov = provenance(prepare(spec.net(), spec.dynamics()), spec.network, spec.solver(), spec.classifier(),
                      spec.n_init, spec.seed, spec.clamp_value)
    return json.loads(json.dumps({**prov, "keep_states": spec.keep_states}))


def differing_keys(a: dict, b: dict) -> list[str]:
    return sorted(k for k in a.keys() | b.keys() if a.get(k) != b.get(k))


def merge_shards(spec: ArmSpec, out: Path) -> tuple[dict, dict, dict]:
    """Merge the shard files into one run. Raises ShardError unless they tile every coalition
    exactly once, each is complete, each was run under exactly the configuration this spec now
    produces (network, model digest, code digest, dynamics, solver, classifier, seeds), and all
    ran on one platform. Returns the arrays, the configuration and the platform."""
    configuration = configuration_of(expected_fingerprint(spec))
    total = 2 ** len(spec.net().node_names)
    files = sorted(f for f in out.glob(f"{spec.stem}.shard_*-*.npz") if not f.name.endswith(".partial.npz"))
    if not files:
        raise ShardError(f"no shard files for {spec.stem} in {out}")
    shards = []
    for f in files:
        with np.load(f) as z:
            fingerprint = json.loads(str(z["fingerprint"]))
            start, end = fingerprint["coalition_range"]
            if int(z["completed"]) != end - start:
                raise ShardError(f"{f.name}: {int(z['completed'])} of {end - start} coalitions complete")
            shards.append((start, end, fingerprint, {k: z[k] for k in z.files if k not in ("completed", "fingerprint")}))
    shards.sort(key=lambda s: s[0])
    expected_start = 0
    for start, end, _, _ in shards:
        if start != expected_start:
            raise ShardError(f"coalitions [{expected_start}, {start}) are missing or overlap for {spec.stem}")
        expected_start = end
    if expected_start != total:
        raise ShardError(f"coalitions [{expected_start}, {total}) are missing for {spec.stem}")
    platform = platform_of(shards[0][2])
    for start, end, fingerprint, _ in shards:
        differing = differing_keys(configuration_of(fingerprint), configuration)
        if differing:
            raise ShardError(f"shard [{start}, {end}) differs from the requested run in {differing}")
        if platform_of(fingerprint) != platform:
            raise ShardError(f"shard [{start}, {end}) was run on a different platform")
    arrays = {k: np.concatenate([s[3][k] for s in shards]) for k in shards[0][3]}
    return arrays, configuration, platform


def replicate_coalitions(n_nodes: int, configuration: dict, arrays: dict) -> tuple[np.ndarray, dict]:
    """A deterministic sample, drawn with a generator seeded from the run's configuration:

      sentinels     every node clamped, no node clamped, and the coalition whose trajectories
                    took the most right-hand-side evaluations in total
      by size       REPLICATE_FRACTION of each coalition size, at least one per size
      retried       coalitions with a retried solver attempt, all of them up to the flag cap and
                    a seeded subsample of the cap beyond it
      extended      coalitions with an extended horizon, likewise
      oscillating   REPLICATE_FRACTION of the coalitions with an oscillating trajectory, at least one

    The flag cap is REPLICATE_FRACTION of all coalitions. A run with any coalition lacking a
    value cannot be scored, so it is not replicated: the sample is empty and the composition
    says why.
    """
    missing = int(np.isnan(arrays["values"]).sum())
    if missing:
        return np.array([], dtype=np.int64), {
            "skipped": f"{missing} of {2**n_nodes} coalitions have no value; the run cannot be scored"}
    digest = hashlib.sha256(json.dumps(configuration, sort_keys=True).encode()).digest()
    rng = np.random.default_rng(int.from_bytes(digest[:8], "little"))
    klass, horizon = arrays["trajectory_class"], arrays["trajectory_horizon_index"]
    cap = max(1, int(round(REPLICATE_FRACTION * 2**n_nodes)))
    coalitions = np.arange(2**n_nodes)
    sizes = np.array([bin(c).count("1") for c in coalitions])
    by_size = []
    for k in range(n_nodes + 1):
        stratum = coalitions[sizes == k]
        by_size.extend(rng.choice(stratum, size=max(1, int(round(REPLICATE_FRACTION * stratum.size))), replace=False))

    def capped(flagged: np.ndarray) -> np.ndarray:
        return rng.choice(flagged, size=cap, replace=False) if flagged.size > cap else flagged

    retried_all = np.flatnonzero((arrays["trajectory_attempts_run"] > horizon + 1).any(axis=1))
    extended_all = np.flatnonzero((horizon > 0).any(axis=1))
    retried, extended = capped(retried_all), capped(extended_all)
    oscillating = np.flatnonzero((klass == TRAJECTORY_CLASSES.index("oscillatory")).any(axis=1))
    osc = (rng.choice(oscillating, size=max(1, int(round(REPLICATE_FRACTION * oscillating.size))), replace=False)
           if oscillating.size else np.array([], dtype=np.int64))
    sentinels = [0, 2**n_nodes - 1, int(np.argmax(arrays["trajectory_nfev"].astype(np.int64).sum(axis=1)))]
    sample = np.unique(np.concatenate([sentinels, by_size, retried, extended, osc]).astype(np.int64))
    composition = {"sentinels": len(set(sentinels)), "by_size": len(set(by_size)), "flag_cap": cap,
                   "retried_total": int(retried_all.size), "retried_sampled": int(retried.size),
                   "extended_total": int(extended_all.size), "extended_sampled": int(extended.size),
                   "oscillating_total": int(oscillating.size), "oscillating_sampled": int(osc.size),
                   "total_distinct": int(sample.size)}
    return sample, composition


def replicate_plan(spec: ArmSpec, out: Path) -> tuple[np.ndarray, dict]:
    arrays, configuration, _ = merge_shards(spec, out)
    return replicate_coalitions(len(spec.net().node_names), configuration, arrays)


def plan_replicates(spec: ArmSpec, out: Path, n_chunks: int) -> list[np.ndarray]:
    sample, _ = replicate_plan(spec, out)
    return [chunk for chunk in np.array_split(sample, n_chunks) if chunk.size]


def replicate_fingerprint(spec: ArmSpec, chunk: int, coalitions: np.ndarray) -> dict:
    return {**expected_fingerprint(spec), "role": "replicate", "chunk": chunk,
            "coalitions_sha256": hashlib.sha256(np.asarray(coalitions, dtype=np.int64).tobytes()).hexdigest()}


def run_replicate(spec: ArmSpec, chunk: int, coalitions: np.ndarray, path: Path) -> None:
    """Recompute the given coalitions from scratch and write their status, class and value,
    under a fingerprint of the run's configuration, platform and this chunk.

    Full tails and attempt logs are kept for a bounded number of coalitions per chunk: the first
    TAILS_IN_SAMPLE_ORDER in sample order, and up to TAILS_FLAGGED more that contain an
    oscillating, extended, retried or failed trajectory. Tails are the raw state; the scored
    output clips it to [0, 1].
    """
    net, solver, classifier = spec.net(), spec.solver(), spec.classifier()
    prep = prepare(net, spec.dynamics())
    n_nodes = len(net.node_names)
    init = initial_states(spec.n_init, n_nodes, spec.seed)
    # A tail accepted at a later horizon has more samples; pad every tail to the longest possible.
    last_horizon = classifier.extension_horizons[-1] if classifier and classifier.extension_horizons else solver.t_max
    width = int(round(solver.n_tail_samples * last_horizon / solver.t_max))
    status = np.zeros((coalitions.size, spec.n_init), dtype=np.int8)
    klass = np.zeros((coalitions.size, spec.n_init), dtype=np.int8)
    values = np.full(coalitions.size, np.nan)
    kept, tails, logs, n_flagged = [], [], {}, 0
    for j, c in enumerate(coalitions):
        r = simulate_coalition(prep, solver, int(c), init, spec.clamp_value, classifier, keep_tails=True)
        status[j], klass[j], values[j] = r.trajectories["status"], r.trajectories["class"], r.value
        flagged = any(t.status != "ok" or t.klass == "oscillatory" or t.horizon_index > 0
                      or len(t.attempt_log) > t.horizon_index + 1 for t in r.results)
        if j < TAILS_IN_SAMPLE_ORDER or (flagged and n_flagged < TAILS_FLAGGED):
            n_flagged += j >= TAILS_IN_SAMPLE_ORDER
            padded = np.full((spec.n_init, n_nodes, width), np.nan, dtype=np.float32)
            for i, t in enumerate(r.results):
                if t.tail is not None:
                    padded[i, :, :t.tail.shape[1]] = t.tail
            kept.append(int(c))
            tails.append(padded)
            logs[str(int(c))] = [t.attempt_log for t in r.results]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".partial.npz")
    np.savez_compressed(tmp, coalitions=coalitions, status=status, klass=klass, values=values,
                        tail_coalitions=np.array(kept, dtype=np.int64),
                        tails=np.stack(tails) if tails else np.zeros((0, spec.n_init, n_nodes, width), np.float32),
                        attempt_logs=np.array(json.dumps(logs)),
                        fingerprint=np.array(json.dumps(replicate_fingerprint(spec, chunk, coalitions), sort_keys=True)))
    tmp.replace(path)


def replicate_is_current(spec: ArmSpec, chunk: int, coalitions: np.ndarray, path: Path) -> bool:
    """True only if the file exists and holds exactly these coalitions under exactly this run's fingerprint."""
    if not path.exists():
        return False
    with np.load(path) as z:
        if "fingerprint" not in z.files or not np.array_equal(z["coalitions"], coalitions):
            return False
        return json.loads(str(z["fingerprint"])) == json.loads(json.dumps(replicate_fingerprint(spec, chunk, coalitions)))


def agree(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    both_nan = np.isnan(x) & np.isnan(y)
    close = np.abs(x - y) <= REPLICATE_ATOL + REPLICATE_RTOL * np.maximum(np.abs(x), np.abs(y))
    return both_nan | close


def compare_replicates(spec: ArmSpec, out: Path, arrays: dict, configuration: dict, platform: dict,
                       planned: np.ndarray) -> dict:
    """Check that every replicate file was produced under this run's configuration and platform,
    for its own chunk, and that together they cover the planned sample exactly; then compare."""
    files = sorted(f for f in out.glob(f"{spec.stem}.replicate_*.npz") if not f.name.endswith(".partial.npz"))
    if not files:
        raise ShardError(f"no replicate files for {spec.stem} in {out}")
    parts = []
    for f in files:
        with np.load(f) as z:
            if "fingerprint" not in z.files:
                raise ShardError(f"{f.name} carries no fingerprint")
            fingerprint = json.loads(str(z["fingerprint"]))
            part = {k: z[k] for k in ("coalitions", "status", "klass", "values")}
        differing = differing_keys(configuration_of(fingerprint), configuration)
        if differing:
            raise ShardError(f"{f.name} differs from the requested run in {differing}")
        if platform_of(fingerprint) != platform:
            raise ShardError(f"{f.name} was run on a different platform from the shards")
        chunk = int(f.name.rsplit(".replicate_", 1)[1].split(".")[0])
        digest = hashlib.sha256(part["coalitions"].astype(np.int64).tobytes()).hexdigest()
        if fingerprint.get("chunk") != chunk or fingerprint.get("coalitions_sha256") != digest:
            raise ShardError(f"{f.name} was written for another chunk or another coalition list")
        parts.append(part)
    coalitions = np.concatenate([p["coalitions"] for p in parts])
    if not np.array_equal(np.sort(coalitions), planned):
        raise ShardError(f"replicate files for {spec.stem} do not cover the planned sample exactly")
    status = np.concatenate([p["status"] for p in parts])
    klass = np.concatenate([p["klass"] for p in parts])
    values = np.concatenate([p["values"] for p in parts])
    same_status = status == arrays["trajectory_status"][coalitions]
    same_class = klass == arrays["trajectory_class"][coalitions]
    same_value = agree(values, arrays["values"][coalitions])
    diff = np.abs(values - arrays["values"][coalitions])
    n_nodes = len(spec.net().node_names)
    return {
        "n_coalitions": int(coalitions.size),
        "sentinels": {"all_clamped": 0, "none_clamped": 2**n_nodes - 1,
                      "most_evaluations": int(np.argmax(arrays["trajectory_nfev"].astype(np.int64).sum(axis=1)))},
        "status_mismatches": int((~same_status).sum()),
        "class_mismatches": int((~same_class).sum()),
        "value_mismatches": int((~same_value).sum()),
        "max_abs_value_diff": float(np.nanmax(diff)) if np.isfinite(diff).any() else None,
        "tolerance": {"atol": REPLICATE_ATOL, "rtol": REPLICATE_RTOL},
        "files": [f.name for f in files],
    }


def classify_delta(delta: float) -> str:
    if delta > NULL_BAND:
        return "creation"
    if delta < -NULL_BAND:
        return "destruction"
    return "null"


def quantiles(x: np.ndarray) -> dict | None:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return None
    q = np.quantile(x, [0.5, 0.9, 0.99, 1.0])
    return {"median": float(q[0]), "q90": float(q[1]), "q99": float(q[2]), "max": float(q[3])}


def summarize(arrays: dict, n_coalitions: int) -> dict:
    status, klass = arrays["trajectory_status"], arrays["trajectory_class"]
    ok = status == STATUSES.index("ok")
    fixed = ok & (klass == TRAJECTORY_CLASSES.index("fixed"))
    return {
        "n_coalitions": n_coalitions,
        "n_coalitions_without_value": int(np.isnan(arrays["values"]).sum()),
        "trajectories": {s: int((status == i).sum()) for i, s in enumerate(STATUSES)},
        "classes_of_accepted": {c: int((klass[ok] == i).sum()) for i, c in enumerate(TRAJECTORY_CLASSES)},
        "accepted_by_attempt": {str(a): int((arrays["trajectory_attempt"][ok] == a).sum())
                                for a in np.unique(arrays["trajectory_attempt"][ok])},
        "accepted_by_horizon_index": {str(h): int((arrays["trajectory_horizon_index"][ok] == h).sum())
                                      for h in np.unique(arrays["trajectory_horizon_index"][ok])},
        "trajectories_with_a_retry": int((arrays["trajectory_attempts_run"] > arrays["trajectory_horizon_index"] + 1).sum()),
        "coalitions_with_an_oscillating_trajectory": int((klass == TRAJECTORY_CLASSES.index("oscillatory")).any(axis=1).sum()),
        "unclassified_with_provisional_output": int(np.isfinite(arrays["trajectory_unclassified_output"]).sum()),
        "nfev": quantiles(arrays["trajectory_nfev"].astype(float)),
        "tail_range_fixed": quantiles(arrays["trajectory_tail_range"][fixed]),
        "final_derivative_fixed": quantiles(arrays["trajectory_final_derivative"][fixed]),
        "tail_range_all_accepted": quantiles(arrays["trajectory_tail_range"][ok]),
    }


def score(network: str, rules: dict, outputs: np.ndarray, n_nodes: int) -> dict:
    """Walsh energies, spectra and Delta_3+ under each registered estimator, from the per-trajectory
    outputs (scripts/walsh_estimators.py). `all_pairs` is primary. Each estimator reports its raw
    order energies, total, and the numerator and denominator of sigma_3+/sigma_2+, and whether its
    spectrum and its ratio are admissible. The local spectrum comes from the rules and has no
    sampling error. The published Boolean value is for reading the record; the registered
    comparisons are made across networks by scripts/analyze_ode_rerun.py.
    """
    local = np.asarray(extract_rule_fourier(rules)["local_energy_spectrum"], dtype=float)
    boolean = json.loads(BOOLEAN_REFERENCE.read_text())["per_network_table"][network]["delta_3plus_pp"] / 100
    estimates = {}
    for name in ESTIMATORS:
        energy = ENERGY[name](outputs, n_nodes)
        entry = {"energy_by_order": energy.tolist(), "total_energy": float(energy.sum()),
                 "energy_3plus": float(energy[3:].sum()), "energy_2plus": float(energy[2:].sum())}
        try:
            spectrum = fractions(energy)
        except EstimatorError as err:
            entry.update(admissible=False, reason=str(err))
        else:
            delta = order3_fraction(spectrum) - order3_fraction(local)
            entry.update(admissible=True, global_spectrum=spectrum.tolist(), global_o3plus=order3_fraction(spectrum),
                         delta_o3plus=delta, **{"class": classify_delta(delta)})
            try:
                entry.update(higher_order_ratio=higher_order_ratio(energy), ratio_admissible=True)
            except EstimatorError as err:
                entry.update(higher_order_ratio=None, ratio_admissible=False, ratio_reason=str(err))
        estimates[name] = entry
    return {
        "local_spectrum": local.tolist(),
        "local_o3plus": order3_fraction(local),
        "estimators": estimates,
        "published_boolean": {"delta_o3plus": boolean, "class": classify_delta(boolean)},
    }


def finalize(spec: ArmSpec, out: Path, launch_manifest: dict | None = None) -> dict:
    """Merge the shards, compare the replicates, and write the record. Returns the record."""
    arrays, configuration, platform = merge_shards(spec, out)
    n_nodes = len(spec.net().node_names)
    planned, composition = replicate_coalitions(n_nodes, configuration, arrays)
    replicate = (compare_replicates(spec, out, arrays, configuration, platform, planned)
                 if planned.size else {})
    record = {
        "network": spec.network,
        "construction": spec.construction,
        "hill_n": spec.hill_n,
        "clamp_value": spec.clamp_value,
        "n_nodes": n_nodes,
        "solver_file": str(spec.solver_path),
        "solver_file_contents": json.loads(spec.solver_path.read_text()),
        "classifier_file": str(spec.classifier_path) if spec.classifier_path else None,
        "classifier_file_contents": json.loads(spec.classifier_path.read_text()) if spec.classifier_path else None,
        "configuration": configuration,
        "configuration_sha256": hashlib.sha256(json.dumps(configuration, sort_keys=True).encode()).hexdigest(),
        "platform": platform,
        "launch_manifest": launch_manifest,
        "launch_manifest_sha256": (hashlib.sha256(json.dumps(launch_manifest, sort_keys=True).encode()).hexdigest()
                                   if launch_manifest else None),
        "summary": summarize(arrays, 2**n_nodes),
        "replicate": {**replicate, "sample_composition": composition},
        "finished_at": datetime.now(timezone.utc).isoformat(),
    }
    values = arrays["values"]
    rules = ALL_MODELS[spec.network]["rules"]
    record["scored"] = bool(not np.isnan(values).any())
    if record["scored"]:
        record["scores"] = score(spec.network, rules, arrays["outputs"], n_nodes)
        if spec.classifier_path is not None and "tail_min_states" in arrays:
            record["basins"] = network_basins(arrays)
    else:
        # If every failure is an unclassified trajectory, the registered sensitivity analysis scores
        # the network with their provisional tail-window outputs. Never the primary score.
        failed = arrays["trajectory_status"] != STATUSES.index("ok")
        if (arrays["trajectory_status"][failed] == STATUSES.index("unclassified")).all():
            provisional = np.where(failed, arrays["trajectory_unclassified_output"], arrays["outputs"])
            record["scores_provisional"] = score(spec.network, rules, provisional, n_nodes)
    json_path = out / f"{spec.stem}.json"
    tmp = json_path.with_name(json_path.stem + ".partial.json")
    tmp.write_text(json.dumps(record, indent=2, default=float))
    tmp.replace(json_path)
    return record


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--network", required=True)
    parser.add_argument("--construction", required=True, choices=CONSTRUCTIONS)
    parser.add_argument("--hill-n", type=float, required=True)
    parser.add_argument("--hill-k", type=float, default=0.5)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--classifier", type=Path)
    parser.add_argument("--n-init", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--clamp-value", type=float, default=0.0)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--keep-states", action="store_true")
    parser.add_argument("--coalitions-per-shard", type=int, default=4096)
    parser.add_argument("--checkpoint-every", type=int, default=64)
    args = parser.parse_args(argv)

    spec = ArmSpec(args.network, args.construction, args.hill_n, args.solver, args.classifier,
                   args.n_init, args.seed, args.clamp_value, args.keep_states, args.hill_k)
    for start, end in plan_shards(len(spec.net().node_names), args.coalitions_per_shard):
        run_shard(spec, start, end, args.out, args.checkpoint_every)
    for i, chunk in enumerate(plan_replicates(spec, args.out, 1)):
        run_replicate(spec, i, chunk, replicate_file(args.out, spec, i))
    record = finalize(spec, args.out)

    s, rep = record["summary"], record["replicate"]
    print(f"wrote {args.out / (spec.stem + '.json')}")
    print(f"  trajectories {s['trajectories']}  classes {s['classes_of_accepted']}  "
          f"coalitions without value {s['n_coalitions_without_value']}/{s['n_coalitions']}")
    if "skipped" in rep["sample_composition"]:
        print(f"  replicate skipped: {rep['sample_composition']['skipped']}")
    else:
        print(f"  replicate: {rep['n_coalitions']} coalitions, mismatches status {rep['status_mismatches']} "
              f"class {rep['class_mismatches']} value {rep['value_mismatches']}")
    # Scores stay in the record and are not printed, so no run shows a number before it is read
    # under its registration.
    print(f"  scored: {record['scored']}")


if __name__ == "__main__":
    main()
