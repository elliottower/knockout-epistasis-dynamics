"""Validate the attractor grouping of scripts/attractors.py before A1 is registered.

    uv run python -m scripts.validate_attractor_grouping

Both parts integrate trajectories under the primary configuration (normalized HillCube,
n_H = 10, K = 0.5, the frozen solver and classifier, 32 initial states, clamp value 0) and
compute no Walsh quantity.

  controls   four systems whose attractors are known, all nodes free, from 32 starts drawn from
             default_rng(SEED): a toggle switch and a 4-node ring (two fixed points each), a
             3-node ring (one limit cycle), and a switch that selects a 3-node or a 5-node ring
             (two limit cycles).
  panel      COALITIONS_PER_NETWORK coalitions of each of the 28 paper networks, drawn with
             default_rng(SEED) in equal numbers from the lower, middle and upper third of
             coalition sizes, from the registered initial states (default_rng(42)). Only pooled
             statistics are written. No value, attractor count or entropy is written or printed
             per network, so the panel cannot inform A1, which is an association across
             networks.

The gate, fixed before the panel was run:

  1. every control has its known number of attractors at every point of the tolerance grid;
  2. the panel holds at least MIN_MULTISTABLE coalitions with more than one fixed-point
     attractor at the primary tolerances, so that criteria 3 and 4 can fail;
  3. in at least STABLE_FRACTION of the panel's coalitions, the partition into attractors is the
     same at every point of the tolerance grid;
  4. in at least STABLE_FRACTION of them, complete and single linkage give the same partition at
     the primary tolerances.

Reported and not gated: the largest within-attractor and smallest between-attractor distances,
the spread of periods within a cycle attractor, the number of coalitions with two or more
oscillating trajectories, and the same stability fractions over those coalitions alone.

Writes results/audit/attractor_grouping_validation.json. Each panel coalition's trajectory
classes, periods and tail states (no outputs or values) are checkpointed under
results/audit/attractor_grouping_validation/, keyed by a digest of the code and configuration,
so that a rerun resumes. Only this script's pooled aggregation reads them.
"""

import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone

import numpy as np

from scripts.attractors import (
    CYCLE_TOL,
    CYCLE_TOL_GRID,
    FIXED,
    FIXED_TOL,
    FIXED_TOL_GRID,
    OSCILLATORY,
    PERIOD_RTOL,
    chebyshev,
    cycle_distance,
    labels_grid,
)
from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.ode_engine import (
    STATUSES,
    DynamicsConfig,
    code_digest,
    compile_model,
    initial_states,
    model_digest,
    prepare,
    simulate_coalition,
)
from scripts.paths import AUDIT, PROJECT_ROOT, RESULTS
from scripts.run_ode_arm import load_classifier, load_solver

SEED = 20260922
N_INIT = 32
INIT_SEED = 42
COALITIONS_PER_NETWORK = 32
MIN_MULTISTABLE = 20
STABLE_FRACTION = 0.99
SOLVER_FILE = PROJECT_ROOT / "experiments/ode_v2/solver_primary.json"
CLASSIFIER_FILE = PROJECT_ROOT / "experiments/ode_v2/classifier_primary.json"
DYNAMICS = DynamicsConfig("hillcube_normalized", hill_n=10.0, hill_k=0.5)
OUT = AUDIT / "attractor_grouping_validation.json"
CHECKPOINTS = AUDIT / "attractor_grouping_validation"
CONTROLS = {
    "toggle": ({"A": "!B", "B": "!A"}, 2),
    "ring4": ({"A": "!D", "B": "!A", "C": "!B", "D": "!C"}, 2),
    "ring3": ({"A": "!C", "B": "!A", "C": "!B"}, 1),
    "switched_rings": ({"S": "S", "A": "S & !C", "B": "S & !A", "C": "S & !B", "D": "!S & !H", "E": "!S & !D",
                        "F": "!S & !E", "G": "!S & !F", "H": "!S & !G"}, 2),
}
STATE_FIELDS = ("tail_means", "tail_mins", "tail_maxs", "cycle_means")


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def settings_digest(network: str) -> str:
    info = ALL_MODELS[network]
    return hashlib.sha256(json.dumps({
        "code_sha256": code_digest(), "model_sha256": model_digest(compile_model(info["rules"], info["output_nodes"])),
        "solver_sha256": hashlib.sha256(SOLVER_FILE.read_bytes()).hexdigest(),
        "classifier_sha256": hashlib.sha256(CLASSIFIER_FILE.read_bytes()).hexdigest(),
        "dynamics": [DYNAMICS.construction, DYNAMICS.hill_n, DYNAMICS.hill_k], "n_init": N_INIT, "seed": INIT_SEED,
    }, sort_keys=True).encode()).hexdigest()


def panel_coalitions(n_nodes: int, rng: np.random.Generator) -> np.ndarray:
    sizes = np.array([bin(c).count("1") for c in range(2**n_nodes)])
    bounds = np.linspace(0, n_nodes + 1, 4)
    per = COALITIONS_PER_NETWORK // 3
    picked = []
    for i, (low, high) in enumerate(zip(bounds[:-1], bounds[1:])):
        stratum = np.flatnonzero((sizes >= low) & (sizes < high))
        k = per + (COALITIONS_PER_NETWORK - 3 * per if i == 2 else 0)
        picked.extend(rng.choice(stratum, size=min(k, stratum.size), replace=False))
    return np.sort(np.array(picked, dtype=np.int64))


def run_coalition(task: tuple[str, int]) -> str:
    """Integrate one panel coalition and checkpoint what grouping needs; values are not kept."""
    network, coalition = task
    path = CHECKPOINTS / network / f"{coalition:07d}.npz"
    digest = settings_digest(network)
    if path.exists():
        with np.load(path) as z:
            if str(z["settings"]) == digest:
                return str(path)
    info = ALL_MODELS[network]
    net = compile_model(info["rules"], info["output_nodes"])
    init = initial_states(N_INIT, len(net.node_names), INIT_SEED)
    r = simulate_coalition(prepare(net, DYNAMICS), load_solver(SOLVER_FILE), coalition, init, 0.0,
                           load_classifier(CLASSIFIER_FILE), keep_states=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".partial.npz")
    np.savez_compressed(tmp, settings=np.array(digest), status=r.trajectories["status"], klass=r.trajectories["class"],
                        period=r.trajectories["period"], **{f: getattr(r, f) for f in STATE_FIELDS})
    tmp.replace(path)
    return str(path)


def canonical(labels: np.ndarray) -> tuple[int, ...]:
    """A partition's labels renumbered by first appearance, so equal partitions compare equal."""
    first: dict[int, int] = {}
    return tuple(first.setdefault(int(label), len(first)) for label in labels)


def distances(klass, mean, low, high, cycle, period, labels) -> dict:
    """Largest within-attractor and smallest between-attractor distance, per class. Cycle
    distances are scaled so that 1 is the primary tolerance."""
    out = {}
    for name in ("fixed", "cycle"):
        idx = np.flatnonzero(klass == (FIXED if name == "fixed" else OSCILLATORY))
        if idx.size < 2:
            continue
        d = (chebyshev(mean[idx].astype(float)) if name == "fixed" else
             cycle_distance(*(a[idx].astype(float) for a in (mean, low, high, cycle)), period[idx].astype(float), CYCLE_TOL))
        same = labels[idx][:, None] == labels[idx][None, :]
        off = ~np.eye(idx.size, dtype=bool)
        out[name] = {"within": float(d[same & off].max()) if (same & off).any() else None,
                     "between": float(d[~same].min()) if (~same).any() else None}
    return out


def group(klass, mean, low, high, cycle, period) -> dict:
    grid = labels_grid(klass, mean, low, high, cycle, period, FIXED_TOL_GRID, CYCLE_TOL_GRID)
    single = labels_grid(klass, mean, low, high, cycle, period, (FIXED_TOL,), (CYCLE_TOL,), "single")[(FIXED_TOL, CYCLE_TOL)]
    primary = grid[(FIXED_TOL, CYCLE_TOL)]
    fixed_groups = len(set(primary[klass == FIXED]))
    cycle_spread = {"cycle_mean": 0.0, "tail_min": 0.0, "tail_max": 0.0}
    for a in set(primary[klass == OSCILLATORY]):
        m = primary == a
        for key, states in (("cycle_mean", cycle), ("tail_min", low), ("tail_max", high)):
            cycle_spread[key] = max(cycle_spread[key], float(chebyshev(states[m].astype(float)).max()))
    periods = [np.ptp(period[primary == a]) / np.max(period[primary == a])
               for a in set(primary[klass == OSCILLATORY]) if np.isfinite(period[primary == a]).all() and (primary == a).sum() > 1]
    return {"n_attractors": {str(k): int(v.max()) + 1 for k, v in grid.items()},
            "same_across_grid": len({canonical(v) for v in grid.values()}) == 1,
            "single_equals_complete": canonical(single) == canonical(primary),
            "fixed_attractors": fixed_groups, "oscillating_trajectories": int((klass == OSCILLATORY).sum()),
            "oscillating_without_period": int(((klass == OSCILLATORY) & np.isnan(period)).sum()),
            "max_relative_period_spread": float(max(periods)) if periods else None,
            "cycle_within_attractor_spread": cycle_spread,
            "distances": distances(klass, mean, low, high, cycle, period, primary)}


def controls() -> dict:
    rng = np.random.default_rng(SEED)
    solver, classifier = load_solver(SOLVER_FILE), load_classifier(CLASSIFIER_FILE)
    out = {}
    for name, (rules, expected) in CONTROLS.items():
        net = compile_model(rules, [next(iter(rules))])
        init = rng.random((N_INIT, len(rules)))
        r = simulate_coalition(prepare(net, DYNAMICS), solver, 2 ** len(rules) - 1, init, 0.0, classifier, keep_states=True)
        if (r.trajectories["status"] != STATUSES.index("ok")).any():
            out[name] = {"expected": expected, "passed": False, "reason": "a trajectory failed"}
            continue
        g = group(r.trajectories["class"], *(getattr(r, f) for f in STATE_FIELDS[:3]), r.cycle_means, r.trajectories["period"])
        out[name] = {"expected": expected, "found": g["n_attractors"],
                     "passed": all(v == expected for v in g["n_attractors"].values())}
    return out


def pooled(rows: list[dict]) -> dict:
    def fraction(key, subset):
        return sum(r[key] for r in subset) / len(subset) if subset else None

    def extreme(cls, which, fn):
        vals = [r["distances"][cls][which] for r in rows if cls in r["distances"] and r["distances"][cls][which] is not None]
        return float(fn(vals)) if vals else None

    oscillating = [r for r in rows if r["oscillating_trajectories"] >= 2]
    spreads = [r["max_relative_period_spread"] for r in rows if r["max_relative_period_spread"] is not None]
    return {
        "n_coalitions": len(rows),
        "n_with_more_than_one_fixed_point_attractor": sum(r["fixed_attractors"] > 1 for r in rows),
        "n_with_two_or_more_oscillating_trajectories": len(oscillating),
        "fraction_same_across_grid": fraction("same_across_grid", rows),
        "fraction_single_equals_complete": fraction("single_equals_complete", rows),
        "oscillating_subset": {"fraction_same_across_grid": fraction("same_across_grid", oscillating),
                               "fraction_single_equals_complete": fraction("single_equals_complete", oscillating)},
        "oscillating_trajectories_without_period": sum(r["oscillating_without_period"] for r in rows),
        "fixed_largest_within_attractor_distance": extreme("fixed", "within", max),
        "fixed_smallest_between_attractor_distance": extreme("fixed", "between", min),
        "cycle_largest_within_attractor_distance_scaled": extreme("cycle", "within", max),
        "cycle_smallest_between_attractor_distance_scaled": extreme("cycle", "between", min),
        "cycle_largest_relative_period_spread": float(max(spreads)) if spreads else None,
        "cycle_largest_within_attractor_spread_by_summary": {
            key: max((r["cycle_within_attractor_spread"][key] for r in rows), default=None)
            for key in ("cycle_mean", "tail_min", "tail_max")},
    }


def main() -> None:
    table = json.loads((RESULTS / "paper_number_reconciliation.json").read_text())["per_network_table"]
    rng = np.random.default_rng(SEED)
    tasks = [(name, int(c)) for name in sorted(table) for c in panel_coalitions(table[name]["n"], rng)]
    print(f"[{now()}] controls", flush=True)
    control = controls()
    print(f"[{now()}] panel: {len(tasks)} coalitions of {len(table)} networks", flush=True)
    with ProcessPoolExecutor() as pool:
        for i, _ in enumerate(pool.map(run_coalition, tasks, chunksize=4), start=1):
            if i % 64 == 0 or i == len(tasks):
                print(f"[{now()}] {i}/{len(tasks)} coalitions", flush=True)
    rows, failed = [], 0
    for network, coalition in tasks:
        with np.load(CHECKPOINTS / network / f"{coalition:07d}.npz") as z:
            if (z["status"] != STATUSES.index("ok")).any():
                failed += 1
                continue
            rows.append(group(z["klass"], z["tail_means"], z["tail_mins"], z["tail_maxs"], z["cycle_means"], z["period"]))
    panel = pooled(rows)
    gate = {
        "controls_have_their_known_attractors": all(c["passed"] for c in control.values()),
        "panel_can_fail": panel["n_with_more_than_one_fixed_point_attractor"] >= MIN_MULTISTABLE,
        "same_partition_across_grid": panel["fraction_same_across_grid"] >= STABLE_FRACTION,
        "complete_equals_single_linkage": panel["fraction_single_equals_complete"] >= STABLE_FRACTION,
    }
    report = {
        "question": "Is the attractor grouping stable across the registered tolerance grid, and does it find known attractors?",
        "settings": {"seed": SEED, "n_init": N_INIT, "init_seed": INIT_SEED, "coalitions_per_network": COALITIONS_PER_NETWORK,
                     "dynamics": [DYNAMICS.construction, DYNAMICS.hill_n, DYNAMICS.hill_k],
                     "solver_file": str(SOLVER_FILE.relative_to(PROJECT_ROOT)),
                     "classifier_file": str(CLASSIFIER_FILE.relative_to(PROJECT_ROOT)),
                     "fixed_tol_grid": FIXED_TOL_GRID, "cycle_tol_grid": CYCLE_TOL_GRID, "primary": [FIXED_TOL, CYCLE_TOL],
                     "period_rtol": PERIOD_RTOL, "code_sha256": code_digest()},
        "gate_criteria": {"min_multistable": MIN_MULTISTABLE, "stable_fraction": STABLE_FRACTION},
        "controls": control,
        "panel": {**panel, "n_networks": len(table), "n_coalitions_with_a_failed_trajectory": failed},
        "gate": gate,
        "passed": all(gate.values()),
        "written_at": now(),
    }
    OUT.write_text(json.dumps(report, indent=2))
    print(f"[{now()}] wrote {OUT.relative_to(PROJECT_ROOT)}: passed {report['passed']}", flush=True)


if __name__ == "__main__":
    main()
