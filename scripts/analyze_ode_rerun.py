"""Evaluate the ODE rerun registrations from their run records.

    uv run python -m scripts.analyze_ode_rerun primary
    uv run python -m scripts.analyze_ode_rerun graded
    uv run python -m scripts.analyze_ode_rerun sensitivity

Each arm reads the finalized records its registration names, copied from the Modal volume into
`results/<run name>/` next to that registration, and the Boolean estimates written there by
`scripts/boolean_estimators.py`. A record is read only if its launch manifest names the commit
the freeze tag points to. Each arm writes `results/analysis.json` beside its registration; the
primary arm also writes the O1 and A1 scatterplots to `results/figures/`.
Frozen with the registrations: experiments/2026-09-21_ode-primary-rerun/PREREG.md,
experiments/2026-09-21_graded-perturbation-rerun/PREREG.md and
experiments/2026-09-21_ode-sensitivity/PREREG.md.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure
from scipy.stats import binomtest, pearsonr, rankdata, spearmanr

from scripts.attractors import CYCLE_TOL_GRID, FIXED_TOL_GRID, grid_key
from scripts.paths import PROJECT_ROOT, RESULTS
from scripts.record_registered_identities import OUT as IDENTITIES
from scripts.verify_freeze import FREEZE_TAG, tag_commit
from scripts.walsh_estimators import MIN_RATIO_FRACTION, MIN_TOTAL_ENERGY

EXPERIMENTS = PROJECT_ROOT / "experiments"
PRIMARY = EXPERIMENTS / "2026-09-21_ode-primary-rerun"
GRADED = EXPERIMENTS / "2026-09-21_graded-perturbation-rerun"
SENSITIVITY = EXPERIMENTS / "2026-09-21_ode-sensitivity"
BOOLEAN = PRIMARY / "results" / "boolean_estimators.json"
NULL_BAND_PP = 0.5
ESTIMATORS = ("all_pairs", "split_half", "plain")  # all_pairs is primary
B1_ESTIMATORS = ("plain", "all_pairs")  # the estimators B1 compares; split-half is reported alongside
TIE = 1e-12
GRADED_LEVELS = (0.0, 0.25, 0.5, 0.75, 1.0)
GATE_ABSOLUTE_O3PLUS = 1e-4       # the v1 amendment's wording, applied only as a sensitivity analysis
MIN_G4_NETWORKS = 10
MAX_NULL_FRACTION_PARTIAL = 0.25
MIN_NETWORKS_WITH_VARIATION = 5
# Published Boolean values carry two decimals, so a recomputation reproduces one if it lies within
# half the last digit.
PUBLISHED_PP_TOL = 0.005
BOOTSTRAP_RESAMPLES = 10_000
BOOTSTRAP_SEED = 0
# The admissibility sensitivity analysis: the registered gates and one decade either side.
TOTAL_ENERGY_GRID = (MIN_TOTAL_ENERGY / 10, MIN_TOTAL_ENERGY, MIN_TOTAL_ENERGY * 10)
RATIO_FRACTION_GRID = (MIN_RATIO_FRACTION / 10, MIN_RATIO_FRACTION, MIN_RATIO_FRACTION * 10)
SETTINGS = {"a": "sensitivity-a", "b": "sensitivity-b", "c": "sensitivity-c", "d": "sensitivity-d", "e": "sensitivity-e"}


class AnalysisError(Exception):
    pass


def classify(delta_pp: float) -> str:
    if delta_pp > NULL_BAND_PP:
        return "creation"
    if delta_pp < -NULL_BAND_PP:
        return "destruction"
    return "null"


def frozen_commit() -> str:
    return tag_commit(PROJECT_ROOT, FREEZE_TAG)


def load_records(run_dir: Path, commit: str) -> dict[str, dict]:
    """Every finalized record in a run directory. Refuses a record not launched from `commit`."""
    records = {}
    for f in sorted(run_dir.glob("*.json")):
        if f.name.endswith(".partial.json") or f.name.endswith(".containers.json"):
            continue
        record = json.loads(f.read_text())
        launched_from = ((record.get("launch_manifest") or {}).get("freeze") or {}).get("commit")
        if launched_from != commit:
            raise AnalysisError(f"{f} was launched from {launched_from}, not from the frozen commit {commit}")
        if record["network"] in records:
            raise AnalysisError(f"two records for {record['network']} in {run_dir}")
        records[record["network"]] = record
    return records


def replicate_mismatches(record: dict) -> int:
    rep = record["replicate"]
    return sum(rep.get(k, 0) for k in ("status_mismatches", "class_mismatches", "value_mismatches"))


def scoreable(record: dict | None) -> tuple[bool, str]:
    """Whether a record can enter any test, whatever the estimator."""
    if record is None:
        return False, "no record"
    if not record["scored"]:
        return False, "cannot be scored"
    if replicate_mismatches(record):
        return False, "withheld for a replicate mismatch"
    return True, ""


def usable(record: dict | None, estimator: str) -> tuple[bool, str]:
    """Whether a record enters the registered tests under an estimator, and why not if it does not."""
    ok, why = scoreable(record)
    if not ok:
        return ok, why
    if not record["scores"]["estimators"][estimator]["admissible"]:
        return False, f"inadmissible under {estimator}"
    return True, ""


def ode_delta_pp(record: dict, estimator: str) -> float:
    return record["scores"]["estimators"][estimator]["delta_o3plus"] * 100


def rescored_pp(estimate: dict, local_o3plus: float, min_total: float) -> float | None:
    """Delta_3+ in pp recomputed from the raw order energies under another total-energy gate."""
    energy = np.asarray(estimate["energy_by_order"], dtype=float)
    total = float(energy.sum())
    return (float(energy[3:].sum()) / total - local_o3plus) * 100 if total >= min_total else None


def ratio_under(estimate: dict, min_ratio: float) -> float | None:
    """sigma_3+/sigma_2+ under another ratio gate, at the registered total-energy gate."""
    total, e3, e2 = estimate["total_energy"], estimate["energy_3plus"], estimate["energy_2plus"]
    if total >= MIN_TOTAL_ENERGY and e3 >= min_ratio * total and e2 >= min_ratio * total:
        return e3 / e2
    return None


def boolean_rows() -> dict[str, dict]:
    return json.loads(BOOLEAN.read_text())["networks"]


def boolean_deltas(estimator: str) -> dict[str, float]:
    """Boolean Delta_3+ in pp under an estimator, for the networks where it is admissible."""
    return {name: row[estimator]["delta_o3plus"] * 100 for name, row in boolean_rows().items()
            if row[estimator]["admissible"]}


def compare(ode: dict[str, dict], boolean: dict[str, float], networks: list[str], estimator: str) -> dict:
    """Class preservation, correlation, median difference and the magnitude sign test."""
    rows, excluded = {}, {}
    for name in networks:
        ok, why = usable(ode.get(name), estimator)
        if not ok or name not in boolean:
            excluded[name] = why or f"Boolean inadmissible under {estimator}"
            continue
        d_ode, d_bool = ode_delta_pp(ode[name], estimator), boolean[name]
        rows[name] = {"boolean_pp": d_bool, "ode_pp": d_ode, "boolean_class": classify(d_bool),
                      "ode_class": classify(d_ode)}
    non_null = [n for n, r in rows.items() if r["boolean_class"] != "null"]
    preserved = [n for n in non_null if rows[n]["ode_class"] == rows[n]["boolean_class"]]
    changed = sorted(set(non_null) - set(preserved))
    b = np.array([rows[n]["boolean_pp"] for n in rows])
    o = np.array([rows[n]["ode_pp"] for n in rows])
    diffs = [abs(rows[n]["ode_pp"] - rows[n]["boolean_pp"]) for n in non_null]
    larger = [n for n in non_null if abs(rows[n]["ode_pp"]) - abs(rows[n]["boolean_pp"]) > TIE]
    smaller = [n for n in non_null if abs(rows[n]["boolean_pp"]) - abs(rows[n]["ode_pp"]) > TIE]
    trials = len(larger) + len(smaller)
    return {
        "rows": rows, "excluded": excluded, "n_scored": len(rows), "n_non_null": len(non_null),
        "n_preserved": len(preserved), "changed_class": changed,
        "pearson_r": float(pearsonr(b, o).statistic) if len(rows) >= 3 else None,
        "spearman_rho": float(spearmanr(b, o).statistic) if len(rows) >= 3 else None,
        "median_abs_diff_pp": float(np.median(diffs)) if diffs else None,
        "n_ode_larger": len(larger), "n_ode_smaller": len(smaller), "n_ties": len(non_null) - trials,
        "sign_test_p": float(binomtest(len(larger), trials, 0.5, alternative="greater").pvalue) if trials else None,
    }


def preservation_under_gate(ode: dict[str, dict], networks: list[str], estimator: str, min_total: float,
                            max_excluded: int) -> dict:
    """Class preservation with both spectra rescored under another total-energy gate."""
    boolean = boolean_rows()
    rows, excluded = {}, {}
    for name in networks:
        ok, why = scoreable(ode.get(name))
        if not ok:
            excluded[name] = why
            continue
        d_ode = rescored_pp(ode[name]["scores"]["estimators"][estimator], ode[name]["scores"]["local_o3plus"], min_total)
        d_bool = (rescored_pp(boolean[name][estimator], boolean[name]["local_o3plus"], min_total)
                  if name in boolean else None)
        if d_ode is None or d_bool is None:
            excluded[name] = f"inadmissible at total energy {min_total:g}"
            continue
        rows[name] = (classify(d_bool), classify(d_ode))
    non_null = [n for n, (b, _) in rows.items() if b != "null"]
    changed = sorted(n for n in non_null if rows[n][0] != rows[n][1])
    void = len(excluded) > max_excluded
    return {"min_total_energy": min_total, "excluded": excluded, "n_non_null": len(non_null),
            "changed_class": changed, "void": void, "all_preserved": None if void else not changed}


def evaluate_primary(ode: dict[str, dict], networks: list[str]) -> dict:
    out = {}
    for estimator in ESTIMATORS:
        c = compare(ode, boolean_deltas(estimator), networks, estimator)
        void = len(c["excluded"]) > 3
        c["hypotheses"] = {
            "void": void,
            "H1": None if void else c["n_preserved"] == c["n_non_null"],
            "H2": None if void or c["pearson_r"] is None else c["pearson_r"] >= 0.90,
            "H3": None if void or c["median_abs_diff_pp"] is None else c["median_abs_diff_pp"] <= 2.0,
            "H4": None if void or c["sign_test_p"] is None else c["sign_test_p"] < 0.05,
        }
        c["admissibility_sensitivity"] = [preservation_under_gate(ode, networks, estimator, t, 3)
                                          for t in TOTAL_ENERGY_GRID]
        out[estimator] = c
    return out


def boolean_estimator_check() -> dict:
    """B1, over the registered 28 networks and their registered published values.

    Void unless every registered network has an admissible plain and all-pairs estimate and the
    plain estimate reproduces its published value to within PUBLISHED_PP_TOL and in class. B1
    holds if no all-pairs class differs from the published class. Split-half is reported where it
    is admissible and decides nothing."""
    registered = {n: r["boolean_delta_3plus_pp_published"] for n, r in json.loads(IDENTITIES.read_text())["networks"].items()}
    data = boolean_rows()
    rows, missing = {}, []
    for name, published in registered.items():
        row = data.get(name)
        if row is None or not all(row[e]["admissible"] for e in B1_ESTIMATORS):
            missing.append(name)
            continue
        values = {e: row[e]["delta_o3plus"] * 100 if row[e]["admissible"] else None for e in ESTIMATORS}
        rows[name] = {"published_pp": published, "published_class": classify(published),
                      **{f"{e}_pp": v for e, v in values.items()},
                      **{f"{e}_class": classify(v) if v is not None else None for e, v in values.items()}}
    not_reproduced = sorted(n for n, r in rows.items()
                            if abs(r["plain_pp"] - r["published_pp"]) > PUBLISHED_PP_TOL or r["plain_class"] != r["published_class"])
    void = bool(missing or not_reproduced)
    moved = {e: sorted(n for n, r in rows.items() if r[f"{e}_class"] is not None and r[f"{e}_class"] != r["published_class"])
             for e in ("all_pairs", "split_half")}
    return {"n_registered": len(registered), "rows": rows, "missing_or_inadmissible": sorted(missing),
            "plain_does_not_reproduce_publication": not_reproduced, "void": void,
            "all_pairs_class_differs_from_published": moved["all_pairs"],
            "split_half_class_differs_from_published": moved["split_half"],
            "split_half_inadmissible": sorted(n for n, r in rows.items() if r["split_half_class"] is None),
            "B1": None if void else not moved["all_pairs"]}


def spearman_interval(x: np.ndarray, y: np.ndarray) -> dict:
    """Percentile interval for Spearman's rho from BOOTSTRAP_RESAMPLES network-level resamples."""
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    idx = rng.integers(0, x.size, size=(BOOTSTRAP_RESAMPLES, x.size))
    xr, yr = rankdata(x[idx], axis=1), rankdata(y[idx], axis=1)
    xr -= xr.mean(axis=1, keepdims=True)
    yr -= yr.mean(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = (xr * yr).sum(axis=1) / np.sqrt((xr ** 2).sum(axis=1) * (yr ** 2).sum(axis=1))
    finite = rho[np.isfinite(rho)]
    return {"low": float(np.percentile(finite, 2.5)) if finite.size else None,
            "high": float(np.percentile(finite, 97.5)) if finite.size else None,
            "n_resamples": BOOTSTRAP_RESAMPLES, "n_undefined": int(rho.size - finite.size), "seed": BOOTSTRAP_SEED}


def association(ode: dict[str, dict], networks: list[str], feature, alternative: str, estimator: str) -> dict:
    """One-sided Spearman test of a per-network feature against ODE Delta_3+ over scored networks.

    A network that cannot be scored, or whose record lacks the feature, is excluded by name, and
    counts toward the void rule. Void if more than 3 networks are excluded, or fewer than
    MIN_NETWORKS_WITH_VARIATION have a nonzero feature."""
    rows, excluded = {}, {}
    for name in networks:
        ok, why = usable(ode.get(name), estimator)
        value = feature(ode[name]) if ok else None
        if not ok:
            excluded[name] = why
        elif value is None:
            excluded[name] = f"no {feature.__name__}"
        else:
            rows[name] = {"feature": value, "ode_pp": ode_delta_pp(ode[name], estimator)}
    x = np.array([r["feature"] for r in rows.values()])
    y = np.array([r["ode_pp"] for r in rows.values()])
    nonzero = int((x > 0).sum())
    void = len(excluded) > 3 or nonzero < MIN_NETWORKS_WITH_VARIATION
    test = spearmanr(x, y, alternative=alternative) if not void else None
    return {"rows": rows, "excluded": excluded, "n": len(rows), "n_nonzero_feature": nonzero,
            "ties": {"feature": int(x.size - np.unique(x).size), "ode_pp": int(y.size - np.unique(y).size)},
            "void": void, "rho": float(test.statistic) if test else None,
            "p_one_sided": float(test.pvalue) if test else None,
            "rho_interval_95": spearman_interval(x, y) if not void else None,
            "holds": None if void else bool(test.pvalue < 0.05)}


def oscillation_fraction(record: dict) -> float | None:
    classes = (record.get("summary") or {}).get("classes_of_accepted")
    if not classes:
        return None
    accepted = sum(classes.values())
    return classes["oscillatory"] / accepted if accepted else None


def mean_basin_entropy(record: dict) -> float | None:
    return record["basins"]["mean_entropy_bits"] if "basins" in record else None


def entropy_at(fixed_tol: float, cycle_tol: float):
    key = grid_key(fixed_tol, cycle_tol)

    def mean_basin_entropy_at(record: dict) -> float | None:
        return record["basins"]["sensitivity"][key]["mean_entropy_bits"] if "basins" in record else None
    return mean_basin_entropy_at


def holm(p: dict[str, float | None]) -> dict[str, float | None]:
    """Holm-adjusted p-values over the tests that are not void."""
    ranked = sorted((v, k) for k, v in p.items() if v is not None)
    out, running = {k: None for k in p}, 0.0
    for i, (v, k) in enumerate(ranked):
        running = max(running, min(1.0, (len(ranked) - i) * v))
        out[k] = running
    return out


def mechanism_tests(ode: dict[str, dict], networks: list[str]) -> dict:
    """O1: oscillation fraction falls as Delta_3+ rises. A1: basin entropy rises with Delta_3+.
    Each is interpreted on its own p-value; the Holm-adjusted pair is reported alongside. A1 is
    also recomputed at every point of the attractor tolerance grid."""
    out = {}
    for estimator in ESTIMATORS:
        o1 = association(ode, networks, oscillation_fraction, "less", estimator)
        a1 = association(ode, networks, mean_basin_entropy, "greater", estimator)
        out[estimator] = {
            "O1": o1, "A1": a1,
            "holm_adjusted_p": holm({"O1": o1["p_one_sided"], "A1": a1["p_one_sided"]}),
            "A1_tolerance_sensitivity": {grid_key(ft, ct): association(ode, networks, entropy_at(ft, ct), "greater", estimator)
                                         for ft in FIXED_TOL_GRID for ct in CYCLE_TOL_GRID},
        }
    return out


def basin_table(ode: dict[str, dict], networks: list[str]) -> dict[str, dict]:
    """Each network's basin statistics at every point of the tolerance grid: mean entropy, the
    fraction of multistable coalitions, the largest attractor count and the count histogram."""
    return {n: ode[n]["basins"]["sensitivity"] for n in networks if n in ode and "basins" in ode[n]}


FEATURE_LABELS = {"O1": "fraction of trajectories classified oscillatory", "A1": "mean basin entropy (bits)"}


def scatterplot(assoc: dict, test: str, estimator: str, path: Path) -> Path | None:
    """A mechanism test's feature against ODE Delta_3+ over its usable networks, whatever the
    outcome, void included. None if no network is usable. The file carries no timestamp or
    version, so the same rows give the same bytes."""
    rows = assoc["rows"]
    if not rows:
        return None
    names = sorted(rows)
    fig = Figure(figsize=(4.5, 3.5))
    ax = fig.subplots()
    ax.axhspan(-NULL_BAND_PP, NULL_BAND_PP, color="0.9", zorder=0)
    ax.scatter([rows[n]["feature"] for n in names], [rows[n]["ode_pp"] for n in names], s=16, color="black", zorder=2)
    ax.set_xlabel(FEATURE_LABELS[test])
    ax.set_ylabel("ODE Δ3+ (pp)")
    outcome = "void" if assoc["void"] else f"ρ = {assoc['rho']:+.2f}, one-sided p = {assoc['p_one_sided']:.2g}"
    ax.set_title(f"{test}, {estimator}: n = {assoc['n']}, {outcome}", fontsize=9)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, metadata={"Software": None})
    return path


def scatterplots(mechanisms: dict, directory: Path) -> dict[str, dict[str, str | None]]:
    """O1 and A1 under every estimator, as `directory/<test>_<estimator>.png`."""
    out = {}
    for estimator, tests in mechanisms.items():
        out[estimator] = {}
        for test in ("O1", "A1"):
            path = scatterplot(tests[test], test, estimator, directory / f"{test}_{estimator}.png")
            out[estimator][test] = str(path.relative_to(PROJECT_ROOT)) if path else None
    return out


def published_ode_pp() -> dict[str, float]:
    """The published ODE arm's Delta_3+, from the files the paper's numbers were reconciled against."""
    out = {}
    for pattern in ("results/ode_full/*_ode.json", "results/grn_v2/ode_full/*_ode.json"):
        for f in sorted(PROJECT_ROOT.glob(pattern)):
            d = json.loads(f.read_text())
            if d.get("ode_delta_o3plus") is not None:
                out.setdefault(d["model"], d["ode_delta_o3plus"] * 100)
    return out


def legacy_audit(legacy: dict[str, dict], networks: list[str]) -> dict:
    published = published_ode_pp()
    rows, unscored = {}, {}
    for name in networks:
        ok, why = usable(legacy.get(name), "plain")
        if not ok:
            unscored[name] = why
            continue
        new = ode_delta_pp(legacy[name], "plain")
        rows[name] = {"published_pp": published.get(name), "without_timeout_pp": new,
                      "change_pp": new - published[name] if name in published else None,
                      "class_changed": name in published and classify(new) != classify(published[name])}
    return {"rows": rows, "unscored": unscored,
            "n_moved_over_0p5pp": sum(1 for r in rows.values() if r["change_pp"] is not None and abs(r["change_pp"]) > 0.5),
            "class_changed": sorted(n for n, r in rows.items() if r["class_changed"])}


def oscillation(records: dict[str, dict]) -> dict:
    out = {}
    for name, r in records.items():
        s = r["summary"]
        accepted = sum(s["classes_of_accepted"].values())
        out[name] = {"oscillatory_trajectory_fraction": s["classes_of_accepted"]["oscillatory"] / accepted if accepted else None,
                     "coalitions_with_an_oscillating_trajectory_fraction":
                     s["coalitions_with_an_oscillating_trajectory"] / s["n_coalitions"]}
    return out


def provisional(records: dict[str, dict]) -> dict:
    return {name: {e: (r["scores_provisional"]["estimators"][e]["delta_o3plus"] * 100
                       if r["scores_provisional"]["estimators"][e]["admissible"] else None) for e in ESTIMATORS}
            for name, r in records.items() if "scores_provisional" in r}


def g4(ratios_f0: dict[str, float], ratios_f075: dict[str, float]) -> dict:
    """G4 over the networks with a ratio at both levels; undecided below MIN_G4_NETWORKS."""
    both = sorted(set(ratios_f0) & set(ratios_f075))
    m0 = float(np.median([ratios_f0[n] for n in both])) if both else None
    m75 = float(np.median([ratios_f075[n] for n in both])) if both else None
    enough = len(both) >= MIN_G4_NETWORKS
    return {"n_admissible_at_both": len(both), "median_ratio_f0": m0, "median_ratio_f075": m75,
            "enough_networks": enough, "holds": abs(m0 - m75) <= 0.05 if enough else None}


def evaluate_graded(by_level: dict[float, dict[str, dict]], networks: list[str]) -> dict:
    out = {}
    partial = GRADED_LEVELS[1:]
    for estimator in ESTIMATORS:
        levels = {}
        for f, records in by_level.items():
            counts = {"creation": 0, "destruction": 0, "null": 0}
            classes, ratios, gated, absolute_ratios, excluded, deltas = {}, {}, {}, {}, {}, []
            for name in networks:
                ok, why = usable(records.get(name), estimator)
                if not ok:
                    excluded[name] = why
                    continue
                est = records[name]["scores"]["estimators"][estimator]
                d = est["delta_o3plus"] * 100
                classes[name] = classify(d)
                counts[classes[name]] += 1
                deltas.append(d)
                parts = {"energy_3plus": est["energy_3plus"], "energy_2plus": est["energy_2plus"],
                         "total_energy": est["total_energy"]}
                if est["ratio_admissible"]:
                    ratios[name] = {"ratio": est["higher_order_ratio"], **parts}
                    if est["energy_3plus"] >= GATE_ABSOLUTE_O3PLUS:
                        absolute_ratios[name] = est["higher_order_ratio"]
                else:
                    gated[name] = {"reason": est["ratio_reason"], **parts}
            n_scored = sum(counts.values())
            levels[f] = {"counts": counts, "n_scored": n_scored, "classes": classes, "ratios": ratios,
                         "gated": gated, "excluded": excluded,
                         "median_delta_pp": float(np.median(deltas)) if deltas else None,
                         "null_fraction": counts["null"] / n_scored if n_scored else None,
                         "void": len(excluded) > 3, "absolute_gate_ratios": absolute_ratios}
        primary_g4 = g4({n: r["ratio"] for n, r in levels[0.0]["ratios"].items()},
                        {n: r["ratio"] for n, r in levels[0.75]["ratios"].items()})
        g4_void = levels[0.0]["void"] or levels[0.75]["void"] or not primary_g4["enough_networks"]
        g5_void = any(levels[f]["void"] for f in partial)
        inverted = [f for f in GRADED_LEVELS if levels[f]["counts"]["destruction"] > levels[f]["counts"]["creation"]]
        changes = {}
        for name in levels[0.0]["classes"]:
            later = [f for f in partial if name in levels[f]["classes"]
                     and levels[f]["classes"][name] != levels[0.0]["classes"][name]]
            changes[name] = later[0] if later else None
        common = sorted(set.intersection(*(set(levels[f]["classes"]) for f in GRADED_LEVELS)))
        common_null = {str(f): sum(levels[f]["classes"][n] == "null" for n in common) for f in partial}
        ratio_gate = {}
        for r in RATIO_FRACTION_GRID:
            at_level = []
            for f in (0.0, 0.75):
                values = {n: ratio_under(by_level[f][n]["scores"]["estimators"][estimator], r)
                          for n in networks if usable(by_level[f].get(n), estimator)[0]}
                at_level.append({n: v for n, v in values.items() if v is not None})
            ratio_gate[f"{r:g}"] = g4(*at_level)
        m0, m75 = primary_g4["median_ratio_f0"], primary_g4["median_ratio_f075"]
        out[estimator] = {
            "levels": {str(f): v for f, v in levels.items()},
            **{k: v for k, v in primary_g4.items() if k != "holds"},
            "v1_factor_f0_over_f075": m0 / m75 if m0 is not None and m75 else None,
            "absolute_gate_sensitivity": g4(levels[0.0]["absolute_gate_ratios"], levels[0.75]["absolute_gate_ratios"]),
            "ratio_gate_sensitivity": ratio_gate,
            "G5_on_networks_scored_at_every_level": {
                "networks": common, "n": len(common), "null_counts": common_null,
                "holds": all(c <= MAX_NULL_FRACTION_PARTIAL * len(common) for c in common_null.values()) if common else None},
            "lowest_level_destruction_outnumbers_creation": inverted[0] if inverted else None,
            "first_class_change_by_network": changes,
            "hypotheses": {
                "G1": None if levels[0.0]["void"] else levels[0.0]["counts"]["creation"] > levels[0.0]["counts"]["destruction"],
                "G2": None if levels[0.5]["void"] else levels[0.5]["counts"]["destruction"] > levels[0.5]["counts"]["creation"],
                "G3": None if levels[0.25]["void"] else levels[0.25]["counts"]["null"] <= 2,
                "G4": None if g4_void else primary_g4["holds"],
                "G4_void_because": None if not g4_void else (
                    f"only {primary_g4['n_admissible_at_both']} networks admissible at both f = 0 and f = 0.75"
                    if not primary_g4["enough_networks"] else "a level has more than 3 networks that cannot be scored"),
                "G5": None if g5_void else all(levels[f]["null_fraction"] <= MAX_NULL_FRACTION_PARTIAL for f in partial),
            },
        }
    return out


def evaluate_sensitivity(by_setting: dict[str, dict[str, dict]], networks: list[str]) -> dict:
    out = {}
    for estimator in ESTIMATORS:
        boolean = boolean_deltas(estimator)
        settings = {}
        for key, records in by_setting.items():
            c = compare(records, boolean, networks, estimator)
            c["void"] = len(c["excluded"]) > 2
            c["admissibility_sensitivity"] = [preservation_under_gate(records, networks, estimator, t, 2)
                                              for t in TOTAL_ENERGY_GRID]
            settings[key] = c
        s1_parts = [None if settings[k]["void"] else settings[k]["n_preserved"] == settings[k]["n_non_null"] for k in "abcd"]
        out[estimator] = {
            "settings": settings,
            "hypotheses": {
                "S1": None if None in s1_parts else all(s1_parts),
                "S1_by_setting": dict(zip("abcd", s1_parts)),
                "S2": None if settings["e"]["void"] else settings["e"]["n_preserved"] == settings["e"]["n_non_null"],
            },
        }
    return out


def paper_networks() -> dict[str, int]:
    table = json.loads((RESULTS / "paper_number_reconciliation.json").read_text())["per_network_table"]
    return {name: table[name]["n"] for name in sorted(table, key=lambda m: (table[m]["n"], m))}


def write(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".partial.json")
    tmp.write_text(json.dumps(record, indent=2, default=float))
    tmp.replace(path)
    print(f"wrote {path.relative_to(PROJECT_ROOT)}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("arm", choices=("primary", "graded", "sensitivity"))
    args = parser.parse_args(argv)
    sizes = paper_networks()
    networks = list(sizes)
    commit = frozen_commit()
    primary = load_records(PRIMARY / "results" / "primary-hillcube-n10", commit)
    if args.arm == "primary":
        mechanisms = mechanism_tests(primary, networks)
        write(PRIMARY / "results" / "analysis.json", {
            "registration": "experiments/2026-09-21_ode-primary-rerun/PREREG.md",
            "frozen_commit": commit,
            "tests": evaluate_primary(primary, networks),
            "boolean_estimator_check": boolean_estimator_check(),
            "mechanism_tests": mechanisms,
            "scatterplots": scatterplots(mechanisms, PRIMARY / "results" / "figures"),
            "basins_by_network_and_tolerance": basin_table(primary, networks),
            "legacy_audit": legacy_audit(load_records(PRIMARY / "results" / "legacy-audit", commit), networks),
            "oscillation": oscillation(primary),
            "provisional_sensitivity": provisional(primary),
        })
    elif args.arm == "graded":
        by_level = {0.0: primary}
        for f in GRADED_LEVELS[1:]:
            by_level[f] = load_records(GRADED / "results" / f"graded-f{f:g}", commit)
        write(GRADED / "results" / "analysis.json", {
            "registration": "experiments/2026-09-21_graded-perturbation-rerun/PREREG.md",
            "frozen_commit": commit,
            "tests": evaluate_graded(by_level, networks),
            "provisional_sensitivity": {str(f): provisional(r) for f, r in by_level.items()},
        })
    else:
        small = [n for n in networks if sizes[n] <= 12]
        by_setting = {k: load_records(SENSITIVITY / "results" / run, commit) for k, run in SETTINGS.items()}
        write(SENSITIVITY / "results" / "analysis.json", {
            "registration": "experiments/2026-09-21_ode-sensitivity/PREREG.md",
            "frozen_commit": commit,
            "networks": small,
            "reference": {e: compare(primary, boolean_deltas(e), small, e) for e in ESTIMATORS},
            "tests": evaluate_sensitivity(by_setting, small),
            "provisional_sensitivity": {k: provisional(r) for k, r in by_setting.items()},
        })


if __name__ == "__main__":
    main()
