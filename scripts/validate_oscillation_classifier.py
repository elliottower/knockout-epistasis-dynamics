"""Does the classifier call sustained oscillations oscillatory, and nothing else?

A panel of trajectories whose fate is settled independently, by integrating each to t = 4000
and reading the window [3600, 4000]:

  sustained    ranges above 1e-3 in the window, and the window's two halves repeat their
               per-node envelopes to within 1% of each node's range
  settles      range below 1e-5 in the window
  undetermined anything else, reported and not scored

Panel: 3- and 5-node repressor rings (sustained at high n_H), a scan of the 3-node ring across
its Hopf point (damped just below it, slow to reach the cycle just above it), a 4-node ring
(bistable, settles), the sustained lambda_phage and arabidopsis_cellcycle coalitions and the two
damped arabidopsis_cellcycle coalitions found in calibration. Four initial states per case.

Each trajectory is then run through the engine under the solver and classifier files, and its
half-envelope disagreement is recorded at every horizon, so the margin between the tolerance
and each class is visible. The pass condition: every sustained trajectory is accepted as
oscillatory, and no trajectory that settles is called oscillatory.

    uv run python -m scripts.validate_oscillation_classifier \
        --solver experiments/ode_v2/solver_primary.json \
        --classifier experiments/ode_v2/classifier_primary.json

Each trajectory is computed in its own process and written to its own file under
results/audit/oscillation_classifier_validation/<settings digest>/ as it finishes; a rerun
under the same solver and classifier skips trajectories already on disk, and other settings
get their own directory. Writes the summary to
results/audit/oscillation_classifier_validation_<settings digest>.json.
"""

import argparse
import dataclasses
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.ode_engine import (
    DynamicsConfig,
    clamp_mask_for,
    compile_model,
    initial_states,
    integrate,
    prepare,
    tail_stats,
)
from scripts.paths import AUDIT
from scripts.run_ode_arm import load_classifier, load_solver

CASES = AUDIT / "oscillation_classifier_validation"
SEED = 20260921
N_INIT = 4
T_TRUTH = 4000.0
TRUTH_WINDOW = 400.0
SUSTAINED_RANGE = 1e-3
SUSTAINED_ENVELOPE = 0.01
SETTLED_RANGE = 1e-5
SAMPLES_PER_UNIT = 5

RING3 = {"A": "!C", "B": "!A", "C": "!B"}
RING4 = {"A": "!D", "B": "!A", "C": "!B", "D": "!C"}
RING5 = {"A": "!E", "B": "!A", "C": "!B", "D": "!C", "E": "!D"}
# (label, rules, outputs, coalition or None for all free, n_H values)
PANEL = [
    ("ring3", RING3, ["A"], None, (4.0, 10.0)),
    ("ring5", RING5, ["A"], None, (4.0, 10.0)),
    ("ring3 near Hopf", RING3, ["A"], None, (2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8)),
    ("ring3 low n_H", RING3, ["A"], None, (2.0,)),
    ("ring4", RING4, ["A"], None, (10.0,)),
    ("lambda_phage 107", "lambda_phage", None, 107, (10.0,)),
    ("arabidopsis_cellcycle 3575", "arabidopsis_cellcycle", None, 3575, (2.0,)),
    ("arabidopsis_cellcycle 6270", "arabidopsis_cellcycle", None, 6270, (2.0,)),
    ("arabidopsis_cellcycle 10998", "arabidopsis_cellcycle", None, 10998, (2.0,)),
]


def ground_truth(prep, solver, mask, y0, static_below):
    long = dataclasses.replace(solver, t_max=T_TRUTH, t_tail=TRUTH_WINDOW,
                               n_tail_samples=int(TRUTH_WINDOW * SAMPLES_PER_UNIT))
    r = integrate(prep, long, mask, 0.0, y0, keep_tail=True)
    if r.status != "ok":
        return {"fate": "undetermined", "reason": r.status}
    stats = tail_stats(r.tail.astype(float), r.stats.final_derivative, 1.0, static_below)
    if stats.tail_range < SETTLED_RANGE:
        fate = "settles"
    elif stats.tail_range > SUSTAINED_RANGE and stats.halves_envelope_diff <= SUSTAINED_ENVELOPE:
        fate = "sustained"
    else:
        fate = "undetermined"
    return {"fate": fate, "window_range": stats.tail_range, "window_envelope_diff": stats.halves_envelope_diff}


def per_horizon(prep, solver, classifier, mask, y0):
    rows = []
    for horizon in (solver.t_max,) + classifier.extension_horizons:
        scale = horizon / solver.t_max
        fixed_horizon = dataclasses.replace(solver, t_max=horizon, t_tail=solver.t_tail * scale,
                                            n_tail_samples=int(round(solver.n_tail_samples * scale)))
        r = integrate(prep, fixed_horizon, mask, 0.0, y0, keep_tail=True)
        if r.status != "ok":
            rows.append({"horizon": horizon, "status": r.status})
            continue
        s = tail_stats(r.tail.astype(float), r.stats.final_derivative, classifier.bounded_margin,
                       classifier.fixed_tail_range_max)
        rows.append({"horizon": horizon, "tail_range": s.tail_range, "final_derivative": s.final_derivative,
                     "halves_envelope_diff": s.halves_envelope_diff})
    return rows


def settings_digest(solver_path: Path, classifier_path: Path) -> str:
    settings = {"solver": json.loads(solver_path.read_text()), "classifier": json.loads(classifier_path.read_text())}
    return hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()[:12]


def run_case(job: tuple) -> dict:
    label, rules, outputs, coalition, hill_n, i, solver_path, classifier_path = job
    path = CASES / settings_digest(solver_path, classifier_path) / f"{label.replace(' ', '_')}__n{hill_n:g}__init{i}.json"
    if path.exists():
        return json.loads(path.read_text())
    solver, classifier = load_solver(solver_path), load_classifier(classifier_path)
    if isinstance(rules, str):
        info = ALL_MODELS[rules]
        net = compile_model(info["rules"], info["output_nodes"])
    else:
        net = compile_model(rules, outputs)
    n = len(net.node_names)
    mask = clamp_mask_for(coalition, n) if coalition is not None else np.zeros(n, dtype=bool)
    prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=hill_n))
    y0 = np.where(mask, 0.0, initial_states(N_INIT, n, SEED)[i])
    truth = ground_truth(prep, solver, mask, y0, classifier.fixed_tail_range_max)
    r = integrate(prep, solver, mask, 0.0, y0, classifier)
    row = {"case": label, "hill_n": hill_n, "initial_state": i, **truth,
           "status": r.status, "class": r.klass, "horizon_index": r.horizon_index,
           "per_horizon": per_horizon(prep, solver, classifier, mask, y0),
           "solver": json.loads(solver_path.read_text()), "classifier": json.loads(classifier_path.read_text())}
    tmp = path.with_name(path.stem + ".partial.json")
    tmp.write_text(json.dumps(row, indent=2))
    tmp.replace(path)
    print(f"{label:<30} n_H={hill_n:<4g} init {i}  truth {truth['fate']:<12} "
          f"engine {r.status}/{r.klass} at horizon {r.horizon_index}", flush=True)
    return row


def summarize(rows: list[dict], classifier) -> dict:
    """Counts, margins and the pass gate: every sustained trajectory accepted as oscillatory, and
    no settling trajectory called oscillatory."""
    def envelopes(fate, moving_only):
        out = []
        for row in rows:
            if row["fate"] != fate:
                continue
            for h in row["per_horizon"]:
                if "halves_envelope_diff" in h and (not moving_only or h["tail_range"] > classifier.oscillatory_tail_range_min):
                    out.append(h["halves_envelope_diff"])
        return out

    sustained_env = envelopes("sustained", moving_only=True)
    settling_env = envelopes("settles", moving_only=True)
    summary = {
        "confusion": {fate: {f"{r['status']}/{r['class']}": sum(1 for x in rows if x["fate"] == fate
                                                                 and (x["status"], x["class"]) == (r["status"], r["class"]))
                             for r in rows if r["fate"] == fate}
                      for fate in ("sustained", "settles", "undetermined")},
        "n_sustained": sum(1 for r in rows if r["fate"] == "sustained"),
        "sustained_called_oscillatory": sum(1 for r in rows if r["fate"] == "sustained"
                                            and (r["status"], r["class"]) == ("ok", "oscillatory")),
        "settling_called_oscillatory": sum(1 for r in rows if r["fate"] == "settles" and r["class"] == "oscillatory"),
        "sustained_called_fixed": sum(1 for r in rows if r["fate"] == "sustained" and r["class"] == "fixed"),
        "sustained_not_accepted": sum(1 for r in rows if r["fate"] == "sustained" and r["status"] != "ok"),
        "envelope_tol": classifier.halves_envelope_tol,
        "largest_envelope_diff_sustained_any_horizon": max(sustained_env) if sustained_env else None,
        "smallest_envelope_diff_settling_while_moving": min(settling_env) if settling_env else None,
    }
    summary["passed"] = (summary["n_sustained"] > 0
                         and summary["sustained_called_oscillatory"] == summary["n_sustained"]
                         and summary["settling_called_oscillatory"] == 0)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--classifier", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    classifier = load_classifier(args.classifier)
    (CASES / settings_digest(args.solver, args.classifier)).mkdir(parents=True, exist_ok=True)
    jobs = [(label, rules, outputs, coalition, hill_n, i, args.solver, args.classifier)
            for label, rules, outputs, coalition, hill_ns in PANEL for hill_n in hill_ns for i in range(N_INIT)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(run_case, jobs))
    stale = [r["case"] for r in rows if r["solver"] != json.loads(args.solver.read_text())
             or r["classifier"] != json.loads(args.classifier.read_text())]
    if stale:
        raise ValueError(f"case files on disk were computed under other settings: {sorted(set(stale))}")

    summary = summarize(rows, classifier)
    out = AUDIT / f"oscillation_classifier_validation_{settings_digest(args.solver, args.classifier)}.json"
    out.write_text(json.dumps({
        "solver_file": str(args.solver), "solver": json.loads(args.solver.read_text()),
        "classifier_file": str(args.classifier), "classifier": json.loads(args.classifier.read_text()),
        "ground_truth": {"t_end": T_TRUTH, "window": TRUTH_WINDOW, "sustained_range": SUSTAINED_RANGE,
                         "sustained_envelope": SUSTAINED_ENVELOPE, "settled_range": SETTLED_RANGE},
        "seed": SEED, "n_init": N_INIT, "summary": summary, "trajectories": rows}, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
