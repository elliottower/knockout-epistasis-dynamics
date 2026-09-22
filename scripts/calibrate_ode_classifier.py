"""Where do the classifier's diagnostics fall on real networks, relative to its thresholds?

Integrates a small deterministic sample of coalitions under a solver file and a classifier
file, and records for every trajectory its status, class, horizon, tail range, final
derivative and half-envelope disagreement. Planning only: no coalition value is kept and
nothing is scored. The question it answers is whether fixed points sit clear of the
solver's noise floor and whether oscillations pass the envelope test.

    uv run python -m scripts.calibrate_ode_classifier \
        --solver experiments/ode_v2/solver_primary.json \
        --classifier experiments/ode_v2/classifier_primary.json

Writes results/audit/classifier_calibration.json.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.ode_engine import (
    STATUSES,
    TRAJECTORY_CLASSES,
    DynamicsConfig,
    compile_model,
    initial_states,
    prepare,
    simulate_coalition,
)
from scripts.paths import AUDIT
from scripts.run_ode_arm import load_classifier, load_solver

OUT = AUDIT / "classifier_calibration.json"
# Small to large, and the two with the most Boolean cycling.
NETWORKS = ("lambda_phage", "davidich_yeast", "faure_cellcycle", "arabidopsis_cellcycle", "fanconi_anemia")
HILL_NS = (10.0, 2.0, 1.0)
N_COALITIONS = 24
N_INIT = 4
SEED = 20260921


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--classifier", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--networks", default=",".join(NETWORKS))
    parser.add_argument("--hill-ns", default=",".join(f"{h:g}" for h in HILL_NS))
    parser.add_argument("--n-coalitions", type=int, default=N_COALITIONS)
    args = parser.parse_args()
    solver, classifier = load_solver(args.solver), load_classifier(args.classifier)
    networks = tuple(args.networks.split(","))
    hill_ns = tuple(float(h) for h in args.hill_ns.split(","))
    rng = np.random.default_rng(SEED)

    rows = []
    for name in networks:
        info = ALL_MODELS[name]
        net = compile_model(info["rules"], info["output_nodes"])
        n = len(net.node_names)
        init = initial_states(N_INIT, n, seed=SEED)
        coalitions = rng.integers(0, 2**n, size=args.n_coalitions)
        for hill_n in hill_ns:
            prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=hill_n))
            for c in coalitions:
                r = simulate_coalition(prep, solver, int(c), init, 0.0, classifier)
                for t in r.results:
                    rows.append({
                        "network": name, "hill_n": hill_n, "coalition": int(c),
                        "status": t.status, "class": t.klass, "horizon_index": t.horizon_index,
                        "attempt": t.attempt, "nfev": t.nfev,
                        "tail_range": t.stats.tail_range if t.stats else None,
                        "final_derivative": t.stats.final_derivative if t.stats else None,
                        "halves_envelope_diff": t.stats.halves_envelope_diff if t.stats else None,
                    })
            done = [x for x in rows if x["network"] == name and x["hill_n"] == hill_n]
            print(f"{name:<24} n_H={hill_n:<4g} " + " ".join(
                f"{s}={sum(x['status'] == s for x in done)}" for s in STATUSES if any(x['status'] == s for x in done))
                + "  " + " ".join(f"{c}={sum(x['class'] == c for x in done)}" for c in TRAJECTORY_CLASSES
                                  if any(x['class'] == c for x in done)))

    def q(values):
        v = np.array([x for x in values if x is not None], dtype=float)
        return None if v.size == 0 else {k: float(np.quantile(v, p)) for k, p in
                                         (("min", 0), ("median", 0.5), ("q90", 0.9), ("q99", 0.99), ("max", 1))}

    summary = {}
    for key in sorted({(x["network"], x["hill_n"]) for x in rows}):
        sel = [x for x in rows if (x["network"], x["hill_n"]) == key]
        summary[f"{key[0]}__n{key[1]:g}"] = {
            "statuses": {s: sum(x["status"] == s for x in sel) for s in STATUSES},
            "classes": {c: sum(x["class"] == c for x in sel) for c in TRAJECTORY_CLASSES},
            "horizon_index": {str(h): sum(x["horizon_index"] == h for x in sel)
                              for h in range(len(classifier.extension_horizons) + 1)},
            "fixed_final_derivative": q([x["final_derivative"] for x in sel if x["class"] == "fixed"]),
            "fixed_tail_range": q([x["tail_range"] for x in sel if x["class"] == "fixed"]),
            "oscillatory_envelope_diff": q([x["halves_envelope_diff"] for x in sel if x["class"] == "oscillatory"]),
            "unclassified_final_derivative": q([x["final_derivative"] for x in sel if x["status"] == "unclassified"]),
            "unclassified_tail_range": q([x["tail_range"] for x in sel if x["status"] == "unclassified"]),
            "unclassified_envelope_diff": q([x["halves_envelope_diff"] for x in sel if x["status"] == "unclassified"]),
        }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "solver_file": str(args.solver), "solver": json.loads(args.solver.read_text()),
        "classifier_file": str(args.classifier), "classifier": json.loads(args.classifier.read_text()),
        "networks": networks, "hill_ns": hill_ns, "n_coalitions": args.n_coalitions, "n_init": N_INIT,
        "seed": SEED, "summary": summary, "trajectories": rows}, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
