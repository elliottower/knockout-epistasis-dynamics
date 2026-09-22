"""How long does one trajectory take, what would the reruns cost, and which stiff solver to fall back on?

Times the engine on random coalitions of three networks of increasing size, under given
solver and classifier files, and extrapolates to every coalition of every paper network.
Then runs Radau and BDF alone on the same trajectories and compares each with a tight RK45
reference, to choose the fallback. Planning only: no coalition value is kept or scored.

    uv run python -m scripts.benchmark_ode_engine \
        --solver experiments/ode_v2/solver_primary.json \
        --classifier experiments/ode_v2/classifier_primary.json

Writes results/audit/engine_timing.json.
"""

import argparse
import dataclasses
import json
import time
from pathlib import Path

import numpy as np

from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.ode_engine import (
    STATUSES,
    DynamicsConfig,
    SolverAttempt,
    clamp_mask_for,
    compile_model,
    initial_states,
    integrate,
    prepare,
    simulate_coalition,
)
from scripts.paths import AUDIT, PROJECT_ROOT, RESULTS
from scripts.run_ode_arm import load_classifier, load_solver

OUT = AUDIT / "engine_timing.json"
NETWORKS = ("lambda_phage", "davidich_yeast", "arabidopsis_cellcycle")
CONFIGS = (("hillcube_normalized", 10.0), ("hillcube_normalized", 2.0), ("operatorwise_normalized", 2.0))
N_COALITIONS = 12
N_INIT = 4
N_STIFF_TRAJECTORIES = 8
STIFF_METHODS = ("Radau", "BDF")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--classifier", type=Path, required=True)
    args = parser.parse_args()
    solver, classifier = load_solver(args.solver), load_classifier(args.classifier)
    rng = np.random.default_rng()

    rows, stiff = [], []
    for name in NETWORKS:
        info = ALL_MODELS[name]
        net = compile_model(info["rules"], info["output_nodes"])
        n = len(net.node_names)
        init = initial_states(N_INIT, n, seed=0)
        coalitions = rng.integers(0, 2**n, size=N_COALITIONS)
        for construction, hill_n in CONFIGS:
            prep = prepare(net, DynamicsConfig(construction, hill_n=hill_n))
            t0 = time.perf_counter()
            results = [simulate_coalition(prep, solver, int(c), init, 0.0, classifier) for c in coalitions]
            elapsed = time.perf_counter() - t0
            traj = {k: np.concatenate([r.trajectories[k] for r in results]) for k in results[0].trajectories}
            rows.append({
                "network": name, "n_nodes": n, "construction": construction, "hill_n": hill_n,
                "coalitions": coalitions.tolist(),
                "statuses": {st: int((traj["status"] == i).sum()) for i, st in enumerate(STATUSES)},
                "seconds_per_trajectory": elapsed / (N_COALITIONS * N_INIT),
                "nfev_max": int(traj["nfev"].max()),
                "trajectories_not_ok": int((traj["status"] != 0).sum()),
                "trajectories_on_fallback": int((traj["attempt"] > 0).sum()),
                "trajectories_extended": int((traj["horizon_index"] > 0).sum()),
            })
            print(f"{name:<24} n={n:<3} {construction:<24} n_H={hill_n:<4g} "
                  f"{rows[-1]['seconds_per_trajectory']*1e3:7.1f} ms/traj  nfev_max {rows[-1]['nfev_max']}  "
                  f"not ok {rows[-1]['trajectories_not_ok']}  extended {rows[-1]['trajectories_extended']}")

        # The stiff fallback, forced, against a tight RK45 reference. No classifier: fixed horizon t_max.
        prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=10.0))
        reference = dataclasses.replace(solver, attempts=(SolverAttempt("RK45", 10_000_000),), rtol=1e-10, atol=1e-12)
        for c, y0 in zip(rng.integers(0, 2**n, size=N_STIFF_TRAJECTORIES), initial_states(N_STIFF_TRAJECTORIES, n, seed=1)):
            mask = clamp_mask_for(int(c), n)
            y0 = np.where(mask, 0.0, y0)
            ref = integrate(prep, reference, mask, 0.0, y0)
            for method in STIFF_METHODS:
                forced = dataclasses.replace(solver, attempts=(SolverAttempt(method, 10_000_000),))
                t0 = time.perf_counter()
                r = integrate(prep, forced, mask, 0.0, y0)
                stiff.append({"network": name, "coalition": int(c), "method": method, "status": r.status, "nfev": r.nfev,
                              "seconds": time.perf_counter() - t0,
                              "abs_error_vs_tight_rk45": abs(r.output - ref.output)})

    stiff_summary = {m: {
        "n": sum(s["method"] == m for s in stiff),
        "not_ok": sum(s["method"] == m and s["status"] != "ok" for s in stiff),
        "median_seconds": float(np.median([s["seconds"] for s in stiff if s["method"] == m])),
        "max_abs_error_vs_tight_rk45": float(np.nanmax([s["abs_error_vs_tight_rk45"] for s in stiff if s["method"] == m])),
    } for m in STIFF_METHODS}

    # Seconds per trajectory grows with network size; fit log(sec) on n per configuration.
    table = json.loads((RESULTS / "paper_number_reconciliation.json").read_text())["per_network_table"]
    sizes = {m: v["n"] for m, v in table.items()}
    extrapolation = {}
    for construction, hill_n in CONFIGS:
        pts = [(r["n_nodes"], r["seconds_per_trajectory"]) for r in rows
               if r["construction"] == construction and r["hill_n"] == hill_n]
        slope, intercept = np.polyfit([p[0] for p in pts], np.log([p[1] for p in pts]), 1)
        per_net = {m: float(np.exp(intercept + slope * n)) * (2**n) * 32 / 3600 for m, n in sizes.items()}
        extrapolation[f"{construction}__n{hill_n:g}"] = {
            "cpu_hours_all_28_networks_32_inits": sum(per_net.values()),
            "cpu_hours_networks_n_le_12": sum(v for m, v in per_net.items() if sizes[m] <= 12),
            "largest": sorted(per_net.items(), key=lambda kv: -kv[1])[:3],
        }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"solver_file": str(args.solver), "solver": json.loads(args.solver.read_text()),
                               "classifier_file": str(args.classifier),
                               "classifier": json.loads(args.classifier.read_text()),
                               "n_coalitions_timed": N_COALITIONS, "n_init_timed": N_INIT, "rows": rows,
                               "stiff_fallback": {"summary": stiff_summary, "trajectories": stiff},
                               "extrapolation_single_core": extrapolation}, indent=2))
    print(f"wrote {OUT.relative_to(PROJECT_ROOT)}")
    for m, s in stiff_summary.items():
        print(f"  {m:<6} not ok {s['not_ok']}/{s['n']}  median {s['median_seconds']*1e3:.0f} ms  "
              f"max |error| vs tight RK45 {s['max_abs_error_vs_tight_rk45']:.2e}")
    for key, e in extrapolation.items():
        print(f"  {key:<34} all 28: {e['cpu_hours_all_28_networks_32_inits']:8.0f} CPU-h   "
              f"n<=12: {e['cpu_hours_networks_n_le_12']:6.1f} CPU-h   largest {e['largest'][0][0]}")


if __name__ == "__main__":
    main()
