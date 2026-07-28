"""Modal launcher for the shapiq budget sweep.

Runs the full pre-registered shapiq sweep on Modal CPU instances,
parallelizing across (circuit, budget, seed) triples.

Usage:
    modal run --detach scripts/modal_shapiq_sweep.py

    # Subset for testing:
    modal run scripts/modal_shapiq_sweep.py --circuits weight_ioi --budgets 0.05 --n-seeds 2
"""
from __future__ import annotations

import os

import modal

REPO = os.path.join(os.path.expanduser("~"), "Documents/GitHub/epistasis-bench")
DATA_DIR = os.path.join(
    os.path.expanduser("~"),
    "Documents/GitHub/weight-circuit-discovery/experiments_batch2/genetics",
)

DATA_FILES = [
    "weight_ioi_zero_v2_coalition_values.npz",
    "ioi_zero_v2_coalition_values.npz",
    "random15_zero_v2_coalition_values.npz",
]

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "numpy==2.2.6",
        "scipy==1.15.3",
        "scikit-learn==1.6.1",
        "shapiq==1.6.0",
        "tqdm==4.67.1",
        "matplotlib==3.10.3",
    )
    .env({"PYTHONPATH": "/root/repo"})
    .add_local_file(os.path.join(REPO, "data_utils.py"), "/root/repo/data_utils.py")
    .add_local_file(os.path.join(REPO, "budget_sweep_harness.py"), "/root/repo/budget_sweep_harness.py")
    .add_local_file(os.path.join(REPO, "shapiq_approximators.py"), "/root/repo/shapiq_approximators.py")
    .add_local_file(os.path.join(REPO, "shapiq_sweep_harness.py"), "/root/repo/shapiq_sweep_harness.py")
    .add_local_file(os.path.join(REPO, "irf_interaction_discovery.py"), "/root/repo/irf_interaction_discovery.py")
    .add_local_file(os.path.join(REPO, "lasso_walsh_oracle.py"), "/root/repo/lasso_walsh_oracle.py")
)

for f in DATA_FILES:
    fpath = os.path.join(DATA_DIR, f)
    if os.path.exists(fpath):
        image = image.add_local_file(fpath, f"/root/data/{f}")

app = modal.App("epistasis-bench-shapiq-sweep", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

CIRCUIT_MAP = {
    "weight_ioi": "weight_ioi_zero_v2_coalition_values.npz",
    "ioi": "ioi_zero_v2_coalition_values.npz",
    "random15": "random15_zero_v2_coalition_values.npz",
}

DEFAULT_BUDGETS = [0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.40]
DEFAULT_CIRCUITS = ["weight_ioi", "ioi", "random15"]
DEFAULT_MAX_ORDERS = [2, 3]
DEFAULT_METHODS = ["lasso_walsh", "irf", "kernelshapiq", "shapiq_mc", "svarmiq"]
MAX_ORDER_4_MIN_BUDGET = 0.20


@app.function(cpu=2, memory=4096, timeout=86400, volumes={"/results": results_vol})
def run_single_trial(circuit_name: str, budget_frac: float, seed: int,
                     max_orders: list[int], methods: list[str]):
    """Run all methods for a single (circuit, budget, seed) triple."""
    import numpy as np
    from data_utils import load_coalition_table, pooled_mean_values
    from shapiq_sweep_harness import run_trial

    data_path = f"/root/data/{CIRCUIT_MAP[circuit_name]}"
    table = load_coalition_table(data_path)
    v = pooled_mean_values(table)
    n = table["n_players"]

    active_orders = [o for o in max_orders if o != 4 or budget_frac >= MAX_ORDER_4_MIN_BUDGET]
    result = run_trial(v, n, budget_frac, seed, active_orders, methods)

    return {
        "circuit": circuit_name,
        "budget_fraction": budget_frac,
        "seed": seed,
        "max_orders": active_orders,
        "trial": result,
    }


@app.function(cpu=4, memory=8192, timeout=86400, volumes={"/results": results_vol})
def orchestrate(circuits: list[str], budgets: list[float], n_seeds: int,
                max_orders: list[int], methods: list[str]):
    """Fan out all trials in parallel, collect and save results."""
    import importlib.metadata
    import json
    from datetime import datetime, timezone

    runtime_versions = {}
    for pkg in ["shapiq", "numpy", "scipy", "scikit-learn"]:
        try:
            runtime_versions[pkg] = importlib.metadata.version(pkg)
        except Exception:
            runtime_versions[pkg] = "NOT_FOUND"

    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] Runtime versions: {runtime_versions}")
    print(f"[{ts}] Circuits: {circuits}, Budgets: {budgets}, Seeds: {n_seeds}")
    print(f"[{ts}] Max orders: {max_orders}, Methods: {methods}")

    jobs = []
    for circuit in circuits:
        for frac in budgets:
            n_seeds_this = n_seeds if frac <= 0.10 else max(n_seeds // 2, 5)
            for seed_offset in range(n_seeds_this):
                seed = 42 + seed_offset
                jobs.append((circuit, frac, seed))

    print(f"[{ts}] Launching {len(jobs)} parallel trials")

    futures = []
    for circuit, frac, seed in jobs:
        futures.append(
            run_single_trial.spawn(circuit, frac, seed, max_orders, methods)
        )

    all_results = {
        "harness_version": "2.1-modal-auroc",
        "prereg": "prereg_shapiq_budget_sweep_v2.md",
        "prereg_sha256": "2caaffff9590e906305eb4fc2c6ee1e8a797fcfec31d4bc86f6b7782b985a976",
        "runtime_versions": runtime_versions,
        "heldout_fraction": 0.10,
        "budgets": budgets,
        "n_seeds": n_seeds,
        "max_orders": max_orders,
        "methods": methods,
        "n_jobs": len(jobs),
        "timestamp_start": ts,
        "trials": [],
    }

    completed = 0
    for future in futures:
        try:
            result = future.get()
            all_results["trials"].append(result)
            completed += 1
            if completed % 50 == 0 or completed == len(futures):
                ts_now = datetime.now(timezone.utc).isoformat()
                print(f"[{ts_now}] {completed}/{len(futures)} trials complete")
                with open("/results/shapiq_sweep_v2_incremental.json", "w") as f:
                    json.dump(all_results, f, indent=2, default=str)
                results_vol.commit()
        except Exception as e:
            print(f"Trial failed: {e}")
            all_results["trials"].append({"error": str(e)})

    all_results["timestamp_end"] = datetime.now(timezone.utc).isoformat()
    all_results["n_completed"] = completed
    all_results["n_failed"] = len(futures) - completed

    with open("/results/shapiq_sweep_v2.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    results_vol.commit()

    print(f"\n[{all_results['timestamp_end']}] Done: {completed}/{len(futures)} trials")
    print(f"Results at /results/shapiq_sweep_v2.json")

    return all_results


@app.local_entrypoint()
def main(
    circuits: str = ",".join(DEFAULT_CIRCUITS),
    budgets: str = ",".join(str(b) for b in DEFAULT_BUDGETS),
    n_seeds: int = 20,
    max_orders: str = "2,3",
    methods: str = ",".join(DEFAULT_METHODS),
):
    circuit_list = circuits.split(",")
    budget_list = [float(b) for b in budgets.split(",")]
    order_list = [int(o) for o in max_orders.split(",")]
    method_list = methods.split(",")

    result = orchestrate.remote(circuit_list, budget_list, n_seeds, order_list, method_list)
    print(f"\nCompleted: {result['n_completed']} trials, {result['n_failed']} failed")
