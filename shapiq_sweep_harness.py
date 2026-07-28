"""Shapiq budget-sweep harness for EpistasisBench.

Runs KernelSHAP-IQ, SHAPIQ MC, and SVARM-IQ alongside LASSO-Walsh
and iRF on the SAME pre-designated held-out set. All five methods
are scored on identical test coalitions.

Pre-designated held-out set: 10% of 2^n coalitions drawn once per
seed before any method runs. LASSO/iRF draw training samples from
the remaining pool. Shapiq methods may query any coalition in the
pool via smart sampling (tracked and logged).

Usage:
    uv run python shapiq_sweep_harness.py \
        --circuits weight_ioi ioi random15 \
        --data-dir /path/to/genetics/ \
        --output results/shapiq_sweep_v1.json

    uv run python shapiq_sweep_harness.py \
        --circuits weight_ioi \
        --data-dir /path/to/genetics/ \
        --budgets 0.05 \
        --n-seeds 1 \
        --max-orders 2 3 \
        --output results/shapiq_timing_probe.json
"""

import argparse
import importlib.metadata
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from tqdm import tqdm

from budget_sweep_harness import (
    CIRCUIT_FILES,
    LASSO_CV_FOLDS,
    LASSO_MAX_ORDER,
    IRF_N_ITERATIONS,
    IRF_N_TREES,
    IRF_LEAF_GRID,
    IRF_DEPTH_GRID,
    build_walsh_design_no_intercept,
    coalition_to_features,
    run_lasso,
    run_irf,
    summarize_trials,
)
from data_utils import (
    compute_recovery_metrics,
    load_coalition_table,
    normalized_wht,
    energy_spectrum,
    compute_k99,
    pooled_mean_values,
)
from shapiq_approximators import (
    APPROXIMATORS,
    make_game_function,
    reconstruct_from_interaction_values,
)

HELDOUT_FRACTION = 0.10
DEFAULT_BUDGETS = [0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.40]
DEFAULT_N_SEEDS = 20
DEFAULT_MAX_ORDERS = [2, 3]
SHAPIQ_INDEX = "k-SII"
MAX_ORDER_4_MIN_BUDGET = 0.20


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def get_runtime_versions():
    packages = ["shapiq", "numpy", "scipy", "scikit-learn"]
    versions = {}
    for pkg in packages:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = "NOT_FOUND"
    return versions


def make_tracked_game(v, n, heldout_set):
    """Game function that tracks which coalitions shapiq queries.

    Returns (game_fn, get_stats) where get_stats() returns a dict
    of query statistics including held-out overlap.
    """
    powers = 1 << np.arange(n)
    queried_indices = []

    def game_fn(coalitions):
        indices = (np.asarray(coalitions, dtype=bool) * powers).sum(axis=1).astype(int)
        queried_indices.extend(indices.tolist())
        return v[indices]

    def get_stats():
        if not queried_indices:
            return {"n_queries": 0, "n_heldout_queries": 0, "heldout_query_frac": 0.0}
        q = set(queried_indices)
        overlap = q & heldout_set
        return {
            "n_queries": len(queried_indices),
            "n_unique_queries": len(q),
            "n_heldout_queries": len(overlap),
            "heldout_query_frac": round(len(overlap) / len(queried_indices), 4),
        }

    return game_fn, get_stats


def draw_heldout_and_pool(n, seed):
    """Draw pre-designated held-out set and available pool.

    Returns (heldout_indices, pool_indices) where heldout is 10%
    of 2^n coalitions and pool is the remaining 90%.
    """
    N = 2 ** n
    n_heldout = int(round(HELDOUT_FRACTION * N))
    rng = np.random.default_rng(seed)
    perm = rng.permutation(N)
    heldout = np.sort(perm[:n_heldout])
    pool = np.sort(perm[n_heldout:])
    return heldout, pool


def run_shapiq_on_heldout(v, n, budget_n_evals, method_name, seed,
                          max_order, heldout_indices, heldout_set):
    """Run a shapiq approximator and score on pre-designated held-out set."""
    N = 2 ** n
    w_true = normalized_wht(v)

    game_fn, get_query_stats = make_tracked_game(v, n, heldout_set)
    approx_cls = APPROXIMATORS[method_name]

    t0 = time.monotonic()
    approx = approx_cls(n=n, max_order=max_order, index=SHAPIQ_INDEX, random_state=seed)
    iv = approx.approximate(budget=budget_n_evals, game=game_fn)
    v_pred = reconstruct_from_interaction_values(iv, n)
    elapsed = time.monotonic() - t0

    w_recovered = normalized_wht(v_pred)

    metrics = compute_recovery_metrics(
        w_recovered, w_true, n,
        v_pred=v_pred, v_true=v,
        heldout_indices=heldout_indices,
    )

    metrics["method"] = method_name
    metrics["wall_seconds"] = round(elapsed, 2)
    metrics["shapiq_index"] = SHAPIQ_INDEX
    metrics["max_order"] = max_order
    metrics["budget_n_evals"] = budget_n_evals
    metrics["baseline_value"] = float(iv.baseline_value)
    metrics["n_interactions"] = len([t for t in iv.dict_values if len(t) > 0])
    metrics["query_stats"] = get_query_stats()

    return metrics


def run_trial(v, n, budget_frac, seed, max_orders, methods):
    """Run all methods on the same pre-designated held-out set."""
    N = 2 ** n
    m = int(round(budget_frac * N))
    assert m == int(round(budget_frac * (2 ** n))), (
        f"Budget must be fraction of 2^n: {m} != {int(round(budget_frac * (2 ** n)))}"
    )
    m = max(m, n + 1)

    heldout_indices, pool_indices = draw_heldout_and_pool(n, seed)
    heldout_set = set(heldout_indices.tolist())

    trial_rng = np.random.default_rng(seed + 10000)
    sample_from_pool = trial_rng.choice(pool_indices, size=min(m, len(pool_indices)), replace=False)
    sample_indices = np.sort(sample_from_pool)

    results = {}

    if "lasso_walsh" in methods:
        lr = run_lasso(v, n, sample_indices, heldout_indices)
        lr["m"] = m
        lr["budget_fraction"] = budget_frac
        lr["seed"] = seed
        results["lasso_walsh"] = lr

    if "irf" in methods:
        ir = run_irf(v, n, sample_indices, heldout_indices, seed=seed)
        ir["m"] = m
        ir["budget_fraction"] = budget_frac
        ir["seed"] = seed
        results["irf"] = ir

    for method_name in ["kernelshapiq", "shapiq_mc", "svarmiq"]:
        if method_name not in methods:
            continue
        for max_order in max_orders:
            if max_order == 4 and budget_frac < MAX_ORDER_4_MIN_BUDGET:
                continue
            key = f"{method_name}_order{max_order}"
            sr = run_shapiq_on_heldout(
                v, n, m, method_name, seed, max_order,
                heldout_indices, heldout_set,
            )
            sr["m"] = m
            sr["budget_fraction"] = budget_frac
            sr["seed"] = seed
            results[key] = sr

    return results


def save_incremental(all_results, output_path):
    """Save results incrementally so crashes don't lose data."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)


def run_circuit(circuit_name, data_path, budgets, n_seeds, max_orders, methods,
                all_results=None, output_path=None):
    print(f"\n{'='*60}")
    print(f"[{timestamp()}] Circuit: {circuit_name}")

    table = load_coalition_table(str(data_path))
    n = table["n_players"]
    v = pooled_mean_values(table)
    N = 2 ** n
    w_true = normalized_wht(v)
    e_true = energy_spectrum(w_true, n)
    k99 = compute_k99(w_true, n)

    print(f"  n={n}, N={N}, k99={k99}")
    print(f"  order-0: {e_true[0]*100:.1f}%, order-1: {e_true[1]*100:.1f}%, "
          f"order-2: {e_true[2]*100:.1f}%, order-3+: {e_true[3:].sum()*100:.1f}%")

    circuit_results = {
        "circuit": circuit_name,
        "n_players": n,
        "n_prompts": table["n_prompts"],
        "k99": k99,
        "energy_spectrum": e_true.tolist(),
        "budgets": [],
    }

    for frac in budgets:
        m = max(int(round(frac * N)), n + 1)
        n_seeds_this = n_seeds if frac <= 0.10 else max(n_seeds // 2, 5)

        active_orders = [o for o in max_orders if o != 4 or frac >= MAX_ORDER_4_MIN_BUDGET]
        n_methods = len([m for m in methods if m in ("lasso_walsh", "irf")])
        n_methods += len([m for m in methods if m in APPROXIMATORS]) * len(active_orders)

        print(f"\n  Budget {frac*100:.0f}% (m={m}), {n_seeds_this} seeds, "
              f"{n_methods} methods, orders {active_orders}")

        all_trials = []
        for seed_offset in tqdm(range(n_seeds_this), desc=f"  {frac*100:.0f}%"):
            seed = 42 + seed_offset
            trial = run_trial(v, n, frac, seed, active_orders, methods)
            all_trials.append(trial)

        method_trials = {}
        for trial in all_trials:
            for key, result in trial.items():
                method_trials.setdefault(key, []).append(result)

        summaries = {}
        for key, trials in method_trials.items():
            summaries[key] = summarize_trials(trials, key)
            r2_med = summaries[key]["heldout_r2"]["median"]
            print(f"    {key}: R²={r2_med:.4f}" if r2_med else f"    {key}: R²=None")

        budget_result = {
            "budget_fraction": frac,
            "m": m,
            "n_seeds": n_seeds_this,
            "max_orders": active_orders,
            "summaries": summaries,
            "trials": method_trials,
        }
        circuit_results["budgets"].append(budget_result)

        if all_results is not None and output_path is not None:
            all_results["last_save"] = timestamp()
            save_incremental(all_results, output_path)
            print(f"  [incremental save at {timestamp()}]")

    return circuit_results


def main():
    parser = argparse.ArgumentParser(description="EpistasisBench shapiq budget sweep")
    parser.add_argument("--circuits", nargs="+", required=True,
                        choices=list(CIRCUIT_FILES.keys()))
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--budgets", type=float, nargs="+", default=DEFAULT_BUDGETS)
    parser.add_argument("--n-seeds", type=int, default=DEFAULT_N_SEEDS)
    parser.add_argument("--max-orders", type=int, nargs="+", default=DEFAULT_MAX_ORDERS)
    parser.add_argument("--methods", nargs="+",
                        default=["lasso_walsh", "irf", "kernelshapiq", "shapiq_mc", "svarmiq"])
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    runtime_versions = get_runtime_versions()
    print(f"[{timestamp()}] Runtime versions: {runtime_versions}")

    all_results = {
        "harness_version": "2.0",
        "prereg": "prereg_shapiq_budget_sweep_v2.md",
        "prereg_sha256": "2caaffff9590e906305eb4fc2c6ee1e8a797fcfec31d4bc86f6b7782b985a976",
        "runtime_versions": runtime_versions,
        "heldout_fraction": HELDOUT_FRACTION,
        "budgets": args.budgets,
        "n_seeds": args.n_seeds,
        "max_orders": args.max_orders,
        "methods": args.methods,
        "lasso_config": {
            "max_order": LASSO_MAX_ORDER,
            "cv_folds": LASSO_CV_FOLDS,
        },
        "irf_config": {
            "n_iterations": IRF_N_ITERATIONS,
            "n_trees": IRF_N_TREES,
            "leaf_grid": IRF_LEAF_GRID,
            "depth_grid": [str(d) for d in IRF_DEPTH_GRID],
        },
        "shapiq_config": {
            "index": SHAPIQ_INDEX,
            "max_order_4_min_budget": MAX_ORDER_4_MIN_BUDGET,
        },
        "timestamp_start": timestamp(),
        "circuits": [],
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    for circuit in args.circuits:
        data_path = data_dir / CIRCUIT_FILES[circuit]
        if not data_path.exists():
            print(f"WARNING: {data_path} not found, skipping {circuit}")
            continue
        result = run_circuit(circuit, data_path, args.budgets, args.n_seeds,
                             args.max_orders, args.methods,
                             all_results=all_results, output_path=args.output)
        all_results["circuits"].append(result)
        save_incremental(all_results, args.output)
        print(f"\n[{timestamp()}] Circuit {circuit} complete, saved to {args.output}")

    all_results["timestamp_end"] = timestamp()
    save_incremental(all_results, args.output)
    print(f"\n[{timestamp()}] All results saved to {args.output}")


if __name__ == "__main__":
    main()
