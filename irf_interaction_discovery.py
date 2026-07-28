"""Iterative Random Forests (iRF) for stable interaction detection.

Basu et al. PNAS 2018. The one genuine genetics import that's plausibly
competitive with Walsh-based methods on higher-order detection.

Scored by predicting on all 2^n coalitions and decomposing predictions
via WHT — same metrics as every other method, in a basis iRF never
touches directly.

min_samples_leaf scales with training set size to avoid handicapping
the method at low budgets.

Usage:
    uv run python irf_interaction_discovery.py \
        --data /path/to/circuit_coalition_values.npz \
        --output results/irf_results.json
"""

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from tqdm import tqdm

from data_utils import (
    compute_recovery_metrics,
    energy_spectrum,
    load_coalition_table,
    normalized_wht,
    pooled_mean_values,
    wht,
)

N_ITERATIONS = 5
N_TREES = 500
MAX_INTERACTION_ORDER = 4


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def coalition_to_features(n):
    N = 2 ** n
    X = np.zeros((N, n), dtype=np.float64)
    for i in range(n):
        X[:, i] = (np.arange(N) >> i) & 1
    return X


def irf_pipeline(X, y, n_iterations=N_ITERATIONS, n_trees=N_TREES,
                 max_order=MAX_INTERACTION_ORDER, seed=42):
    n_features = X.shape[1]
    n_samples = X.shape[0]
    weights = np.ones(n_features)

    min_leaf = max(1, n_samples // 100)

    importance_history = []

    for it in range(n_iterations):
        print(f"  [iRF iteration {it+1}/{n_iterations}] {n_trees} trees, "
              f"min_leaf={min_leaf}, m={n_samples}")

        n_virtual = max(n_features * 5, 50)
        w = weights / weights.sum()
        rng = np.random.default_rng(seed + it + 1000)
        virtual_idx = rng.choice(n_features, size=n_virtual, p=w)
        X_weighted = X[:, virtual_idx]

        rf = RandomForestRegressor(
            n_estimators=n_trees,
            max_features=None,
            min_samples_leaf=min_leaf,
            random_state=seed + it,
            n_jobs=-1,
        )
        rf.fit(X_weighted, y)

        raw_imp = rf.feature_importances_
        importances = np.zeros(n_features)
        for vi, orig in enumerate(virtual_idx):
            importances[orig] += raw_imp[vi]
        total = importances.sum()
        if total > 0:
            importances /= total
        importance_history.append(importances.copy())

        r2 = rf.score(X_weighted, y)
        print(f"    R² = {r2:.4f}")

        weights = importances + 1e-10

    return {
        "importance_history": importance_history,
        "final_importances": importance_history[-1],
        "final_rf": rf,
        "final_virtual_idx": virtual_idx,
    }


def irf_recover(v, n, m, seed=42, n_iterations=N_ITERATIONS, n_trees=N_TREES,
                max_order=MAX_INTERACTION_ORDER):
    """iRF recovery with the same train/heldout split protocol as LASSO."""
    N = 2 ** n
    rng = np.random.default_rng(seed)
    w_true = normalized_wht(v)

    all_indices = rng.permutation(N)
    sample_indices = np.sort(all_indices[:m])
    heldout_indices = np.sort(all_indices[m:m + min(m, N - m)])

    X_full = coalition_to_features(n)
    X_train = X_full[sample_indices]
    y_train = v[sample_indices]

    irf_out = irf_pipeline(
        X_train, y_train,
        n_iterations=n_iterations,
        n_trees=n_trees,
        max_order=max_order,
        seed=seed,
    )

    rf = irf_out["final_rf"]
    vidx = irf_out["final_virtual_idx"]
    v_pred = rf.predict(X_full[:, vidx])

    w_recovered = normalized_wht(v_pred)

    metrics = compute_recovery_metrics(
        w_recovered, w_true, n,
        v_pred=v_pred, v_true=v,
        heldout_indices=heldout_indices,
    )
    metrics["m"] = m
    metrics["sample_fraction"] = m / N
    metrics["seed"] = seed
    metrics["feature_importances"] = irf_out["final_importances"].tolist()
    metrics["train_r2"] = float(rf.score(X_train[:, vidx], y_train))

    return metrics


def run_budget_sweep(v, n, budgets, n_trials=3, n_iterations=N_ITERATIONS,
                     n_trees=N_TREES, max_order=MAX_INTERACTION_ORDER):
    N = 2 ** n
    results = []

    for frac in budgets:
        m = max(int(frac * N), n + 1)
        print(f"\n[{timestamp()}] Budget {frac*100:.0f}% — m={m} samples, {n_trials} trials")

        trial_results = []
        for trial in tqdm(range(n_trials), desc=f"Budget {frac*100:.0f}%"):
            r = irf_recover(
                v, n, m,
                seed=42 + trial,
                n_iterations=n_iterations,
                n_trees=n_trees,
                max_order=max_order,
            )
            r["trial"] = trial
            trial_results.append(r)

        def mean_std(key):
            vals = [r[key] for r in trial_results if r.get(key) is not None]
            if not vals:
                return float('nan'), float('nan')
            return float(np.mean(vals)), float(np.std(vals))

        nmse_m, nmse_s = mean_std("spectrum_nmse_no_intercept")
        held_m, held_s = mean_std("heldout_r2")
        pair_m, pair_s = mean_std("pairwise_spearman_rho")

        summary = {
            "budget_fraction": frac,
            "m": m,
            "n_trials": n_trials,
            "spectrum_nmse_no_intercept_mean": nmse_m,
            "spectrum_nmse_no_intercept_std": nmse_s,
            "heldout_r2_mean": held_m,
            "heldout_r2_std": held_s,
            "pairwise_rho_mean": pair_m,
            "pairwise_rho_std": pair_s,
            "trials": trial_results,
        }

        print(f"  NMSE (no intercept): {nmse_m:.4f} ± {nmse_s:.4f}")
        print(f"  Heldout R²:          {held_m:.4f} ± {held_s:.4f}")
        print(f"  Pairwise Spearman:   {pair_m:.4f} ± {pair_s:.4f}")

        results.append(summary)

    return results


def main():
    parser = argparse.ArgumentParser(description="iRF interaction discovery")
    parser.add_argument("--data", required=True, help="Path to coalition table npz")
    parser.add_argument("--output", required=True, help="Path to save results JSON")
    parser.add_argument("--circuit", default="unknown", help="Circuit name for labeling")
    parser.add_argument("--n-iterations", type=int, default=N_ITERATIONS)
    parser.add_argument("--n-trees", type=int, default=N_TREES)
    parser.add_argument("--max-order", type=int, default=MAX_INTERACTION_ORDER)
    parser.add_argument("--budget-sweep", action="store_true")
    parser.add_argument("--budgets", type=float, nargs="+",
                        default=[0.01, 0.03, 0.05, 0.10, 0.20, 1.0])
    args = parser.parse_args()

    print(f"[{timestamp()}] Loading coalition table: {args.data}")
    table = load_coalition_table(args.data)
    n = table["n_players"]
    v = pooled_mean_values(table)
    N = 2 ** n
    print(f"  n_players={n}, n_prompts={table['n_prompts']}, N={N}")

    w_true = normalized_wht(v)
    e_true = energy_spectrum(w_true, n)
    print(f"  Ground truth energy spectrum:")
    for k in range(min(6, n + 1)):
        print(f"    order-{k}: {e_true[k]*100:.2f}%")
    print(f"    order-3+: {e_true[3:].sum()*100:.2f}%")

    if args.budget_sweep:
        sweep_results = run_budget_sweep(
            v, n,
            budgets=args.budgets,
            n_iterations=args.n_iterations,
            n_trees=args.n_trees,
            max_order=args.max_order,
        )
        output = {
            "method": "irf_budget_sweep",
            "circuit": args.circuit,
            "data_path": str(args.data),
            "n_players": n,
            "n_prompts": table["n_prompts"],
            "ground_truth_energy_spectrum": e_true.tolist(),
            "timestamp": timestamp(),
            "budget_sweep": sweep_results,
        }
    else:
        print(f"\n[{timestamp()}] Running iRF on full coalition table...")
        metrics = irf_recover(
            v, n, m=N,
            n_iterations=args.n_iterations,
            n_trees=args.n_trees,
            max_order=args.max_order,
        )
        output = {
            "method": "irf_interaction_discovery",
            "circuit": args.circuit,
            "data_path": str(args.data),
            "n_players": n,
            "n_prompts": table["n_prompts"],
            "ground_truth_energy_spectrum": e_true.tolist(),
            "timestamp": timestamp(),
            "scores": metrics,
        }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\n[{timestamp()}] Results saved to {args.output}")


if __name__ == "__main__":
    main()
