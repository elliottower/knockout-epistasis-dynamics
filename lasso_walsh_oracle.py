"""LASSO sparse recovery on the Walsh-Hadamard basis (oracle ceiling).

This is NOT a benchmark competitor — it is the theoretical ceiling that
all other methods are measured against. On Pareto plots, this appears as
the dashed line at the top.

The T=0 (intercept) column is excluded from the penalized design matrix
and recovered via fit_intercept=True, so the largest coefficient is not
L1-shrunk.

Usage:
    uv run python lasso_walsh_oracle.py \
        --data /path/to/circuit_coalition_values.npz \
        --output results/lasso_walsh_results.json
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from sklearn.linear_model import LassoCV
from tqdm import tqdm

from data_utils import (
    compute_recovery_metrics,
    energy_spectrum,
    load_coalition_table,
    normalized_wht,
    pooled_mean_values,
    popcount_array,
    wht,
)

BUDGET_FRACTIONS = [0.01, 0.03, 0.05, 0.10, 0.15, 0.20]
N_TRIALS = 10
CV_FOLDS = 5


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def build_walsh_design_matrix_no_intercept(sample_indices, n, max_order=None):
    """Build Walsh design matrix EXCLUDING T=0 (intercept column).

    The intercept is handled by sklearn's fit_intercept=True, which
    does not penalize it.
    """
    N = 2 ** n
    pc_all = popcount_array(n)

    if max_order is not None:
        col_indices = np.where((pc_all <= max_order) & (pc_all > 0))[0]
    else:
        col_indices = np.where(pc_all > 0)[0]

    m = len(sample_indices)
    p = len(col_indices)
    Phi = np.empty((m, p), dtype=np.float64)

    for col_idx, T in enumerate(tqdm(col_indices, desc="Building design matrix", leave=False)):
        bits = sample_indices & T
        pc = np.zeros(m, dtype=np.int32)
        tmp = bits.copy()
        while tmp.any():
            pc += (tmp & 1).astype(np.int32)
            tmp >>= 1
        Phi[:, col_idx] = (-1.0) ** pc

    return Phi, col_indices


def lasso_recover(v, n, m, seed=42, max_order=None):
    N = 2 ** n
    rng = np.random.default_rng(seed)

    w_true = normalized_wht(v)

    all_indices = rng.permutation(N)
    sample_indices = np.sort(all_indices[:m])
    heldout_indices = np.sort(all_indices[m:m + min(m, N - m)])
    y_sample = v[sample_indices]

    Phi, col_indices = build_walsh_design_matrix_no_intercept(
        sample_indices, n, max_order
    )

    lasso = LassoCV(
        cv=CV_FOLDS,
        max_iter=10000,
        n_jobs=-1,
        fit_intercept=True,
    )
    lasso.fit(Phi, y_sample)

    w_recovered = np.zeros(N, dtype=np.float64)
    w_recovered[0] = lasso.intercept_
    w_recovered[col_indices] = lasso.coef_

    v_pred = np.zeros(N, dtype=np.float64)
    Phi_full, _ = build_walsh_design_matrix_no_intercept(
        np.arange(N, dtype=np.int64), n, max_order
    )
    v_pred = Phi_full @ lasso.coef_ + lasso.intercept_

    metrics = compute_recovery_metrics(
        w_recovered, w_true, n,
        v_pred=v_pred, v_true=v,
        heldout_indices=heldout_indices,
    )
    metrics["alpha"] = float(lasso.alpha_)
    metrics["n_nonzero"] = int(np.sum(np.abs(lasso.coef_) > 1e-10))
    metrics["m"] = m
    metrics["sample_fraction"] = m / N
    metrics["seed"] = seed

    return metrics


def run_budget_sweep(v, n, budgets=None, n_trials=N_TRIALS, max_order=None):
    if budgets is None:
        budgets = BUDGET_FRACTIONS

    N = 2 ** n
    results = []

    for frac in budgets:
        m = max(int(frac * N), n + 1)
        print(f"\n[{timestamp()}] Budget {frac*100:.0f}% — m={m} samples, {n_trials} trials")

        trial_results = []
        for trial in tqdm(range(n_trials), desc=f"Budget {frac*100:.0f}%"):
            seed = 42 + trial
            r = lasso_recover(v, n, m, seed=seed, max_order=max_order)
            r["trial"] = trial
            trial_results.append(r)

        def mean_std(key):
            vals = [r[key] for r in trial_results if r[key] is not None]
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
    parser = argparse.ArgumentParser(description="LASSO Walsh oracle ceiling")
    parser.add_argument("--data", required=True, help="Path to coalition table npz")
    parser.add_argument("--output", required=True, help="Path to save results JSON")
    parser.add_argument("--circuit", default="unknown", help="Circuit name for labeling")
    parser.add_argument("--max-order", type=int, default=None,
                        help="Max interaction order in Walsh basis (None = full)")
    parser.add_argument("--n-trials", type=int, default=N_TRIALS)
    parser.add_argument("--budgets", type=float, nargs="+", default=None)
    args = parser.parse_args()

    print(f"[{timestamp()}] Loading coalition table: {args.data}")
    table = load_coalition_table(args.data)
    n = table["n_players"]
    v = pooled_mean_values(table)
    print(f"  n_players={n}, n_prompts={table['n_prompts']}, N={2**n}")

    w_true = normalized_wht(v)
    e_true = energy_spectrum(w_true, n)
    print(f"  Ground truth energy spectrum:")
    for k in range(min(6, n + 1)):
        print(f"    order-{k}: {e_true[k]*100:.2f}%")
    print(f"    order-3+: {e_true[3:].sum()*100:.2f}%")

    budgets = args.budgets if args.budgets else BUDGET_FRACTIONS
    results = run_budget_sweep(v, n, budgets=budgets, n_trials=args.n_trials,
                                max_order=args.max_order)

    output = {
        "method": "lasso_walsh_oracle",
        "circuit": args.circuit,
        "data_path": str(args.data),
        "n_players": n,
        "n_prompts": table["n_prompts"],
        "max_order": args.max_order,
        "ground_truth_energy_spectrum": e_true.tolist(),
        "timestamp": timestamp(),
        "budget_sweep": results,
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n[{timestamp()}] Results saved to {args.output}")


if __name__ == "__main__":
    main()
