"""Paired budget-sweep harness for EpistasisBench.

Runs LASSO-Walsh and iRF on the SAME train/heldout splits across
multiple budget levels, seeds, and circuits. This is the experimental
unit of the benchmark — every future method plugs into this harness.

Paired seeds: at each (circuit, budget, seed), both methods receive
exactly the same sampled and held-out coalitions.

Usage:
    uv run python budget_sweep_harness.py \
        --circuits weight_ioi ioi random15 \
        --data-dir /path/to/genetics/ \
        --output results/budget_sweep.json

    uv run python budget_sweep_harness.py \
        --circuits weight_ioi \
        --data-dir /path/to/genetics/ \
        --budgets 0.01 0.05 0.10 \
        --n-seeds 5 \
        --output results/pilot_sweep.json
"""

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from tqdm import tqdm

from data_utils import (
    compute_k99,
    compute_recovery_metrics,
    energy_spectrum,
    load_coalition_table,
    normalized_wht,
    pooled_mean_values,
    popcount_array,
    wht,
)

CIRCUIT_FILES = {
    "weight_ioi": "weight_ioi_zero_v2_coalition_values.npz",
    "ioi": "ioi_zero_v2_coalition_values.npz",
    "random15": "random15_zero_v2_coalition_values.npz",
}

CIRCUIT_FILES_MEAN = {
    "weight_ioi": "weight_ioi_mean_v2_coalition_values.npz",
    "ioi": "ioi_mean_v2_coalition_values.npz",
    "random15": "random15_mean_v2_coalition_values.npz",
}

DEFAULT_BUDGETS = [0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.40]
DEFAULT_N_SEEDS = 20

LASSO_MAX_ORDER = 4
LASSO_CV_FOLDS = 5
IRF_N_ITERATIONS = 5
IRF_N_TREES = 500


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def coalition_to_features(n):
    N = 2 ** n
    X = np.zeros((N, n), dtype=np.float64)
    for i in range(n):
        X[:, i] = (np.arange(N) >> i) & 1
    return X


def build_walsh_design_no_intercept(sample_indices, n, max_order=None):
    N = 2 ** n
    pc_all = popcount_array(n)
    if max_order is not None:
        col_indices = np.where((pc_all <= max_order) & (pc_all > 0))[0]
    else:
        col_indices = np.where(pc_all > 0)[0]

    m = len(sample_indices)
    Phi = np.empty((m, len(col_indices)), dtype=np.float64)
    for col_idx, T in enumerate(col_indices):
        bits = sample_indices & T
        pc = np.zeros(m, dtype=np.int32)
        tmp = bits.copy()
        while tmp.any():
            pc += (tmp & 1).astype(np.int32)
            tmp >>= 1
        Phi[:, col_idx] = (-1.0) ** pc

    return Phi, col_indices


def run_lasso(v, n, sample_indices, heldout_indices, max_order=LASSO_MAX_ORDER):
    from sklearn.linear_model import LassoCV

    N = 2 ** n
    w_true = normalized_wht(v)
    y_sample = v[sample_indices]

    t0 = time.monotonic()

    Phi, col_indices = build_walsh_design_no_intercept(
        sample_indices, n, max_order
    )

    lasso = LassoCV(
        cv=LASSO_CV_FOLDS,
        max_iter=10000,
        n_jobs=-1,
        fit_intercept=True,
    )
    lasso.fit(Phi, y_sample)

    w_recovered = np.zeros(N, dtype=np.float64)
    w_recovered[0] = lasso.intercept_
    w_recovered[col_indices] = lasso.coef_

    Phi_full, _ = build_walsh_design_no_intercept(
        np.arange(N, dtype=np.int64), n, max_order
    )
    v_pred = Phi_full @ lasso.coef_ + lasso.intercept_

    elapsed = time.monotonic() - t0

    metrics = compute_recovery_metrics(
        w_recovered, w_true, n,
        v_pred=v_pred, v_true=v,
        heldout_indices=heldout_indices,
    )
    metrics["method"] = "lasso_walsh"
    metrics["wall_seconds"] = round(elapsed, 2)
    metrics["alpha"] = float(lasso.alpha_)
    metrics["n_nonzero"] = int(np.sum(np.abs(lasso.coef_) > 1e-10))

    alpha_path = lasso.alphas_
    metrics["alpha_at_boundary"] = bool(
        lasso.alpha_ == alpha_path[0] or lasso.alpha_ == alpha_path[-1]
    )
    metrics["alpha_path_min"] = float(alpha_path[-1])
    metrics["alpha_path_max"] = float(alpha_path[0])

    return metrics


IRF_LEAF_GRID = [1, 2, 5, 10, 20]
IRF_DEPTH_GRID = [None, 6, 10]
IRF_CV_FOLDS = 3


def tune_irf_hyperparams(X_train, y_train, n_trees, seed):
    """Nested CV to select min_samples_leaf and max_depth.

    Uses only the training data — never touches held-out coalitions.
    Returns (best_leaf, best_depth, cv_scores_dict).
    """
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import KFold

    kf = KFold(n_splits=IRF_CV_FOLDS, shuffle=True, random_state=seed)
    best_score = -np.inf
    best_leaf = 1
    best_depth = None
    cv_scores = {}

    for leaf in IRF_LEAF_GRID:
        for depth in IRF_DEPTH_GRID:
            scores = []
            for train_idx, val_idx in kf.split(X_train):
                rf = RandomForestRegressor(
                    n_estimators=min(n_trees, 100),
                    max_features=None,
                    min_samples_leaf=leaf,
                    max_depth=depth,
                    random_state=seed,
                    n_jobs=-1,
                )
                rf.fit(X_train[train_idx], y_train[train_idx])
                scores.append(rf.score(X_train[val_idx], y_train[val_idx]))
            mean_score = np.mean(scores)
            cv_scores[f"leaf={leaf}_depth={depth}"] = round(mean_score, 4)
            if mean_score > best_score:
                best_score = mean_score
                best_leaf = leaf
                best_depth = depth

    return best_leaf, best_depth, cv_scores


def run_irf(v, n, sample_indices, heldout_indices, n_iterations=IRF_N_ITERATIONS,
            n_trees=IRF_N_TREES, seed=42):
    from sklearn.ensemble import RandomForestRegressor

    N = 2 ** n
    w_true = normalized_wht(v)
    X_full = coalition_to_features(n)
    X_train = X_full[sample_indices]
    y_train = v[sample_indices]
    n_features = n

    t0 = time.monotonic()

    best_leaf, best_depth, cv_scores = tune_irf_hyperparams(
        X_train, y_train, n_trees, seed
    )

    weights = np.ones(n_features)

    for it in range(n_iterations):
        n_virtual = max(n_features * 5, 50)
        w = weights / weights.sum()
        rng = np.random.default_rng(seed + it + 1000)
        virtual_idx = rng.choice(n_features, size=n_virtual, p=w)

        rf = RandomForestRegressor(
            n_estimators=n_trees,
            max_features=None,
            min_samples_leaf=best_leaf,
            max_depth=best_depth,
            random_state=seed + it,
            n_jobs=-1,
        )
        rf.fit(X_train[:, virtual_idx], y_train)

        raw_imp = rf.feature_importances_
        importances = np.zeros(n_features)
        for vi, orig in enumerate(virtual_idx):
            importances[orig] += raw_imp[vi]
        total = importances.sum()
        if total > 0:
            importances /= total
        weights = importances + 1e-10

    v_pred = rf.predict(X_full[:, virtual_idx])
    w_recovered = normalized_wht(v_pred)

    elapsed = time.monotonic() - t0

    metrics = compute_recovery_metrics(
        w_recovered, w_true, n,
        v_pred=v_pred, v_true=v,
        heldout_indices=heldout_indices,
    )
    metrics["method"] = "irf"
    metrics["wall_seconds"] = round(elapsed, 2)
    metrics["train_r2"] = float(rf.score(X_train[:, virtual_idx], y_train))
    metrics["feature_importances"] = importances.tolist()
    metrics["tuned_min_leaf"] = best_leaf
    metrics["tuned_max_depth"] = best_depth
    metrics["cv_scores"] = cv_scores

    return metrics


def run_paired_trial(v, n, budget_frac, seed):
    """Run both methods on the same train/heldout split."""
    N = 2 ** n
    m = max(int(budget_frac * N), n + 1)
    rng = np.random.default_rng(seed)

    all_indices = rng.permutation(N)
    sample_indices = np.sort(all_indices[:m])
    heldout_indices = np.sort(all_indices[m:m + min(m, N - m)])

    lasso_result = run_lasso(v, n, sample_indices, heldout_indices)
    lasso_result["m"] = m
    lasso_result["budget_fraction"] = budget_frac
    lasso_result["seed"] = seed

    irf_result = run_irf(v, n, sample_indices, heldout_indices, seed=seed)
    irf_result["m"] = m
    irf_result["budget_fraction"] = budget_frac
    irf_result["seed"] = seed

    return lasso_result, irf_result


def summarize_trials(trials, method_name):
    """Compute median and 95% interval for key metrics."""
    def stats(key):
        vals = [t[key] for t in trials if t.get(key) is not None and np.isfinite(t[key])]
        if not vals:
            return {"median": None, "p2.5": None, "p97.5": None}
        vals = sorted(vals)
        return {
            "median": float(np.median(vals)),
            "p2.5": float(np.percentile(vals, 2.5)),
            "p97.5": float(np.percentile(vals, 97.5)),
            "mean": float(np.mean(vals)),
            "std": float(np.std(vals)),
        }

    return {
        "method": method_name,
        "n_trials": len(trials),
        "heldout_r2": stats("heldout_r2"),
        "spectrum_nmse_no_intercept": stats("spectrum_nmse_no_intercept"),
        "pairwise_spearman_rho": stats("pairwise_spearman_rho"),
        "marginal_spearman_rho": stats("marginal_spearman_rho"),
        "energy_spectrum_corr": stats("energy_spectrum_corr"),
        "wall_seconds": stats("wall_seconds"),
    }


def run_circuit(circuit_name, data_path, budgets, n_seeds):
    print(f"\n{'='*60}")
    print(f"[{timestamp()}] Circuit: {circuit_name}")
    print(f"  Data: {data_path}")

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
        m = max(int(frac * N), n + 1)
        n_seeds_this = n_seeds if frac <= 0.10 else max(n_seeds // 2, 5)

        print(f"\n  Budget {frac*100:.0f}% (m={m}), {n_seeds_this} seeds")

        lasso_trials = []
        irf_trials = []

        for seed in tqdm(range(n_seeds_this), desc=f"  {frac*100:.0f}%"):
            lr, ir = run_paired_trial(v, n, frac, seed=42 + seed)
            lasso_trials.append(lr)
            irf_trials.append(ir)

        lasso_summary = summarize_trials(lasso_trials, "lasso_walsh")
        irf_summary = summarize_trials(irf_trials, "irf")

        budget_result = {
            "budget_fraction": frac,
            "m": m,
            "n_seeds": n_seeds_this,
            "summaries": {
                "lasso_walsh": lasso_summary,
                "irf": irf_summary,
            },
            "trials": {
                "lasso_walsh": lasso_trials,
                "irf": irf_trials,
            },
        }

        lr_r2 = lasso_summary["heldout_r2"]["median"]
        ir_r2 = irf_summary["heldout_r2"]["median"]
        lr_nmse = lasso_summary["spectrum_nmse_no_intercept"]["median"]
        ir_nmse = irf_summary["spectrum_nmse_no_intercept"]["median"]

        print(f"    LASSO: R²={lr_r2:.4f}, NMSE={lr_nmse:.4f}")
        print(f"    iRF:   R²={ir_r2:.4f}, NMSE={ir_nmse:.4f}")

        circuit_results["budgets"].append(budget_result)

    return circuit_results


def main():
    parser = argparse.ArgumentParser(description="EpistasisBench paired budget sweep")
    parser.add_argument("--circuits", nargs="+", required=True,
                        choices=list(CIRCUIT_FILES.keys()))
    parser.add_argument("--data-dir", required=True,
                        help="Directory containing coalition table npz files")
    parser.add_argument("--output", required=True)
    parser.add_argument("--budgets", type=float, nargs="+", default=DEFAULT_BUDGETS)
    parser.add_argument("--n-seeds", type=int, default=DEFAULT_N_SEEDS)
    parser.add_argument("--ablation", choices=["zero", "mean"], default="zero")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    file_map = CIRCUIT_FILES_MEAN if args.ablation == "mean" else CIRCUIT_FILES

    all_results = {
        "harness_version": "1.0",
        "ablation": args.ablation,
        "budgets": args.budgets,
        "n_seeds": args.n_seeds,
        "lasso_config": {
            "max_order": LASSO_MAX_ORDER,
            "cv_folds": LASSO_CV_FOLDS,
        },
        "irf_config": {
            "n_iterations": IRF_N_ITERATIONS,
            "n_trees": IRF_N_TREES,
            "tuning": "nested 3-fold CV on training set",
            "leaf_grid": IRF_LEAF_GRID,
            "depth_grid": [str(d) for d in IRF_DEPTH_GRID],
            "tuning_n_trees": 100,
        },
        "timestamp_start": timestamp(),
        "circuits": [],
    }

    for circuit in args.circuits:
        data_path = data_dir / file_map[circuit]
        if not data_path.exists():
            print(f"WARNING: {data_path} not found, skipping {circuit}")
            continue
        result = run_circuit(circuit, data_path, args.budgets, args.n_seeds)
        all_results["circuits"].append(result)

    all_results["timestamp_end"] = timestamp()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n[{timestamp()}] All results saved to {args.output}")


if __name__ == "__main__":
    main()
