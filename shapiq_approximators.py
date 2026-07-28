"""Shapiq approximator wrappers for EpistasisBench.

Wraps KernelSHAP-IQ, SHAPIQ (Monte Carlo), and SVARM-IQ as budget-limited
interaction estimators. Each approximator chooses its own coalitions (smart
sampling), unlike LASSO/iRF which receive a random sample.

Scored by reconstructing v(S) for all 2^n coalitions via plug-in prediction:
  v_hat(S) = baseline + sum_{T subseteq S} phi_T

This reconstruction is exact at max_order=n and approximate at lower orders.
"""

import time
from itertools import combinations

import numpy as np
import shapiq

from data_utils import (
    compute_k99,
    compute_recovery_metrics,
    normalized_wht,
    popcount_array,
)

SHAPIQ_MAX_ORDER = 4
SHAPIQ_INDEX = "k-SII"

APPROXIMATORS = {
    "kernelshapiq": shapiq.KernelSHAPIQ,
    "shapiq_mc": shapiq.SHAPIQ,
    "svarmiq": shapiq.SVARMIQ,
}


def make_game_function(v, n):
    """Create a callable game from the full coalition value vector.

    The shapiq API calls game(coalitions) where coalitions is a
    boolean array of shape (batch, n_players). Returns the values.
    """
    powers = 1 << np.arange(n)

    def game_fn(coalitions):
        indices = (np.asarray(coalitions, dtype=bool) * powers).sum(axis=1).astype(int)
        return v[indices]

    return game_fn


def reconstruct_from_interaction_values(iv, n):
    """Reconstruct all 2^n coalition values from InteractionValues.

    Uses the plug-in predictor: v_hat(S) = baseline + sum_{T⊆S, |T|>0} phi_T
    """
    N = 2 ** n
    v_pred = np.full(N, iv.baseline_value, dtype=np.float64)

    subset_vals = {}
    for T_tuple, val in iv.dict_values.items():
        if len(T_tuple) == 0:
            continue
        bitmask = 0
        for j in T_tuple:
            bitmask |= (1 << j)
        subset_vals[bitmask] = val

    for s_int in range(N):
        for T_mask, val in subset_vals.items():
            if (s_int & T_mask) == T_mask:
                v_pred[s_int] += val

    return v_pred


def run_shapiq_approximator(v, n, budget_n_evals, method_name, seed,
                            max_order=SHAPIQ_MAX_ORDER, index=SHAPIQ_INDEX):
    """Run a single shapiq approximator and score its reconstruction.

    Args:
        v: Full coalition value vector (2^n,)
        n: Number of players
        budget_n_evals: Number of game evaluations allowed
        method_name: Key into APPROXIMATORS dict
        seed: Random seed
        max_order: Maximum interaction order
        index: Shapley interaction index type

    Returns:
        Dict of metrics compatible with compute_recovery_metrics
    """
    N = 2 ** n
    w_true = normalized_wht(v)
    game_fn = make_game_function(v, n)

    approx_cls = APPROXIMATORS[method_name]

    t0 = time.monotonic()

    approx = approx_cls(
        n=n,
        max_order=max_order,
        index=index,
        random_state=seed,
    )
    iv = approx.approximate(budget=budget_n_evals, game=game_fn)

    v_pred = reconstruct_from_interaction_values(iv, n)
    w_recovered = normalized_wht(v_pred)

    elapsed = time.monotonic() - t0

    all_indices = np.arange(N, dtype=np.int64)
    metrics = compute_recovery_metrics(
        w_recovered, w_true, n,
        v_pred=v_pred, v_true=v,
        heldout_indices=all_indices,
    )

    metrics["method"] = method_name
    metrics["wall_seconds"] = round(elapsed, 2)
    metrics["shapiq_index"] = index
    metrics["max_order"] = max_order
    metrics["budget_n_evals"] = budget_n_evals
    metrics["estimation_budget"] = int(iv.estimation_budget) if hasattr(iv, 'estimation_budget') else None
    metrics["baseline_value"] = float(iv.baseline_value)
    metrics["n_interactions"] = len([t for t in iv.dict_values if len(t) > 0])

    return metrics
