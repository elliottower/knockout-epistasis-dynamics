"""Shared data loading and Walsh-Hadamard utilities for EpistasisBench.

Coalition tables are npz files with keys:
  target_logits: (2^n, n_prompts) float64
  foil_logits:   (2^n, n_prompts) float64
  coalition_indices: (2^n,) int64  — bitmask index of each coalition
  circuit_heads: list of (layer, head) tuples
  n_players: int
  n_prompts: int

Value function: v(S) = target_logit(S) - foil_logit(S) (logit-diff).
"""

import numpy as np


def load_coalition_table(path):
    data = np.load(path, allow_pickle=True)
    target = data["target_logits"]
    foil = data["foil_logits"]
    n_players = int(data["n_players"])
    n_prompts = int(data["n_prompts"])
    assert target.shape[0] == 2 ** n_players, (
        f"Expected {2**n_players} coalitions, got {target.shape[0]}. "
        f"Possibly loading stub file — check path."
    )
    return {
        "target_logits": target,
        "foil_logits": foil,
        "values": target - foil,
        "circuit_heads": data["circuit_heads"],
        "n_players": n_players,
        "n_prompts": n_prompts,
    }


def pooled_mean_values(table):
    return table["values"].mean(axis=1)


def wht(values):
    """Walsh-Hadamard Transform (unnormalized). O(n * 2^n) butterfly."""
    n = int(np.log2(len(values)))
    assert len(values) == 2 ** n
    a = values.astype(np.float64).copy()
    h = 1
    for _ in range(n):
        a_view = a.reshape(-1, 2 * h)
        lo = a_view[:, :h].copy()
        hi = a_view[:, h:]
        a_view[:, :h] = lo + hi
        a_view[:, h:] = lo - hi
        h *= 2
    return a


def normalized_wht(values):
    """Normalized Walsh coefficients: w_T = (1/2^n) * WHT(v)[T]."""
    N = len(values)
    return wht(values) / N


def popcount_array(n):
    idx = np.arange(2 ** n, dtype=np.uint32)
    pc = np.zeros(2 ** n, dtype=np.int32)
    while idx.any():
        pc += (idx & 1).astype(np.int32)
        idx >>= 1
    return pc


def energy_by_order(coeffs, n):
    pc = popcount_array(n)
    c2 = coeffs.astype(np.float64) ** 2
    energy = np.zeros(n + 1)
    for k in range(n + 1):
        energy[k] = c2[pc == k].sum()
    return energy


def energy_spectrum(coeffs, n):
    e = energy_by_order(coeffs, n)
    total = e.sum()
    if total == 0:
        return e
    return e / total


def walsh_design_row(sample_index, n):
    """Compute one row of the Walsh design matrix.

    For coalition S (encoded as sample_index), the Walsh basis
    function at index T is (-1)^{popcount(S & T)}.
    Returns a (2^n,) vector.
    """
    N = 2 ** n
    T = np.arange(N, dtype=np.int64)
    bits = sample_index & T
    pc = np.zeros(N, dtype=np.int32)
    tmp = bits.copy()
    while tmp.any():
        pc += (tmp & 1).astype(np.int32)
        tmp >>= 1
    return (-1.0) ** pc


def walsh_design_matrix(sample_indices, n):
    """Build the (m, 2^n) Walsh design matrix for sampled coalitions."""
    N = 2 ** n
    m = len(sample_indices)
    T = np.arange(N, dtype=np.int64)
    Phi = np.empty((m, N), dtype=np.float64)
    for row, s in enumerate(sample_indices):
        bits = s & T
        pc = np.zeros(N, dtype=np.int32)
        tmp = bits.copy()
        while tmp.any():
            pc += (tmp & 1).astype(np.int32)
            tmp >>= 1
        Phi[row] = (-1.0) ** pc
    return Phi


def index_to_heads(index, n):
    return [b for b in range(n) if index & (1 << b)]


def compute_k99(w_true, n, exclude_order_zero=True):
    """Number of Walsh coefficients capturing 99% of energy."""
    coeffs = w_true.copy()
    if exclude_order_zero:
        coeffs[0] = 0.0
    magnitudes = np.abs(coeffs)
    sorted_mag = np.sort(magnitudes)[::-1]
    total = np.sum(sorted_mag ** 2)
    if total == 0:
        return 0
    cumulative = np.cumsum(sorted_mag ** 2) / total
    return int(np.searchsorted(cumulative, 0.99) + 1)


def compute_recovery_metrics(w_recovered, w_true, n, v_pred=None, v_true=None,
                             heldout_indices=None):
    """Shared scoring for all methods. Fair comparison requires metrics
    that no single method is natively optimizing.

    Primary metrics (method-neutral):
      heldout_r2:       R² on held-out coalitions (fair to all methods)
      heldout_mse:      MSE on held-out coalitions

    Spectrum metrics (order-0 excluded to avoid intercept domination):
      spectrum_nmse_no_intercept: NMSE excluding order-0 coefficient
      per_order_nmse:   NMSE within each interaction order separately

    Interaction ranking metrics:
      pairwise_spearman:  Spearman rho of pairwise Walsh coefficient rankings
      marginal_spearman:  Spearman rho of order-1 coefficient rankings

    Energy metrics:
      energy_spectrum_corr: correlation of energy-by-order profiles
    """
    from itertools import combinations
    from scipy.stats import spearmanr

    pc = popcount_array(n)
    results = {}

    # --- Held-out prediction (method-neutral primary metric) ---
    if v_pred is not None and v_true is not None and heldout_indices is not None:
        y_true = v_true[heldout_indices]
        y_pred = v_pred[heldout_indices]
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - y_true.mean()) ** 2)
        results["heldout_r2"] = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0
        results["heldout_mse"] = float(ss_res / len(y_true))

    # --- Spectrum NMSE with and without intercept ---
    results["spectrum_nmse"] = float(
        np.sum((w_recovered - w_true) ** 2) / np.sum(w_true ** 2)
    )

    mask_no_intercept = pc > 0
    denom_no_intercept = np.sum(w_true[mask_no_intercept] ** 2)
    if denom_no_intercept > 0:
        results["spectrum_nmse_no_intercept"] = float(
            np.sum((w_recovered[mask_no_intercept] - w_true[mask_no_intercept]) ** 2)
            / denom_no_intercept
        )
    else:
        results["spectrum_nmse_no_intercept"] = float('inf')

    # --- Per-order NMSE ---
    per_order = {}
    for k in range(n + 1):
        mask_k = pc == k
        denom_k = np.sum(w_true[mask_k] ** 2)
        if denom_k > 0:
            per_order[str(k)] = float(
                np.sum((w_recovered[mask_k] - w_true[mask_k]) ** 2) / denom_k
            )
        else:
            per_order[str(k)] = None
    results["per_order_nmse"] = per_order

    # --- Support recovery with fixed k ---
    k99 = compute_k99(w_true, n, exclude_order_zero=True)
    if k99 > 0:
        true_top = set(np.argsort(np.abs(w_true))[-k99:])
        recov_top = set(np.argsort(np.abs(w_recovered))[-k99:])
        overlap = len(true_top & recov_top)
        results["support_jaccard_k99"] = float(overlap / len(true_top | recov_top))
        results["support_overlap_k99"] = overlap
        results["k99"] = k99

    # --- Energy spectrum ---
    e_true = energy_spectrum(w_true, n)
    e_recovered = energy_spectrum(w_recovered, n)
    results["energy_true"] = e_true.tolist()
    results["energy_recovered"] = e_recovered.tolist()
    if not np.all(e_recovered == 0):
        results["energy_spectrum_corr"] = float(np.corrcoef(e_true, e_recovered)[0, 1])
    else:
        results["energy_spectrum_corr"] = 0.0

    # --- Pairwise ranking ---
    pred_pairs = []
    true_pairs = []
    for i, j in combinations(range(n), 2):
        idx = (1 << i) | (1 << j)
        pred_pairs.append(abs(w_recovered[idx]))
        true_pairs.append(abs(w_true[idx]))
    pred_pairs = np.array(pred_pairs)
    true_pairs = np.array(true_pairs)
    rho, pval = spearmanr(pred_pairs, true_pairs)
    results["pairwise_spearman_rho"] = float(rho)
    results["pairwise_spearman_pval"] = float(pval)

    # --- Pairwise detection AUROC ---
    from sklearn.metrics import roc_auc_score
    labels = (true_pairs >= np.median(true_pairs)).astype(int)
    if labels.sum() > 0 and labels.sum() < len(labels):
        results["pairwise_auroc"] = float(roc_auc_score(labels, pred_pairs))
    else:
        results["pairwise_auroc"] = None

    # --- Marginal ranking ---
    pred_order1 = [abs(w_recovered[1 << i]) for i in range(n)]
    true_order1 = [abs(w_true[1 << i]) for i in range(n)]
    m_rho, _ = spearmanr(pred_order1, true_order1)
    results["marginal_spearman_rho"] = float(m_rho)

    return results
