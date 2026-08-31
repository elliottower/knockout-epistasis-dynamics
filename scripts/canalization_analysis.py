"""Canalization analysis of Boolean GRN rules.

Computes per-network:
  - Canalization depth (max canalizing inputs per rule)
  - Nested canalizing depth (max nested canalizing inputs per rule)
  - Mean rule bias (fraction of truth table that is 1)
  - Derrida parameter (sensitivity to single-bit perturbations)

Then correlates each feature with Δ₃₊ across all 27 networks
using Spearman rank correlation with Bonferroni correction.

All computation is from truth tables only — no simulation.

Usage:
    uv run python scripts/canalization_analysis.py
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))

from grn_coalition_sweep import BUILTIN_MODELS, compile_network
from scripts.run_batch2_blind_sweep import WEB_MODELS
from scripts.run_batch2b_extra_models import EXTRA_MODELS


RESULTS_DIR = Path(__file__).parent.parent / "results" / "canalization"


def is_canalizing(truth_table, k, input_idx):
    """Check if input_idx is a canalizing input for the truth table.

    A Boolean function f is canalizing in input x_i if there exists
    a canalizing value a and canalizing output b such that
    f(x_1, ..., x_i=a, ..., x_k) = b for all other inputs.

    Returns (is_canal, canal_value, canal_output) or (False, None, None).
    """
    n = 2**k
    for canal_val in [0, 1]:
        outputs_when_canal = []
        for row in range(n):
            bit = (row >> input_idx) & 1
            if bit == canal_val:
                outputs_when_canal.append(truth_table[row])
        if len(set(outputs_when_canal)) == 1:
            return True, canal_val, outputs_when_canal[0]
    return False, None, None


def canalizing_depth(truth_table, k):
    """Count how many inputs are canalizing for this truth table."""
    count = 0
    for i in range(k):
        is_canal, _, _ = is_canalizing(truth_table, k, i)
        if is_canal:
            count += 1
    return count


def nested_canalizing_depth(truth_table, k):
    """Compute nested canalizing depth.

    Repeatedly find a canalizing input, fix it to the canalizing value,
    reduce the truth table, and repeat. The depth is the number of
    layers before no more canalizing inputs exist.
    """
    if k == 0:
        return 0

    remaining_inputs = list(range(k))
    current_tt = np.array(truth_table, dtype=np.int8)
    current_k = k
    depth = 0

    while remaining_inputs:
        found = False
        for idx_pos, input_idx in enumerate(remaining_inputs):
            is_canal, canal_val, canal_out = is_canalizing(
                current_tt, current_k, idx_pos
            )
            if is_canal:
                non_canal_val = 1 - canal_val
                reduced_tt = []
                for row in range(2**current_k):
                    bit = (row >> idx_pos) & 1
                    if bit == non_canal_val:
                        reduced_tt.append(current_tt[row])
                current_tt = np.array(reduced_tt, dtype=np.int8)
                remaining_inputs.pop(idx_pos)
                current_k -= 1
                depth += 1
                found = True
                break
        if not found:
            break

    return depth


def rule_bias(truth_table):
    """Fraction of truth table entries that are 1."""
    if len(truth_table) == 0:
        return 0.5
    return float(np.mean(truth_table))


def derrida_parameter(compiled):
    """Compute the Derrida parameter for the network.

    The Derrida parameter measures sensitivity to single-bit
    perturbations. For each node's truth table, compute the probability
    that flipping one input changes the output. Average over all nodes
    and all inputs.

    Lambda > 1 indicates chaotic regime; lambda < 1 indicates ordered.
    """
    n = len(compiled)
    if n == 0:
        return 0.0

    total_sensitivity = 0.0
    total_inputs = 0

    for node_name, reg_indices, truth_table in compiled:
        k = len(reg_indices)
        if k == 0:
            continue
        for i in range(k):
            flips = 0
            for row in range(2**k):
                flipped_row = row ^ (1 << i)
                if truth_table[row] != truth_table[flipped_row]:
                    flips += 1
            sensitivity = flips / (2**k)
            total_sensitivity += sensitivity
            total_inputs += 1

    if total_inputs == 0:
        return 0.0
    return total_sensitivity / total_inputs * n


def analyze_one_network(name, rules, output_nodes):
    """Compute all canalization features for one network."""
    compiled, node_names = compile_network(rules)
    n = len(node_names)

    canal_depths = []
    nested_depths = []
    biases = []

    for node_name, reg_indices, truth_table in compiled:
        k = len(reg_indices)
        if k == 0:
            biases.append(float(truth_table[0]) if len(truth_table) > 0 else 0.5)
            canal_depths.append(0)
            nested_depths.append(0)
            continue

        canal_depths.append(canalizing_depth(truth_table, k))
        nested_depths.append(nested_canalizing_depth(truth_table, k))
        biases.append(rule_bias(truth_table))

    derrida = derrida_parameter(compiled)

    mean_canal = float(np.mean(canal_depths)) if canal_depths else 0.0
    max_canal = max(canal_depths) if canal_depths else 0
    frac_any_canal = float(np.mean([d > 0 for d in canal_depths])) if canal_depths else 0.0

    mean_nested = float(np.mean(nested_depths)) if nested_depths else 0.0
    max_nested = max(nested_depths) if nested_depths else 0
    frac_fully_canal = float(np.mean([
        nested_depths[i] == len(compiled[i][1])
        for i in range(len(compiled))
        if len(compiled[i][1]) > 0
    ])) if any(len(c[1]) > 0 for c in compiled) else 0.0

    mean_bias = float(np.mean(biases)) if biases else 0.5
    bias_variance = float(np.var(biases)) if biases else 0.0
    mean_bias_deviation = float(np.mean([abs(b - 0.5) for b in biases])) if biases else 0.0

    return {
        "model": name,
        "n": n,
        "mean_canalizing_depth": mean_canal,
        "max_canalizing_depth": max_canal,
        "frac_any_canalizing": frac_any_canal,
        "mean_nested_depth": mean_nested,
        "max_nested_depth": max_nested,
        "frac_fully_canalizing": frac_fully_canal,
        "mean_bias": mean_bias,
        "bias_variance": bias_variance,
        "mean_bias_deviation": mean_bias_deviation,
        "derrida_parameter": derrida,
        "per_node_canal_depth": canal_depths,
        "per_node_nested_depth": nested_depths,
        "per_node_bias": biases,
    }


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_models = {}
    for src in [BUILTIN_MODELS, WEB_MODELS, EXTRA_MODELS]:
        all_models.update(src)

    merged_path = Path(__file__).parent.parent / "results" / "grn_v2" / "merged_all_27_analysis.json"
    with open(merged_path) as f:
        merged = json.load(f)
    real_deltas = {m["model"]: m.get("delta_o3plus", 0) * 100 for m in merged["all_models"]}

    results = {}
    for name in sorted(all_models.keys()):
        info = all_models[name]
        n = len(info["rules"])
        print(f"  {name} (n={n})...", end=" ", flush=True)
        result = analyze_one_network(name, info["rules"], info["output_nodes"])
        result["real_delta_o3plus"] = real_deltas.get(name, None)
        results[name] = result
        print(f"canal={result['mean_canalizing_depth']:.2f}, "
              f"nested={result['mean_nested_depth']:.2f}, "
              f"derrida={result['derrida_parameter']:.2f}, "
              f"bias_dev={result['mean_bias_deviation']:.3f}")

    with open(RESULTS_DIR / "canalization_features.json", "w") as f:
        json.dump(results, f, indent=2)

    features = [
        ("mean_canalizing_depth", "Mean canalizing depth"),
        ("max_canalizing_depth", "Max canalizing depth"),
        ("frac_any_canalizing", "Fraction with any canalizing input"),
        ("mean_nested_depth", "Mean nested canalizing depth"),
        ("max_nested_depth", "Max nested canalizing depth"),
        ("frac_fully_canalizing", "Fraction fully canalizing"),
        ("mean_bias", "Mean rule bias"),
        ("bias_variance", "Bias variance"),
        ("mean_bias_deviation", "Mean |bias - 0.5|"),
        ("derrida_parameter", "Derrida parameter"),
        ("n", "Network size"),
    ]

    models_with_delta = [
        name for name, r in results.items()
        if r.get("real_delta_o3plus") is not None
    ]

    print(f"\n{'='*70}")
    print(f"Spearman correlations with Δ₃₊ (n={len(models_with_delta)} networks)")
    print(f"{'='*70}")
    print(f"{'Feature':40s} {'rho':>8s} {'p-value':>10s} {'p*k':>10s} {'sig':>5s}")
    print("-" * 75)

    n_tests = len(features)
    deltas = np.array([results[m]["real_delta_o3plus"] for m in models_with_delta])

    corr_results = []
    for feat_key, feat_name in features:
        feat_vals = np.array([results[m][feat_key] for m in models_with_delta])
        if np.std(feat_vals) < 1e-10:
            print(f"{feat_name:40s} {'N/A':>8s} {'N/A':>10s} {'N/A':>10s} {'':>5s}")
            continue
        rho, p = stats.spearmanr(feat_vals, deltas)
        p_bonf = min(p * n_tests, 1.0)
        sig = "*" if p_bonf < 0.05 else ""
        print(f"{feat_name:40s} {rho:>+8.3f} {p:>10.4f} {p_bonf:>10.4f} {sig:>5s}")
        corr_results.append({
            "feature": feat_key,
            "feature_name": feat_name,
            "rho": float(rho),
            "p_value": float(p),
            "p_bonferroni": float(p_bonf),
            "significant": p_bonf < 0.05,
        })

    summary = {
        "n_models": len(models_with_delta),
        "n_tests": n_tests,
        "correction": "bonferroni",
        "alpha": 0.05,
        "correlations": corr_results,
        "models": models_with_delta,
    }
    with open(RESULTS_DIR / "canalization_correlations.json", "w") as f:
        json.dump(summary, f, indent=2, default=lambda o: bool(o) if isinstance(o, np.bool_) else float(o))

    print(f"\nResults saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
