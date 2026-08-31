"""Robustness check: max-aggregation vs sum-aggregation for local pairwise magnitudes.

The default scorer sums |local_w| across all genes for each pair.
This script tests whether using max instead changes the Spearman rho results.
"""
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_utils import normalized_wht


RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "grn_v2"

ALL_MODELS = [
    "faure_cellcycle", "tournier_apoptosis", "davidich_yeast",
    "drosophila_cellcycle", "fanconi_anemia", "arabidopsis_cellcycle",
    "lambda_phage", "arellano_rootstem", "asymmetric_cell_division",
    "cell_cycle_transcription", "remy_p53_mdm2", "albert_segment_polarity",
    "blood_stem_cell", "calzone_cellfate_reduced", "li_budding_yeast",
    "myeloid_progenitors", "pair_rule_module", "emt_switch",
    "morphogenetic_checkpoint", "zanudo_tlgl", "lac_operon",
    "mendoza_thelper", "saadatpour_guardcell", "fumia_cellcycle",
    "hematopoiesis_aging", "irons_cardiac", "calzone_cell_fate",
]


def local_pairwise_max(rule_fourier):
    """Local pairwise magnitudes using MAX instead of SUM across genes."""
    strengths = defaultdict(float)
    for gene, info in rule_fourier["per_gene"].items():
        for interaction in info.get("interactions", []):
            if interaction["order"] == 2:
                key = frozenset(interaction["regulators"])
                strengths[key] = max(strengths[key], abs(interaction["coefficient"]))
    return dict(strengths)


def local_pairwise_sum(rule_fourier):
    """Local pairwise magnitudes using SUM (the default method)."""
    strengths = defaultdict(float)
    for gene, info in rule_fourier["per_gene"].items():
        for interaction in info.get("interactions", []):
            if interaction["order"] == 2:
                key = frozenset(interaction["regulators"])
                strengths[key] += abs(interaction["coefficient"])
    return dict(strengths)


def global_pairwise_coefficients(v_mean, n, node_names):
    """Extract global Walsh pairwise coefficients."""
    w = normalized_wht(v_mean)
    pairs = {}
    for i, j in combinations(range(n), 2):
        T = (1 << i) | (1 << j)
        pairs[(node_names[i], node_names[j])] = float(w[T])
    return pairs


def run_model(model_name):
    wiring_path = RESULTS_DIR / f"{model_name}_wiring_blind.json"
    coalition_path = RESULTS_DIR / f"{model_name}_coalition_blind.npz"

    if not wiring_path.exists() or not coalition_path.exists():
        return None

    with open(wiring_path) as f:
        wiring = json.load(f)

    npz = np.load(coalition_path, allow_pickle=True)
    n = int(npz["n_players"])
    node_names = list(npz["circuit_heads"])
    v_mean = npz["target_logits"].mean(axis=1)

    rule_fourier = wiring["rule_fourier"]

    local_sum = local_pairwise_sum(rule_fourier)
    local_max = local_pairwise_max(rule_fourier)
    global_pw = global_pairwise_coefficients(v_mean, n, node_names)

    all_pairs = list(combinations(node_names, 2))
    global_vec = np.array([abs(global_pw[(a, b)]) for a, b in all_pairs])

    sum_vec = np.array([local_sum.get(frozenset([a, b]), 0.0) for a, b in all_pairs])
    max_vec = np.array([local_max.get(frozenset([a, b]), 0.0) for a, b in all_pairs])

    rho_sum, p_sum = spearmanr(sum_vec, global_vec)
    rho_max, p_max = spearmanr(max_vec, global_vec)

    return {
        "model": model_name,
        "n": n,
        "rho_sum": float(rho_sum),
        "p_sum": float(p_sum),
        "rho_max": float(rho_max),
        "p_max": float(p_max),
        "rho_diff": float(rho_max - rho_sum),
        "sign_flip": bool((rho_sum > 0) != (rho_max > 0)),
    }


if __name__ == "__main__":
    results = []
    for model in ALL_MODELS:
        r = run_model(model)
        if r is None:
            print(f"SKIP {model}: missing files")
            continue
        results.append(r)
        sign = "FLIP" if r["sign_flip"] else ""
        print(f"{model:35s}  sum={r['rho_sum']:+.3f}  max={r['rho_max']:+.3f}  diff={r['rho_diff']:+.3f}  {sign}")

    rho_sums = [r["rho_sum"] for r in results]
    rho_maxs = [r["rho_max"] for r in results]
    n_flips = sum(1 for r in results if r["sign_flip"])

    print(f"\n{'='*60}")
    print(f"SUMMARY (N={len(results)})")
    print(f"{'='*60}")
    print(f"Median rho (sum): {np.median(rho_sums):+.3f}")
    print(f"Median rho (max): {np.median(rho_maxs):+.3f}")
    print(f"Mean rho (sum):   {np.mean(rho_sums):+.3f}")
    print(f"Mean rho (max):   {np.mean(rho_maxs):+.3f}")
    print(f"Sign flips:       {n_flips}/{len(results)}")
    print(f"Median diff:      {np.median([r['rho_diff'] for r in results]):+.3f}")

    n_sig_sum = sum(1 for r in results if r["p_sum"] < 0.05)
    n_sig_max = sum(1 for r in results if r["p_max"] < 0.05)
    print(f"Significant (sum): {n_sig_sum}/{len(results)}")
    print(f"Significant (max): {n_sig_max}/{len(results)}")

    output = {
        "description": "Max-aggregation vs sum-aggregation robustness check",
        "summary": {
            "n_models": len(results),
            "median_rho_sum": float(np.median(rho_sums)),
            "median_rho_max": float(np.median(rho_maxs)),
            "mean_rho_sum": float(np.mean(rho_sums)),
            "mean_rho_max": float(np.mean(rho_maxs)),
            "n_sign_flips": n_flips,
            "n_significant_sum": n_sig_sum,
            "n_significant_max": n_sig_max,
            "median_diff": float(np.median([r["rho_diff"] for r in results])),
        },
        "per_model": results,
    }

    out_path = RESULTS_DIR / "robustness_max_aggregation.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {out_path}")
