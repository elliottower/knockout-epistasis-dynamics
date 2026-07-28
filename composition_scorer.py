"""Score how much local rule Fourier structure survives dynamical composition.

The core question: do pairs with strong local rule interactions also have
strong attractor-level pairwise epistasis? And does composition create
higher-order structure that doesn't exist in any single rule?

Inputs:
  - Coalition table (.npz) from grn_coalition_sweep.py
  - Wiring JSON (contains rule_fourier with per-gene local Walsh coefficients)

Outputs:
  - Spearman correlation: local |w_{i,j}| vs global |w_{i,j}| for all pairs
  - Destruction rate: fraction of strong local interactions that vanish globally
  - Creation rate: fraction of strong global interactions with no local source
  - Per-order energy comparison: local rule spectrum vs attractor spectrum

Usage:
    uv run python composition_scorer.py \
        --coalition results/faure_cellcycle_coalition.npz \
        --wiring results/faure_cellcycle_coalition_wiring.json \
        --output results/faure_cellcycle_composition.json
"""

import argparse
import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from data_utils import energy_by_order, normalized_wht, popcount_array, pooled_mean_values


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def bootstrap_spearman_ci(x, y, n_boot=2000, alpha=0.05, seed=42):
    """Bootstrap confidence interval for Spearman rho.

    Returns (rho, ci_lo, ci_hi) where CI is a (1-alpha)*100% interval.
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    rho_point, _ = spearmanr(x, y)
    rhos = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rhos[b], _ = spearmanr(x[idx], y[idx])
    lo = np.nanpercentile(rhos, 100 * alpha / 2)
    hi = np.nanpercentile(rhos, 100 * (1 - alpha / 2))
    return float(rho_point), float(lo), float(hi)


def spectrum_gate(v_mean, n, min_abs_energy_3plus=1e-4):
    """Check whether attractor-level order-3+ structure exceeds the noise floor.

    Returns (passes, abs_energy_3plus, frac_energy_3plus).
    The gate requires absolute order-3+ energy (sum of squared coefficients)
    above a floor, not just a percentage of total variance.
    """
    w = normalized_wht(v_mean)
    pc = popcount_array(n)
    energy_3plus = float(np.sum(w[pc >= 3] ** 2))
    total_energy = float(np.sum(w[pc > 0] ** 2))
    frac = energy_3plus / total_energy if total_energy > 0 else 0.0
    return energy_3plus >= min_abs_energy_3plus, energy_3plus, frac


def global_pairwise_coefficients(v_mean, n, node_names):
    """Extract global Walsh pairwise coefficients for all C(n,2) pairs.

    Returns dict mapping (name_i, name_j) -> float coefficient.
    """
    w = normalized_wht(v_mean)
    pairs = {}
    for i, j in combinations(range(n), 2):
        T = (1 << i) | (1 << j)
        pairs[(node_names[i], node_names[j])] = float(w[T])
    return pairs


def global_triple_coefficients(v_mean, n, node_names):
    """Extract global Walsh triple coefficients for all C(n,3) triples."""
    w = normalized_wht(v_mean)
    triples = {}
    for i, j, k in combinations(range(n), 3):
        T = (1 << i) | (1 << j) | (1 << k)
        triples[(node_names[i], node_names[j], node_names[k])] = float(w[T])
    return triples


def local_pairwise_strengths(rule_fourier):
    """Extract local pairwise strengths from rule_fourier output.

    Returns dict mapping frozenset({name_i, name_j}) -> float (sum of |local w|).
    """
    strengths = {}
    for pair_list in rule_fourier["pairwise"]:
        key = frozenset(pair_list)
        strengths[key] = 0.0

    for gene, info in rule_fourier["per_gene"].items():
        for interaction in info.get("interactions", []):
            if interaction["order"] == 2:
                key = frozenset(interaction["regulators"])
                strengths[key] = strengths.get(key, 0.0) + abs(interaction["coefficient"])

    return strengths


def local_triple_strengths(rule_fourier):
    """Extract local triple strengths from rule_fourier output."""
    strengths = {}
    for gene, info in rule_fourier["per_gene"].items():
        for interaction in info.get("interactions", []):
            if interaction["order"] == 3:
                key = frozenset(interaction["regulators"])
                strengths[key] = strengths.get(key, 0.0) + abs(interaction["coefficient"])
    return strengths


def local_energy_spectrum(rule_fourier, n):
    """Get aggregate local rule energy spectrum (all orders including 0 and 1).

    Uses the pre-computed spectrum from extract_rule_fourier, which sums
    squared Walsh coefficients across all genes' truth tables at each order.
    Pads or truncates to match the global spectrum's n+1 bins.
    """
    raw = rule_fourier.get("local_energy_spectrum", [])
    spectrum = np.zeros(n + 1, dtype=np.float64)
    for k in range(min(len(raw), n + 1)):
        spectrum[k] = raw[k]
    total = spectrum.sum()
    if total > 0:
        spectrum /= total
    return spectrum


def score_composition(v_mean, n, node_names, rule_fourier, n_boot=2000):
    """Score how much local rule structure survives into attractor-level epistasis.

    Returns a dict with all metrics including bootstrap CIs and spectrum gate.
    """
    global_pw = global_pairwise_coefficients(v_mean, n, node_names)
    local_pw = local_pairwise_strengths(rule_fourier)

    all_pairs = list(combinations(node_names, 2))

    global_vec = np.array([abs(global_pw[(a, b)]) for a, b in all_pairs])
    local_vec = np.array([local_pw.get(frozenset([a, b]), 0.0) for a, b in all_pairs])

    n_pairs = len(all_pairs)
    n_local_nonzero = int(np.sum(local_vec > 0))
    n_global_nonzero = int(np.sum(global_vec > 1e-10))

    rho, pvalue = spearmanr(local_vec, global_vec)
    _, ci_lo, ci_hi = bootstrap_spearman_ci(local_vec, global_vec, n_boot=n_boot)

    top_k = max(1, n_local_nonzero)
    global_top_idx = np.argsort(global_vec)[-top_k:]
    local_top_idx = np.argsort(local_vec)[-top_k:]

    created_count = int(np.sum(local_vec[global_top_idx] == 0.0))
    destroyed_count = int(np.sum(global_vec[local_top_idx] < 1e-10))

    pair_details = []
    for idx, (a, b) in enumerate(all_pairs):
        pair_details.append({
            "pair": [a, b],
            "global_magnitude": float(global_vec[idx]),
            "local_magnitude": float(local_vec[idx]),
            "global_coefficient": float(global_pw[(a, b)]),
        })
    pair_details.sort(key=lambda x: x["global_magnitude"], reverse=True)

    global_tri = global_triple_coefficients(v_mean, n, node_names)
    local_tri = local_triple_strengths(rule_fourier)
    all_triples = list(combinations(node_names, 3))

    global_tri_vec = np.array([abs(global_tri[t]) for t in all_triples])
    local_tri_vec = np.array([local_tri.get(frozenset(t), 0.0) for t in all_triples])

    if len(all_triples) > 2:
        tri_rho, tri_pvalue = spearmanr(local_tri_vec, global_tri_vec)
        _, tri_ci_lo, tri_ci_hi = bootstrap_spearman_ci(
            local_tri_vec, global_tri_vec, n_boot=n_boot
        )
    else:
        tri_rho, tri_pvalue = float("nan"), float("nan")
        tri_ci_lo, tri_ci_hi = float("nan"), float("nan")

    n_global_tri_nonzero = int(np.sum(global_tri_vec > 1e-10))
    n_local_tri_nonzero = int(np.sum(local_tri_vec > 0))

    w_global = normalized_wht(v_mean)
    global_energy = energy_by_order(w_global, n)
    global_spectrum = global_energy / global_energy.sum() if global_energy.sum() > 0 else global_energy

    local_spectrum = local_energy_spectrum(rule_fourier, n)

    gate_passes, abs_e3plus, frac_e3plus = spectrum_gate(v_mean, n)

    return {
        "pairwise": {
            "n_pairs": n_pairs,
            "n_local_nonzero": n_local_nonzero,
            "n_global_nonzero": n_global_nonzero,
            "local_positive_rate": n_local_nonzero / n_pairs,
            "spearman_rho": float(rho),
            "spearman_pvalue": float(pvalue),
            "spearman_ci_95": [ci_lo, ci_hi],
            "creation_rate": created_count / top_k,
            "destruction_rate": destroyed_count / top_k,
            "created_count": created_count,
            "destroyed_count": destroyed_count,
            "top_k": top_k,
        },
        "triples": {
            "n_triples": len(all_triples),
            "n_local_nonzero": n_local_tri_nonzero,
            "n_global_nonzero": n_global_tri_nonzero,
            "spearman_rho": float(tri_rho) if not np.isnan(tri_rho) else None,
            "spearman_pvalue": float(tri_pvalue) if not np.isnan(tri_pvalue) else None,
            "spearman_ci_95": [tri_ci_lo, tri_ci_hi] if not np.isnan(tri_ci_lo) else None,
        },
        "spectrum_gate": {
            "passes": gate_passes,
            "abs_energy_3plus": abs_e3plus,
            "frac_energy_3plus": frac_e3plus,
        },
        "energy_spectrum": {
            "global": global_spectrum.tolist(),
            "local_rules": local_spectrum.tolist(),
            "order_labels": [f"order-{k}" for k in range(n + 1)],
        },
        "pair_details": pair_details,
    }


def main():
    parser = argparse.ArgumentParser(description="Score local-to-global Fourier composition")
    parser.add_argument("--coalition", required=True, help="Coalition table npz")
    parser.add_argument("--wiring", required=True, help="Wiring JSON with rule_fourier")
    parser.add_argument("--output", required=True, help="Output JSON path")
    args = parser.parse_args()

    print(f"[{timestamp()}] Loading coalition table: {args.coalition}")
    data = np.load(args.coalition, allow_pickle=True)
    n = int(data["n_players"])
    n_init = int(data["n_prompts"])
    node_names = list(data["circuit_heads"])
    values = data["target_logits"]
    v_mean = values.mean(axis=1)

    print(f"  n={n}, n_init={n_init}, nodes={node_names}")

    print(f"[{timestamp()}] Loading rule Fourier data: {args.wiring}")
    with open(args.wiring) as f:
        wiring_data = json.load(f)

    if "rule_fourier" not in wiring_data:
        print("  ERROR: wiring JSON has no rule_fourier — regenerate with updated grn_coalition_sweep.py")
        return

    rule_fourier = wiring_data["rule_fourier"]
    print(f"  Local: {rule_fourier['n_pairwise']} pairwise, {rule_fourier['n_triples']} triples")

    input_nodes = wiring_data.get("input_nodes", [])
    input_config = wiring_data.get("input_config")
    if input_nodes:
        print(f"  Input nodes: {input_nodes} (config={input_config})")

    print(f"\n[{timestamp()}] Scoring composition...")
    result = score_composition(v_mean, n, node_names, rule_fourier)

    pw = result["pairwise"]
    ci = pw["spearman_ci_95"]
    print(f"\n  Pairwise composition:")
    print(f"    Spearman rho:      {pw['spearman_rho']:.4f}  95% CI [{ci[0]:.4f}, {ci[1]:.4f}]  (p={pw['spearman_pvalue']:.2e})")
    print(f"    Local nonzero:     {pw['n_local_nonzero']}/{pw['n_pairs']} ({pw['local_positive_rate']:.0%})")
    print(f"    Global nonzero:    {pw['n_global_nonzero']}/{pw['n_pairs']}")
    print(f"    Creation rate:     {pw['creation_rate']:.0%} ({pw['created_count']}/{pw['top_k']} top global have no local source)")
    print(f"    Destruction rate:  {pw['destruction_rate']:.0%} ({pw['destroyed_count']}/{pw['top_k']} top local vanish globally)")

    tri = result["triples"]
    if tri["spearman_rho"] is not None:
        tri_ci = tri["spearman_ci_95"]
        ci_str = f"  95% CI [{tri_ci[0]:.4f}, {tri_ci[1]:.4f}]" if tri_ci else ""
        print(f"\n  Triple composition:")
        print(f"    Spearman rho:      {tri['spearman_rho']:.4f}{ci_str}  (p={tri['spearman_pvalue']:.2e})")
        print(f"    Local nonzero:     {tri['n_local_nonzero']}/{tri['n_triples']}")
        print(f"    Global nonzero:    {tri['n_global_nonzero']}/{tri['n_triples']}")

    gate = result["spectrum_gate"]
    status = "PASS" if gate["passes"] else "FAIL"
    print(f"\n  Spectrum gate ({status}):")
    print(f"    Absolute order-3+ energy:    {gate['abs_energy_3plus']:.6e}")
    print(f"    Fraction order-3+ of total:  {gate['frac_energy_3plus']:.4f}")

    spec = result["energy_spectrum"]
    print(f"\n  Energy spectrum (global attractor vs local rules):")
    for k in range(min(5, n + 1)):
        g = spec["global"][k]
        l = spec["local_rules"][k]
        delta = g - l
        arrow = "+" if delta > 0 else ""
        print(f"    order-{k}: global {g:.1%}  local {l:.1%}  ({arrow}{delta:.1%})")
    if n >= 3:
        g3p = sum(spec["global"][3:])
        l3p = sum(spec["local_rules"][3:])
        delta = g3p - l3p
        arrow = "+" if delta > 0 else ""
        print(f"    order-3+: global {g3p:.1%}  local {l3p:.1%}  ({arrow}{delta:.1%})")

    output = {
        "model": wiring_data.get("model", "unknown"),
        "n_players": n,
        "n_init": n_init,
        "node_names": node_names,
        "input_nodes": input_nodes,
        "input_config": input_config,
        "timestamp": timestamp(),
        **result,
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n[{timestamp()}] Saved: {args.output}")


if __name__ == "__main__":
    main()
