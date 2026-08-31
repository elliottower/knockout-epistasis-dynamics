"""Finalize grieco_bladder integration: merge Boolean shards, recompute correlations, generate paper text.

Usage:
    uv run python scripts/finalize_grieco.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))
from composition_scorer import normalized_wht, energy_by_order
from grn_coalition_sweep import compile_network, extract_rule_fourier
from scripts.run_batch2_blind_sweep import WEB_MODELS

MODEL_NAME = "grieco_bladder"
SHARD_DIR = Path("results/boolean_full/grieco_shards_v2")
ODE_PATH = Path("results/ode_full/grieco_bladder_ode.json")
MERGED_27 = Path("results/grn_v2/merged_all_27_analysis.json")
CANAL_FEATURES = Path("results/canalization/canalization_features.json")
N_SHARDS = 8
N_PLAYERS = 18
N_INIT = 512


def merge_boolean_shards():
    """Merge Boolean shards and return analysis dict."""
    arrays = []
    total_fixed = 0
    total_cycling = 0

    for i in range(N_SHARDS):
        shard_path = SHARD_DIR / f"shard_{i:03d}.npy"
        meta_path = SHARD_DIR / f"shard_{i:03d}_meta.json"
        fix_shard = SHARD_DIR / f"shard_{i:03d}_fix.npy"
        fix_meta = SHARD_DIR / f"shard_{i:03d}_fix_meta.json"

        if fix_meta.exists():
            fix_m = json.loads(fix_meta.read_text())
            if fix_m.get("status") == "complete":
                shard_path = fix_shard
                meta_path = fix_meta

        if not shard_path.exists():
            print(f"Missing shard {i}: {shard_path}")
            return None

        meta = json.loads(meta_path.read_text())
        if meta.get("status") != "complete":
            print(f"Shard {i} not complete: {meta.get('n_completed', '?')}/{meta.get('n_total', meta.get('n_coalitions', '?'))}")
            return None

        arr = np.load(shard_path)
        arrays.append(arr)
        total_fixed += meta.get("total_fixed", 0)
        total_cycling += meta.get("total_cycling", 0)

    values = np.concatenate(arrays, axis=0)
    expected = 2**N_PLAYERS
    assert values.shape[0] == expected, f"Got {values.shape[0]}, expected {expected}"

    v_mean = values.mean(axis=1)

    model_info = WEB_MODELS[MODEL_NAME]
    compiled, node_names = compile_network(model_info["rules"])
    rf = extract_rule_fourier(model_info["rules"])

    from composition_scorer import score_composition
    scores = score_composition(v_mean, N_PLAYERS, node_names, rf)

    pw = scores["pairwise"]
    spec = scores["energy_spectrum"]
    global_o3plus = sum(spec["global"][3:])
    local_o3plus = sum(spec["local_rules"][3:])
    delta = global_o3plus - local_o3plus

    total_samples = values.shape[0] * values.shape[1]
    cycling_frac = total_cycling / total_samples if total_samples > 0 else 0

    return {
        "model": MODEL_NAME,
        "description": model_info.get("description", ""),
        "n_players": N_PLAYERS,
        "n_init": N_INIT,
        "global_spectrum": spec["global"],
        "global_o3plus": global_o3plus,
        "local_spectrum": spec["local_rules"],
        "local_o3plus": local_o3plus,
        "delta_o3plus": delta,
        "n_fixed_point": total_fixed,
        "n_cycling": total_cycling,
        "cycling_fraction": cycling_frac,
        "spearman_rho": pw["spearman_rho"],
        "spearman_p": pw["spearman_pvalue"],
        "spearman_ci": pw["spearman_ci_95"],
    }


def recompute_canalization_correlations(bool_delta):
    """Recompute Spearman correlations with 28 networks."""
    features = json.loads(CANAL_FEATURES.read_text())

    features[MODEL_NAME]["real_delta_o3plus"] = bool_delta * 100

    feat_keys = [
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
        name for name, r in features.items()
        if r.get("real_delta_o3plus") is not None
    ]

    n_tests = len(feat_keys)
    deltas = np.array([features[m]["real_delta_o3plus"] for m in models_with_delta])

    results = {}
    for key, name in feat_keys:
        vals = np.array([features[m][key] for m in models_with_delta])
        if np.std(vals) < 1e-10:
            continue
        rho, p = stats.spearmanr(vals, deltas)
        p_adj = min(p * n_tests, 1.0)
        results[key] = {"name": name, "rho": rho, "p_raw": p, "p_adjusted": p_adj, "significant": p_adj < 0.05}

    return results, len(models_with_delta)


def main():
    print("=" * 60)
    print("FINALIZE GRIECO_BLADDER INTEGRATION")
    print("=" * 60)

    # 1. Merge Boolean shards
    print("\n1. Merging Boolean shards...")
    result = merge_boolean_shards()
    if result is None:
        print("SHARDS NOT READY. Run again when complete.")
        return

    bool_delta = result["delta_o3plus"]
    cycling = result["cycling_fraction"]
    rho = result["spearman_rho"]
    p_val = result["spearman_p"]

    print(f"   Boolean delta o3+: {bool_delta:+.4f} ({bool_delta*100:+.1f} pp)")
    print(f"   Cycling fraction: {cycling:.1%}")
    print(f"   Pairwise rho: {rho:.2f} (p={p_val:.4f})")

    # Save individual analysis
    out_path = Path("results/grn_v2/grieco_bladder_analysis.json")
    out_path.write_text(json.dumps(result, indent=2))

    # Save merged 28-network analysis
    merged = json.loads(MERGED_27.read_text())
    existing = {m["model"] for m in merged["all_models"]}
    if MODEL_NAME not in existing:
        merged["all_models"].append(result)
    else:
        merged["all_models"] = [m if m["model"] != MODEL_NAME else result for m in merged["all_models"]]
    merged["n_models"] = len(merged["all_models"])
    merged_28_path = Path("results/grn_v2/merged_all_28_analysis.json")
    merged_28_path.write_text(json.dumps(merged, indent=2))
    print(f"   Saved: {merged_28_path}")

    # 2. ODE comparison
    print("\n2. ODE comparison...")
    if ODE_PATH.exists():
        ode = json.loads(ODE_PATH.read_text())
        ode_delta = ode["ode_delta_o3plus"]
        bool_sign = "creation" if bool_delta > 0.005 else ("destruction" if bool_delta < -0.005 else "null")
        ode_sign = "creation" if ode_delta > 0.005 else ("destruction" if ode_delta < -0.005 else "null")
        sign_preserved = (bool_sign == ode_sign)
        print(f"   ODE delta: {ode_delta:+.4f} ({ode_delta*100:+.1f} pp)")
        print(f"   Boolean sign: {bool_sign}, ODE sign: {ode_sign}")
        print(f"   Sign preserved: {sign_preserved}")

        ode["boolean_delta_o3plus"] = bool_delta
        ode["sign_preserved"] = sign_preserved
        ODE_PATH.write_text(json.dumps(ode, indent=2))
    else:
        ode_delta = None
        sign_preserved = None
        print("   ODE result not found!")

    # 3. Recompute canalization correlations
    print("\n3. Recomputing canalization correlations with 28 networks...")
    corr, n_models = recompute_canalization_correlations(bool_delta)
    print(f"   {n_models} networks with delta values")
    smallest_adj_p = min(c["p_adjusted"] for c in corr.values())
    for key, c in sorted(corr.items(), key=lambda x: abs(x[1]["rho"]), reverse=True):
        sig = "*" if c["significant"] else " "
        print(f"   {c['name']:35s} rho={c['rho']:+.3f}  p_adj={c['p_adjusted']:.4f} {sig}")

    # Save updated correlations
    corr_out = {
        "n_models": n_models,
        "n_tests": len(corr),
        "correction": "bonferroni",
        "alpha": 0.05,
        "correlations": [{**v, "feature": k} for k, v in corr.items()],
        "models": sorted([name for name, r in json.loads(CANAL_FEATURES.read_text()).items()
                         if r.get("real_delta_o3plus") is not None]),
    }
    corr_path = Path("results/canalization/canalization_correlations_28.json")
    corr_path.write_text(json.dumps(corr_out, indent=2))
    print(f"   Saved: {corr_path}")

    # 4. Generate paper text updates
    print("\n" + "=" * 60)
    print("PAPER TEXT UPDATES FOR composition_gap_v10.tex -> v11")
    print("=" * 60)

    # Determine classification
    if abs(bool_delta * 100) < 0.5:
        classification = "null"
    elif bool_delta > 0:
        classification = "creation"
    else:
        classification = "destruction"

    # Count classification totals (from existing 27 + grieco)
    all_models = merged["all_models"]
    n_creation = sum(1 for m in all_models if m.get("delta_o3plus", 0) > 0.005)
    n_destruction = sum(1 for m in all_models if m.get("delta_o3plus", 0) < -0.005)
    n_null = len(all_models) - n_creation - n_destruction
    n_nonnull = n_creation + n_destruction

    # Compute summary stats
    deltas_pp = [m.get("delta_o3plus", 0) * 100 for m in all_models]
    median_delta = np.median(deltas_pp)
    mean_delta = np.mean(deltas_pp)

    # Format p-value for table
    def fmt_p(p):
        if p < 1e-12:
            return "$< 10^{-12}$"
        elif p < 1e-7:
            return "$< 10^{-7}$"
        elif p < 1e-4:
            return "$< 10^{-4}$"
        elif p < 1e-3:
            return "$< 10^{-3}$"
        else:
            return f"{p:.3f}"

    print(f"\n--- Grieco row for main table (sorted by cycling fraction: {cycling:.1%}) ---")
    null_marker = "$^\\dagger$" if classification == "null" else ""
    delta_str = f"${'+'if bool_delta>0 else '-'}${abs(bool_delta*100):.1f}"
    rho_str = f"{rho:.2f}" if rho >= 0 else f"$-${abs(rho):.2f}"
    print(f"Grieco bladder{null_marker} & {N_PLAYERS} & {rho_str} & {fmt_p(p_val)} & {delta_str} & {cycling:.1%} \\\\")

    print(f"\n--- Grieco row for ODE table ---")
    if ode_delta is not None:
        ode_str = f"${'+'if ode_delta>0 else '-'}${abs(ode_delta*100):.1f}"
        print(f"Grieco bladder & {N_PLAYERS} & {delta_str} & {ode_str} & {'preserved' if sign_preserved else 'CHANGED'} \\\\")

    print(f"\n--- Updated counts ---")
    print(f"Total networks: 28 (was 27)")
    print(f"Gene range: 7--18 (was 7--17)")
    print(f"Classification: {n_creation} creation, {n_destruction} destruction, {n_null} null (was 18/6/3)")
    print(f"Non-null: {n_nonnull} (was 24)")
    print(f"Median delta: {median_delta:+.1f} pp (was +8.5)")
    print(f"Mean delta: {mean_delta:+.0f} pp (was +20)")

    print(f"\n--- Canalization text (28 networks) ---")
    key_feats = ["mean_canalizing_depth", "mean_nested_depth", "frac_fully_canalizing", "derrida_parameter", "mean_bias_deviation"]
    for k in key_feats:
        if k in corr:
            c = corr[k]
            print(f"  {c['name']}: rho={c['rho']:+.2f}")
    print(f"  Smallest adjusted p: {smallest_adj_p:.3f}")

    print(f"\n--- Null model text ---")
    print(f"  grieco is n=18, too large for null model sweeps (limit n<=15)")
    print(f"  Update '26 testable networks' to keep at 26 (grieco excluded)")

    print(f"\n--- ODE table ---")
    n_ode = 26  # existing
    if ode_delta is not None:
        n_ode += 1  # grieco
    print(f"  ODE networks: {n_ode} of 28 (arabidopsis still pending)")


if __name__ == "__main__":
    main()
