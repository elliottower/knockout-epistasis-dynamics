"""Finalize grieco_bladder on Modal: merge Boolean shards, compute WHT + pairwise rho,
recompute canalization correlations for 28 networks. Downloads only small JSON results.

Reads shards from Modal volume, writes results back to volume + prints JSON to stdout.

Usage:
    modal run scripts/modal_finalize_grieco.py
"""
from __future__ import annotations

import os

import modal

REPO = os.path.join(os.path.expanduser("~"), "Documents/GitHub/epistasis-bench")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "numpy==2.2.6",
        "scipy==1.15.3",
        "scikit-learn==1.6.1",
        "tqdm==4.67.1",
        "matplotlib==3.10.3",
    )
    .env({"PYTHONPATH": "/root/repo:/root/repo/scripts"})
    .add_local_file(os.path.join(REPO, "data_utils.py"), "/root/repo/data_utils.py", copy=True)
    .add_local_file(os.path.join(REPO, "grn_coalition_sweep.py"), "/root/repo/grn_coalition_sweep.py", copy=True)
    .add_local_file(os.path.join(REPO, "composition_scorer.py"), "/root/repo/composition_scorer.py", copy=True)
    .add_local_file(
        os.path.join(REPO, "scripts/run_batch2_blind_sweep.py"),
        "/root/repo/scripts/run_batch2_blind_sweep.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "scripts/structural_analysis.py"),
        "/root/repo/scripts/structural_analysis.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "results/canalization/canalization_features.json"),
        "/root/repo/canalization_features.json",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "results/ode_full/grieco_bladder_ode.json"),
        "/root/repo/grieco_bladder_ode.json",
        copy=True,
    )
    .run_commands("touch /root/repo/scripts/__init__.py")
)

app = modal.App("epistasis-bench-finalize-grieco", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

MODEL_NAME = "grieco_bladder"
N_SHARDS = 8
N_PLAYERS = 18
N_INIT = 512


@app.function(cpu=4, memory=16384, timeout=3600, volumes={"/results": results_vol})
def finalize():
    """Merge shards, compute WHT, pairwise rho, return JSON results."""
    import json

    import numpy as np
    from scipy import stats

    from composition_scorer import score_composition
    from grn_coalition_sweep import compile_network, extract_rule_fourier
    from scripts.run_batch2_blind_sweep import WEB_MODELS

    shard_dir = "/results/boolean_full/grieco_shards_v2"

    # 1. Merge Boolean shards
    print("=== Merging Boolean shards ===")
    arrays = []
    total_fixed = 0
    total_cycling = 0

    for i in range(N_SHARDS):
        shard_path = f"{shard_dir}/shard_{i:03d}.npy"
        meta_path = f"{shard_dir}/shard_{i:03d}_meta.json"
        fix_shard = f"{shard_dir}/shard_{i:03d}_fix.npy"
        fix_meta = f"{shard_dir}/shard_{i:03d}_fix_meta.json"

        if os.path.exists(fix_meta):
            with open(fix_meta) as f:
                fix_m = json.load(f)
            if fix_m.get("status") == "complete":
                shard_path = fix_shard
                meta_path = fix_meta
                print(f"  Shard {i}: using _fix version")

        if not os.path.exists(shard_path):
            raise FileNotFoundError(f"Missing shard {i}: {shard_path}")

        with open(meta_path) as f:
            meta = json.load(f)
        if meta.get("status") != "complete":
            raise RuntimeError(f"Shard {i} not complete: {meta.get('n_completed', '?')}/{meta.get('n_total', meta.get('n_coalitions', '?'))}")

        arr = np.load(shard_path)
        print(f"  Shard {i}: {arr.shape}")
        arrays.append(arr)
        total_fixed += meta.get("total_fixed", 0)
        total_cycling += meta.get("total_cycling", 0)

    values = np.concatenate(arrays, axis=0)
    expected = 2**N_PLAYERS
    assert values.shape[0] == expected, f"Got {values.shape[0]}, expected {expected}"
    print(f"  Merged: {values.shape} ({expected} coalitions x {values.shape[1]} inits)")

    # 2. Compute WHT + pairwise rho via score_composition
    print("\n=== Computing WHT + pairwise rho ===")
    v_mean = values.mean(axis=1)

    model_info = WEB_MODELS[MODEL_NAME]
    compiled, node_names = compile_network(model_info["rules"])
    rf = extract_rule_fourier(model_info["rules"])

    scores = score_composition(v_mean, N_PLAYERS, node_names, rf)

    pw = scores["pairwise"]
    spec = scores["energy_spectrum"]
    global_o3plus = sum(spec["global"][3:])
    local_o3plus = sum(spec["local_rules"][3:])
    delta = global_o3plus - local_o3plus

    total_samples = values.shape[0] * values.shape[1]
    cycling_frac = total_cycling / total_samples if total_samples > 0 else 0

    grieco_result = {
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
        "pairwise": pw,
        "triples": scores["triples"],
        "spectrum_gate": scores["spectrum_gate"],
    }

    print(f"  Boolean delta o3+: {delta:+.4f} ({delta*100:+.1f} pp)")
    print(f"  Cycling fraction:  {cycling_frac:.1%}")
    print(f"  Pairwise rho:      {pw['spearman_rho']:.2f} (p={pw['spearman_pvalue']:.4f})")

    # 3. ODE comparison
    print("\n=== ODE comparison ===")
    ode_path = "/root/repo/grieco_bladder_ode.json"
    if os.path.exists(ode_path):
        with open(ode_path) as f:
            ode = json.load(f)
        ode_delta = ode["ode_delta_o3plus"]
        bool_sign = "creation" if delta > 0.005 else ("destruction" if delta < -0.005 else "null")
        ode_sign = "creation" if ode_delta > 0.005 else ("destruction" if ode_delta < -0.005 else "null")
        sign_preserved = (bool_sign == ode_sign)
        print(f"  ODE delta: {ode_delta:+.4f} ({ode_delta*100:+.1f} pp)")
        print(f"  Boolean sign: {bool_sign}, ODE sign: {ode_sign}")
        print(f"  Sign preserved: {sign_preserved}")
        grieco_result["ode_delta"] = ode_delta
        grieco_result["ode_sign_preserved"] = sign_preserved
    else:
        ode_delta = None
        sign_preserved = None
        print("  ODE result not found!")

    # 4. Canalization correlations with 28 networks
    print("\n=== Canalization correlations (28 networks) ===")
    with open("/root/repo/canalization_features.json") as f:
        features = json.load(f)

    features[MODEL_NAME]["real_delta_o3plus"] = delta * 100

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

    canal_results = {}
    for key, name in feat_keys:
        vals = np.array([features[m][key] for m in models_with_delta])
        if np.std(vals) < 1e-10:
            continue
        rho, p = stats.spearmanr(vals, deltas)
        p_adj = min(p * n_tests, 1.0)
        canal_results[key] = {
            "name": name,
            "rho": float(rho),
            "p_raw": float(p),
            "p_adjusted": float(p_adj),
            "significant": bool(p_adj < 0.05),
        }

    for key, c in sorted(canal_results.items(), key=lambda x: abs(x[1]["rho"]), reverse=True):
        sig = "*" if c["significant"] else " "
        print(f"  {c['name']:35s} rho={c['rho']:+.3f}  p_adj={c['p_adjusted']:.4f} {sig}")

    canal_out = {
        "n_models": len(models_with_delta),
        "n_tests": len(canal_results),
        "correction": "bonferroni",
        "alpha": 0.05,
        "correlations": [{**v, "feature": k} for k, v in canal_results.items()],
        "models": sorted(models_with_delta),
    }

    # 5. Save results to volume
    print("\n=== Saving results to volume ===")
    out_dir = "/results/grn_v2"
    os.makedirs(out_dir, exist_ok=True)

    with open(f"{out_dir}/grieco_bladder_analysis.json", "w") as f:
        json.dump(grieco_result, f, indent=2)
    print(f"  Saved grieco_bladder_analysis.json")

    canal_dir = "/results/canalization"
    os.makedirs(canal_dir, exist_ok=True)
    with open(f"{canal_dir}/canalization_correlations_28.json", "w") as f:
        json.dump(canal_out, f, indent=2)
    print(f"  Saved canalization_correlations_28.json")

    results_vol.commit()

    # Return the results so they print to local stdout
    return {
        "grieco_result": grieco_result,
        "canalization_correlations": canal_out,
    }


@app.local_entrypoint()
def main():
    import json
    from pathlib import Path

    print("Running finalize on Modal (no local download needed)...")
    result = finalize.remote()

    # Save the small JSON results locally
    out_dir = Path("results/grn_v2")
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "grieco_bladder_analysis.json", "w") as f:
        json.dump(result["grieco_result"], f, indent=2)
    print(f"\nSaved locally: {out_dir / 'grieco_bladder_analysis.json'}")

    canal_dir = Path("results/canalization")
    canal_dir.mkdir(parents=True, exist_ok=True)
    with open(canal_dir / "canalization_correlations_28.json", "w") as f:
        json.dump(result["canalization_correlations"], f, indent=2)
    print(f"Saved locally: {canal_dir / 'canalization_correlations_28.json'}")

    # Print summary
    gr = result["grieco_result"]
    print(f"\n{'='*60}")
    print(f"GRIECO BOOLEAN RESULTS")
    print(f"{'='*60}")
    print(f"  Delta o3+:        {gr['delta_o3plus']:+.4f} ({gr['delta_o3plus']*100:+.1f} pp)")
    print(f"  Cycling fraction: {gr['cycling_fraction']:.1%}")
    print(f"  Pairwise rho:     {gr['spearman_rho']:.2f} (p={gr['spearman_p']:.4f})")
    if gr.get("ode_sign_preserved") is not None:
        print(f"  ODE sign match:   {gr['ode_sign_preserved']}")
