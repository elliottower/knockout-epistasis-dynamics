"""Robustness check: exclude constitutive input nodes from the player set.

Input nodes (f(x)=x) persist indefinitely once set. They act as external
parameters, not dynamically regulated genes. This script re-runs coalition
sweeps for the 4 models with input nodes, excluding those nodes from the
player set, and compares delta and rho to the original analysis.

Only 4 of 27 models have input nodes:
  - davidich_yeast (Start)
  - drosophila_cellcycle (Ago, CycD, Notch)
  - faure_cellcycle (CycD)
  - tournier_apoptosis (TNF)
"""
import json
import numpy as np
from pathlib import Path
from datetime import datetime, timezone

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from grn_coalition_sweep import (
    BUILTIN_MODELS, identify_input_nodes, sweep_coalitions,
    extract_rule_fourier,
)
from composition_scorer import score_composition


RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "grn_v2"

MODELS_WITH_INPUTS = [
    "davidich_yeast",
    "drosophila_cellcycle",
    "faure_cellcycle",
    "tournier_apoptosis",
]


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def run_model(model_name, input_config=0):
    info = BUILTIN_MODELS[model_name]
    rules = info["rules"]
    output_nodes = info["output_nodes"]
    input_nodes = identify_input_nodes(rules)

    print(f"  Input nodes: {input_nodes} (fixed to {input_config})")

    rf = extract_rule_fourier(rules)

    print(f"  [{timestamp()}] Starting coalition sweep (exclude_inputs=True, input_config={input_config})...")
    result = sweep_coalitions(
        rules, output_nodes, n_init=512, max_steps=200, seed=42,
        exclude_inputs=True, input_config=input_config,
    )

    player_names = result["player_names"]
    n = result["n_players"]
    v_mean = result["values"].mean(axis=1)
    conv = result["convergence"]

    print(f"  Players: {n} (was {len(rules)} with inputs)")
    print(f"  Cycling: {conv.get('cycling_fraction', 0):.1%}")

    npz_path = RESULTS_DIR / f"{model_name}_coalition_excl_inputs_cfg{input_config}.npz"
    np.savez_compressed(
        npz_path,
        target_logits=result["values"],
        foil_logits=np.zeros_like(result["values"]),
        coalition_indices=np.arange(2**n, dtype=np.int64),
        circuit_heads=np.array(player_names, dtype=object),
        n_players=np.int64(n),
        n_prompts=np.int64(512),
        model_name=model_name,
        output_nodes=np.array(output_nodes, dtype=object),
    )
    print(f"  Saved: {npz_path}")

    print(f"  [{timestamp()}] Scoring composition...")
    scores = score_composition(v_mean, n, player_names, rf)
    pw = scores["pairwise"]
    spec = scores["energy_spectrum"]
    g3p = sum(spec["global"][3:])
    l3p = sum(spec["local_rules"][3:])

    print(f"  Spearman rho={pw['spearman_rho']:.4f} (p={pw['spearman_pvalue']:.2e})")
    print(f"  Global o3+={g3p:.1%}, Local o3+={l3p:.1%}, Delta={g3p-l3p:+.1%}")

    comp_path = RESULTS_DIR / f"{model_name}_composition_excl_inputs_cfg{input_config}.json"
    comp_output = {
        "model": model_name,
        "n_players": n,
        "n_original": len(rules),
        "excluded_inputs": input_nodes,
        "input_config": input_config,
        "player_names": player_names,
        "timestamp": timestamp(),
        **scores,
    }
    with open(comp_path, "w") as f:
        json.dump(comp_output, f, indent=2)
    print(f"  Saved: {comp_path}")

    return {
        "model": model_name,
        "n_original": len(rules),
        "n_players": n,
        "excluded": input_nodes,
        "input_config": input_config,
        "rho": pw["spearman_rho"],
        "p": pw["spearman_pvalue"],
        "global_o3plus": g3p,
        "local_o3plus": l3p,
        "delta_o3plus": g3p - l3p,
        "label": "creation" if g3p > l3p + 0.005 else ("destruction" if l3p > g3p + 0.005 else "null"),
        "cycling_fraction": conv.get("cycling_fraction", 0),
    }


if __name__ == "__main__":
    clamp0_data = {}
    for summary_file in ["blind_experiment_summary.json", "blind_batch2_summary.json"]:
        path = RESULTS_DIR / summary_file
        if path.exists():
            with open(path) as f:
                data = json.load(f)
                items = data if isinstance(data, list) else data.get("results", [])
                for m in items:
                    clamp0_data[m["model"]] = m

    all_results = []
    for model_name in MODELS_WITH_INPUTS:
        for input_cfg in [0, 1]:
            n = len(BUILTIN_MODELS[model_name]["rules"])
            print(f"\n{'='*60}")
            print(f"[{timestamp()}] === {model_name} (n={n}, input_config={input_cfg}) ===")
            print(f"{'='*60}")
            result = run_model(model_name, input_config=input_cfg)
            all_results.append(result)

    print(f"\n{'='*80}")
    print(f"COMPARISON: all-players vs exclude-inputs")
    print(f"{'='*80}")
    print(f"{'Model':25s} {'cfg':>3s} {'n_orig':>6s} {'n_excl':>6s} {'rho_orig':>9s} {'rho_excl':>9s} {'d_orig':>9s} {'d_excl':>9s} {'same?':>6s}")
    print("-" * 100)

    NULL_THRESHOLD = 0.005
    for r in all_results:
        orig = clamp0_data.get(r["model"], {})
        rho_orig = orig.get("spearman_rho", float("nan"))
        d_orig = orig.get("delta_o3plus", float("nan"))
        d_excl = r["delta_o3plus"]

        null_o = abs(d_orig) < NULL_THRESHOLD if not np.isnan(d_orig) else True
        null_e = abs(d_excl) < NULL_THRESHOLD

        if null_o and null_e:
            same = "null"
        elif null_o or null_e:
            same = "~"
        else:
            same = "YES" if (d_orig > 0) == (d_excl > 0) else "NO"

        print(f"{r['model']:25s} {r['input_config']:3d} {r['n_original']:6d} {r['n_players']:6d} "
              f"{rho_orig:+8.3f}  {r['rho']:+8.3f}  {d_orig*100:+8.2f}pp {d_excl*100:+8.2f}pp {same:>6s}")

    summary_path = RESULTS_DIR / "robustness_exclude_inputs_summary.json"
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n[{timestamp()}] Done. Summary saved to {summary_path}")
