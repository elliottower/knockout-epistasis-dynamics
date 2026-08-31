"""Run blind composition experiment on the 2 remaining models:
fanconi_anemia (32768 coalitions) and arabidopsis_cellcycle (16384 coalitions).

Appends to existing blind_experiment_summary.json.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from composition_scorer import score_composition
from grn_coalition_sweep import (
    BUILTIN_MODELS,
    extract_interaction_graph,
    extract_rule_fourier,
    sweep_coalitions,
)


def timestamp():
    return datetime.now(timezone.utc).isoformat()


REMAINING_MODELS = [
    "fanconi_anemia",
    "arabidopsis_cellcycle",
]

N_INIT = 512
SEED = 42
MAX_STEPS = 200
RESULTS_DIR = Path(__file__).parent.parent / "results" / "grn_v2"


def run_model(model_name):
    model_info = BUILTIN_MODELS[model_name]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]
    n = len(rules)

    print(f"\n{'='*60}")
    print(f"[{timestamp()}] {model_name} -- {model_info['description']}")
    print(f"  n={n}, coalitions={2**n}, output={output_nodes}, n_init={N_INIT}")
    print(f"{'='*60}")

    # Extract local Fourier structure
    rf = extract_rule_fourier(rules)
    print(f"  Local Fourier: {rf['n_pairwise']} pairwise, {rf['n_triples']} triples")
    local_spec = rf["local_energy_spectrum"]
    local_o3p = sum(local_spec[3:]) if len(local_spec) > 3 else 0.0
    print(f"  Local order-3+ energy: {local_o3p:.4f}")

    # Extract wiring
    wiring = extract_interaction_graph(rules)
    n_edges = sum(len(edges) for edges in wiring.values())

    # Run coalition sweep
    print(f"\n  [Running coalition sweep...]")
    result = sweep_coalitions(
        rules, output_nodes,
        n_init=N_INIT, max_steps=MAX_STEPS, seed=SEED,
        update_scheme="sync", clamp_value=0,
    )

    node_names = result["player_names"]
    n_players = result["n_players"]
    v_mean = result["values"].mean(axis=1)
    n_unique = len(np.unique(np.round(v_mean, 6)))
    total_var = float(np.var(v_mean))
    conv = result["convergence"]

    print(f"\n  Value function: unique={n_unique}, var={total_var:.4e}")
    print(f"  v(empty)={v_mean[0]:.4f}, v(full)={v_mean[-1]:.4f}")
    if conv["update_scheme"] == "sync":
        print(f"  Fixed: {conv['total_fixed_point']}, Cycling: {conv['total_cycling']} ({conv['cycling_fraction']:.1%})")

    # Save coalition table
    npz_path = RESULTS_DIR / f"{model_name}_coalition_blind.npz"
    np.savez_compressed(
        npz_path,
        target_logits=result["values"],
        foil_logits=np.zeros_like(result["values"]),
        coalition_indices=np.arange(2**n_players, dtype=np.int64),
        circuit_heads=np.array(node_names, dtype=object),
        n_players=np.int64(n_players),
        n_prompts=np.int64(N_INIT),
        model_name=model_name,
        output_nodes=np.array(output_nodes, dtype=object),
        update_scheme="sync",
        clamp_value=np.int64(0),
        n_unique_values=np.int64(n_unique),
        total_variance=np.float64(total_var),
    )
    print(f"  Saved coalition table: {npz_path}")

    # Save wiring
    wiring_path = RESULTS_DIR / f"{model_name}_wiring_blind.json"
    wiring_output = {
        "model": model_name,
        "citation": model_info["citation"],
        "description": model_info["description"],
        "n_nodes": n,
        "node_names": list(rules.keys()),
        "player_names": node_names,
        "output_nodes": output_nodes,
        "interaction_graph": wiring,
        "n_edges": n_edges,
        "rule_fourier": rf,
        "convergence": conv,
        "n_unique_values": n_unique,
        "total_variance": total_var,
        "timestamp": timestamp(),
    }
    with open(wiring_path, "w") as f:
        json.dump(wiring_output, f, indent=2)
    print(f"  Saved wiring: {wiring_path}")

    # Score composition
    print(f"\n  [Scoring composition...]")
    scores = score_composition(v_mean, n_players, node_names, rf)

    pw = scores["pairwise"]
    spec = scores["energy_spectrum"]
    g3p = sum(spec["global"][3:])
    l3p = sum(spec["local_rules"][3:])

    print(f"\n  RESULTS:")
    print(f"    Spearman rho = {pw['spearman_rho']:.4f}")
    ci = pw["spearman_ci_95"]
    print(f"    95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]")
    print(f"    p-value: {pw['spearman_pvalue']:.2e}")
    print(f"    Created: {pw['created_count']}/{pw['top_k']}  Destroyed: {pw['destroyed_count']}/{pw['top_k']}")

    print(f"\n    Energy spectrum:")
    for k in range(min(5, n_players + 1)):
        g = spec["global"][k]
        l = spec["local_rules"][k]
        delta = g - l
        arrow = "+" if delta > 0 else ""
        print(f"      order-{k}: global {g:.4f}  local {l:.4f}  ({arrow}{delta:.4f})")
    print(f"      order-3+: global {g3p:.4f}  local {l3p:.4f}  (delta={g3p-l3p:+.4f})")
    print(f"      Direction: {'CREATION' if g3p > l3p else 'DESTRUCTION'}")

    # Triple composition
    tri = scores["triples"]
    if tri["spearman_rho"] is not None:
        print(f"\n    Triple Spearman rho = {tri['spearman_rho']:.4f}")

    # Save scores
    scores_path = RESULTS_DIR / f"{model_name}_composition_blind.json"
    scores_output = {
        "model": model_name,
        "n_players": n_players,
        "n_init": N_INIT,
        "node_names": node_names,
        "timestamp": timestamp(),
        **scores,
    }
    with open(scores_path, "w") as f:
        json.dump(scores_output, f, indent=2)
    print(f"  Saved composition scores: {scores_path}")

    return {
        "model": model_name,
        "n": n_players,
        "spearman_rho": pw["spearman_rho"],
        "spearman_ci": pw["spearman_ci_95"],
        "spearman_p": pw["spearman_pvalue"],
        "global_spectrum": spec["global"],
        "local_spectrum": spec["local_rules"],
        "global_o3plus": g3p,
        "local_o3plus": l3p,
        "delta_o3plus": g3p - l3p,
        "creation_or_destruction": "creation" if g3p > l3p else "destruction",
        "creation_rate": pw["creation_rate"],
        "destruction_rate": pw["destruction_rate"],
        "cycling_fraction": conv.get("cycling_fraction", 0),
        "triple_rho": tri["spearman_rho"],
    }


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load existing summaries
    summary_path = RESULTS_DIR / "blind_experiment_summary.json"
    with open(summary_path) as f:
        summaries = json.load(f)
    existing_models = {s["model"] for s in summaries}
    print(f"[{timestamp()}] Resuming blind experiment. Existing: {sorted(existing_models)}")
    print(f"Remaining: {REMAINING_MODELS}")

    for model_name in REMAINING_MODELS:
        if model_name in existing_models:
            print(f"  Skipping {model_name} (already done)")
            continue

        summary = run_model(model_name)
        summaries.append(summary)

        # Save intermediate summary
        with open(summary_path, "w") as f:
            json.dump(summaries, f, indent=2, default=str)
        print(f"  [Intermediate summary saved: {summary_path}]")

    # Final summary table
    print(f"\n\n{'='*80}")
    print(f"FINAL SUMMARY ({len(summaries)} models)")
    print(f"{'='*80}")
    header = f"{'Model':<25} {'rho':>6} {'95% CI':>18} {'p':>10} {'g_o3+':>8} {'l_o3+':>8} {'gap':>8} {'Dir':>10}"
    print(header)
    print("-" * len(header))
    for s in summaries:
        ci = f"[{s['spearman_ci'][0]:.3f}, {s['spearman_ci'][1]:.3f}]"
        print(f"{s['model']:<25} {s['spearman_rho']:>6.3f} {ci:>18} {s['spearman_p']:>10.2e} {s['global_o3plus']:>8.4f} {s['local_o3plus']:>8.4f} {s['delta_o3plus']:>+8.4f} {s['creation_or_destruction']:>10}")

    # Save final summary
    with open(summary_path, "w") as f:
        json.dump(summaries, f, indent=2, default=str)
    print(f"\n[{timestamp()}] All done. Summary: {summary_path}")


if __name__ == "__main__":
    main()
