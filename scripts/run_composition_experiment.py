"""Run the composition experiment: n_init=512, all 3 models, score each.

Saves coalition tables and composition scores to results/grn_v2/.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from composition_scorer import score_composition
from grn_coalition_sweep import (
    BUILTIN_MODELS,
    compile_network,
    extract_interaction_graph,
    extract_rule_fourier,
    sweep_coalitions,
)


def timestamp():
    return datetime.now(timezone.utc).isoformat()


MODELS = ["faure_cellcycle", "tournier_apoptosis", "davidich_yeast"]
N_INIT = 512
SEED = 42
MAX_STEPS = 200
RESULTS_DIR = Path(__file__).parent.parent / "results" / "grn_v2"


def run_model(model_name):
    model_info = BUILTIN_MODELS[model_name]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]

    print(f"\n{'='*60}")
    print(f"[{timestamp()}] {model_name} — {model_info['description']}")
    print(f"  Output: {output_nodes}, n_init={N_INIT}, seed={SEED}")
    print(f"{'='*60}")

    rf = extract_rule_fourier(rules)
    print(f"  Local Fourier: {rf['n_pairwise']} pairwise, {rf['n_triples']} triples")

    wiring = extract_interaction_graph(rules)
    n_edges = sum(len(edges) for edges in wiring.values())

    result = sweep_coalitions(
        rules, output_nodes,
        n_init=N_INIT, max_steps=MAX_STEPS, seed=SEED,
        update_scheme="sync", clamp_value=0,
    )

    node_names = result["node_names"]
    n = result["n_players"]
    v_mean = result["values"].mean(axis=1)
    n_unique = len(np.unique(np.round(v_mean, 6)))
    total_var = float(np.var(v_mean))
    conv = result["convergence"]

    print(f"\n  Value function: unique={n_unique}, var={total_var:.4e}")
    if conv["update_scheme"] == "sync":
        print(f"  Fixed: {conv['total_fixed_point']}, Cycling: {conv['total_cycling']} ({conv['cycling_fraction']:.1%})")

    npz_path = RESULTS_DIR / f"{model_name}_coalition.npz"
    np.savez_compressed(
        npz_path,
        target_logits=result["values"],
        foil_logits=np.zeros_like(result["values"]),
        coalition_indices=np.arange(2**n, dtype=np.int64),
        circuit_heads=np.array(node_names, dtype=object),
        n_players=np.int64(n),
        n_prompts=np.int64(N_INIT),
        model_name=model_name,
        output_nodes=np.array(output_nodes, dtype=object),
        update_scheme="sync",
        clamp_value=np.int64(0),
        n_unique_values=np.int64(n_unique),
        total_variance=np.float64(total_var),
    )
    print(f"  Saved: {npz_path}")

    wiring_path = RESULTS_DIR / f"{model_name}_wiring.json"
    wiring_output = {
        "model": model_name,
        "citation": model_info["citation"],
        "description": model_info["description"],
        "n_nodes": n,
        "node_names": node_names,
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
    print(f"  Saved: {wiring_path}")

    print(f"\n  [Scoring composition...]")
    scores = score_composition(v_mean, n, node_names, rf)

    pw = scores["pairwise"]
    spec = scores["energy_spectrum"]
    g3p = sum(spec["global"][3:])
    l3p = sum(spec["local_rules"][3:])

    print(f"  Spearman rho={pw['spearman_rho']:.4f} (p={pw['spearman_pvalue']:.2e})")
    ci = pw["spearman_ci_95"]
    print(f"    95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]")
    print(f"  Created: {pw['created_count']}/{pw['top_k']}  Destroyed: {pw['destroyed_count']}/{pw['top_k']}")
    print(f"  Global spectrum: o1={spec['global'][1]:.1%} o2={spec['global'][2]:.1%} o3+={g3p:.1%}")
    print(f"  Local spectrum:  o1={spec['local_rules'][1]:.1%} o2={spec['local_rules'][2]:.1%} o3+={l3p:.1%}")
    print(f"  Composition delta o3+: {g3p - l3p:+.1%} ({'creates' if g3p > l3p else 'destroys'})")

    scores_path = RESULTS_DIR / f"{model_name}_composition.json"
    scores_output = {
        "model": model_name,
        "n_players": n,
        "n_init": N_INIT,
        "node_names": node_names,
        "timestamp": timestamp(),
        **scores,
    }
    with open(scores_path, "w") as f:
        json.dump(scores_output, f, indent=2)
    print(f"  Saved: {scores_path}")

    return {
        "model": model_name,
        "n": n,
        "spearman_rho": pw["spearman_rho"],
        "spearman_ci": (pw["spearman_ci_95"][0], pw["spearman_ci_95"][1]),
        "spearman_p": pw["spearman_pvalue"],
        "global_o3plus": g3p,
        "local_o3plus": l3p,
        "delta_o3plus": g3p - l3p,
        "creation_rate": pw["creation_rate"],
        "destruction_rate": pw["destruction_rate"],
        "cycling_fraction": conv.get("cycling_fraction", 0),
    }


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[{timestamp()}] Composition experiment: {len(MODELS)} models, n_init={N_INIT}")

    summaries = []
    for model_name in MODELS:
        summary = run_model(model_name)
        summaries.append(summary)

    print(f"\n\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"{'Model':<25} {'rho':>6} {'95% CI':>16} {'p':>10} {'o3+ gap':>10} {'Dir':>10}")
    print("-" * 80)
    for s in summaries:
        ci = f"[{s['spearman_ci'][0]:.3f}, {s['spearman_ci'][1]:.3f}]"
        gap = f"{s['delta_o3plus']:+.1%}"
        direction = "creates" if s["delta_o3plus"] > 0 else "destroys"
        print(f"{s['model']:<25} {s['spearman_rho']:>6.3f} {ci:>16} {s['spearman_p']:>10.2e} {gap:>10} {direction:>10}")

    summary_path = RESULTS_DIR / "composition_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summaries, f, indent=2)
    print(f"\n[{timestamp()}] All done. Summary: {summary_path}")


if __name__ == "__main__":
    main()
