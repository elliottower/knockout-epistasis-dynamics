"""Run remaining 3 blind models with multiprocessing parallelism.

Drosophila (16k coalitions), fanconi (32k), arabidopsis (16k) at n_init=512.
Each coalition is independent, so we parallelize across CPU cores.
"""
import json
import multiprocessing as mp
import numpy as np
from datetime import datetime, timezone
from functools import partial
from tqdm import tqdm

import sys
sys.path.insert(0, "/Users/elliottower/Documents/GitHub/epistasis-bench")

from grn_coalition_sweep import (
    BUILTIN_MODELS, compile_network, identify_input_nodes,
    simulate_sync_output, extract_rule_fourier, extract_interaction_graph,
)
from composition_scorer import score_composition


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def _run_one_coalition(args):
    coalition, player_indices, init_states, compiled, clamp_value, output_indices, max_steps, n_total = args
    clamp_mask = np.zeros(n_total, dtype=bool)
    for bit_pos, node_idx in enumerate(player_indices):
        if not (coalition & (1 << bit_pos)):
            clamp_mask[node_idx] = True
    output, info = simulate_sync_output(
        init_states, compiled, clamp_mask, clamp_value,
        output_indices, max_steps=max_steps,
    )
    return coalition, output, info["n_fixed_point"], info["n_cycling"]


def sweep_coalitions_parallel(rules, output_nodes, n_init=512, max_steps=200,
                               seed=42, clamp_value=0, n_workers=10):
    compiled, node_names = compile_network(rules)
    n_total = len(node_names)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    player_indices = list(range(n_total))
    player_names = node_names
    n_players = len(player_indices)
    N = 2**n_players

    rng = np.random.default_rng(seed)
    init_states = rng.integers(0, 2, size=(n_init, n_total), dtype=np.int8)

    values = np.zeros((N, n_init), dtype=np.float64)
    total_cycling = 0
    total_fixed = 0

    task_args = [
        (c, player_indices, init_states, compiled, clamp_value, output_indices, max_steps, n_total)
        for c in range(N)
    ]

    with mp.Pool(n_workers) as pool:
        for coalition, output, n_fixed, n_cyc in tqdm(
            pool.imap_unordered(_run_one_coalition, task_args, chunksize=64),
            total=N,
            desc=f"Coalitions (parallel, {n_workers} workers)",
        ):
            values[coalition] = output
            total_fixed += n_fixed
            total_cycling += n_cyc

    convergence = {
        "update_scheme": "sync",
        "clamp_value": clamp_value,
        "total_simulations": N * n_init,
        "total_fixed_point": total_fixed,
        "total_cycling": total_cycling,
        "cycling_fraction": total_cycling / (N * n_init),
    }

    return {
        "values": values,
        "node_names": node_names,
        "player_names": player_names,
        "output_nodes": output_nodes,
        "n_players": n_players,
        "n_total_nodes": n_total,
        "n_init_states": n_init,
        "convergence": convergence,
    }


def run_model(model_name, n_init=512, max_steps=200, seed=42, n_workers=10):
    results_dir = "/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2"

    info = BUILTIN_MODELS[model_name]
    rules = info["rules"]
    output_nodes = info["output_nodes"]

    rf = extract_rule_fourier(rules)
    print(f"  Local Fourier: {rf['n_pairwise']} pairwise, {rf['n_triples']} triples")

    wiring = extract_interaction_graph(rules)

    print(f"  [{timestamp()}] Starting coalition sweep...")
    result = sweep_coalitions_parallel(
        rules, output_nodes, n_init=n_init, max_steps=max_steps,
        seed=seed, clamp_value=0, n_workers=n_workers,
    )

    node_names = result["node_names"]
    n = result["n_players"]
    v_mean = result["values"].mean(axis=1)
    n_unique = len(np.unique(np.round(v_mean, 6)))
    total_var = float(np.var(v_mean))
    conv = result["convergence"]

    print(f"  Value function: unique={n_unique}, var={total_var:.4e}")
    print(f"  Fixed: {conv['total_fixed_point']}, Cycling: {conv['total_cycling']} ({conv['cycling_fraction']:.1%})")

    npz_path = f"{results_dir}/{model_name}_coalition_blind.npz"
    np.savez_compressed(
        npz_path,
        target_logits=result["values"],
        foil_logits=np.zeros_like(result["values"]),
        coalition_indices=np.arange(2**n, dtype=np.int64),
        circuit_heads=np.array(node_names, dtype=object),
        n_players=np.int64(n),
        n_prompts=np.int64(n_init),
        model_name=model_name,
        output_nodes=np.array(output_nodes, dtype=object),
    )
    print(f"  Saved: {npz_path}")

    wiring_path = f"{results_dir}/{model_name}_wiring_blind.json"
    wiring_output = {
        "model": model_name,
        "n_nodes": n,
        "node_names": node_names,
        "output_nodes": output_nodes,
        "interaction_graph": wiring,
        "rule_fourier": rf,
        "convergence": conv,
        "n_unique_values": n_unique,
        "total_variance": total_var,
        "timestamp": timestamp(),
    }
    with open(wiring_path, "w") as f:
        json.dump(wiring_output, f, indent=2)
    print(f"  Saved: {wiring_path}")

    print(f"  [{timestamp()}] Scoring composition...")
    scores = score_composition(v_mean, n, node_names, rf)
    pw = scores["pairwise"]
    spec = scores["energy_spectrum"]
    g3p = sum(spec["global"][3:])
    l3p = sum(spec["local_rules"][3:])
    tri = scores.get("triples", {})

    print(f"  Spearman rho={pw['spearman_rho']:.4f} (p={pw['spearman_pvalue']:.2e})")
    print(f"  CI: [{pw['spearman_ci_95'][0]:.4f}, {pw['spearman_ci_95'][1]:.4f}]")
    print(f"  Global o3+={g3p:.1%}, Local o3+={l3p:.1%}, Delta={g3p-l3p:+.1%} ({'creates' if g3p > l3p else 'destroys'})")

    comp_path = f"{results_dir}/{model_name}_composition_blind.json"
    comp_output = {
        "model": model_name,
        "n_players": n,
        "n_init": n_init,
        "node_names": node_names,
        "timestamp": timestamp(),
        **scores,
    }
    with open(comp_path, "w") as f:
        json.dump(comp_output, f, indent=2)
    print(f"  Saved: {comp_path}")

    tri_rho = tri.get("spearman_rho", None)
    summary_line = (
        f"  SUMMARY: rho={pw['spearman_rho']:.3f} "
        f"[{pw['spearman_ci_95'][0]:.3f},{pw['spearman_ci_95'][1]:.3f}] "
        f"p={pw['spearman_pvalue']:.2e} o3+gap={g3p-l3p:+.1%} triple_rho={tri_rho}"
    )
    print(summary_line)
    return {
        "model": model_name,
        "n": n,
        "spearman_rho": pw["spearman_rho"],
        "spearman_ci": pw["spearman_ci_95"],
        "spearman_p": pw["spearman_pvalue"],
        "global_o3plus": g3p,
        "local_o3plus": l3p,
        "delta_o3plus": g3p - l3p,
        "creation_or_destruction": "creation" if g3p > l3p else "destruction",
        "cycling_fraction": conv["cycling_fraction"],
        "triple_rho": tri_rho,
        "global_spectrum": spec["global"],
        "local_spectrum": spec["local_rules"],
    }


if __name__ == "__main__":
    MODELS = ["fanconi_anemia", "arabidopsis_cellcycle"]
    N_WORKERS = 10

    all_results = []
    for model_name in MODELS:
        print(f"\n{'='*60}")
        print(f"[{timestamp()}] === {model_name} ===")
        print(f"{'='*60}")
        result = run_model(model_name, n_workers=N_WORKERS)
        all_results.append(result)

    summary_path = "/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2/blind_remaining_3_summary.json"
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n[{timestamp()}] All 3 models complete. Summary: {summary_path}")
