"""Run blind composition experiment on remaining models with multiprocessing.

Uses 6 worker processes to parallelize the coalition sweep.
Saves checkpoints every 2000 coalitions so it can resume after interruption.
"""

import json
import sys
from datetime import datetime, timezone
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from composition_scorer import score_composition
from grn_coalition_sweep import (
    BUILTIN_MODELS,
    compile_network,
    extract_interaction_graph,
    extract_rule_fourier,
    identify_input_nodes,
    simulate_sync_output,
)


def timestamp():
    return datetime.now(timezone.utc).isoformat()


N_INIT = 512
SEED = 42
MAX_STEPS = 200
N_WORKERS = 6
CHUNK_SIZE = 2000
RESULTS_DIR = Path(__file__).parent.parent / "results" / "grn_v2"


def process_coalition_batch(args):
    """Process a batch of coalitions. Runs in a worker process."""
    (coalition_start, coalition_end, compiled_data, init_states,
     player_indices, n_total, clamp_value, output_indices, max_steps) = args

    compiled = compiled_data
    results = []
    for coalition in range(coalition_start, coalition_end):
        clamp_mask = np.zeros(n_total, dtype=bool)
        for bit_pos, node_idx in enumerate(player_indices):
            if not (coalition & (1 << bit_pos)):
                clamp_mask[node_idx] = True

        output, info = simulate_sync_output(
            init_states, compiled, clamp_mask, clamp_value,
            output_indices, max_steps=max_steps,
        )
        results.append((coalition, output, info["n_fixed_point"], info["n_cycling"]))

    return results


def sweep_coalitions_parallel(rules, output_nodes, n_init, max_steps, seed,
                              clamp_value, checkpoint_path=None):
    """Parallelized coalition sweep with checkpointing."""
    compiled, node_names = compile_network(rules)
    n_total = len(node_names)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    player_indices = list(range(n_total))
    player_names = [node_names[i] for i in player_indices]
    n_players = len(player_indices)
    N = 2 ** n_players

    rng = np.random.default_rng(seed)
    init_states = rng.integers(0, 2, size=(n_init, n_total), dtype=np.int8)

    values = np.zeros((N, n_init), dtype=np.float64)
    total_cycling = 0
    total_fixed = 0
    start_coalition = 0

    # Resume from checkpoint if available
    if checkpoint_path and checkpoint_path.exists():
        ckpt = np.load(checkpoint_path)
        start_coalition = int(ckpt["completed_coalitions"])
        values[:start_coalition] = ckpt["values"][:start_coalition]
        total_cycling = int(ckpt["total_cycling"])
        total_fixed = int(ckpt["total_fixed"])
        print(f"  Resuming from checkpoint: {start_coalition}/{N} coalitions done")

    remaining = N - start_coalition
    if remaining == 0:
        print(f"  Already complete!")
    else:
        # Process in chunks for checkpointing, within each chunk use multiprocessing
        chunk_start = start_coalition
        while chunk_start < N:
            chunk_end = min(chunk_start + CHUNK_SIZE, N)
            batch_size = (chunk_end - chunk_start + N_WORKERS - 1) // N_WORKERS

            batches = []
            for w in range(N_WORKERS):
                b_start = chunk_start + w * batch_size
                b_end = min(b_start + batch_size, chunk_end)
                if b_start >= chunk_end:
                    break
                batches.append((
                    b_start, b_end, compiled, init_states,
                    player_indices, n_total, clamp_value, output_indices, max_steps
                ))

            with Pool(N_WORKERS) as pool:
                batch_results = pool.map(process_coalition_batch, batches)

            for batch in batch_results:
                for coalition, output, n_fp, n_cyc in batch:
                    values[coalition] = output
                    total_fixed += n_fp
                    total_cycling += n_cyc

            chunk_start = chunk_end
            pct = chunk_start / N * 100
            print(f"  [{timestamp()}] {chunk_start}/{N} ({pct:.1f}%) coalitions complete")
            sys.stdout.flush()

            # Save checkpoint
            if checkpoint_path and chunk_start < N:
                np.savez_compressed(
                    checkpoint_path,
                    values=values,
                    completed_coalitions=np.int64(chunk_start),
                    total_cycling=np.int64(total_cycling),
                    total_fixed=np.int64(total_fixed),
                )

    cycling_frac = total_cycling / (N * n_init)
    convergence = {
        "update_scheme": "sync",
        "clamp_value": clamp_value,
        "total_simulations": N * n_init,
        "total_fixed_point": total_fixed,
        "total_cycling": total_cycling,
        "cycling_fraction": cycling_frac,
    }

    return {
        "values": values,
        "node_names": node_names,
        "player_names": player_names,
        "output_nodes": output_nodes,
        "n_players": n_players,
        "n_total_nodes": n_total,
        "n_init_states": n_init,
        "input_nodes": [],
        "input_config": None,
        "convergence": convergence,
    }


def run_model(model_name):
    model_info = BUILTIN_MODELS[model_name]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]
    n = len(rules)

    print(f"\n{'='*60}")
    print(f"[{timestamp()}] {model_name} -- {model_info['description']}")
    print(f"  n={n}, coalitions={2**n}, output={output_nodes}, n_init={N_INIT}")
    print(f"  Workers={N_WORKERS}, chunk_size={CHUNK_SIZE}")
    print(f"{'='*60}")

    rf = extract_rule_fourier(rules)
    print(f"  Local Fourier: {rf['n_pairwise']} pairwise, {rf['n_triples']} triples")
    local_spec = rf["local_energy_spectrum"]
    local_o3p = sum(local_spec[3:]) if len(local_spec) > 3 else 0.0
    print(f"  Local order-3+ energy: {local_o3p:.4f}")

    wiring = extract_interaction_graph(rules)
    n_edges = sum(len(edges) for edges in wiring.values())

    checkpoint_path = RESULTS_DIR / f"{model_name}_checkpoint.npz"

    print(f"\n  [Running parallel coalition sweep...]")
    result = sweep_coalitions_parallel(
        rules, output_nodes,
        n_init=N_INIT, max_steps=MAX_STEPS, seed=SEED,
        clamp_value=0, checkpoint_path=checkpoint_path,
    )

    node_names = result["player_names"]
    n_players = result["n_players"]
    v_mean = result["values"].mean(axis=1)
    n_unique = len(np.unique(np.round(v_mean, 6)))
    total_var = float(np.var(v_mean))
    conv = result["convergence"]

    print(f"\n  Value function: unique={n_unique}, var={total_var:.4e}")
    print(f"  v(empty)={v_mean[0]:.4f}, v(full)={v_mean[-1]:.4f}")
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

    tri = scores["triples"]
    if tri["spearman_rho"] is not None:
        print(f"\n    Triple Spearman rho = {tri['spearman_rho']:.4f}")

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

    # Clean up checkpoint
    if checkpoint_path.exists():
        checkpoint_path.unlink()
        print(f"  Removed checkpoint: {checkpoint_path}")

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

    summary_path = RESULTS_DIR / "blind_experiment_summary.json"
    with open(summary_path) as f:
        summaries = json.load(f)
    existing_models = {s["model"] for s in summaries}

    remaining = [m for m in ["fanconi_anemia", "arabidopsis_cellcycle"]
                 if m not in existing_models]

    print(f"[{timestamp()}] Parallel blind experiment. Existing: {sorted(existing_models)}")
    print(f"Remaining: {remaining}")
    print(f"Workers: {N_WORKERS}, Chunk: {CHUNK_SIZE}")

    for model_name in remaining:
        summary = run_model(model_name)
        summaries.append(summary)

        with open(summary_path, "w") as f:
            json.dump(summaries, f, indent=2, default=str)
        print(f"  [Summary saved: {summary_path}]")

    # Final summary
    print(f"\n\n{'='*80}")
    print(f"FINAL SUMMARY ({len(summaries)} models)")
    print(f"{'='*80}")
    header = f"{'Model':<25} {'rho':>6} {'95% CI':>18} {'p':>10} {'g_o3+':>8} {'l_o3+':>8} {'gap':>8} {'Dir':>10}"
    print(header)
    print("-" * len(header))
    for s in summaries:
        ci = f"[{s['spearman_ci'][0]:.3f}, {s['spearman_ci'][1]:.3f}]"
        print(f"{s['model']:<25} {s['spearman_rho']:>6.3f} {ci:>18} {s['spearman_p']:>10.2e} {s['global_o3plus']:>8.4f} {s['local_o3plus']:>8.4f} {s['delta_o3plus']:>+8.4f} {s['creation_or_destruction']:>10}")

    with open(summary_path, "w") as f:
        json.dump(summaries, f, indent=2, default=str)
    print(f"\n[{timestamp()}] All done. Summary: {summary_path}")


if __name__ == "__main__":
    main()
