"""Robustness check: re-run all 27 coalition sweeps with clamp_value=1.

The original analysis uses clamp-to-0 (knockout). This script runs the
same sweeps with clamp-to-1 (constitutive activation) to test whether
the creation dominance result is baseline-dependent.

Boolean AND/OR asymmetry makes the value function baseline-dependent,
so this is a meaningful sensitivity analysis.
"""
import json
import multiprocessing as mp
import numpy as np
from datetime import datetime, timezone
from tqdm import tqdm

import sys
sys.path.insert(0, "/Users/elliottower/Documents/GitHub/epistasis-bench")

from grn_coalition_sweep import (
    BUILTIN_MODELS, compile_network, identify_input_nodes,
    simulate_sync_output, extract_rule_fourier, extract_interaction_graph,
)
from composition_scorer import score_composition
from scripts.run_batch2b_extra_models import EXTRA_MODELS
from scripts.run_batch2_blind_sweep import WEB_MODELS

ALL_KNOWN_MODELS = {**BUILTIN_MODELS, **EXTRA_MODELS, **WEB_MODELS}


RESULTS_DIR = "/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2"
CLAMP_VALUE = 1

ALL_MODELS_SORTED = [
    "lambda_phage",
    "arellano_rootstem",
    "asymmetric_cell_division",
    "cell_cycle_transcription",
    "faure_cellcycle",
    "davidich_yeast",
    "remy_p53_mdm2",
    "albert_segment_polarity",
    "blood_stem_cell",
    "calzone_cellfate_reduced",
    "li_budding_yeast",
    "myeloid_progenitors",
    "pair_rule_module",
    "tournier_apoptosis",
    "emt_switch",
    "morphogenetic_checkpoint",
    "zanudo_tlgl",
    "lac_operon",
    "mendoza_thelper",
    "saadatpour_guardcell",
    "drosophila_cellcycle",
    "arabidopsis_cellcycle",
    "fumia_cellcycle",
    "fanconi_anemia",
    "hematopoiesis_aging",
    "irons_cardiac",
    "calzone_cell_fate",
]


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
                               seed=42, clamp_value=1, n_workers=10):
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
            desc=f"Coalitions (clamp={clamp_value}, {n_workers} workers)",
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
    info = ALL_KNOWN_MODELS[model_name]
    rules = info["rules"]
    output_nodes = info["output_nodes"]

    rf = extract_rule_fourier(rules)
    print(f"  Local Fourier: {rf['n_pairwise']} pairwise, {rf['n_triples']} triples")

    print(f"  [{timestamp()}] Starting coalition sweep (clamp=1)...")
    result = sweep_coalitions_parallel(
        rules, output_nodes, n_init=n_init, max_steps=max_steps,
        seed=seed, clamp_value=CLAMP_VALUE, n_workers=n_workers,
    )

    node_names = result["node_names"]
    n = result["n_players"]
    v_mean = result["values"].mean(axis=1)
    n_unique = len(np.unique(np.round(v_mean, 6)))
    total_var = float(np.var(v_mean))
    conv = result["convergence"]

    print(f"  Value function: unique={n_unique}, var={total_var:.4e}")
    print(f"  Fixed: {conv['total_fixed_point']}, Cycling: {conv['total_cycling']} ({conv['cycling_fraction']:.1%})")

    npz_path = f"{RESULTS_DIR}/{model_name}_coalition_clamp1.npz"
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

    wiring_path = f"{RESULTS_DIR}/{model_name}_wiring_clamp1.json"
    wiring_output = {
        "model": model_name,
        "n_nodes": n,
        "node_names": node_names,
        "output_nodes": output_nodes,
        "rule_fourier": rf,
        "convergence": conv,
        "n_unique_values": n_unique,
        "total_variance": total_var,
        "clamp_value": CLAMP_VALUE,
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

    print(f"  Spearman rho={pw['spearman_rho']:.4f} (p={pw['spearman_pvalue']:.2e})")
    print(f"  CI: [{pw['spearman_ci_95'][0]:.4f}, {pw['spearman_ci_95'][1]:.4f}]")
    print(f"  Global o3+={g3p:.1%}, Local o3+={l3p:.1%}, Delta={g3p-l3p:+.1%} ({'creates' if g3p > l3p else 'destroys'})")

    comp_path = f"{RESULTS_DIR}/{model_name}_composition_clamp1.json"
    comp_output = {
        "model": model_name,
        "n_players": n,
        "n_init": n_init,
        "node_names": node_names,
        "clamp_value": CLAMP_VALUE,
        "timestamp": timestamp(),
        **scores,
    }
    with open(comp_path, "w") as f:
        json.dump(comp_output, f, indent=2)
    print(f"  Saved: {comp_path}")

    tri_rho = scores.get("triples", {}).get("spearman_rho", None)
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
    N_WORKERS = 10

    # Load clamp-0 results for comparison
    clamp0_deltas = {}
    with open(f"{RESULTS_DIR}/blind_experiment_summary.json") as f:
        for m in json.load(f):
            clamp0_deltas[m["model"]] = m["delta_o3plus"]
    with open(f"{RESULTS_DIR}/blind_batch2_summary.json") as f:
        for m in json.load(f)["results"]:
            clamp0_deltas[m["model"]] = m["delta_o3plus"]

    all_results = []
    for model_name in ALL_MODELS_SORTED:
        n = len(ALL_KNOWN_MODELS[model_name]["rules"])
        print(f"\n{'='*60}")
        print(f"[{timestamp()}] === {model_name} (n={n}) ===")
        print(f"{'='*60}")
        result = run_model(model_name, n_workers=N_WORKERS)
        all_results.append(result)

        # Intermediate save after each model
        summary_path = f"{RESULTS_DIR}/robustness_clamp1_summary.json"
        with open(summary_path, "w") as f:
            json.dump(all_results, f, indent=2)

    # Print comparison table
    print(f"\n{'='*80}")
    print(f"COMPARISON: clamp-to-0 vs clamp-to-1")
    print(f"{'='*80}")
    print(f"{'Model':35s} {'n':>3s} {'clamp0':>9s} {'clamp1':>9s} {'same?':>6s}")
    print("-" * 80)

    n_same = 0
    n_total_non_null = 0
    NULL_THRESHOLD = 0.005

    for r in all_results:
        name = r["model"]
        d0 = clamp0_deltas.get(name, 0)
        d1 = r["delta_o3plus"]
        d0_pp = d0 * 100
        d1_pp = d1 * 100

        null_0 = abs(d0) < NULL_THRESHOLD
        null_1 = abs(d1) < NULL_THRESHOLD

        if null_0 and null_1:
            same = "null"
        elif null_0 or null_1:
            same = "~"
            n_total_non_null += 1
        else:
            same_dir = (d0 > 0) == (d1 > 0)
            same = "YES" if same_dir else "NO"
            n_total_non_null += 1
            if same_dir:
                n_same += 1

        print(f"{name:35s} {r['n']:3d} {d0_pp:+8.2f}pp {d1_pp:+8.2f}pp {same:>6s}")

    print(f"\nSame direction: {n_same}/{n_total_non_null} non-null models")

    # Count creation/destruction for clamp-1
    labels_1 = []
    for r in all_results:
        d = r["delta_o3plus"]
        if abs(d) < NULL_THRESHOLD:
            labels_1.append("null")
        elif d > 0:
            labels_1.append("creation")
        else:
            labels_1.append("destruction")

    n_c1 = labels_1.count("creation")
    n_d1 = labels_1.count("destruction")
    n_n1 = labels_1.count("null")
    print(f"Clamp-0: 18 creation, 6 destruction, 3 null")
    print(f"Clamp-1: {n_c1} creation, {n_d1} destruction, {n_n1} null")

    print(f"\n[{timestamp()}] All {len(all_results)} models complete.")
    print(f"Summary saved to: {summary_path}")
