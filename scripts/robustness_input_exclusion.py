"""Robustness check: exclude input nodes (f(x)=x) from the player set.

Input nodes are external signals (e.g., GF, TNF, Stimuli) whose Boolean rule
is the identity. They are never influenced by other nodes, so knocking them
out tests an experimental intervention that is impossible in vivo.

This script re-runs coalition sweeps with exclude_inputs=True, fixing input
nodes to 0. The reduced player set yields a smaller coalition table (2^(n-k)
instead of 2^n for k input nodes), so we recompute both global and local
Walsh spectra over the reduced set and compare delta_o3plus with the original.
"""
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from data_utils import energy_by_order, normalized_wht, popcount_array
from grn_coalition_sweep import (
    BUILTIN_MODELS,
    extract_rule_fourier,
    identify_input_nodes,
    sweep_coalitions,
)
from scripts.run_batch2_blind_sweep import WEB_MODELS
from scripts.run_batch2b_extra_models import EXTRA_MODELS

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "grn_v2"


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def collect_all_models():
    all_models = {}
    for src in (BUILTIN_MODELS, WEB_MODELS, EXTRA_MODELS):
        for name, info in src.items():
            if name == "grieco_bladder":
                continue
            all_models[name] = info
    return all_models


def local_energy_spectrum_filtered(rule_fourier, player_names, n_players):
    """Recompute local energy spectrum excluding interactions with input nodes.

    Only counts Walsh coefficients whose regulators are ALL in player_names.
    """
    player_set = set(player_names)
    spectrum = np.zeros(n_players + 1, dtype=np.float64)

    for gene, info in rule_fourier["per_gene"].items():
        if gene not in player_set:
            continue
        for interaction in info.get("interactions", []):
            regs = interaction["regulators"]
            if all(r in player_set for r in regs):
                order = interaction["order"]
                if order <= n_players:
                    spectrum[order] += interaction["coefficient"] ** 2

    total = spectrum.sum()
    if total > 0:
        spectrum /= total
    return spectrum


def compute_delta_o3plus(v_mean, n_players, rule_fourier, player_names):
    """Compute delta_o3plus (global - local order-3+ energy fraction)."""
    w = normalized_wht(v_mean)
    global_energy = energy_by_order(w, n_players)
    total_global = global_energy.sum()
    if total_global > 0:
        global_spectrum = global_energy / total_global
    else:
        global_spectrum = global_energy

    local_spectrum = local_energy_spectrum_filtered(
        rule_fourier, player_names, n_players
    )

    g3p = float(sum(global_spectrum[3:])) if len(global_spectrum) > 3 else 0.0
    l3p = float(sum(local_spectrum[3:])) if len(local_spectrum) > 3 else 0.0
    return g3p - l3p, g3p, l3p, global_spectrum.tolist(), local_spectrum.tolist()


def load_original_delta(model_name):
    """Load original delta_o3plus from merged analysis."""
    merged_path = RESULTS_DIR / "merged_all_27_analysis.json"
    with open(merged_path) as f:
        merged = json.load(f)
    for m in merged["all_models"]:
        if m["model"] == model_name:
            return m["delta_o3plus"], m["label"]
    return None, None


def run_model_exclude_inputs(model_name, rules, output_nodes, input_config=0):
    """Run sweep with input nodes excluded and compute delta_o3plus."""
    input_nodes = identify_input_nodes(rules)
    if not input_nodes:
        return None

    n_total = len(rules)
    n_inputs = len(input_nodes)
    n_players = n_total - n_inputs

    print(f"  {model_name}: n={n_total}, inputs={input_nodes} -> {n_players} players")
    print(f"  Coalitions: 2^{n_players} = {2**n_players}")

    result = sweep_coalitions(
        rules, output_nodes,
        n_init=512, max_steps=200, seed=42,
        update_scheme="sync", clamp_value=0,
        exclude_inputs=True, input_config=input_config,
    )

    v_mean = result["values"].mean(axis=1)
    player_names = result["player_names"]
    rf = extract_rule_fourier(rules)

    delta, g3p, l3p, g_spec, l_spec = compute_delta_o3plus(
        v_mean, n_players, rf, player_names
    )

    conv = result["convergence"]
    cycling = conv.get("cycling_fraction", 0.0)

    return {
        "model": model_name,
        "n_total": n_total,
        "n_players": n_players,
        "input_nodes": input_nodes,
        "input_config": input_config,
        "delta_o3plus_excl": delta,
        "global_o3plus_excl": g3p,
        "local_o3plus_excl": l3p,
        "cycling_fraction_excl": cycling,
        "global_spectrum_excl": g_spec,
        "local_spectrum_excl": l_spec,
    }


if __name__ == "__main__":
    print(f"[{timestamp()}] Input-node exclusion robustness check")
    all_models = collect_all_models()

    models_with_inputs = {}
    models_without_inputs = []
    for name, info in sorted(all_models.items()):
        inputs = identify_input_nodes(info["rules"])
        if inputs:
            models_with_inputs[name] = info
        else:
            models_without_inputs.append(name)

    print(f"  Models with inputs: {len(models_with_inputs)}")
    print(f"  Models without inputs: {len(models_without_inputs)} (results unchanged)")

    results = []
    for model_name in sorted(models_with_inputs.keys()):
        info = models_with_inputs[model_name]
        orig_delta, orig_label = load_original_delta(model_name)

        excl = run_model_exclude_inputs(
            model_name, info["rules"], info["output_nodes"], input_config=0
        )

        if excl is not None and orig_delta is not None:
            excl["delta_o3plus_orig"] = orig_delta
            excl["label_orig"] = orig_label
            excl_label = "creation" if excl["delta_o3plus_excl"] > 0.005 else (
                "destruction" if excl["delta_o3plus_excl"] < -0.005 else "null"
            )
            excl["label_excl"] = excl_label
            excl["sign_flip"] = (orig_label in ("creation", "destruction") and
                                 excl_label in ("creation", "destruction") and
                                 orig_label != excl_label)
            results.append(excl)

            flip_str = " FLIP!" if excl["sign_flip"] else ""
            print(f"    orig={orig_delta:+.4f} ({orig_label})  "
                  f"excl={excl['delta_o3plus_excl']:+.4f} ({excl_label}){flip_str}")

    for name in models_without_inputs:
        orig_delta, orig_label = load_original_delta(name)
        if orig_delta is not None:
            results.append({
                "model": name,
                "n_total": len(all_models[name]["rules"]),
                "n_players": len(all_models[name]["rules"]),
                "input_nodes": [],
                "input_config": None,
                "delta_o3plus_excl": orig_delta,
                "delta_o3plus_orig": orig_delta,
                "label_orig": orig_label,
                "label_excl": orig_label,
                "sign_flip": False,
                "note": "no input nodes, results identical",
            })

    results.sort(key=lambda r: r["model"])

    deltas_orig = [r["delta_o3plus_orig"] for r in results]
    deltas_excl = [r["delta_o3plus_excl"] for r in results]
    n_flips = sum(1 for r in results if r.get("sign_flip", False))
    n_with_inputs = sum(1 for r in results if r["input_nodes"])

    orig_creation = sum(1 for r in results if r["label_orig"] == "creation")
    excl_creation = sum(1 for r in results if r["label_excl"] == "creation")
    orig_destruction = sum(1 for r in results if r["label_orig"] == "destruction")
    excl_destruction = sum(1 for r in results if r["label_excl"] == "destruction")

    rho_deltas, p_deltas = spearmanr(deltas_orig, deltas_excl)

    print(f"\n{'='*60}")
    print(f"SUMMARY (N={len(results)}, {n_with_inputs} with inputs)")
    print(f"{'='*60}")
    print(f"Sign flips (creation<->destruction): {n_flips}/{n_with_inputs}")
    print(f"Creation count: orig={orig_creation}, excl={excl_creation}")
    print(f"Destruction count: orig={orig_destruction}, excl={excl_destruction}")
    print(f"Delta correlation (orig vs excl): rho={rho_deltas:.3f} (p={p_deltas:.2e})")
    print(f"Median delta orig: {np.median(deltas_orig):+.4f}")
    print(f"Median delta excl: {np.median(deltas_excl):+.4f}")

    output = {
        "description": "Input-node exclusion robustness check",
        "input_config": 0,
        "timestamp": timestamp(),
        "summary": {
            "n_models": len(results),
            "n_with_inputs": n_with_inputs,
            "n_without_inputs": len(models_without_inputs),
            "n_sign_flips": n_flips,
            "creation_count_orig": orig_creation,
            "creation_count_excl": excl_creation,
            "destruction_count_orig": orig_destruction,
            "destruction_count_excl": excl_destruction,
            "delta_correlation_rho": float(rho_deltas),
            "delta_correlation_p": float(p_deltas),
            "median_delta_orig": float(np.median(deltas_orig)),
            "median_delta_excl": float(np.median(deltas_excl)),
        },
        "per_model": results,
    }

    out_path = RESULTS_DIR / "robustness_input_exclusion.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {out_path}")
