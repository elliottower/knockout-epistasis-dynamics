"""Sensitivity analysis: convergence of delta_3+ with number of initial conditions.

Runs batch-1 networks at N = 64, 128, 256, 512, 1024 and reports
delta_3+ at each N to verify convergence.

Usage:
    uv run python scripts/sensitivity_n_init.py
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from data_utils import normalized_wht, energy_by_order
from grn_coalition_sweep import BUILTIN_MODELS, extract_rule_fourier, sweep_coalitions

BATCH1 = [
    "faure_cellcycle",
    "davidich_yeast",
    "tournier_apoptosis",
    "drosophila_cellcycle",
    "arabidopsis_cellcycle",
    "fanconi_anemia",
]

N_VALUES = [64, 128, 256, 512, 1024]

from scripts.run_batch2_blind_sweep import WEB_MODELS
from scripts.run_batch2b_extra_models import EXTRA_MODELS

ALL_MODELS = {**BUILTIN_MODELS, **WEB_MODELS, **EXTRA_MODELS}


def run_sensitivity(model_name, n_values):
    model_info = ALL_MODELS[model_name]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]
    n_nodes = len(rules)

    rf = extract_rule_fourier(rules)
    local_spectrum = rf.get("local_energy_spectrum", [])
    local_o3plus = sum(local_spectrum[3:]) if len(local_spectrum) > 3 else 0.0

    results = {}
    for n_init in n_values:
        ts = datetime.now(timezone.utc).strftime("%H:%M:%S")
        print(f"  [{ts}] {model_name} n_init={n_init} (2^{n_nodes}={2**n_nodes} coalitions)")

        sweep = sweep_coalitions(
            rules, output_nodes,
            n_init=n_init, seed=42, max_steps=200,
        )
        values = sweep["values"]
        v_mean = values.mean(axis=1)
        w = normalized_wht(v_mean)
        energy = energy_by_order(w, sweep["n_players"])
        total = energy.sum()
        if total > 0:
            spectrum = (energy / total).tolist()
        else:
            spectrum = energy.tolist()
        global_o3plus = sum(spectrum[3:]) if len(spectrum) > 3 else 0.0
        delta = global_o3plus - local_o3plus

        results[n_init] = {
            "global_o3plus": global_o3plus,
            "delta_o3plus": delta,
            "n_unique": len(np.unique(np.round(v_mean, 6))),
        }
        print(f"         delta_o3+ = {delta:+.4f}")

    return results


def main():
    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] Sensitivity analysis: N convergence on batch-1 networks")
    print(f"N values: {N_VALUES}")
    print()

    all_results = {}
    for model_name in BATCH1:
        n_nodes = len(ALL_MODELS[model_name]["rules"])
        print(f"\n{model_name} (n={n_nodes}):")
        all_results[model_name] = run_sensitivity(model_name, N_VALUES)

    out_path = Path("results/sensitivity/n_init_convergence.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({
            "batch1_networks": BATCH1,
            "n_values": N_VALUES,
            "results": all_results,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }, f, indent=2)
    print(f"\nSaved to {out_path}")

    print("\n\nConvergence table (delta_3+ in pp):")
    header = f"{'Model':<25}"
    for n in N_VALUES:
        header += f"  N={n:>4}"
    header += "  max_diff"
    print(header)
    print("-" * len(header))

    for model_name in BATCH1:
        r = all_results[model_name]
        row = f"{model_name:<25}"
        deltas = []
        for n in N_VALUES:
            d = r[n]["delta_o3plus"] * 100
            deltas.append(d)
            row += f"  {d:>+6.1f}"
        max_diff = max(deltas) - min(deltas)
        row += f"  {max_diff:>8.2f}"
        print(row)


if __name__ == "__main__":
    main()
