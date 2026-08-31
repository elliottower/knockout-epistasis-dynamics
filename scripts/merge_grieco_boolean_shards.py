"""Merge grieco Boolean shards into a single result and add to merged analysis.

Usage:
    uv run python scripts/merge_grieco_boolean_shards.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from composition_scorer import normalized_wht, energy_by_order
from grn_coalition_sweep import compile_network, extract_rule_fourier
from scripts.run_batch2_blind_sweep import WEB_MODELS

MODEL_NAME = "grieco_bladder"
SHARD_DIR = Path("results/boolean_full/grieco_shards_v2")
MERGED_PATH = Path("results/grn_v2/merged_all_27_analysis.json")
OUT_ANALYSIS = Path("results/grn_v2/grieco_bladder_analysis.json")
N_SHARDS = 8
N_PLAYERS = 18
N_INIT = 512


def main():
    arrays = []
    total_fixed = 0
    total_cycling = 0

    for i in range(N_SHARDS):
        shard_path = SHARD_DIR / f"shard_{i:03d}.npy"
        meta_path = SHARD_DIR / f"shard_{i:03d}_meta.json"
        if not shard_path.exists():
            print(f"Missing shard {i}: {shard_path}")
            return

        meta = json.loads(meta_path.read_text())
        if meta.get("status") != "complete":
            print(f"Shard {i} not complete: {meta}")
            return

        arr = np.load(shard_path)
        print(f"Shard {i}: {arr.shape} [{meta['start_coalition']}, {meta['end_coalition']})")
        arrays.append(arr)
        total_fixed += meta.get("total_fixed", 0)
        total_cycling += meta.get("total_cycling", 0)

    values = np.concatenate(arrays, axis=0)
    expected = 2**N_PLAYERS
    print(f"Merged: {values.shape} (expected {expected} coalitions)")
    assert values.shape[0] == expected, f"Got {values.shape[0]}, expected {expected}"

    v_mean = values.mean(axis=1)
    w = normalized_wht(v_mean)
    energy = energy_by_order(w, N_PLAYERS)
    total_energy = energy.sum()
    spectrum = (energy / total_energy).tolist() if total_energy > 0 else energy.tolist()
    global_o3plus = sum(spectrum[3:]) if len(spectrum) > 3 else 0.0

    model_info = WEB_MODELS[MODEL_NAME]
    rf = extract_rule_fourier(model_info["rules"])
    local_spectrum_raw = rf.get("local_energy_spectrum", [])
    local_o3plus = sum(local_spectrum_raw[3:]) if len(local_spectrum_raw) > 3 else 0.0
    delta = global_o3plus - local_o3plus

    result = {
        "model": MODEL_NAME,
        "description": model_info.get("description", "Grieco bladder cancer (n=18)"),
        "n_players": N_PLAYERS,
        "n_init": N_INIT,
        "global_spectrum": spectrum,
        "global_o3plus": global_o3plus,
        "local_spectrum": rf.get("local_energy_spectrum", []),
        "local_o3plus": local_o3plus,
        "delta_o3plus": delta,
        "n_fixed_point": total_fixed,
        "n_cycling": total_cycling,
    }

    OUT_ANALYSIS.write_text(json.dumps(result, indent=2))
    print(f"\nSaved: {OUT_ANALYSIS}")
    print(f"  Boolean delta o3+: {delta:+.4f} ({delta*100:+.1f} pp)")
    print(f"  local o3+: {local_o3plus:.4f}")
    print(f"  global o3+: {global_o3plus:.4f}")

    if MERGED_PATH.exists():
        merged = json.loads(MERGED_PATH.read_text())
        existing_names = {m["model"] for m in merged["all_models"]}
        if MODEL_NAME not in existing_names:
            merged["all_models"].append(result)
            merged["n_models"] = len(merged["all_models"])
            merged_28_path = Path("results/grn_v2/merged_all_28_analysis.json")
            merged_28_path.write_text(json.dumps(merged, indent=2))
            print(f"  Added to merged analysis: {merged_28_path}")
        else:
            print(f"  {MODEL_NAME} already in merged analysis")


if __name__ == "__main__":
    main()
