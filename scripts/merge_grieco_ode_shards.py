"""Merge grieco ODE shards into a single result JSON.

Usage:
    uv run python scripts/merge_grieco_ode_shards.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from grn_coalition_sweep import extract_rule_fourier
from scripts.ode_coalition_sweep import ALL_MODELS, score_ode_result

MODEL_NAME = "grieco_bladder"
SHARD_DIR = Path("results/ode_full/grieco_shards/grieco_shards")
OUT_PATH = Path("results/ode_full/grieco_bladder_ode.json")
BOOL_PATH = Path("results/grn_v2/merged_all_27_analysis.json")
N_SHARDS = 8
N_PLAYERS = 18

ODE_PARAMS = {
    "t_max": 30.0,
    "t_tail": 10.0,
    "tau": 1.0,
    "hill_n": 10.0,
    "hill_k": 0.5,
}


def main():
    arrays = []
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

    values = np.concatenate(arrays, axis=0)
    expected = 2**N_PLAYERS
    print(f"Merged: {values.shape} (expected {expected} coalitions)")
    assert values.shape[0] == expected, f"Got {values.shape[0]}, expected {expected}"

    model_info = ALL_MODELS[MODEL_NAME]
    node_names = list(model_info["rules"].keys())
    scores = score_ode_result(values, N_PLAYERS, node_names)

    rf = extract_rule_fourier(model_info["rules"])
    local_spectrum_raw = rf.get("local_energy_spectrum", [])
    local_o3plus = sum(local_spectrum_raw[3:]) if len(local_spectrum_raw) > 3 else 0.0
    ode_delta = scores["global_o3plus"] - local_o3plus

    bool_delta = None
    if BOOL_PATH.exists():
        merged = json.loads(BOOL_PATH.read_text())
        for m in merged["all_models"]:
            if m["model"] == MODEL_NAME:
                bool_delta = m.get("delta_o3plus")
                break

    sign_preserved = None
    if bool_delta is not None:
        bool_sign = "creation" if bool_delta > 0.005 else ("destruction" if bool_delta < -0.005 else "null")
        ode_sign = "creation" if ode_delta > 0.005 else ("destruction" if ode_delta < -0.005 else "null")
        sign_preserved = (bool_sign == ode_sign)

    output = {
        "model": MODEL_NAME,
        "method": "ode_hill",
        "description": model_info["description"],
        "n_players": N_PLAYERS,
        "n_init": 32,
        "ode_params": ODE_PARAMS,
        "local_o3plus": local_o3plus,
        "ode_global_o3plus": scores["global_o3plus"],
        "ode_delta_o3plus": ode_delta,
        "boolean_delta_o3plus": bool_delta,
        "sign_preserved": sign_preserved,
        "ode_global_spectrum": scores["global_spectrum"],
        "n_unique_values": scores["n_unique_values"],
        "v_mean_range": scores["v_mean_range"],
        "v_mean_std": scores["v_mean_std"],
    }

    OUT_PATH.write_text(json.dumps(output, indent=2))
    print(f"\nSaved: {OUT_PATH}")
    print(f"  ODE delta o3+: {ode_delta:+.4f}")
    print(f"  local o3+:     {local_o3plus:.4f}")
    print(f"  ODE global o3+: {scores['global_o3plus']:.4f}")
    if bool_delta is not None:
        print(f"  Boolean delta: {bool_delta:+.4f}")
        print(f"  Sign preserved: {sign_preserved}")


if __name__ == "__main__":
    main()
