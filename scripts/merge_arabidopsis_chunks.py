"""Merge 8 parallel ODE chunks for arabidopsis_cellcycle into final NPZ + JSON."""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.ode_coalition_sweep import ALL_MODELS, score_ode_result
from grn_coalition_sweep import extract_rule_fourier

MODEL_NAME = "arabidopsis_cellcycle"
N_CHUNKS = 8
TOTAL_COALITIONS = 16384
N_INIT = 32

CHUNK_DIR = Path("results/grn_v2/ode_chunks")
OUT_DIR = Path("results/grn_v2/ode_coalitions")
JSON_DIR = Path("results/grn_v2/ode_full")

ODE_PARAMS = {
    "n_init": N_INIT,
    "t_max": 10.0,
    "t_tail": 5.0,
    "tau": 1.0,
    "hill_n": 10.0,
    "hill_k": 0.5,
}

print(f"Merging {N_CHUNKS} chunks for {MODEL_NAME}...")

values = np.zeros((TOTAL_COALITIONS, N_INIT))
for i in range(N_CHUNKS):
    chunk_path = CHUNK_DIR / f"{MODEL_NAME}_chunk_{i}.npz"
    chunk = np.load(chunk_path)
    start = int(chunk["range_start"])
    end = int(chunk["range_end"])
    values[start:end] = chunk["values"]
    print(f"  Chunk {i}: coalitions {start}-{end}, shape={chunk['values'].shape}")

model_info = ALL_MODELS[MODEL_NAME]
node_names = list(model_info["rules"].keys())
n_players = len(node_names)
player_names = node_names

v_mean = values.mean(axis=1)
OUT_DIR.mkdir(parents=True, exist_ok=True)
npz_path = OUT_DIR / f"{MODEL_NAME}_ode_coalition.npz"
np.savez_compressed(
    npz_path,
    v_mean=v_mean,
    n_players=n_players,
    player_names=player_names,
    node_names=node_names,
    model=MODEL_NAME,
)
print(f"Saved merged NPZ: {npz_path}")

scores = score_ode_result(values, n_players, player_names)

rf = extract_rule_fourier(model_info["rules"])
local_spectrum_raw = rf.get("local_energy_spectrum", [])
local_o3plus = sum(local_spectrum_raw[3:]) if len(local_spectrum_raw) > 3 else 0.0
ode_delta = scores["global_o3plus"] - local_o3plus

merged_path = Path("results/grn_v2/merged_all_27_analysis.json")
bool_delta = None
if merged_path.exists():
    with open(merged_path) as f:
        merged = json.load(f)
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
    "n_players": n_players,
    "n_init": N_INIT,
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
    "npz_path": str(npz_path),
    "timestamp": time.strftime('%Y-%m-%dT%H:%M:%S+00:00'),
}

JSON_DIR.mkdir(parents=True, exist_ok=True)
json_path = JSON_DIR / f"{MODEL_NAME}_ode.json"
with open(json_path, "w") as f:
    json.dump(output, f, indent=2)

print(f"\nResult:")
print(f"  local_o3plus:  {local_o3plus:.4f}")
print(f"  ode_global:    {scores['global_o3plus']:.4f}")
print(f"  ode_delta:     {ode_delta:+.4f}")
print(f"  bool_delta:    {bool_delta:+.4f}" if bool_delta is not None else "  bool_delta:    N/A")
print(f"  sign_preserved: {sign_preserved}")
print(f"  spectrum:      {scores['global_spectrum']}")
print(f"\nSaved JSON: {json_path}")
