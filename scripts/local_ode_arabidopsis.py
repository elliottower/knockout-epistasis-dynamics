"""Run arabidopsis ODE coalition sweep locally, saving v_mean NPZ + JSON.

Usage:
    nohup uv run python scripts/local_ode_arabidopsis.py &
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.ode_coalition_sweep import ALL_MODELS, sweep_coalitions_ode, score_ode_result
from grn_coalition_sweep import extract_rule_fourier

MODEL_NAME = "arabidopsis_cellcycle"
ODE_PARAMS = {"n_init": 32, "t_max": 30.0, "t_tail": 10.0, "tau": 1.0, "hill_n": 10.0, "hill_k": 0.5}

OUT_DIR = Path("results/grn_v2/ode_coalitions")
OUT_DIR.mkdir(parents=True, exist_ok=True)
CKPT_PATH = OUT_DIR / f"{MODEL_NAME}_ode_checkpoint.npz"

model_info = ALL_MODELS[MODEL_NAME]
rules = model_info["rules"]
output_nodes = model_info["output_nodes"]
n_nodes = len(rules)

print(f"[{time.strftime('%H:%M:%S')}] Starting ODE sweep: {MODEL_NAME}")
print(f"  {model_info['description']}")
print(f"  n={n_nodes}, 2^n={2**n_nodes} coalitions, n_init={ODE_PARAMS['n_init']}")

result = sweep_coalitions_ode(
    rules, output_nodes,
    n_init=ODE_PARAMS["n_init"], seed=42,
    t_max=ODE_PARAMS["t_max"], t_tail=ODE_PARAMS["t_tail"],
    tau=ODE_PARAMS["tau"], hill_n=ODE_PARAMS["hill_n"], hill_k=ODE_PARAMS["hill_k"],
    checkpoint_path=str(CKPT_PATH), checkpoint_every=1000,
)

v_mean = result["values"].mean(axis=1)
npz_path = OUT_DIR / f"{MODEL_NAME}_ode_coalition.npz"
np.savez_compressed(
    npz_path,
    v_mean=v_mean,
    n_players=result["n_players"],
    player_names=result["player_names"],
    node_names=result["node_names"],
    model=MODEL_NAME,
)
print(f"  Saved NPZ: {npz_path} ({len(v_mean)} values)")

if CKPT_PATH.exists():
    CKPT_PATH.unlink()

scores = score_ode_result(result["values"], result["n_players"], result["player_names"])
rf = extract_rule_fourier(rules)
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
    "model": MODEL_NAME, "method": "ode_hill",
    "description": model_info["description"],
    "n_players": result["n_players"], "n_init": ODE_PARAMS["n_init"],
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
    "timestamp": time.strftime('%Y-%m-%dT%H:%M:%S'),
}

json_path = OUT_DIR / f"{MODEL_NAME}_ode.json"
with open(json_path, "w") as f:
    json.dump(output, f, indent=2)

print(f"[{time.strftime('%H:%M:%S')}] Done: {MODEL_NAME} delta={ode_delta:+.4f} sign_ok={sign_preserved}")
