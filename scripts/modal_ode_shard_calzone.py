"""Sharded ODE sweep for calzone_cell_fate (n=17, 131072 coalitions).

Splits the coalition space across N_SHARDS containers (~8K coalitions
each, ~2.5h per shard instead of ~44h sequential).

Usage:
    modal run --detach scripts/modal_ode_shard_calzone.py
"""
from __future__ import annotations

import os

import modal

REPO = os.path.join(os.path.expanduser("~"), "Documents/GitHub/epistasis-bench")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(
        "numpy==2.2.6",
        "scipy==1.15.3",
        "scikit-learn==1.6.1",
        "tqdm==4.67.1",
        "matplotlib==3.10.3",
    )
    .env({"PYTHONPATH": "/root/repo:/root/repo/scripts"})
    .add_local_file(os.path.join(REPO, "data_utils.py"), "/root/repo/data_utils.py", copy=True)
    .add_local_file(os.path.join(REPO, "grn_coalition_sweep.py"), "/root/repo/grn_coalition_sweep.py", copy=True)
    .add_local_file(os.path.join(REPO, "composition_scorer.py"), "/root/repo/composition_scorer.py", copy=True)
    .add_local_file(
        os.path.join(REPO, "scripts/ode_coalition_sweep.py"),
        "/root/repo/scripts/ode_coalition_sweep.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "scripts/run_batch2b_extra_models.py"),
        "/root/repo/scripts/run_batch2b_extra_models.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "scripts/run_batch2_blind_sweep.py"),
        "/root/repo/scripts/run_batch2_blind_sweep.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "scripts/structural_analysis.py"),
        "/root/repo/scripts/structural_analysis.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "results/grn_v2/merged_all_27_analysis.json"),
        "/root/repo/results/grn_v2/merged_all_27_analysis.json",
        copy=True,
    )
    .run_commands("touch /root/repo/scripts/__init__.py")
)

app = modal.App("epistasis-bench-ode-shard-calzone", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

MODEL_NAME = "calzone_cell_fate"
N_SHARDS = 4

ODE_PARAMS = {
    "n_init": 32,
    "t_max": 30.0,
    "t_tail": 10.0,
    "tau": 1.0,
    "hill_n": 10.0,
    "hill_k": 0.5,
}


@app.function(cpu=4, memory=8192, timeout=86400, volumes={"/results": results_vol})
def run_shard(shard_id: int, start_coalition: int, end_coalition: int):
    """Run ODE sweep for a range of coalitions."""
    import json
    from datetime import datetime, timezone

    import numpy as np
    from tqdm import tqdm

    from scripts.ode_coalition_sweep import (
        ALL_MODELS,
        simulate_ode_output,
    )
    from grn_coalition_sweep import identify_input_nodes

    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] Shard {shard_id}: coalitions [{start_coalition}, {end_coalition})")

    model_info = ALL_MODELS[MODEL_NAME]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]
    node_names = list(rules.keys())
    n_total = len(node_names)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    player_indices = list(range(n_total))
    n_players = len(player_indices)

    rng = np.random.default_rng(42)
    init_states = rng.random(size=(ODE_PARAMS["n_init"], n_total))

    n_coalitions = end_coalition - start_coalition
    values = np.zeros((n_coalitions, ODE_PARAMS["n_init"]))

    for i, coalition in enumerate(tqdm(
        range(start_coalition, end_coalition),
        desc=f"Shard {shard_id} (n={n_players})",
    )):
        clamp_mask = np.zeros(n_total, dtype=bool)
        for bit_pos, node_idx in enumerate(player_indices):
            if not (coalition & (1 << bit_pos)):
                clamp_mask[node_idx] = True

        outputs = simulate_ode_output(
            rules, node_names, clamp_mask, 0,
            output_indices, init_states,
            t_max=ODE_PARAMS["t_max"],
            t_tail=ODE_PARAMS["t_tail"],
            tau=ODE_PARAMS["tau"],
            hill_n=ODE_PARAMS["hill_n"],
            hill_k=ODE_PARAMS["hill_k"],
        )
        values[i] = outputs

    out_dir = "/results/ode_full/calzone_shards"
    os.makedirs(out_dir, exist_ok=True)
    shard_path = f"{out_dir}/shard_{shard_id:03d}.npy"
    np.save(shard_path, values)

    meta = {
        "shard_id": shard_id,
        "start_coalition": start_coalition,
        "end_coalition": end_coalition,
        "n_coalitions": n_coalitions,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    with open(f"{out_dir}/shard_{shard_id:03d}_meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    results_vol.commit()

    ts_end = datetime.now(timezone.utc).isoformat()
    print(f"[{ts_end}] Shard {shard_id} done: {n_coalitions} coalitions")
    return meta


@app.function(cpu=2, memory=8192, timeout=86400, volumes={"/results": results_vol})
def merge_shards():
    """Merge all shard results into final calzone_cell_fate_ode.json."""
    import json
    from datetime import datetime, timezone
    from pathlib import Path

    import numpy as np

    from scripts.ode_coalition_sweep import ALL_MODELS, score_ode_result
    from grn_coalition_sweep import extract_rule_fourier

    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] Merging calzone_cell_fate shards")

    model_info = ALL_MODELS[MODEL_NAME]
    rules = model_info["rules"]
    n_players = len(rules)
    N = 2**n_players

    shard_dir = Path("/results/ode_full/calzone_shards")
    results_vol.reload()
    shard_files = sorted(shard_dir.glob("shard_*.npy"))
    print(f"  Found {len(shard_files)} shard files")

    all_values = []
    for sf in shard_files:
        arr = np.load(sf)
        all_values.append(arr)
        print(f"  {sf.name}: {arr.shape}")

    values = np.concatenate(all_values, axis=0)
    print(f"  Total values shape: {values.shape} (expected ({N}, {ODE_PARAMS['n_init']}))")
    assert values.shape[0] == N, f"Expected {N} coalitions, got {values.shape[0]}"

    scores = score_ode_result(values, n_players, list(rules.keys()))

    rf = extract_rule_fourier(rules)
    local_spectrum_raw = rf.get("local_energy_spectrum", [])
    local_o3plus = sum(local_spectrum_raw[3:]) if len(local_spectrum_raw) > 3 else 0.0
    ode_delta = scores["global_o3plus"] - local_o3plus

    bool_path = Path("/root/repo/results/grn_v2/merged_all_27_analysis.json")
    bool_delta = None
    if bool_path.exists():
        with open(bool_path) as f:
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
        "n_init": ODE_PARAMS["n_init"],
        "n_shards": len(shard_files),
        "ode_params": {
            "t_max": ODE_PARAMS["t_max"],
            "t_tail": ODE_PARAMS["t_tail"],
            "tau": ODE_PARAMS["tau"],
            "hill_n": ODE_PARAMS["hill_n"],
            "hill_k": ODE_PARAMS["hill_k"],
        },
        "local_o3plus": local_o3plus,
        "ode_global_o3plus": scores["global_o3plus"],
        "ode_delta_o3plus": ode_delta,
        "boolean_delta_o3plus": bool_delta,
        "sign_preserved": sign_preserved,
        "ode_global_spectrum": scores["global_spectrum"],
        "n_unique_values": scores["n_unique_values"],
        "v_mean_range": scores["v_mean_range"],
        "v_mean_std": scores["v_mean_std"],
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = "/results/ode_full/calzone_cell_fate_ode.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    results_vol.commit()

    ts_end = datetime.now(timezone.utc).isoformat()
    print(f"[{ts_end}] Merged: ODE delta={ode_delta:+.4f}, sign_preserved={sign_preserved}")
    return output


@app.local_entrypoint()
def main():
    import json

    merged_path = os.path.join(REPO, "results/grn_v2/merged_all_27_analysis.json")
    with open(merged_path) as f:
        merged = json.load(f)

    n_players = None
    for m in merged["all_models"]:
        if m["model"] == MODEL_NAME:
            n_players = m.get("n_players", m.get("n_nodes", m.get("n")))
            break

    if n_players is None:
        n_players = 17

    N = 2**n_players
    chunk = (N + N_SHARDS - 1) // N_SHARDS

    print(f"Sharding {MODEL_NAME}: n={n_players}, N={N}, {N_SHARDS} shards of ~{chunk}")

    futures = []
    for shard_id in range(N_SHARDS):
        start = shard_id * chunk
        end = min(start + chunk, N)
        if start >= N:
            break
        futures.append(run_shard.spawn(shard_id, start, end))

    for future in futures:
        result = future.get()
        print(f"  Shard {result['shard_id']} done: [{result['start_coalition']}, {result['end_coalition']})")

    print("\nAll shards done. Merging...")
    final = merge_shards.remote()
    print(f"\nResult: ODE delta={final['ode_delta_o3plus']:+.4f}, sign_preserved={final['sign_preserved']}")
