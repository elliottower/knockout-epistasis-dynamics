"""Sharded ODE sweep for grieco_bladder (n=18, 262144 coalitions).

Splits the coalition space across 8 shards (~32K coalitions each,
~12h per shard instead of ~100h sequential).

Usage:
    modal run --detach scripts/modal_ode_shard_grieco.py
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

app = modal.App("epistasis-bench-ode-shard-grieco", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

MODEL_NAME = "grieco_bladder"
N_SHARDS = 8

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

    save_interval = 1000
    out_dir = "/results/ode_full/grieco_shards"
    os.makedirs(out_dir, exist_ok=True)

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

        if (i + 1) % save_interval == 0:
            shard_path = f"{out_dir}/shard_{shard_id:03d}.npy"
            np.save(shard_path, values[:i + 1])
            meta = {
                "shard_id": shard_id,
                "start_coalition": start_coalition,
                "end_coalition": end_coalition,
                "n_completed": i + 1,
                "n_total": n_coalitions,
                "timestamp": datetime.now(timezone.utc).isoformat(),
            }
            with open(f"{out_dir}/shard_{shard_id:03d}_meta.json", "w") as f:
                json.dump(meta, f, indent=2)
            results_vol.commit()

    shard_path = f"{out_dir}/shard_{shard_id:03d}.npy"
    np.save(shard_path, values)

    meta = {
        "shard_id": shard_id,
        "start_coalition": start_coalition,
        "end_coalition": end_coalition,
        "n_coalitions": n_coalitions,
        "n_completed": n_coalitions,
        "status": "complete",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    with open(f"{out_dir}/shard_{shard_id:03d}_meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    results_vol.commit()

    ts_end = datetime.now(timezone.utc).isoformat()
    print(f"[{ts_end}] Shard {shard_id} done: {n_coalitions} coalitions")
    return meta


@app.local_entrypoint()
def main():
    n_players = 18
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

    print("\nAll shards done.")
