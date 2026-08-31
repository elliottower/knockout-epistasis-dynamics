"""Replacement shard 5 for grieco Boolean sweep.

Coalition range [163840, 196608) — the other 7 shards are complete.

Usage:
    modal run --detach scripts/modal_boolean_grieco_shard5.py
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
        os.path.join(REPO, "scripts/run_batch2_blind_sweep.py"),
        "/root/repo/scripts/run_batch2_blind_sweep.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "scripts/structural_analysis.py"),
        "/root/repo/scripts/structural_analysis.py",
        copy=True,
    )
    .run_commands("touch /root/repo/scripts/__init__.py")
)

app = modal.App("epistasis-bench-boolean-grieco-shard5-fix", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

MODEL_NAME = "grieco_bladder"
SHARD_ID = 5
START = 163840
END = 196608
N_INIT = 512
MAX_STEPS = 200
SEED = 42


@app.function(cpu=4, memory=8192, timeout=86400, volumes={"/results": results_vol})
def run_shard5():
    """Run Boolean coalition sweep for shard 5 coalition range."""
    import json
    from datetime import datetime, timezone

    import numpy as np
    from tqdm import tqdm

    from grn_coalition_sweep import compile_network, simulate_sync_output
    from scripts.run_batch2_blind_sweep import WEB_MODELS

    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] Boolean shard {SHARD_ID} FIX: coalitions [{START}, {END})")

    model_info = WEB_MODELS[MODEL_NAME]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]

    compiled, node_names = compile_network(rules)
    n_total = len(node_names)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    player_indices = list(range(n_total))
    n_players = n_total

    rng = np.random.default_rng(SEED)
    init_states = rng.integers(0, 2, size=(N_INIT, n_total), dtype=np.int8)

    n_coalitions = END - START
    values = np.zeros((n_coalitions, N_INIT), dtype=np.float64)
    total_fixed = 0
    total_cycling = 0

    save_interval = 5000
    out_dir = "/results/boolean_full/grieco_shards_v2"
    os.makedirs(out_dir, exist_ok=True)

    for i, coalition in enumerate(tqdm(
        range(START, END),
        desc=f"Bool shard {SHARD_ID} FIX (n={n_players})",
    )):
        clamp_mask = np.zeros(n_total, dtype=bool)
        for bit_pos, node_idx in enumerate(player_indices):
            if not (coalition & (1 << bit_pos)):
                clamp_mask[node_idx] = True

        output, info = simulate_sync_output(
            init_states, compiled, clamp_mask, 0,
            output_indices, max_steps=MAX_STEPS,
        )
        values[i] = output
        total_fixed += info["n_fixed_point"]
        total_cycling += info["n_cycling"]

        if (i + 1) % save_interval == 0:
            np.save(f"{out_dir}/shard_{SHARD_ID:03d}_fix.npy", values[:i + 1])
            meta = {
                "shard_id": SHARD_ID,
                "start_coalition": START,
                "end_coalition": END,
                "n_completed": i + 1,
                "n_total": n_coalitions,
                "total_fixed": total_fixed,
                "total_cycling": total_cycling,
                "timestamp": datetime.now(timezone.utc).isoformat(),
            }
            with open(f"{out_dir}/shard_{SHARD_ID:03d}_fix_meta.json", "w") as f:
                json.dump(meta, f, indent=2)
            results_vol.commit()

    np.save(f"{out_dir}/shard_{SHARD_ID:03d}_fix.npy", values)
    meta = {
        "shard_id": SHARD_ID,
        "start_coalition": START,
        "end_coalition": END,
        "n_coalitions": n_coalitions,
        "n_completed": n_coalitions,
        "total_fixed": total_fixed,
        "total_cycling": total_cycling,
        "status": "complete",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    with open(f"{out_dir}/shard_{SHARD_ID:03d}_fix_meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    results_vol.commit()

    ts_end = datetime.now(timezone.utc).isoformat()
    print(f"[{ts_end}] Shard {SHARD_ID} FIX done: {n_coalitions} coalitions, "
          f"cycling={total_cycling}/{n_coalitions * N_INIT}")
    return meta


@app.local_entrypoint()
def main():
    print(f"Launching shard {SHARD_ID} fix: [{START}, {END})")
    run_shard5.spawn()
    print("Spawned — will run independently after detach.")
