"""Parallel ODE sweep for arabidopsis_cellcycle — 8 workers, each handles a chunk.

Usage:
    modal run --detach scripts/modal_ode_arabidopsis_parallel.py
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

app = modal.App("epistasis-ode-arabidopsis-parallel", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

MODEL_NAME = "arabidopsis_cellcycle"
N_CHUNKS = 8
TOTAL_COALITIONS = 16384

ODE_PARAMS = {
    "n_init": 32,
    "t_max": 10.0,
    "t_tail": 5.0,
    "tau": 1.0,
    "hill_n": 10.0,
    "hill_k": 0.5,
}


@app.function(cpu=4, memory=8192, timeout=86400, volumes={"/results": results_vol})
def run_chunk(chunk_id):
    import time

    import numpy as np
    from scripts.ode_coalition_sweep import ALL_MODELS, sweep_coalitions_ode

    chunk_size = (TOTAL_COALITIONS + N_CHUNKS - 1) // N_CHUNKS
    start = chunk_id * chunk_size
    end = min(start + chunk_size, TOTAL_COALITIONS)

    print(f"[{time.strftime('%H:%M:%S')}] Chunk {chunk_id}/{N_CHUNKS}: "
          f"coalitions {start}-{end} ({end - start} total)")

    model_info = ALL_MODELS[MODEL_NAME]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]

    ckpt_dir = "/results/ode_checkpoints"
    os.makedirs(ckpt_dir, exist_ok=True)
    ckpt_path = f"{ckpt_dir}/{MODEL_NAME}_chunk_{chunk_id}.npz"

    result = sweep_coalitions_ode(
        rules, output_nodes,
        n_init=ODE_PARAMS["n_init"],
        seed=42,
        t_max=ODE_PARAMS["t_max"],
        t_tail=ODE_PARAMS["t_tail"],
        tau=ODE_PARAMS["tau"],
        hill_n=ODE_PARAMS["hill_n"],
        hill_k=ODE_PARAMS["hill_k"],
        checkpoint_path=ckpt_path,
        checkpoint_every=100,
        on_checkpoint=lambda: results_vol.commit(),
        coalition_range=(start, end),
    )

    chunk_dir = "/results/ode_chunks"
    os.makedirs(chunk_dir, exist_ok=True)
    chunk_path = f"{chunk_dir}/{MODEL_NAME}_chunk_{chunk_id}.npz"
    np.savez_compressed(
        chunk_path,
        values=result["values"],
        range_start=start,
        range_end=end,
        n_players=result["n_players"],
        player_names=result["player_names"],
        node_names=result["node_names"],
    )
    results_vol.commit()

    print(f"[{time.strftime('%H:%M:%S')}] Chunk {chunk_id} done: saved {chunk_path}")
    return {"chunk_id": chunk_id, "start": start, "end": end, "shape": result["values"].shape}


@app.function(cpu=2, memory=8192, timeout=3600, volumes={"/results": results_vol})
def merge_chunks():
    import json
    import time
    from pathlib import Path

    import numpy as np
    from scripts.ode_coalition_sweep import ALL_MODELS, score_ode_result
    from grn_coalition_sweep import extract_rule_fourier

    print(f"[{time.strftime('%H:%M:%S')}] Merging {N_CHUNKS} chunks...")

    n_init = ODE_PARAMS["n_init"]
    values = np.zeros((TOTAL_COALITIONS, n_init))

    for i in range(N_CHUNKS):
        chunk_path = f"/results/ode_chunks/{MODEL_NAME}_chunk_{i}.npz"
        chunk = np.load(chunk_path)
        start = int(chunk["range_start"])
        end = int(chunk["range_end"])
        values[start:end] = chunk["values"]
        print(f"  Loaded chunk {i}: coalitions {start}-{end}")

    model_info = ALL_MODELS[MODEL_NAME]
    node_names = list(model_info["rules"].keys())
    n_players = len(node_names)
    player_names = node_names

    v_mean = values.mean(axis=1)
    npz_dir = "/results/ode_coalitions"
    os.makedirs(npz_dir, exist_ok=True)
    npz_path = f"{npz_dir}/{MODEL_NAME}_ode_coalition.npz"
    np.savez_compressed(
        npz_path,
        v_mean=v_mean,
        n_players=n_players,
        player_names=player_names,
        node_names=node_names,
        model=MODEL_NAME,
    )
    results_vol.commit()
    print(f"  Saved merged NPZ: {npz_path}")

    scores = score_ode_result(values, n_players, player_names)

    rf = extract_rule_fourier(model_info["rules"])
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
        "npz_path": npz_path,
        "timestamp": time.strftime('%Y-%m-%dT%H:%M:%S+00:00'),
    }

    json_path = f"/results/ode_full/{MODEL_NAME}_ode.json"
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, "w") as f:
        json.dump(output, f, indent=2)
    results_vol.commit()

    for i in range(N_CHUNKS):
        ckpt = f"/results/ode_checkpoints/{MODEL_NAME}_chunk_{i}.npz"
        if Path(ckpt).exists():
            Path(ckpt).unlink()
    results_vol.commit()

    print(f"[{time.strftime('%H:%M:%S')}] Done: delta={ode_delta:+.4f} sign_ok={sign_preserved}")
    return output


@app.local_entrypoint()
def main():
    print(f"Launching {N_CHUNKS} parallel workers for {MODEL_NAME} "
          f"({TOTAL_COALITIONS} coalitions)")

    handles = [run_chunk.spawn(i) for i in range(N_CHUNKS)]
    results = [h.get() for h in handles]

    for r in results:
        print(f"  Chunk {r['chunk_id']}: coalitions {r['start']}-{r['end']} "
              f"shape={r['shape']}")

    print("All chunks done. Merging...")
    final = merge_chunks.remote()
    print(f"Final result: delta={final['ode_delta_o3plus']:+.4f} "
          f"sign_preserved={final['sign_preserved']}")
