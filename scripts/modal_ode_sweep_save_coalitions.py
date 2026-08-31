"""Re-run ODE coalition sweep for 26 networks (skip n=17,18), saving raw coalition values.

Same ODE integration as modal_ode_sweep.py, but saves the v_mean vector as NPZ
so the sample-budget recovery experiment can run on continuous-dynamics data.

Usage:
    modal run -d scripts/modal_ode_sweep_save_coalitions.py
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

app = modal.App("ode-sweep-save-coalitions", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

SKIP_MODELS = {"calzone_cell_fate", "grieco_bladder"}

ODE_PARAMS = {
    "n_init": 32,
    "t_max": 30.0,
    "t_tail": 10.0,
    "tau": 1.0,
    "hill_n": 10.0,
    "hill_k": 0.5,
}


@app.function(cpu=4, memory=16384, timeout=86400, volumes={"/results": results_vol})
def run_single_model(model_name: str):
    import json
    import time
    from pathlib import Path

    import numpy as np
    from scripts.ode_coalition_sweep import (
        ALL_MODELS,
        sweep_coalitions_ode,
        score_ode_result,
    )
    from grn_coalition_sweep import extract_rule_fourier

    print(f"[{time.strftime('%H:%M:%S')}] Starting ODE sweep: {model_name}")

    if model_name not in ALL_MODELS:
        print(f"  Unknown model: {model_name}")
        return {"model": model_name, "error": "unknown model"}

    model_info = ALL_MODELS[model_name]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]
    n_nodes = len(rules)

    print(f"  {model_info['description']}")
    print(f"  n={n_nodes}, 2^n={2**n_nodes} coalitions, n_init={ODE_PARAMS['n_init']}")

    ckpt_dir = "/results/ode_checkpoints"
    os.makedirs(ckpt_dir, exist_ok=True)
    ckpt_path = f"{ckpt_dir}/{model_name}_ode_checkpoint.npz"

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
        checkpoint_every=1000,
        on_checkpoint=lambda: results_vol.commit(),
    )

    v_mean = result["values"].mean(axis=1)

    npz_dir = "/results/ode_coalitions"
    os.makedirs(npz_dir, exist_ok=True)
    npz_path = f"{npz_dir}/{model_name}_ode_coalition.npz"
    np.savez_compressed(
        npz_path,
        v_mean=v_mean,
        n_players=result["n_players"],
        player_names=result["player_names"],
        node_names=result["node_names"],
        model=model_name,
    )
    results_vol.commit()
    print(f"  Saved coalition NPZ: {npz_path} ({len(v_mean)} values)")

    if Path(ckpt_path).exists():
        Path(ckpt_path).unlink()
        results_vol.commit()

    scores = score_ode_result(
        result["values"], result["n_players"], result["player_names"]
    )

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
            if m["model"] == model_name:
                bool_delta = m.get("delta_o3plus")
                break

    sign_preserved = None
    if bool_delta is not None:
        bool_sign = "creation" if bool_delta > 0.005 else ("destruction" if bool_delta < -0.005 else "null")
        ode_sign = "creation" if ode_delta > 0.005 else ("destruction" if ode_delta < -0.005 else "null")
        sign_preserved = (bool_sign == ode_sign)

    output = {
        "model": model_name,
        "method": "ode_hill",
        "description": model_info["description"],
        "n_players": result["n_players"],
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

    json_path = f"/results/ode_full/{model_name}_ode.json"
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, "w") as f:
        json.dump(output, f, indent=2)
    results_vol.commit()

    print(f"[{time.strftime('%H:%M:%S')}] Done: {model_name} "
          f"delta={ode_delta:+.4f} sign_ok={sign_preserved}")
    return output


@app.local_entrypoint()
def main():
    import json
    import time

    merged_path = os.path.join(REPO, "results/grn_v2/merged_all_27_analysis.json")
    with open(merged_path) as f:
        merged = json.load(f)
    all_models = sorted(m["model"] for m in merged["all_models"])
    all_models.append("grieco_bladder")

    model_list = [m for m in all_models if m not in SKIP_MODELS]
    skipped = [m for m in all_models if m in SKIP_MODELS]

    print(f"ODE coalition sweep (saving raw values)")
    print(f"  Running: {len(model_list)} models")
    print(f"  Skipping: {skipped} (too large)")
    print(f"  ODE params: {ODE_PARAMS}")

    t0 = time.time()
    results = list(run_single_model.map(model_list))
    elapsed = time.time() - t0

    completed = [r for r in results if "error" not in r]
    failed = [r for r in results if "error" in r]

    print(f"\nDone in {elapsed:.0f}s")
    print(f"  Completed: {len(completed)}/{len(model_list)}")
    if failed:
        print(f"  Failed: {[r.get('model', '?') for r in failed]}")

    sign_ok = sum(1 for r in completed if r.get("sign_preserved") is True)
    sign_total = sum(1 for r in completed if r.get("sign_preserved") is not None)
    print(f"  Sign preserved: {sign_ok}/{sign_total}")

    for r in sorted(completed, key=lambda x: x.get("n_players", 0)):
        print(f"    n={r['n_players']:2d} {r['model']:40s} "
              f"delta={r['ode_delta_o3plus']:+.4f} sign_ok={r.get('sign_preserved')}")
