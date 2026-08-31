"""Modal launcher for ODE coalition sweep across all 27 GRN models.

Converts Boolean rules to Hill-function ODEs and runs exhaustive
coalition sweeps on Modal CPU instances, one model per container.

Usage:
    modal run --detach scripts/modal_ode_sweep.py

    # Subset for testing:
    modal run scripts/modal_ode_sweep.py --models lambda_phage,faure_cellcycle
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

app = modal.App("epistasis-bench-ode-sweep", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

ODE_PARAMS = {
    "n_init": 32,
    "t_max": 30.0,
    "t_tail": 10.0,
    "tau": 1.0,
    "hill_n": 10.0,
    "hill_k": 0.5,
}


@app.function(cpu=4, memory=8192, timeout=86400, volumes={"/results": results_vol})
def run_single_model(model_name: str):
    """Run ODE coalition sweep for a single model."""
    import json
    from datetime import datetime, timezone
    from pathlib import Path

    from scripts.ode_coalition_sweep import (
        ALL_MODELS,
        sweep_coalitions_ode,
        score_ode_result,
    )
    from grn_coalition_sweep import extract_rule_fourier

    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] Starting ODE sweep: {model_name}")

    if model_name not in ALL_MODELS:
        print(f"  Unknown model: {model_name}")
        return {"model": model_name, "error": "unknown model"}

    model_info = ALL_MODELS[model_name]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]
    n_nodes = len(rules)

    print(f"  {model_info['description']}")
    print(f"  n={n_nodes}, 2^n={2**n_nodes} coalitions, n_init={ODE_PARAMS['n_init']}")

    result = sweep_coalitions_ode(
        rules, output_nodes,
        n_init=ODE_PARAMS["n_init"],
        seed=42,
        t_max=ODE_PARAMS["t_max"],
        t_tail=ODE_PARAMS["t_tail"],
        tau=ODE_PARAMS["tau"],
        hill_n=ODE_PARAMS["hill_n"],
        hill_k=ODE_PARAMS["hill_k"],
    )

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

    out_path = f"/results/ode_full/{model_name}_ode.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    results_vol.commit()

    ts_end = datetime.now(timezone.utc).isoformat()
    print(f"[{ts_end}] Done: {model_name}")
    print(f"  ODE delta o3+: {ode_delta:+.4f}")
    if bool_delta is not None:
        print(f"  Boolean delta: {bool_delta:+.4f}")
        print(f"  Sign preserved: {sign_preserved}")

    return output


@app.function(cpu=2, memory=4096, timeout=86400, volumes={"/results": results_vol})
def orchestrate(model_names: list[str]):
    """Fan out all models in parallel, collect results."""
    import json
    from datetime import datetime, timezone

    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] ODE sweep: launching {len(model_names)} models in parallel")
    print(f"  Models: {model_names}")
    print(f"  ODE params: {ODE_PARAMS}")

    futures = []
    for name in model_names:
        futures.append(run_single_model.spawn(name))

    all_results = []
    completed = 0
    failed = 0
    for future in futures:
        try:
            result = future.get()
            all_results.append(result)
            completed += 1
            ts_now = datetime.now(timezone.utc).isoformat()
            model = result.get("model", "?")
            sign = result.get("sign_preserved", "?")
            delta = result.get("ode_delta_o3plus", 0)
            print(f"[{ts_now}] {completed}/{len(futures)}: {model} delta={delta:+.4f} sign_preserved={sign}")
        except Exception as e:
            failed += 1
            print(f"  FAILED: {e}")
            all_results.append({"error": str(e)})

    summary = {
        "n_models": len(model_names),
        "n_completed": completed,
        "n_failed": failed,
        "ode_params": ODE_PARAMS,
        "timestamp_start": ts,
        "timestamp_end": datetime.now(timezone.utc).isoformat(),
        "models": all_results,
    }

    os.makedirs("/results/ode_full", exist_ok=True)
    with open("/results/ode_full/ode_sweep_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    results_vol.commit()

    sign_preserved_count = sum(
        1 for r in all_results
        if isinstance(r, dict) and r.get("sign_preserved") is True
    )
    sign_total = sum(
        1 for r in all_results
        if isinstance(r, dict) and r.get("sign_preserved") is not None
    )
    print(f"\n[{summary['timestamp_end']}] Done: {completed}/{len(futures)} models")
    print(f"  Sign preserved: {sign_preserved_count}/{sign_total}")

    return summary


ALREADY_DONE = {
    "lambda_phage",
    "arellano_rootstem",
    "faure_cellcycle",
    "davidich_yeast",
    "albert_segment_polarity",
    "asymmetric_cell_division",
    "blood_stem_cell",
    "calzone_cellfate_reduced",
    "cell_cycle_transcription",
    "li_budding_yeast",
    "myeloid_progenitors",
    "remy_p53_mdm2",
    "calzone_cell_fate",
}


@app.local_entrypoint()
def main(models: str = "", include_done: bool = False):
    import json

    if models:
        model_list = models.split(",")
    else:
        merged_path = os.path.join(REPO, "results/grn_v2/merged_all_27_analysis.json")
        with open(merged_path) as f:
            merged = json.load(f)
        model_list = sorted(m["model"] for m in merged["all_models"])

    if not include_done:
        skipped = [m for m in model_list if m in ALREADY_DONE]
        model_list = [m for m in model_list if m not in ALREADY_DONE]
        if skipped:
            print(f"Skipping {len(skipped)} already-done pilots: {skipped}")

    print(f"Running ODE sweep on {len(model_list)} models")
    result = orchestrate.remote(model_list)
    print(f"\nCompleted: {result['n_completed']}, Failed: {result['n_failed']}")
