"""Graded perturbation sweep: is the composition gap about composition, or about
complete removal?

Registered in prereg_graded_perturbation_v1.md (commit 3095900), with
Amendment 1 (98d9bea) and Amendment 2 (2afd478).

Every published result clamps knocked-out genes to zero. This sweeps the clamp
VALUE f over {0, 0.25, 0.5, 0.75, 1}, so f=0 is loss of function, f=1 is gain of
function, and the intermediate levels are the new measurement.

Runs on the ODE path only. The Boolean update indexes a truth table via
states.astype(np.int64), so a fractional clamp would silently cast to 0 and
return the knockout result at every level (Amendment 1.1).

One container per (model, level). Results commit to the volume after each unit.
"""

from __future__ import annotations

import os

import modal

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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
    .add_local_file(os.path.join(REPO, "data_utils.py"),
                    "/root/repo/data_utils.py", copy=True)
    .add_local_file(os.path.join(REPO, "grn_coalition_sweep.py"),
                    "/root/repo/grn_coalition_sweep.py", copy=True)
    .add_local_file(os.path.join(REPO, "composition_scorer.py"),
                    "/root/repo/composition_scorer.py", copy=True)
    .add_local_file(os.path.join(REPO, "scripts/ode_coalition_sweep.py"),
                    "/root/repo/scripts/ode_coalition_sweep.py", copy=True)
    .add_local_file(os.path.join(REPO, "scripts/run_batch2b_extra_models.py"),
                    "/root/repo/scripts/run_batch2b_extra_models.py", copy=True)
    .add_local_file(os.path.join(REPO, "scripts/run_batch2_blind_sweep.py"),
                    "/root/repo/scripts/run_batch2_blind_sweep.py", copy=True)
    .add_local_file(os.path.join(REPO, "scripts/structural_analysis.py"),
                    "/root/repo/scripts/structural_analysis.py", copy=True)
    .run_commands("touch /root/repo/scripts/__init__.py")
)

app = modal.App("graded-perturbation", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results",
                                     create_if_missing=True)

# Registered in Amendment 1.2. f is the clamp value, not a fraction of
# wild-type: the Hill sigmoid bounds states to [0,1].
LEVELS = [0.0, 0.25, 0.5, 0.75, 1.0]

# Identical to the published ODE arm, so f=0 is comparable to it (gate, A2.1).
ODE_PARAMS = {
    "n_init": 32,
    "t_max": 30.0,
    "t_tail": 10.0,
    "tau": 1.0,
    "hill_n": 10.0,
    "hill_k": 0.5,
}
SEED = 42

OUT_DIR = "/results/graded_perturbation"


@app.function(cpu=8, memory=16384, timeout=86400,
              volumes={"/results": results_vol})
def run_unit(model_name: str, level: float):
    """One (model, clamp level) sweep. Writes its own JSON and commits."""
    import json
    import time

    import numpy as np
    from scripts.ode_coalition_sweep import (
        ALL_MODELS,
        sweep_coalitions_ode,
        score_ode_result,
    )

    tag = f"{model_name}_f{level:.2f}"
    out_path = f"{OUT_DIR}/{tag}.json"
    os.makedirs(OUT_DIR, exist_ok=True)

    if os.path.exists(out_path):
        print(f"[{time.strftime('%H:%M:%S')}] {tag}: already done, skipping")
        return {"model": model_name, "level": level, "status": "cached"}

    if model_name not in ALL_MODELS:
        return {"model": model_name, "level": level, "error": "unknown model"}

    info = ALL_MODELS[model_name]
    t0 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] {tag}: starting "
          f"({info['description']})", flush=True)

    result = sweep_coalitions_ode(
        info["rules"], info["output_nodes"],
        n_init=ODE_PARAMS["n_init"],
        seed=SEED,
        clamp_value=level,
        t_max=ODE_PARAMS["t_max"],
        t_tail=ODE_PARAMS["t_tail"],
        tau=ODE_PARAMS["tau"],
        hill_n=ODE_PARAMS["hill_n"],
        hill_k=ODE_PARAMS["hill_k"],
    )
    scores = score_ode_result(result["values"], result["n_players"],
                              result["player_names"])
    v_mean = result["values"].mean(axis=1)

    # Spectrum gate from Amendment 1.4: as the perturbation weakens the value
    # function flattens, so the headline ratio becomes a quotient of small
    # numbers. Gated networks are reported as gated, never as zero.
    spectrum = scores.get("global_spectrum") or scores.get("ode_global_spectrum")
    energy_3plus = float(sum(spectrum[3:])) if spectrum else float("nan")
    gated = bool(energy_3plus < 1e-4)

    o2plus = float(sum(spectrum[2:])) if spectrum else float("nan")
    ratio = float(energy_3plus / o2plus) if o2plus > 0 else float("nan")

    payload = {
        "prereg": "prereg_graded_perturbation_v1.md",
        "prereg_commit": "3095900",
        "amendments": ["98d9bea", "2afd478"],
        "model": model_name,
        "clamp_level": level,
        "n_players": int(result["n_players"]),
        "n_coalitions": int(result["values"].shape[0]),
        "ode_params": ODE_PARAMS,
        "seed": SEED,
        "global_spectrum": [float(x) for x in spectrum] if spectrum else None,
        "energy_3plus": energy_3plus,
        "energy_2plus": o2plus,
        "normalised_higher_order": ratio,
        "spectrum_gated": gated,
        "v_mean_std": float(np.std(v_mean)),
        "v_mean_range": [float(v_mean.min()), float(v_mean.max())],
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    results_vol.commit()

    print(f"[{time.strftime('%H:%M:%S')}] {tag}: E3+={energy_3plus:.5f} "
          f"ratio={ratio:.4f} gated={gated} "
          f"({payload['elapsed_s']}s)", flush=True)
    return {"model": model_name, "level": level, "status": "done",
            "energy_3plus": energy_3plus}


@app.local_entrypoint()
def main(models: str = "", levels: str = ""):
    """Fan out over (model, level). Smallest networks first so failures surface
    fast rather than after the expensive ones."""
    import time

    # The ODE subset, per Amendment 2.2: these are the tractable networks and
    # the only ones where the f=0 gate can be evaluated.
    ode_models = sorted(
        f[:-len("_ode.json")]
        for f in os.listdir(os.path.join(REPO, "results/ode_full"))
        if f.endswith("_ode.json")
    )
    if models:
        wanted = {m.strip() for m in models.split(",") if m.strip()}
        ode_models = [m for m in ode_models if m in wanted]

    lv = [float(x) for x in levels.split(",")] if levels else LEVELS

    units = [(m, f) for m in ode_models for f in lv]
    print(f"[{time.strftime('%H:%M:%S')}] {len(ode_models)} models x "
          f"{len(lv)} levels = {len(units)} units")
    print(f"  levels: {lv}")
    print(f"  models: {', '.join(ode_models)}")

    done = 0
    for res in run_unit.starmap(units, order_outputs=False):
        done += 1
        if res.get("status") == "done":
            print(f"  [{done}/{len(units)}] {res['model']} "
                  f"f={res['level']:.2f} E3+={res.get('energy_3plus'):.5f}")
        elif res.get("error"):
            print(f"  [{done}/{len(units)}] {res['model']} "
                  f"f={res['level']:.2f} ERROR: {res['error']}")

    print(f"\nResults: {OUT_DIR}/<model>_f<level>.json on volume "
          f"epistasis-bench-results")
