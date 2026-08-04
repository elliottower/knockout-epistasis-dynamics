"""Graded perturbation sweep: is the composition gap about composition, or about
complete removal?

Registered in prereg_graded_perturbation_v1.md (3095900), amended 98d9bea and
2afd478.

Every published result clamps knocked-out genes to zero. This sweeps the clamp
VALUE f, so f=0 is loss of function, f=1 is gain of function, and the
intermediate levels are the new measurement.

ODE path only: the Boolean update indexes a truth table via
states.astype(np.int64), so a fractional clamp would silently cast to 0 and
return the knockout result at every level (Amendment 1.1).

TWO THINGS THIS VERSION FIXES, both of which cost a full run:

1. CHECKPOINTS INSIDE THE COALITION LOOP. The previous version wrote one JSON
   after all 2^n coalitions finished, so a unit killed at hour 16 of 17 lost
   everything. It writes partial results every CHECKPOINT_EVERY coalitions and
   resumes from them.

2. REUSES f=0. The saved ODE coalitions in ode_coalitions/*.npz were computed
   with clamp_value defaulting to 0, same seed and ODE params, so they ARE the
   f=0 level. Recomputing them was pure waste. Only the graded levels are new.
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

# f is the clamp value, not a fraction of wild-type (Amendment 1.2): the Hill
# sigmoid bounds states to [0,1]. f=0 is already on the volume, so only the
# graded levels are computed.
GRADED_LEVELS = [0.25, 0.5, 0.75, 1.0]

# Identical to the published ODE arm so f=0 is directly comparable (gate, A2.1).
ODE_PARAMS = {"n_init": 32, "t_max": 30.0, "t_tail": 10.0,
              "tau": 1.0, "hill_n": 10.0, "hill_k": 0.5}
SEED = 42

OUT_DIR = "/results/graded_perturbation"
CKPT_DIR = "/results/graded_perturbation/checkpoints"
EXISTING_F0 = "/results/ode_coalitions"

# Write partial progress this often. Sized so a killed container loses minutes,
# not hours, while keeping commit overhead negligible.
CHECKPOINT_EVERY = 256


@app.function(cpu=8, memory=16384, timeout=86400,
              volumes={"/results": results_vol})
def run_unit(model_name: str, level: float):
    """One (model, clamp level) sweep, resumable at CHECKPOINT_EVERY coalitions."""
    import json
    import time

    import numpy as np
    from scripts.ode_coalition_sweep import (
        ALL_MODELS, simulate_ode_output, score_ode_result,
    )
    from grn_coalition_sweep import identify_input_nodes

    tag = f"{model_name}_f{level:.2f}"
    out_path = f"{OUT_DIR}/{tag}.json"
    ckpt_path = f"{CKPT_DIR}/{tag}.npz"
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(CKPT_DIR, exist_ok=True)

    if os.path.exists(out_path):
        print(f"{tag}: complete, skipping", flush=True)
        return {"model": model_name, "level": level, "status": "cached"}

    info = ALL_MODELS[model_name]
    rules, output_nodes = info["rules"], info["output_nodes"]
    node_names = list(rules.keys())
    n_total = len(node_names)
    name_to_idx = {n: i for i, n in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]
    player_indices = list(range(n_total))
    n_players = len(player_indices)
    N = 2 ** n_players

    rng = np.random.default_rng(SEED)
    init_states = rng.random(size=(ODE_PARAMS["n_init"], n_total))

    # --- resume from checkpoint if one exists -------------------------------
    values = np.zeros((N, ODE_PARAMS["n_init"]))
    start = 0
    if os.path.exists(ckpt_path):
        try:
            ck = np.load(ckpt_path)
            if int(ck["n_players"]) == n_players:
                values = ck["values"]
                start = int(ck["next_index"])
                print(f"{tag}: resuming at coalition {start}/{N}", flush=True)
        except Exception as e:
            print(f"{tag}: checkpoint unreadable ({e}), restarting", flush=True)

    t0 = time.time()
    print(f"{tag}: {N} coalitions, starting at {start}", flush=True)

    for idx in range(start, N):
        clamp_mask = np.zeros(n_total, dtype=bool)
        for bit_pos, node_idx in enumerate(player_indices):
            if not (idx >> bit_pos) & 1:
                clamp_mask[node_idx] = True
        values[idx] = simulate_ode_output(
            rules, node_names, clamp_mask, level,
            output_indices, init_states,
            t_max=ODE_PARAMS["t_max"], t_tail=ODE_PARAMS["t_tail"],
            tau=ODE_PARAMS["tau"], hill_n=ODE_PARAMS["hill_n"],
            hill_k=ODE_PARAMS["hill_k"],
        )

        if (idx + 1) % CHECKPOINT_EVERY == 0 or idx == N - 1:
            np.savez_compressed(ckpt_path, values=values,
                                next_index=idx + 1, n_players=n_players)
            results_vol.commit()
            done = idx + 1
            rate = (time.time() - t0) / max(1, done - start)
            eta = rate * (N - done) / 3600
            print(f"{tag}: {done}/{N} ({100*done/N:.1f}%) "
                  f"eta {eta:.1f}h", flush=True)

    scores = score_ode_result(values, n_players, node_names)
    spectrum = scores["global_spectrum"]
    e3 = float(sum(spectrum[3:]))
    e2 = float(sum(spectrum[2:]))

    payload = {
        "prereg": "prereg_graded_perturbation_v1.md",
        "prereg_commit": "3095900",
        "amendments": ["98d9bea", "2afd478"],
        "model": model_name, "clamp_level": level,
        "n_players": n_players, "n_coalitions": N,
        "ode_params": ODE_PARAMS, "seed": SEED,
        "global_spectrum": [float(x) for x in spectrum],
        "energy_3plus": e3, "energy_2plus": e2,
        "normalised_higher_order": float(e3 / e2) if e2 > 0 else float("nan"),
        "spectrum_gated": bool(e3 < 1e-4),
        "elapsed_s": round(time.time() - t0, 1),
        "source": "computed",
    }
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    if os.path.exists(ckpt_path):
        os.remove(ckpt_path)
    results_vol.commit()
    print(f"{tag}: DONE E3+={e3:.5f}", flush=True)
    return {"model": model_name, "level": level, "status": "done",
            "energy_3plus": e3}


@app.function(cpu=2, memory=8192, timeout=3600,
              volumes={"/results": results_vol})
def adopt_f0(model_name: str):
    """f=0 already exists as an ODE coalition sweep. Score it, don't recompute.

    ode_coalitions/*.npz were produced by modal_ode_sweep_save_coalitions.py,
    which never passes clamp_value, so it used the default 0 with the same seed
    and ODE parameters used here.
    """
    import json

    import numpy as np
    from scripts.ode_coalition_sweep import ALL_MODELS, score_ode_result

    tag = f"{model_name}_f0.00"
    out_path = f"{OUT_DIR}/{tag}.json"
    npz = f"{EXISTING_F0}/{model_name}_ode_coalition.npz"
    os.makedirs(OUT_DIR, exist_ok=True)

    if os.path.exists(out_path):
        return {"model": model_name, "level": 0.0, "status": "cached"}
    if not os.path.exists(npz):
        return {"model": model_name, "level": 0.0, "status": "no_existing_f0"}

    d = np.load(npz, allow_pickle=True)
    v_mean = d["v_mean"]
    n_players = int(d["n_players"])
    node_names = list(ALL_MODELS[model_name]["rules"].keys())

    # score_ode_result expects (N, n_init); the saved artifact is the mean over
    # initial conditions, which is the quantity the spectrum is computed from.
    scores = score_ode_result(v_mean[:, None], n_players, node_names)
    spectrum = scores["global_spectrum"]
    e3 = float(sum(spectrum[3:]))
    e2 = float(sum(spectrum[2:]))

    payload = {
        "prereg": "prereg_graded_perturbation_v1.md",
        "prereg_commit": "3095900",
        "amendments": ["98d9bea", "2afd478"],
        "model": model_name, "clamp_level": 0.0,
        "n_players": n_players, "n_coalitions": int(v_mean.shape[0]),
        "ode_params": ODE_PARAMS, "seed": SEED,
        "global_spectrum": [float(x) for x in spectrum],
        "energy_3plus": e3, "energy_2plus": e2,
        "normalised_higher_order": float(e3 / e2) if e2 > 0 else float("nan"),
        "spectrum_gated": bool(e3 < 1e-4),
        "source": f"adopted from {npz}",
    }
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    results_vol.commit()
    print(f"{tag}: adopted from existing NPZ, E3+={e3:.5f}", flush=True)
    return {"model": model_name, "level": 0.0, "status": "adopted",
            "energy_3plus": e3}


@app.local_entrypoint()
def main(models: str = "", levels: str = "", skip_f0: str = ""):
    ode_models = sorted(
        f[:-len("_ode.json")]
        for f in os.listdir(os.path.join(REPO, "results/ode_full"))
        if f.endswith("_ode.json")
    )
    if models:
        want = {m.strip() for m in models.split(",") if m.strip()}
        ode_models = [m for m in ode_models if m in want]

    lv = [float(x) for x in levels.split(",")] if levels else GRADED_LEVELS

    print(f"{len(ode_models)} models")

    # SPAWN, do not map. A .map/.starmap driver is consumed by THIS process, so
    # if the launcher is backgrounded and later reaped, the remaining units are
    # orphaned and the app shuts down -- which killed two earlier runs even with
    # --detach. Spawned calls are server-side and outlive the launcher.
    handles = []
    if not skip_f0:
        for m in ode_models:
            handles.append(("adopt_f0", m, 0.0, adopt_f0.spawn(m)))
        print(f"  spawned {len(ode_models)} f=0 adoptions")

    units = [(m, f) for m in ode_models for f in lv]
    for m, f in units:
        handles.append(("run_unit", m, f, run_unit.spawn(m, f)))
    print(f"  spawned {len(units)} graded units (levels {lv})")

    print(f"\n{len(handles)} calls submitted. The launcher can now exit safely.")
    print(f"Results: {OUT_DIR}/ on volume epistasis-bench-results")
    print(f"Checkpoints: {CKPT_DIR}/ -- a killed unit resumes from there.")
