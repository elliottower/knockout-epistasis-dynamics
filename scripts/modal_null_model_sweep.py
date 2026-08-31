"""Modal sweep: null model comparison for all 26 Boolean networks (n=7-15).

Runs all three null types (rule-preserving, degree-preserving, Kauffman NK).
Each (model, null_type) pair runs in its own container. Skips calzone_cell_fate
(n=17) — 131K coalitions per null is too expensive.

Stratified replicates by size:
  n <= 11: 100 nulls
  n = 12-14: 30 nulls
  n = 15: 10 nulls

Usage:
    modal run --detach scripts/modal_null_model_sweep.py
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
        "tqdm==4.67.1",
        "matplotlib==3.10.3",
    )
    .env({"PYTHONPATH": "/root/repo:/root/repo/scripts"})
    .add_local_file(os.path.join(REPO, "data_utils.py"), "/root/repo/data_utils.py", copy=True)
    .add_local_file(os.path.join(REPO, "grn_coalition_sweep.py"), "/root/repo/grn_coalition_sweep.py", copy=True)
    .add_local_file(os.path.join(REPO, "composition_scorer.py"), "/root/repo/composition_scorer.py", copy=True)
    .add_local_file(
        os.path.join(REPO, "scripts/null_model_generators.py"),
        "/root/repo/scripts/null_model_generators.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "scripts/run_batch2_blind_sweep.py"),
        "/root/repo/scripts/run_batch2_blind_sweep.py",
        copy=True,
    )
    .add_local_file(
        os.path.join(REPO, "scripts/run_batch2b_extra_models.py"),
        "/root/repo/scripts/run_batch2b_extra_models.py",
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

app = modal.App("epistasis-bench-null-model-sweep", image=image)
results_vol = modal.Volume.from_name("epistasis-bench-results", create_if_missing=True)

NULL_TYPES = ["rule_preserving", "degree_preserving", "kauffman_nk"]

def n_nulls_for_size(n):
    """Stratified replicate count by network size."""
    if n <= 11:
        return 100
    if n <= 14:
        return 30
    return 10  # n=15


@app.function(cpu=4, memory=8192, timeout=86400, volumes={"/results": results_vol})
def run_null_for_model(model_name: str, rules: dict, output_nodes: list,
                       null_type: str, n_nulls: int, n_init: int = 128,
                       seed_base: int = 42):
    """Run one null type for one model on Modal."""
    import json
    from datetime import datetime, timezone

    import numpy as np
    from tqdm import tqdm

    from data_utils import wht
    from grn_coalition_sweep import compile_network, extract_rule_fourier
    from scripts.null_model_generators import (
        compute_delta_3plus,
        degree_preserving_rewire,
        kauffman_nk,
        quality_gate,
        rule_preserving_rewire,
        sweep_compiled,
    )

    generators = {
        "rule_preserving": rule_preserving_rewire,
        "degree_preserving": degree_preserving_rewire,
        "kauffman_nk": kauffman_nk,
    }
    generator = generators[null_type]

    ts = datetime.now(timezone.utc).isoformat()
    n = len(rules)
    print(f"[{ts}] {model_name} (n={n}), null={null_type}, {n_nulls} replicates")

    compiled, node_names = compile_network(rules)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    local_fourier = extract_rule_fourier(rules)
    local_spectrum = local_fourier.get("local_energy_spectrum", [])
    local_o3plus = sum(local_spectrum[3:]) if len(local_spectrum) > 3 else 0.0

    def extract_local_compiled(comp, names):
        max_k = max((len(regs) for _, regs, _ in comp), default=0)
        local_energy = np.zeros(max_k + 1, dtype=np.float64)
        for _, reg_indices, truth_table in comp:
            k = len(reg_indices)
            if k == 0:
                continue
            w = wht(truth_table.astype(np.float64)) / (2**k)
            for T in range(2**k):
                order = bin(T).count("1")
                if order <= max_k:
                    local_energy[order] += w[T] ** 2
        total = local_energy.sum()
        if total < 1e-15:
            return local_energy.tolist()
        return (local_energy / total).tolist()

    all_deltas = []
    gated_deltas = []
    n_degenerate = 0
    n_zero_energy = 0

    for i in tqdm(range(n_nulls), desc=f"{model_name}:{null_type}"):
        rng = np.random.default_rng(seed_base + i * 1000 + hash(null_type) % 10000)
        null_compiled = generator(compiled, rng)

        v_mean, n_players, cycling_frac = sweep_compiled(
            null_compiled, output_indices,
            n_init=n_init, seed=seed_base + i,
        )

        passes, n_unique, max_frac = quality_gate(v_mean)

        delta = compute_delta_3plus(v_mean, n_players)
        if delta is None:
            n_zero_energy += 1
            continue

        null_local = extract_local_compiled(null_compiled, node_names)
        null_local_o3plus = sum(null_local[3:]) if len(null_local) > 3 else 0.0
        gap = delta - 100.0 * null_local_o3plus
        all_deltas.append(gap)

        if passes:
            gated_deltas.append(gap)
        else:
            n_degenerate += 1

        if (i + 1) % 5 == 0:
            out_dir = "/results/null_model"
            os.makedirs(out_dir, exist_ok=True)
            _save_interim(out_dir, model_name, null_type, all_deltas, gated_deltas,
                          n_degenerate, n_zero_energy, n_nulls, i + 1)
            results_vol.commit()

    out_dir = "/results/null_model"
    os.makedirs(out_dir, exist_ok=True)
    result = _save_interim(out_dir, model_name, null_type, all_deltas, gated_deltas,
                           n_degenerate, n_zero_energy, n_nulls, n_nulls)
    results_vol.commit()

    ts_end = datetime.now(timezone.utc).isoformat()
    print(f"[{ts_end}] Done: {model_name}:{null_type} — "
          f"mean={result.get('all_mean', 'N/A')}, "
          f"n_valid={len(all_deltas)}, degen={n_degenerate}")
    return result


def _save_interim(out_dir, model_name, null_type, all_deltas, gated_deltas,
                  n_degenerate, n_zero_energy, n_total_requested, n_completed):
    import json
    from datetime import datetime, timezone

    import numpy as np

    result = {
        "model": model_name,
        "null_type": null_type,
        "n_completed": n_completed,
        "n_requested": n_total_requested,
        "n_valid": len(all_deltas),
        "n_gated": len(gated_deltas),
        "n_degenerate": n_degenerate,
        "n_zero_energy": n_zero_energy,
        "degenerate_rate": n_degenerate / n_completed if n_completed > 0 else 0.0,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    for label, deltas in [("all", all_deltas), ("gated", gated_deltas)]:
        if deltas:
            arr = np.array(deltas)
            result[f"{label}_mean"] = float(np.mean(arr))
            result[f"{label}_std"] = float(np.std(arr))
            result[f"{label}_median"] = float(np.median(arr))
            result[f"{label}_creation_frac"] = float(np.mean(arr > 0.5))
            result[f"{label}_destruction_frac"] = float(np.mean(arr < -0.5))
            result[f"{label}_deltas"] = [float(d) for d in deltas]

    path = f"{out_dir}/{model_name}_{null_type}_modal.json"
    with open(path, "w") as f:
        json.dump(result, f, indent=2)
    return result


@app.local_entrypoint()
def main():
    import json

    eligible_path = os.path.join(REPO, "results/null_model/eligible_models.json")
    with open(eligible_path) as f:
        eligible_raw = json.load(f)

    merged_path = os.path.join(REPO, "results/grn_v2/merged_all_27_analysis.json")
    with open(merged_path) as f:
        merged = json.load(f)

    real_deltas = {m["model"]: m.get("delta_o3plus", 0) * 100
                   for m in merged["all_models"]}

    print(f"Null model sweep: {len(eligible_raw)} networks (n=7-15)")
    for name, info in sorted(eligible_raw.items(), key=lambda x: x[1]["n"]):
        n = info["n"]
        nn = n_nulls_for_size(n)
        real_d = real_deltas.get(name, 0)
        print(f"  {name}: n={n}, {nn} nulls, real Δ₃₊={real_d:+.1f}pp")

    futures = []
    for model_name, info in eligible_raw.items():
        n = info["n"]
        nn = n_nulls_for_size(n)
        for null_type in NULL_TYPES:
            seed_base = 42 + hash(model_name) % 100000
            futures.append(
                run_null_for_model.spawn(
                    model_name, info["rules"], info["output_nodes"],
                    null_type, nn, n_init=128, seed_base=seed_base,
                )
            )

    print(f"\nLaunched {len(futures)} jobs ({len(eligible_raw)} models × {len(NULL_TYPES)} null types)")

    results = []
    for future in futures:
        result = future.get()
        results.append(result)
        print(f"  Done: {result['model']}:{result['null_type']} — "
              f"valid={result['n_valid']}, degen_rate={result['degenerate_rate']:.0%}")

    print("\nAll jobs complete.")
