"""Extension B: Async update dynamics sweep (pre-registered).

Runs batch-1 networks under asynchronous random-order update and compares
the composition gap with synchronous results.

Pre-reg predictions (prereg_extensions_v2.md, SHA 55cba8a):
  B1. Gap sign preserved in >= 4/6 networks
  B2. Median |delta_async| < median |delta_sync|

Protocol: 100 random update-order sequences per initial condition.

Usage:
    uv run python scripts/async_update_sweep.py                    # all 6
    uv run python scripts/async_update_sweep.py --small-only       # n<=12 only
"""
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from data_utils import normalized_wht, energy_by_order
from grn_coalition_sweep import BUILTIN_MODELS, extract_rule_fourier, sweep_coalitions

from scripts.run_batch2_blind_sweep import WEB_MODELS
from scripts.run_batch2b_extra_models import EXTRA_MODELS

ALL_MODELS = {**BUILTIN_MODELS, **WEB_MODELS, **EXTRA_MODELS}

BATCH1 = [
    "faure_cellcycle",
    "davidich_yeast",
    "tournier_apoptosis",
    "drosophila_cellcycle",
    "arabidopsis_cellcycle",
    "fanconi_anemia",
]

SMALL = [m for m in BATCH1 if len(ALL_MODELS[m]["rules"]) <= 12]
LARGE = [m for m in BATCH1 if len(ALL_MODELS[m]["rules"]) > 12]

ASYNC_REPLICATES = 100
N_INIT = 512


def run_one(model_name, out_dir):
    model_info = ALL_MODELS[model_name]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]
    n_nodes = len(rules)

    ts = datetime.now(timezone.utc).strftime("%H:%M:%S")
    print(f"[{ts}] {model_name} (n={n_nodes}, async_replicates={ASYNC_REPLICATES})")

    rf = extract_rule_fourier(rules)
    local_spectrum = rf.get("local_energy_spectrum", [])
    local_o3plus = sum(local_spectrum[3:]) if len(local_spectrum) > 3 else 0.0

    sweep = sweep_coalitions(
        rules, output_nodes,
        n_init=N_INIT, seed=42, max_steps=200,
        update_scheme="async",
        async_replicates=ASYNC_REPLICATES,
    )
    values = sweep["values"]
    v_mean = values.mean(axis=1)

    w = normalized_wht(v_mean)
    energy = energy_by_order(w, sweep["n_players"])
    total = energy.sum()
    if total > 0:
        spectrum = (energy / total).tolist()
    else:
        spectrum = energy.tolist()

    global_o3plus = sum(spectrum[3:]) if len(spectrum) > 3 else 0.0
    delta = global_o3plus - local_o3plus

    result = {
        "model": model_name,
        "n_nodes": n_nodes,
        "update_scheme": "async",
        "async_replicates": ASYNC_REPLICATES,
        "n_init": N_INIT,
        "local_o3plus": local_o3plus,
        "global_o3plus": global_o3plus,
        "delta_o3plus": delta,
        "gap_sign": "creation" if delta > 0.001 else ("destruction" if delta < -0.001 else "null"),
        "spectrum": spectrum,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = out_dir / f"{model_name}_async.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  delta_3+ = {delta:+.4f} ({result['gap_sign']}), saved to {out_path}")

    return result


def load_sync_deltas():
    """Load synchronous composition gap results for comparison."""
    sync = {}
    for model_name in BATCH1:
        f = Path(f"results/grn_v2/{model_name}_composition_blind.json")
        if not f.exists():
            continue
        with open(f) as fh:
            d = json.load(fh)
        es = d["energy_spectrum"]
        g = es["global"]
        l = es["local_rules"]
        global_o3plus = sum(g[3:])
        local_o3plus = sum(l[3:])
        sync[model_name] = {
            "global_o3plus": global_o3plus,
            "local_o3plus": local_o3plus,
            "delta_o3plus": global_o3plus - local_o3plus,
        }
    return sync


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-only", action="store_true",
                        help="Only run n<=12 networks (faure, davidich, tournier)")
    args = parser.parse_args()

    models = SMALL if args.small_only else BATCH1

    out_dir = Path("results/async_update")
    out_dir.mkdir(parents=True, exist_ok=True)

    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] Extension B: Async update sweep")
    print(f"Models: {models}")
    print(f"Async replicates: {ASYNC_REPLICATES}, N_init: {N_INIT}")
    print()

    sync_deltas = load_sync_deltas()

    async_results = {}
    for model_name in models:
        result = run_one(model_name, out_dir)
        async_results[model_name] = result

    print("\n\nComparison table (delta_3+ in pp):")
    print(f"{'Model':<25}  {'Sync':>8}  {'Async':>8}  {'Sign match':>10}")
    print("-" * 60)

    sign_matches = 0
    for model_name in models:
        a = async_results[model_name]
        s = sync_deltas.get(model_name, {})
        sync_d = s.get("delta_o3plus", float("nan")) * 100
        async_d = a["delta_o3plus"] * 100
        sync_sign = "+" if sync_d > 0.1 else ("-" if sync_d < -0.1 else "0")
        async_sign = "+" if async_d > 0.1 else ("-" if async_d < -0.1 else "0")
        match = sync_sign == async_sign
        if match:
            sign_matches += 1
        print(f"{model_name:<25}  {sync_d:>+7.2f}%  {async_d:>+7.2f}%  {'YES' if match else 'NO':>10}")

    n = len(models)
    print(f"\nSign preserved: {sign_matches}/{n}")
    print(f"Pre-reg prediction B1: >= 4/6 sign preserved → {'PASS' if sign_matches >= 4 else 'FAIL'}")

    sync_abs = [abs(sync_deltas[m]["delta_o3plus"]) for m in models if m in sync_deltas]
    async_abs = [abs(async_results[m]["delta_o3plus"]) for m in models]
    if sync_abs and async_abs:
        med_sync = np.median(sync_abs)
        med_async = np.median(async_abs)
        print(f"Median |delta| sync: {med_sync:.4f}, async: {med_async:.4f}")
        print(f"Pre-reg prediction B2: async < sync → {'PASS' if med_async < med_sync else 'FAIL'}")

    summary = {
        "models_run": models,
        "async_replicates": ASYNC_REPLICATES,
        "n_init": N_INIT,
        "sign_matches": sign_matches,
        "results": {m: async_results[m] for m in models},
        "sync_comparison": {m: sync_deltas.get(m) for m in models},
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    summary_path = out_dir / "async_sweep_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved summary to {summary_path}")


if __name__ == "__main__":
    main()
