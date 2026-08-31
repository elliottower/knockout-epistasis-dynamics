"""Run null model comparison on n<=11 networks locally.

Tests whether the composition gap in real Boolean networks differs from
randomly wired networks, using three nested nulls that isolate different
aspects of network structure.

For each real network (n<=11, 13 networks total):
  - Compute real Δ₃₊ from existing results
  - Generate 100 null networks per null type
  - Run coalition sweep on each null
  - Report: real Δ₃₊ vs null distribution (mean, std, percentile rank)

Usage:
    uv run python scripts/run_null_model_local.py
    uv run python scripts/run_null_model_local.py --n-nulls 30 --max-n 9
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from grn_coalition_sweep import BUILTIN_MODELS, compile_network, extract_rule_fourier
from scripts.run_batch2_blind_sweep import WEB_MODELS
from scripts.run_batch2b_extra_models import EXTRA_MODELS
from scripts.null_model_generators import (
    compute_delta_3plus,
    degree_preserving_rewire,
    kauffman_nk,
    quality_gate,
    rule_preserving_rewire,
    sweep_compiled,
)


RESULTS_DIR = Path(__file__).parent.parent / "results" / "null_model"


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def collect_all_models():
    """Merge all model sources into one dict."""
    all_models = {}
    for name, info in BUILTIN_MODELS.items():
        all_models[name] = info
    for name, info in WEB_MODELS.items():
        all_models[name] = info
    for name, info in EXTRA_MODELS.items():
        all_models[name] = info
    return all_models


def load_real_results():
    """Load real Δ₃₊ values from the merged analysis."""
    path = Path(__file__).parent.parent / "results" / "grn_v2" / "merged_all_27_analysis.json"
    with open(path) as f:
        data = json.load(f)
    return {m["model"]: m for m in data["all_models"]}


def run_one_model(model_name, rules, output_nodes, n_nulls, n_init, seed_base,
                   null_type_names=None):
    """Run null model comparison for a single real network."""
    compiled, node_names = compile_network(rules)
    n = len(node_names)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    local_fourier = extract_rule_fourier(rules)
    local_spectrum = local_fourier.get("local_energy_spectrum", [])
    local_o3plus = sum(local_spectrum[3:]) if len(local_spectrum) > 3 else 0.0

    all_null_types = {
        "rule_preserving": rule_preserving_rewire,
        "degree_preserving": degree_preserving_rewire,
        "kauffman_nk": kauffman_nk,
    }
    if null_type_names:
        null_types = {k: v for k, v in all_null_types.items() if k in null_type_names}
    else:
        null_types = all_null_types

    results = {}
    for null_name, generator in null_types.items():
        all_deltas = []
        gated_deltas = []
        n_degenerate = 0
        n_zero_energy = 0

        for i in tqdm(range(n_nulls), desc=f"  {null_name}", leave=False):
            rng = np.random.default_rng(seed_base + i * 1000 + hash(null_name) % 10000)
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

            null_local = extract_rule_fourier_compiled(null_compiled, node_names)
            null_local_o3plus = sum(null_local[3:]) if len(null_local) > 3 else 0.0
            gap = delta - 100.0 * null_local_o3plus
            all_deltas.append(gap)

            if passes:
                gated_deltas.append(gap)
            else:
                n_degenerate += 1

        results[null_name] = {
            "all_deltas": all_deltas,
            "gated_deltas": gated_deltas,
            "n_total": len(all_deltas),
            "n_gated": len(gated_deltas),
            "n_degenerate": n_degenerate,
            "n_zero_energy": n_zero_energy,
            "degenerate_rate": n_degenerate / n_nulls if n_nulls > 0 else 0.0,
        }

        for label, deltas in [("all", all_deltas), ("gated", gated_deltas)]:
            if deltas:
                arr = np.array(deltas)
                results[null_name][f"{label}_mean"] = float(np.mean(arr))
                results[null_name][f"{label}_std"] = float(np.std(arr))
                results[null_name][f"{label}_median"] = float(np.median(arr))
                results[null_name][f"{label}_creation_frac"] = float(np.mean(arr > 0.5))
                results[null_name][f"{label}_destruction_frac"] = float(np.mean(arr < -0.5))

    return results


def extract_rule_fourier_compiled(compiled, node_names):
    """Extract local energy spectrum directly from compiled truth tables."""
    from data_utils import wht

    n = len(compiled)
    max_k = max((len(regs) for _, regs, _ in compiled), default=0)
    local_energy = np.zeros(max_k + 1, dtype=np.float64)

    for node_name, reg_indices, truth_table in compiled:
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


def main():
    parser = argparse.ArgumentParser(description="Run null model comparison locally")
    parser.add_argument("--max-n", type=int, default=11,
                        help="Maximum network size to run (default: 11)")
    parser.add_argument("--n-nulls", type=int, default=100,
                        help="Number of null networks per type (default: 100)")
    parser.add_argument("--n-init", type=int, default=128,
                        help="Initial states per coalition (default: 128)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--models", nargs="*", default=None,
                        help="Specific model names to run (default: all n<=max_n)")
    parser.add_argument("--null-types", nargs="*",
                        default=["rule_preserving", "degree_preserving", "kauffman_nk"],
                        choices=["rule_preserving", "degree_preserving", "kauffman_nk"],
                        help="Which null types to run (default: all three)")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_models = collect_all_models()
    real_results = load_real_results()

    eligible = {}
    for name, info in all_models.items():
        n = len(info["rules"])
        if n <= args.max_n:
            if args.models is None or name in args.models:
                eligible[name] = info

    print(f"[{timestamp()}] Null model comparison")
    print(f"  Networks: {len(eligible)} (n <= {args.max_n})")
    print(f"  Nulls per type: {args.n_nulls}")
    print(f"  Initial states: {args.n_init}")
    print(f"  Models: {sorted(eligible.keys())}")
    print()

    all_output = {}
    for model_name in sorted(eligible.keys(), key=lambda k: len(eligible[k]["rules"])):
        info = eligible[model_name]
        n = len(info["rules"])
        real_data = real_results.get(model_name, {})
        real_delta_frac = real_data.get("delta_o3plus", 0.0)
        real_delta = real_delta_frac * 100.0  # convert fraction to pp

        print(f"[{timestamp()}] {model_name} (n={n}, real Δ₃₊={real_delta:+.1f}pp)")

        null_results = run_one_model(
            model_name, info["rules"], info["output_nodes"],
            n_nulls=args.n_nulls, n_init=args.n_init,
            seed_base=args.seed + hash(model_name) % 100000,
            null_type_names=args.null_types,
        )

        model_output = {
            "model": model_name,
            "n": n,
            "real_delta_o3plus": real_delta,
            "null_types": {},
        }

        for null_name, nr in null_results.items():
            entry = {
                "n_total": nr["n_total"],
                "n_gated": nr["n_gated"],
                "n_degenerate": nr["n_degenerate"],
                "n_zero_energy": nr["n_zero_energy"],
                "degenerate_rate": nr["degenerate_rate"],
            }

            for label in ["all", "gated"]:
                deltas = nr[f"{label}_deltas"]
                if deltas:
                    percentile = float(np.mean(np.array(deltas) <= real_delta)) * 100
                    entry[f"{label}_mean"] = nr[f"{label}_mean"]
                    entry[f"{label}_std"] = nr[f"{label}_std"]
                    entry[f"{label}_creation_frac"] = nr[f"{label}_creation_frac"]
                    entry[f"{label}_real_percentile"] = percentile

            model_output["null_types"][null_name] = entry

            if nr["n_total"] > 0:
                pctl = float(np.mean(np.array(nr["all_deltas"]) <= real_delta)) * 100
                print(f"    {null_name}: "
                      f"mean={nr['all_mean']:+.1f}pp, "
                      f"std={nr['all_std']:.1f}, "
                      f"real at {pctl:.0f}th pctl, "
                      f"creation={nr['all_creation_frac']:.0%}, "
                      f"degen={nr['degenerate_rate']:.0%} "
                      f"({nr['n_gated']}/{nr['n_total']} pass gate)")
            else:
                print(f"    {null_name}: all {args.n_nulls} nulls had zero energy")

        all_output[model_name] = model_output

        per_model_path = RESULTS_DIR / f"{model_name}_null_comparison.json"
        with open(per_model_path, "w") as f:
            json.dump(model_output, f, indent=2, default=float)
        print()

    summary_path = RESULTS_DIR / "null_model_summary.json"
    summary = {
        "timestamp": timestamp(),
        "max_n": args.max_n,
        "n_nulls": args.n_nulls,
        "n_init": args.n_init,
        "n_models": len(all_output),
        "models": all_output,
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=float)

    print(f"\n[{timestamp()}] Summary saved to {summary_path}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    n_real_creation = sum(1 for mo in all_output.values() if (mo["real_delta_o3plus"] or 0) > 0.5)
    print(f"Real networks: {n_real_creation}/{len(all_output)} creation")

    for null_type in ["rule_preserving", "degree_preserving", "kauffman_nk"]:
        creation_fracs = []
        degen_rates = []
        for model_name, mo in all_output.items():
            nt = mo["null_types"].get(null_type, {})
            if "all_creation_frac" in nt:
                creation_fracs.append(nt["all_creation_frac"])
            if "degenerate_rate" in nt:
                degen_rates.append(nt["degenerate_rate"])
        if creation_fracs:
            print(f"\n{null_type}:")
            print(f"  Mean null creation fraction: {np.mean(creation_fracs):.1%}")
            print(f"  Mean degeneracy rate: {np.mean(degen_rates):.1%}")


if __name__ == "__main__":
    main()
