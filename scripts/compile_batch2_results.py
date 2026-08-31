"""Compile all batch 2 + 2b blind experiment results into a combined summary and scorecard.

Reads:
  - blind_predictions_batch2.json (frozen before sweeps)
  - blind_predictions_batch2b.json (frozen before sweeps)
  - *_composition_blind.json (actual results for each model)

Outputs:
  - blind_batch2_summary.json (combined summary)
  - blind_batch2_scorecard.md (predictions vs actuals)
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

RESULTS_DIR = Path(__file__).parent.parent / "results" / "grn_v2"


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def main():
    # Load all predictions
    all_predictions = {}
    pred_shas = {}

    for pred_file in ["blind_predictions_batch2.json", "blind_predictions_batch2b.json"]:
        pred_path = RESULTS_DIR / pred_file
        if pred_path.exists():
            with open(pred_path) as f:
                data = json.load(f)
            pred_shas[pred_file] = data.get("sha256", "unknown")
            for model_name, pred_data in data["models"].items():
                all_predictions[model_name] = pred_data
            print(f"Loaded {len(data['models'])} predictions from {pred_file}")
            print(f"  SHA-256: {pred_shas[pred_file]}")

    # The 6 batch-1 models (already swept, skip from scorecard)
    batch1_models = {
        "faure_cellcycle", "tournier_apoptosis", "davidich_yeast",
        "drosophila_cellcycle", "fanconi_anemia", "arabidopsis_cellcycle",
    }

    # Load all composition results
    all_results = []
    for comp_file in sorted(RESULTS_DIR.glob("*_composition_blind.json")):
        model_name = comp_file.stem.replace("_composition_blind", "")
        if model_name in batch1_models:
            continue  # Skip batch 1 models (not in our predictions)

        with open(comp_file) as f:
            comp_data = json.load(f)

        pw = comp_data["pairwise"]
        spec = comp_data["energy_spectrum"]
        tri = comp_data.get("triples", {})
        n = comp_data["n_players"]
        g3p = sum(spec["global"][3:])
        l3p = sum(spec["local_rules"][3:])
        conv_path = RESULTS_DIR / f"{model_name}_wiring_blind.json"
        cycling_frac = 0.0
        if conv_path.exists():
            with open(conv_path) as f:
                wiring = json.load(f)
            cycling_frac = wiring.get("convergence", {}).get("cycling_fraction", 0.0)

        result = {
            "model": model_name,
            "n": n,
            "spearman_rho": pw["spearman_rho"],
            "spearman_ci": pw["spearman_ci_95"],
            "spearman_p": pw["spearman_pvalue"],
            "global_o3plus": g3p,
            "local_o3plus": l3p,
            "delta_o3plus": g3p - l3p,
            "creation_or_destruction": "creation" if g3p > l3p else "destruction",
            "cycling_fraction": cycling_frac,
            "triple_rho": tri.get("spearman_rho"),
            "global_spectrum": spec["global"],
            "local_spectrum": spec["local_rules"],
        }
        all_results.append(result)

    print(f"\nLoaded {len(all_results)} composition results")

    # Save combined summary
    summary = {
        "experiment": "batch2_combined_blind_sweep",
        "timestamp": timestamp(),
        "prediction_shas": pred_shas,
        "n_models": len(all_results),
        "parameters": {
            "n_init": 512,
            "max_steps": 200,
            "seed": 42,
            "n_workers": 10,
            "clamp_value": 0,
        },
        "results": sorted(all_results, key=lambda x: x["n"]),
    }

    summary_path = RESULTS_DIR / "blind_batch2_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Summary saved: {summary_path}")

    # Build scorecard
    scorecard_lines = [
        "# Blind Batch 2 Combined Scorecard",
        "",
        "## Prediction Provenance",
        "",
    ]
    for fname, sha in pred_shas.items():
        scorecard_lines.append(f"- `{fname}`: SHA-256 `{sha[:16]}...`")
    scorecard_lines.extend([
        "",
        "## Predictions vs Actuals",
        "",
        "| Model | n | Pred Dir | Actual Dir | Pred rho Range | Actual rho | Pred o3+ | Actual o3+ | rho Hit | Dir Hit |",
        "|-------|---|---------|-----------|---------------|------------|----------|-----------|---------|---------|",
    ])

    n_dir_hits = 0
    n_rho_hits = 0
    n_scored = 0
    rho_errors = []

    for result in sorted(all_results, key=lambda x: x["n"]):
        model_name = result["model"]
        if model_name not in all_predictions:
            scorecard_lines.append(
                f"| {model_name} | {result['n']} | ? | {result['creation_or_destruction']} | ? | {result['spearman_rho']:.3f} | ? | {result['global_o3plus']:.4f} | ? | ? |"
            )
            continue

        n_scored += 1
        pred = all_predictions[model_name]["prediction"]
        actual_dir = result["creation_or_destruction"]
        pred_dir = pred["creation_or_destruction"]
        dir_hit = actual_dir == pred_dir
        if dir_hit:
            n_dir_hits += 1

        actual_rho = result["spearman_rho"]
        rho_lo, rho_hi = pred["spearman_rho_range"]
        rho_hit = rho_lo <= actual_rho <= rho_hi
        if rho_hit:
            n_rho_hits += 1

        # Track rho prediction error
        rho_mid = (rho_lo + rho_hi) / 2
        rho_errors.append(abs(actual_rho - rho_mid))

        actual_o3 = result["global_o3plus"]
        pred_o3 = pred["global_o3plus_estimate"]

        scorecard_lines.append(
            f"| {model_name} | {result['n']} | {pred_dir} | {actual_dir} | "
            f"[{rho_lo:.2f}, {rho_hi:.2f}] | {actual_rho:.3f} | "
            f"{pred_o3:.4f} | {actual_o3:.4f} | "
            f"{'YES' if rho_hit else 'NO'} | {'YES' if dir_hit else 'NO'} |"
        )

    import numpy as np
    mean_rho_error = float(np.mean(rho_errors)) if rho_errors else 0.0

    scorecard_lines.extend([
        "",
        "## Summary Statistics",
        "",
        f"- **Models scored**: {n_scored}",
        f"- **Direction hits**: {n_dir_hits}/{n_scored} ({n_dir_hits/max(1,n_scored):.0%})",
        f"- **Rho in range**: {n_rho_hits}/{n_scored} ({n_rho_hits/max(1,n_scored):.0%})",
        f"- **Mean rho prediction error**: {mean_rho_error:.3f}",
        "",
        "## All Results (sorted by n)",
        "",
        "| Model | n | rho | 95% CI | p-value | Global o3+ | Local o3+ | Delta | Direction | Cycling |",
        "|-------|---|-----|--------|---------|-----------|----------|-------|-----------|---------|",
    ])

    for result in sorted(all_results, key=lambda x: x["n"]):
        ci = result["spearman_ci"]
        scorecard_lines.append(
            f"| {result['model']} | {result['n']} | {result['spearman_rho']:.3f} | "
            f"[{ci[0]:.3f}, {ci[1]:.3f}] | {result['spearman_p']:.2e} | "
            f"{result['global_o3plus']:.4f} | {result['local_o3plus']:.4f} | "
            f"{result['delta_o3plus']:+.4f} | {result['creation_or_destruction']} | "
            f"{result['cycling_fraction']:.1%} |"
        )

    scorecard_path = RESULTS_DIR / "blind_batch2_scorecard.md"
    with open(scorecard_path, "w") as f:
        f.write("\n".join(scorecard_lines) + "\n")
    print(f"Scorecard saved: {scorecard_path}")

    # Print summary
    print(f"\n{'='*70}")
    print(f"COMBINED RESULTS ({n_scored} models scored)")
    print(f"{'='*70}")
    print(f"  Direction accuracy: {n_dir_hits}/{n_scored} ({n_dir_hits/max(1,n_scored):.0%})")
    print(f"  Rho in range:       {n_rho_hits}/{n_scored} ({n_rho_hits/max(1,n_scored):.0%})")
    print(f"  Mean rho error:     {mean_rho_error:.3f}")
    print()
    for r in sorted(all_results, key=lambda x: x["n"]):
        print(f"  {r['model']:35s} n={r['n']:2d}  rho={r['spearman_rho']:+.3f}  "
              f"o3+={r['global_o3plus']:.4f}  {r['creation_or_destruction']:12s}  "
              f"cycling={r['cycling_fraction']:.1%}")


if __name__ == "__main__":
    main()
