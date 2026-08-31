"""D2: Threshold sensitivity analysis for composition gap classification.

Tests whether the creation/destruction/null classification is stable across
different null-threshold values.

Usage:
    uv run python scripts/threshold_sensitivity.py
"""
import json
import os
from pathlib import Path

import numpy as np


def classify(delta_pp, threshold_pp):
    if delta_pp > threshold_pp:
        return "creation"
    elif delta_pp < -threshold_pp:
        return "destruction"
    else:
        return "null"


def main():
    results_dir = Path("results/grn_v2")

    models = {}
    for f in sorted(results_dir.glob("*_composition_blind.json")):
        with open(f) as fh:
            d = json.load(fh)
        es = d["energy_spectrum"]
        g = sum(es["global"][3:])
        l = sum(es["local_rules"][3:])
        delta_pp = (g - l) * 100
        models[d["model"]] = delta_pp

    thresholds = [0.1, 0.25, 0.5, 1.0, 2.0]

    print(f"{'Model':<35}", end="")
    print(f"{'delta':>7}", end="")
    for t in thresholds:
        print(f"  {t:.2f}pp", end="")
    print()
    print("-" * (35 + 7 + len(thresholds) * 8))

    counts = {t: {"creation": 0, "destruction": 0, "null": 0} for t in thresholds}
    sensitive_models = []

    for model, delta in sorted(models.items()):
        row = f"{model:<35}{delta:>+7.2f}"
        classifications = {}
        for t in thresholds:
            c = classify(delta, t)
            classifications[t] = c
            symbol = {"creation": "+", "destruction": "-", "null": "0"}[c]
            row += f"  {symbol:>6}"
            counts[t][c] += 1

        # Check if classification changes across thresholds
        unique = set(classifications.values())
        if len(unique) > 1:
            sensitive_models.append(model)
            row += "  *"
        print(row)

    print()
    print("Classification counts:")
    print(f"{'':35}{'':>7}", end="")
    for t in thresholds:
        print(f"  {t:.2f}pp", end="")
    print()
    for label in ["creation", "destruction", "null"]:
        row = f"{label:<35}{'':>7}"
        for t in thresholds:
            row += f"  {counts[t][label]:>6}"
        print(row)

    print(f"\nThreshold-sensitive networks ({len(sensitive_models)}):")
    for m in sensitive_models:
        print(f"  {m}: delta = {models[m]:+.2f} pp")

    n_total = len(models)
    stable = n_total - len(sensitive_models)
    print(f"\nStable across all thresholds: {stable}/{n_total} ({100*stable/n_total:.0f}%)")

    out_dir = Path("results/sensitivity")
    out_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "thresholds_pp": thresholds,
        "models": {m: {"delta_pp": d, "classifications": {str(t): classify(d, t) for t in thresholds}} for m, d in models.items()},
        "counts": {str(t): counts[t] for t in thresholds},
        "sensitive_models": sensitive_models,
        "n_total": n_total,
        "n_stable": stable,
    }
    out_path = out_dir / "threshold_sensitivity.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
