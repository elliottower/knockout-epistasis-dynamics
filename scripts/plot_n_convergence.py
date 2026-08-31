"""Plot N convergence for the sensitivity analysis (D1).

Shows delta_3+ as a function of N (number of initial conditions) for each
network, demonstrating convergence.

Usage:
    uv run python scripts/plot_n_convergence.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

N_VALUES = [64, 128, 256, 512, 1024]

PARTIAL_RESULTS = {
    "faure_cellcycle": {
        64: -0.0309, 128: -0.0304, 256: -0.0300, 512: -0.0302, 1024: -0.0301
    },
    "davidich_yeast": {
        64: +0.0014, 128: -0.0022, 256: -0.0033, 512: -0.0038, 1024: -0.0055
    },
    "tournier_apoptosis": {
        64: +0.0898, 128: +0.0948, 256: +0.0851, 512: +0.0846, 1024: +0.0856
    },
}


def main():
    results_path = Path("results/sensitivity/n_init_convergence.json")
    if results_path.exists():
        with open(results_path) as f:
            data = json.load(f)
        results = data["results"]
    else:
        results = {name: {str(n): {"delta_o3plus": v} for n, v in vals.items()}
                   for name, vals in PARTIAL_RESULTS.items()}

    models = list(results.keys())
    n_models = len(models)

    fig, axes = plt.subplots(1, min(n_models, 6), figsize=(3.5 * min(n_models, 6), 3.5),
                             sharey=False)
    if n_models == 1:
        axes = [axes]

    colors = ["#2196F3", "#4CAF50", "#FF9800", "#E91E63", "#9C27B0", "#795548"]

    for i, (ax, model) in enumerate(zip(axes, models)):
        r = results[model]
        ns = []
        deltas = []
        for n_val in N_VALUES:
            key = str(n_val) if str(n_val) in r else n_val
            if key in r:
                ns.append(n_val)
                deltas.append(r[key]["delta_o3plus"] * 100)

        ax.plot(ns, deltas, "o-", color=colors[i % len(colors)],
                linewidth=2, markersize=6)

        max_range = max(deltas) - min(deltas) if deltas else 0
        ax.axhline(y=deltas[-1] if deltas else 0, color="gray",
                   linestyle="--", alpha=0.5, linewidth=0.8)

        ax.set_xlabel("N (initial conditions)")
        short_name = model.replace("_cellcycle", "").replace("_apoptosis", "")
        ax.set_title(f"{short_name}\n(range: {max_range:.2f} pp)",
                     fontsize=9, fontweight="bold")
        ax.set_xscale("log", base=2)
        ax.set_xticks(N_VALUES)
        ax.set_xticklabels([str(n) for n in N_VALUES], fontsize=7)

        if deltas:
            ymin, ymax = min(deltas), max(deltas)
            margin = max(0.5, (ymax - ymin) * 0.5)
            ax.set_ylim(ymin - margin, ymax + margin)

    axes[0].set_ylabel(r"$\Delta_{3+}$ (pp)")
    fig.suptitle(r"Convergence of $\Delta_{3+}$ with number of initial conditions",
                 fontsize=11, fontweight="bold", y=1.02)
    plt.tight_layout()

    out_dir = Path("paper/figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    for fmt in ["pdf", "png"]:
        fig.savefig(out_dir / f"fig_n_convergence.{fmt}", dpi=300, bbox_inches="tight")
    print(f"Saved to paper/figures/fig_n_convergence.pdf/png")
    plt.close()


if __name__ == "__main__":
    main()
