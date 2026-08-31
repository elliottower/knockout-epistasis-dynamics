"""Plot Walsh spectra for all empirical fitness landscapes including Franke.

Generates a multi-panel figure showing the Walsh energy spectrum for each
empirical landscape analyzed.

Usage:
    uv run python scripts/plot_empirical_spectra_v2.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_spectra():
    results_dir = Path("results/empirical")

    # Weinreich TEM-1
    with open(results_dir / "weinreich_2006_walsh.json") as f:
        w = json.load(f)
    weinreich_log = {
        "label": r"Weinreich TEM-1 ($\log$ MIC)",
        "organism": "E. coli",
        "n": 5,
        "spectrum": w["log_mic"]["spectrum"],
    }
    weinreich_raw = {
        "label": "Weinreich TEM-1 (raw MIC)",
        "organism": "E. coli",
        "n": 5,
        "spectrum": w["raw_mic"]["spectrum"],
    }

    # Hall 2010
    with open(results_dir / "hall_2010_walsh.json") as f:
        h = json.load(f)
    hall_growth = {
        "label": "Hall haploid growth",
        "organism": "S. cerevisiae",
        "n": 6,
        "spectrum": h["results_by_component"]["haploid_growth_rate"]["spectrum"],
    }
    hall_mating = {
        "label": "Hall mating efficiency",
        "organism": "S. cerevisiae",
        "n": 6,
        "spectrum": h["results_by_component"]["mating_efficiency"]["spectrum"],
    }

    # Franke 2011
    with open(results_dir / "franke2011_spectrum.json") as f:
        fk = json.load(f)
    franke = {
        "label": "Franke growth rate",
        "organism": "A. niger",
        "n": 8,
        "spectrum": fk["spectrum"],
        "note": "186/256 measured",
    }

    return [weinreich_log, weinreich_raw, hall_growth, hall_mating, franke]


def main():
    datasets = load_spectra()
    n_panels = len(datasets)

    fig, axes = plt.subplots(1, n_panels, figsize=(3.0 * n_panels, 3.5),
                             sharey=False)

    colors = ["#2196F3", "#4CAF50", "#FF9800", "#E91E63", "#9C27B0"]

    for i, (ax, ds) in enumerate(zip(axes, datasets)):
        spec = ds["spectrum"]
        n = ds["n"]
        orders = list(range(len(spec)))

        # Skip order 0 (grand mean) for cleaner visualization
        orders_plot = orders[1:]
        spec_plot = spec[1:]
        total_epistatic = sum(spec_plot)

        bars = ax.bar(orders_plot, [s * 100 for s in spec_plot],
                      color=colors[i], alpha=0.85, edgecolor="white", linewidth=0.5)

        o3plus = sum(spec[3:]) * 100
        ax.axhline(y=0, color="black", linewidth=0.5)

        ax.set_xlabel("Walsh order")
        ax.set_title(f"{ds['label']}\n({ds['organism']}, n={n})",
                     fontsize=9, fontweight="bold")
        ax.set_xlim(0.3, len(spec) - 0.3)

        ax.text(0.95, 0.95, f"order 3+:\n{o3plus:.1f}%",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="gray", alpha=0.9))

    axes[0].set_ylabel("Walsh energy (% of total)")

    fig.suptitle("Walsh spectra of empirical fitness landscapes",
                 fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()

    out_dir = Path("paper/figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    for fmt in ["pdf", "png"]:
        out_path = out_dir / f"fig_empirical_spectra_v2.{fmt}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved to paper/figures/fig_empirical_spectra_v2.pdf/png")
    plt.close()


if __name__ == "__main__":
    main()
