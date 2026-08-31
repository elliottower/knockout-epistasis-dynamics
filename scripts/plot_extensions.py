"""Generate figures for composition gap paper extensions (ODE + empirical).

Saves PDF + PNG to paper/figures/.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = "/Users/elliottower/Documents/GitHub/epistasis-bench"
FIGDIR = os.path.join(BASE, "paper", "figures")
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

BLUE = "#2166ac"
ORANGE = "#e66101"
GRAY = "#999999"
TEAL = "#1b7837"
PURPLE = "#762a83"


def fig_ode_comparison():
    """Side-by-side bars: Boolean vs ODE delta for 4 pilot networks (broken y-axis)."""
    pilot_files = [
        ("lambda_phage", "Lambda phage\n(n=7)"),
        ("arellano_rootstem", "Arellano root stem\n(n=9)"),
        ("faure_cellcycle", "Faure cell cycle\n(n=10)"),
        ("davidich_yeast", "Davidich yeast\n(n=10)"),
    ]

    bool_deltas = []
    ode_deltas = []
    labels = []
    for fname, label in pilot_files:
        path = os.path.join(BASE, "results", "grn_v2", "ode_pilot", f"{fname}_ode.json")
        with open(path) as f:
            d = json.load(f)
        bool_deltas.append(d["boolean_delta_o3plus"] * 100)
        ode_deltas.append(d["ode_delta_o3plus"] * 100)
        labels.append(label)

    x = np.arange(len(labels))
    width = 0.35

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(5.5, 3.8), sharex=True,
        gridspec_kw={"height_ratios": [1, 1.2], "hspace": 0.08},
    )

    for ax in (ax_top, ax_bot):
        ax.bar(x - width / 2, bool_deltas, width, label="Boolean",
               color=BLUE, edgecolor="white", linewidth=0.5)
        ax.bar(x + width / 2, ode_deltas, width, label="Hill ODE",
               color=ORANGE, edgecolor="white", linewidth=0.5)
        ax.axhline(0, color="black", linewidth=0.4, zorder=0)

    ax_top.set_ylim(42, 55)
    ax_bot.set_ylim(-5, 5)

    ax_bot.axhline(0.5, color=GRAY, linewidth=0.3, linestyle="--", zorder=0)
    ax_bot.axhline(-0.5, color=GRAY, linewidth=0.3, linestyle="--", zorder=0)

    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    ax_top.tick_params(bottom=False)

    d_mark = 0.015
    kwargs = dict(transform=ax_top.transAxes, color="black",
                  clip_on=False, linewidth=0.6)
    ax_top.plot((-d_mark, d_mark), (-d_mark, d_mark), **kwargs)
    ax_top.plot((1 - d_mark, 1 + d_mark), (-d_mark, d_mark), **kwargs)
    kwargs["transform"] = ax_bot.transAxes
    ax_bot.plot((-d_mark, d_mark), (1 - d_mark, 1 + d_mark), **kwargs)
    ax_bot.plot((1 - d_mark, 1 + d_mark), (1 - d_mark, 1 + d_mark), **kwargs)

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(labels, fontsize=8)
    fig.text(0.02, 0.5, r"$\Delta_{3+}$ (pp)", va="center",
             rotation="vertical", fontsize=9)
    ax_top.legend(frameon=False, fontsize=8, loc="upper left")

    ax_bot.text(3.4, 1.0, "null\nthreshold", fontsize=6.5, color=GRAY,
                ha="right", va="bottom")

    for fmt in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"fig_ode_comparison.{fmt}"))
    plt.close(fig)
    print(f"  Saved fig_ode_comparison.pdf/png")


def fig_empirical_spectra():
    """Energy-by-order for Weinreich TEM-1 and Hall yeast (3 panels)."""
    with open(os.path.join(BASE, "results", "empirical", "weinreich_2006_walsh.json")) as f:
        wein = json.load(f)
    with open(os.path.join(BASE, "results", "empirical", "hall_2010_walsh.json")) as f:
        hall = json.load(f)

    panels = [
        ("Weinreich TEM-1\nlog(MIC), n=5",
         wein["log_mic"]["spectrum"],
         BLUE),
        ("Hall yeast\nhaploid growth, n=6",
         hall["results_by_component"]["haploid_growth_rate"]["spectrum"],
         TEAL),
        ("Hall yeast\nmating efficiency, n=6",
         hall["results_by_component"]["mating_efficiency"]["spectrum"],
         PURPLE),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.8), sharey=False)

    for ax, (title, spectrum, color) in zip(axes, panels):
        orders = list(range(len(spectrum)))
        bars = ax.bar(orders, [s * 100 for s in spectrum], color=color,
                      edgecolor="white", linewidth=0.5, width=0.7)

        o3plus = sum(spectrum[3:]) * 100
        ax.set_xlabel("Walsh order")
        ax.set_title(title, fontsize=8.5, pad=6)
        ax.set_xticks(orders)

        for i, bar in enumerate(bars):
            if i >= 3 and spectrum[i] * 100 > 0.3:
                ax.annotate(f"{spectrum[i]*100:.1f}%",
                            xy=(bar.get_x() + bar.get_width() / 2,
                                bar.get_height()),
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=6,
                            color=color)

        if o3plus > 0.05:
            ax.text(0.97, 0.95, f"$E_{{3+}}$ = {o3plus:.1f}%",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=7.5, color=color,
                    bbox=dict(boxstyle="round,pad=0.2",
                              facecolor="white", edgecolor=color,
                              alpha=0.8, linewidth=0.5))

    axes[0].set_ylabel("Energy fraction (%)")

    # Truncate order-0 bar for panels with dominant intercept
    for idx in [1, 2]:
        ax = axes[idx]
        spectrum = panels[idx][1]
        if spectrum[0] > 0.5:
            max_non0 = max(spectrum[1:]) * 100
            cap = max(max_non0 * 2.5, 15)
            ax.set_ylim(0, cap)
            bar0 = ax.patches[0]
            bar0.set_height(cap * 0.95)
            ax.annotate(f"{spectrum[0]*100:.0f}%",
                        xy=(bar0.get_x() + bar0.get_width() / 2,
                            cap * 0.92),
                        ha="center", va="top", fontsize=6.5,
                        color="white", fontweight="bold")
            ax.plot([bar0.get_x(), bar0.get_x() + bar0.get_width()],
                    [cap * 0.97, cap * 0.97], color="white",
                    linewidth=1.5, zorder=5)
            ax.plot([bar0.get_x(), bar0.get_x() + bar0.get_width()],
                    [cap * 0.99, cap * 0.99], color="white",
                    linewidth=1.5, zorder=5)

    fig.tight_layout(w_pad=1.5)

    for fmt in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"fig_empirical_spectra.{fmt}"))
    plt.close(fig)
    print(f"  Saved fig_empirical_spectra.pdf/png")


if __name__ == "__main__":
    print("Generating extension figures...")
    fig_ode_comparison()
    fig_empirical_spectra()
    print("Done.")
