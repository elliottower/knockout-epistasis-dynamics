"""Generate publication-quality figures for composition_gap paper.

Reads data from results/grn_v2/ JSON files.
Saves PDF + PNG to paper/figures/.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

BASE = "/Users/elliottower/Documents/GitHub/epistasis-bench"
RESULTS = os.path.join(BASE, "results", "grn_v2")
FIGDIR = os.path.join(BASE, "paper", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# -- Style --
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
RED = "#b2182b"
GRAY = "#999999"

COLOR_MAP = {"creation": BLUE, "destruction": RED, "null": GRAY}

DISPLAY_NAMES = {
    "faure_cellcycle": "Faure cell cycle",
    "tournier_apoptosis": "Tournier apoptosis",
    "davidich_yeast": "Davidich yeast",
    "drosophila_cellcycle": "Drosophila cell cycle",
    "fanconi_anemia": "Fanconi anemia",
    "arabidopsis_cellcycle": "Arabidopsis cell cycle",
    "lambda_phage": "Lambda phage",
    "arellano_rootstem": "Arellano root stem",
    "asymmetric_cell_division": "Asymmetric cell div.",
    "cell_cycle_transcription": "Cell cycle transcription",
    "remy_p53_mdm2": "Remy p53-Mdm2",
    "albert_segment_polarity": "Albert segment polarity",
    "blood_stem_cell": "Blood stem cell",
    "calzone_cellfate_reduced": "Calzone cell fate (red.)",
    "li_budding_yeast": "Li budding yeast",
    "myeloid_progenitors": "Myeloid progenitors",
    "pair_rule_module": "Pair-rule module",
    "emt_switch": "EMT switch",
    "morphogenetic_checkpoint": "Morphogenetic checkpoint",
    "zanudo_tlgl": "Zanudo T-LGL",
    "lac_operon": "Lac operon",
    "mendoza_thelper": "Mendoza T-helper",
    "saadatpour_guardcell": "Saadatpour guard cell",
    "fumia_cellcycle": "Fumia cell cycle",
    "hematopoiesis_aging": "Hematopoiesis aging",
    "irons_cardiac": "Irons cardiac",
    "calzone_cell_fate": "Calzone cell fate",
}


def load_data():
    with open(os.path.join(RESULTS, "merged_all_27_analysis.json")) as f:
        merged = json.load(f)
    with open(os.path.join(RESULTS, "blind_experiment_summary.json")) as f:
        batch1_raw = json.load(f)
    with open(os.path.join(RESULTS, "blind_batch2_summary.json")) as f:
        batch2_raw = json.load(f)
    spectrum_lookup = {}
    for m in batch1_raw:
        spectrum_lookup[m["model"]] = {
            "global": m["global_spectrum"],
            "local": m["local_spectrum"],
        }
    for m in batch2_raw["results"]:
        spectrum_lookup[m["model"]] = {
            "global": m["global_spectrum"],
            "local": m["local_spectrum"],
        }
    return merged, spectrum_lookup


def fig1_ranked_bars(merged):
    models = sorted(merged["all_models"], key=lambda m: m["delta_o3plus"])
    names = [DISPLAY_NAMES.get(m["model"], m["model"]) for m in models]
    deltas_pp = [m["delta_o3plus"] * 100 for m in models]
    colors = [COLOR_MAP[m["label"]] for m in models]
    batches = [m["batch"] for m in models]

    fig, ax = plt.subplots(figsize=(5, 7))
    y = np.arange(len(models))
    bars = ax.barh(y, deltas_pp, height=0.7, color=colors, edgecolor="none", zorder=2)

    for i, (bar, batch) in enumerate(zip(bars, batches)):
        if batch == "batch1":
            bar.set_edgecolor("black")
            bar.set_linewidth(1.2)
        else:
            bar.set_edgecolor("#555555")
            bar.set_linewidth(0.3)

    ax.axvline(0, color="black", linewidth=0.5, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=7.5)
    ax.set_xlabel(r"$\Delta_{3+}$ (percentage points)")
    ax.set_xlim(-12, 90)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=BLUE, lw=6, label="Creation"),
        Line2D([0], [0], color=RED, lw=6, label="Destruction"),
        Line2D([0], [0], color=GRAY, lw=6, label="Null"),
        Line2D([0], [0], marker="s", color="w", markeredgecolor="black",
               markeredgewidth=1.2, markersize=7, label="Batch 1"),
        Line2D([0], [0], marker="s", color="w", markeredgecolor="#555555",
               markeredgewidth=0.3, markersize=7, label="Batch 2"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=7,
              frameon=True, framealpha=0.9, edgecolor="#cccccc")

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig1_delta_ranked.pdf"))
    fig.savefig(os.path.join(FIGDIR, "fig1_delta_ranked.png"))
    plt.close(fig)
    print("Saved fig1_delta_ranked")


def fig2_scatter(merged):
    models = merged["all_models"]
    fig, ax = plt.subplots(figsize=(5, 4))

    for m in models:
        x = m["cycling_fraction"] * 100
        y = m["delta_o3plus"] * 100
        c = COLOR_MAP[m["label"]]
        if m["batch"] == "batch1":
            ax.scatter(x, y, c=c, s=40, edgecolors="black", linewidths=0.8, zorder=3)
        else:
            ax.scatter(x, y, c="white", s=40, edgecolors=c, linewidths=1.2, zorder=3)

    label_models = {
        "emt_switch": (-15, 10),
        "blood_stem_cell": (-10, -12),
    }
    for m in models:
        if m["model"] in label_models:
            dx, dy = label_models[m["model"]]
            ax.annotate(
                DISPLAY_NAMES[m["model"]],
                (m["cycling_fraction"] * 100, m["delta_o3plus"] * 100),
                xytext=(dx, dy), textcoords="offset points",
                fontsize=7, fontstyle="italic",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#666666"),
            )

    ax.axhline(0, color="black", linewidth=0.4, linestyle="--", zorder=1)

    ax.text(0.97, 0.03,
            "Batch 2 only:\n" + r"$\rho = -0.37,\ p = 0.10$",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.5, color="#444444",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#cccccc", alpha=0.9))

    ax.set_xlabel("Cycling fraction (%)")
    ax.set_ylabel(r"$\Delta_{3+}$ (percentage points)")

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE,
               markeredgecolor="black", markeredgewidth=0.8, markersize=7,
               label="Creation (batch 1)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
               markeredgecolor=BLUE, markeredgewidth=1.2, markersize=7,
               label="Creation (batch 2)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=RED,
               markeredgecolor="black", markeredgewidth=0.8, markersize=7,
               label="Destruction (batch 1)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
               markeredgecolor=RED, markeredgewidth=1.2, markersize=7,
               label="Destruction (batch 2)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=6.5,
              frameon=True, framealpha=0.9, edgecolor="#cccccc")

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig2_delta_vs_cycling.pdf"))
    fig.savefig(os.path.join(FIGDIR, "fig2_delta_vs_cycling.png"))
    plt.close(fig)
    print("Saved fig2_delta_vs_cycling")


def fig3_spectra(spectrum_lookup, merged):
    panels = [
        ("cell_cycle_transcription", +83.25),
        ("fanconi_anemia", -5.48),
        ("davidich_yeast", -0.38),
    ]
    panel_labels = ["(a)", "(b)", "(c)"]
    panel_titles = ["Creation", "Destruction", "Null"]

    fig, axes = plt.subplots(1, 3, figsize=(7, 2.5), sharey=False)

    for ax, (model, delta_pp), plabel, ptitle in zip(axes, panels, panel_labels, panel_titles):
        spec = spectrum_lookup[model]
        g = np.array(spec["global"])
        l = np.array(spec["local"])

        n_orders = max(len(g), len(l))
        if len(g) < n_orders:
            g = np.pad(g, (0, n_orders - len(g)))
        if len(l) < n_orders:
            l = np.pad(l, (0, n_orders - len(l)))

        # Group into orders 0, 1, 2, 3, 4, 5+
        max_shown = 6
        g_grouped = list(g[:5]) + [g[5:].sum()]
        l_grouped = list(l[:5]) + [l[5:].sum()]

        x = np.arange(max_shown)
        w = 0.35

        ax.bar(x - w/2, l_grouped, w, color=BLUE, alpha=0.35, edgecolor=BLUE,
               linewidth=0.5, label="Local", hatch="///")
        ax.bar(x + w/2, g_grouped, w, color=BLUE, alpha=0.85, edgecolor=BLUE,
               linewidth=0.5, label="Global")

        ax.set_xticks(x)
        ax.set_xticklabels(["0", "1", "2", "3", "4", "5+"], fontsize=8)
        ax.set_xlabel("Interaction order", fontsize=8)
        if ax == axes[0]:
            ax.set_ylabel("Energy fraction", fontsize=8)

        display = DISPLAY_NAMES.get(model, model)
        ax.set_title(f"{plabel} {display}\n" + r"$\Delta_{3+}$" + f" = {delta_pp:+.1f}pp",
                     fontsize=8, pad=4)

    axes[0].legend(fontsize=7, loc="upper right", frameon=True,
                   framealpha=0.9, edgecolor="#cccccc")

    fig.tight_layout(w_pad=1.5)
    fig.savefig(os.path.join(FIGDIR, "fig3_spectrum_comparison.pdf"))
    fig.savefig(os.path.join(FIGDIR, "fig3_spectrum_comparison.png"))
    plt.close(fig)
    print("Saved fig3_spectrum_comparison")


if __name__ == "__main__":
    merged, spectrum_lookup = load_data()
    fig1_ranked_bars(merged)
    fig2_scatter(merged)
    fig3_spectra(spectrum_lookup, merged)
    print(f"\nAll figures saved to {FIGDIR}")
