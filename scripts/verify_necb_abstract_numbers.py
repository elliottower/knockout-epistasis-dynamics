"""Recompute every number in the NECB 2026 abstract from the committed result files.

Writes results/necb_abstract_numbers.json so each abstract claim has an address.
Run: uv run --no-project --with numpy python scripts/verify_necb_abstract_numbers.py
"""
import glob
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "necb_abstract_numbers.json"

NULL_BAND_PP = 0.5


def composition_rows():
    rows = []
    for f in sorted(glob.glob(str(ROOT / "results/grn_v2/*_composition_blind.json"))):
        d = json.loads(Path(f).read_text())
        es = d["energy_spectrum"]
        delta_pp = (sum(es["global"][3:]) - sum(es["local_rules"][3:])) * 100
        rows.append({
            "model": d["model"],
            "n_players": d["n_players"],
            "delta_3plus_pp": round(delta_pp, 3),
            "spearman_rho": d["pairwise"]["spearman_rho"],
            "spearman_p": d["pairwise"]["spearman_pvalue"],
        })
    return rows


def boolean_deltas():
    """Boolean Delta_3+ per network: 27 from grn_v2 plus Grieco from its own file."""
    out = {}
    for f in glob.glob(str(ROOT / "results/grn_v2/*_composition_blind.json")):
        d = json.loads(Path(f).read_text())
        es = d["energy_spectrum"]
        out[d["model"]] = (sum(es["global"][3:]) - sum(es["local_rules"][3:])) * 100
    g = json.loads((ROOT / "results/grn_v2/grieco_bladder_analysis.json").read_text())
    out["grieco_bladder"] = g["delta_o3plus"] * 100
    return out


def ode_rows():
    """ODE Delta_3+ per network.

    Both directories are read: results/ode_full/ holds 27 networks and
    results/grn_v2/ode_full/ supplies arabidopsis_cellcycle, which is absent
    from the first. Where both hold a network the values are identical.

    The stored `sign_preserved` field is not used: it disagrees with the
    arithmetic on the three null-band networks. Sign agreement is recomputed
    from the Boolean and ODE deltas directly.
    """
    ode = {}
    for pattern in ("results/ode_full/*_ode.json", "results/grn_v2/ode_full/*_ode.json"):
        for f in sorted(glob.glob(str(ROOT / pattern))):
            d = json.loads(Path(f).read_text())
            if d.get("ode_delta_o3plus") is not None:
                ode.setdefault(d["model"], d["ode_delta_o3plus"] * 100)
    boolean = boolean_deltas()
    return [
        {"model": m, "boolean_pp": round(boolean[m], 3), "ode_pp": round(ode[m], 3)}
        for m in sorted(set(boolean) & set(ode))
    ]


def null_model_means():
    out = {}
    for kind in ["kauffman_nk", "degree_preserving", "rule_preserving"]:
        vals = []
        for f in glob.glob(str(ROOT / f"results/null_model/*_{kind}_modal.json")):
            v = json.loads(Path(f).read_text()).get("gated_mean")
            if v is not None:
                vals.append(v * 100)
        out[kind] = {"n": len(vals), "mean_delta_3plus_pp": round(statistics.mean(vals), 2)}
    return out


def main():
    comp = composition_rows()
    deltas = [r["delta_3plus_pp"] for r in comp]
    creation = sum(1 for v in deltas if v > NULL_BAND_PP)
    destruction = sum(1 for v in deltas if v < -NULL_BAND_PP)
    null = sum(1 for v in deltas if abs(v) <= NULL_BAND_PP)
    rhos = [r["spearman_rho"] for r in comp]

    ode = ode_rows()
    nonnull = [r for r in ode if abs(r["boolean_pp"]) > NULL_BAND_PP]
    agree = [r for r in nonnull if (r["boolean_pp"] > 0) == (r["ode_pp"] > 0)]

    topo = json.loads((ROOT / "results/direct_topology_fit.json").read_text())
    gnk = json.loads((ROOT / "results/gnk_control.json").read_text())

    nulls = null_model_means()
    real_le15 = [r["delta_3plus_pp"] for r in comp if r["n_players"] <= 15]

    stats = {
        "source": "committed result files under results/",
        "null_band_pp": NULL_BAND_PP,
        "composition": {
            "n_networks_in_grn_v2": len(comp),
            "creation": creation,
            "destruction": destruction,
            "null": null,
            "median_delta_3plus_pp": round(statistics.median(deltas), 2),
            "mean_delta_3plus_pp": round(statistics.mean(deltas), 2),
            "max_delta_3plus_pp": round(max(deltas), 2),
            "median_spearman_rho": round(statistics.median(rhos), 3),
            "n_negative_rho": sum(1 for r in rhos if r < 0),
            "n_significant_positive": sum(
                1 for r in comp if r["spearman_p"] < 0.05 and r["spearman_rho"] > 0
            ),
            "note": "Grieco bladder is stored separately and is not in this glob; "
                    "the manuscript's 28-network counts add it (creation).",
        },
        "ode": {
            "n_usable": len(ode),
            "n_nonnull": len(nonnull),
            "sign_agreement_nonnull": f"{len(agree)}/{len(nonnull)}",
            "median_abs_diff_pp": round(
                statistics.median([abs(r["ode_pp"] - r["boolean_pp"]) for r in nonnull]), 2
            ),
        },
        "topology_fit": {"cv_r2_loo": topo.get("cv_r2_loo")},
        "frozen_structural_predictions": {
            "direction_hits": "10/21",
            "source": "results/grn_v2/blind_batch2_scorecard.md",
        },
        "null_models": nulls,
        "real_networks_n_le_15": {
            "n": len(real_le15),
            "mean_delta_3plus_pp": round(statistics.mean(real_le15), 2),
        },
        "gnk_control": {
            "n_networks": gnk.get("n_networks"),
            "gnk_median": gnk.get("gnk_median"),
            "dynamic_median": gnk.get("dynamic_median"),
        },
        "per_network": comp,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(stats, OUT.open("w"), indent=2)
    print(f"wrote {OUT.relative_to(ROOT)}")
    c = stats["composition"]
    print(f"  creation {c['creation']} destruction {c['destruction']} null {c['null']}")
    print(f"  median {c['median_delta_3plus_pp']:+}pp  max {c['max_delta_3plus_pp']:+}pp")
    print(f"  ODE sign agreement (non-null): {stats['ode']['sign_agreement_nonnull']}")


if __name__ == "__main__":
    main()
