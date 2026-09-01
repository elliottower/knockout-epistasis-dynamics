"""Reconcile the numbers claimed in composition_gap_v18.tex against the committed results.

Writes results/paper_number_reconciliation.json. Every value the manuscript states in
its headline claims is recomputed here from the result files and compared.

Run: uv run --no-project python scripts/reconcile_paper_numbers.py
"""
import glob
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "paper_number_reconciliation.json"
NULL_BAND_PP = 0.5


def boolean_deltas():
    """Delta_3+ per network, in percentage points. 27 from grn_v2 plus Grieco."""
    out = {}
    for f in glob.glob(str(ROOT / "results/grn_v2/*_composition_blind.json")):
        d = json.loads(Path(f).read_text())
        es = d["energy_spectrum"]
        out[d["model"]] = ((sum(es["global"][3:]) - sum(es["local_rules"][3:])) * 100,
                           d["n_players"])
    g = json.loads((ROOT / "results/grn_v2/grieco_bladder_analysis.json").read_text())
    out["grieco_bladder"] = (g["delta_o3plus"] * 100, g.get("n_players", 18))
    return out


def ode_deltas():
    """ODE Delta_3+ in pp. Both directories: grn_v2/ode_full supplies arabidopsis."""
    out = {}
    for pattern in ("results/ode_full/*_ode.json", "results/grn_v2/ode_full/*_ode.json"):
        for f in sorted(glob.glob(str(ROOT / pattern))):
            d = json.loads(Path(f).read_text())
            if d.get("ode_delta_o3plus") is not None:
                out.setdefault(d["model"], d["ode_delta_o3plus"] * 100)
    return out


def null_means():
    out = {}
    for kind in ["kauffman_nk", "degree_preserving", "rule_preserving"]:
        vals = [json.loads(Path(f).read_text()).get("gated_mean")
                for f in glob.glob(str(ROOT / f"results/null_model/*_{kind}_modal.json"))]
        vals = [v for v in vals if v is not None]
        out[kind] = {"n": len(vals), "mean_pp": round(statistics.mean(vals), 3)}
    return out



def secondary_results():
    """Numbers the manuscript reports outside the composition-gap headline."""
    out = {}

    topo = json.loads((ROOT / "results/direct_topology_fit.json").read_text())
    out["topology_cv_r2_loo"] = topo.get("cv_r2_loo")

    gnk = json.loads((ROOT / "results/gnk_control.json").read_text())
    out["gnk"] = {k: gnk.get(k) for k in
                  ("n_networks", "gnk_median", "dynamic_median", "gap", "verdicts")}

    scorecard = (ROOT / "results/grn_v2/blind_batch2_scorecard.md").read_text()
    hits = [ln.strip() for ln in scorecard.splitlines() if "Direction hits" in ln]
    out["blind_direction_hits"] = hits[0] if hits else None

    # local-global pairwise correlation across the composition panel
    rhos, sig, neg = [], 0, 0
    for f in glob.glob(str(ROOT / "results/grn_v2/*_composition_blind.json")):
        d = json.loads(Path(f).read_text())["pairwise"]
        rhos.append(d["spearman_rho"])
        if d["spearman_pvalue"] < 0.05 and d["spearman_rho"] > 0:
            sig += 1
        if d["spearman_rho"] < 0:
            neg += 1
    g = json.loads((ROOT / "results/grn_v2/grieco_bladder_analysis.json").read_text())
    rhos.append(g["spearman_rho"])
    if g.get("spearman_p", 1) < 0.05 and g["spearman_rho"] > 0:
        sig += 1
    if g["spearman_rho"] < 0:
        neg += 1
    out["local_global_rho"] = {
        "n": len(rhos),
        "median": round(statistics.median(rhos), 3),
        "mean": round(statistics.mean(rhos), 3),
        "n_significant_positive": sig,
        "n_negative": neg,
    }

    emp = ROOT / "results/empirical/consolidated_comparison.json"
    if emp.exists():
        out["empirical"] = json.loads(emp.read_text())
    return out


def main():
    bo = boolean_deltas()
    deltas = [v for v, _ in bo.values()]
    creation = sum(1 for v in deltas if v > NULL_BAND_PP)
    destruction = sum(1 for v in deltas if v < -NULL_BAND_PP)
    null = sum(1 for v in deltas if abs(v) <= NULL_BAND_PP)
    sub = [v for v, n in bo.values() if n <= 15]

    ode = ode_deltas()
    common = sorted(set(bo) & set(ode))
    nonnull = [m for m in common if abs(bo[m][0]) > NULL_BAND_PP]
    agree_all = [m for m in common if (bo[m][0] > 0) == (ode[m] > 0)]
    agree_nn = [m for m in nonnull if (bo[m][0] > 0) == (ode[m] > 0)]
    flips = [(m, round(bo[m][0], 2), round(ode[m], 2))
             for m in common if m not in agree_all]

    nulls = null_means()

    checks = [
        ("creation / destruction, non-null", "19 / 6", f"{creation} / {destruction}",
         creation == 19 and destruction == 6),
        ("null-band networks", "3", str(null), null == 3),
        ("median Delta_3+ (pp)", "+7.4", f"{statistics.median(deltas):+.2f}",
         abs(statistics.median(deltas) - 7.4) < 0.05),
        ("mean Delta_3+ (pp)", "+19", f"{statistics.mean(deltas):+.2f}",
         abs(statistics.mean(deltas) - 19) < 0.5),
        ("real mean, n<=15 (pp)", "+18.6", f"{statistics.mean(sub):+.2f}",
         abs(statistics.mean(sub) - 18.6) < 0.05),
        ("Kauffman NK null (pp)", "+3.5", f"{nulls['kauffman_nk']['mean_pp']:+.2f}",
         abs(nulls["kauffman_nk"]["mean_pp"] - 3.5) < 0.05),
        ("degree-preserving null (pp)", "+3.8",
         f"{nulls['degree_preserving']['mean_pp']:+.2f}",
         abs(nulls["degree_preserving"]["mean_pp"] - 3.8) < 0.05),
        ("rule-preserving null (pp)", "+10.3",
         f"{nulls['rule_preserving']['mean_pp']:+.2f}",
         abs(nulls["rule_preserving"]["mean_pp"] - 10.3) < 0.05),
        ("ODE sign preserved, non-null", "25/25",
         f"{len(agree_nn)}/{len(nonnull)}",
         len(agree_nn) == len(nonnull) == 25),
        ("ODE raw sign agreement, all 28", "(null band excluded)",
         f"{len(agree_all)}/{len(common)}", True),
    ]

    stats = {
        "manuscript": "paper/composition_gap_v18.tex",
        "null_band_pp": NULL_BAND_PP,
        "checks": [
            {"quantity": q, "manuscript": m, "computed": c, "agrees": bool(a)}
            for q, m, c, a in checks
        ],
        "ode_sign_flips": flips,
        "null_models": nulls,
        "n_boolean_networks": len(bo),
        "n_with_ode": len(common),
        "secondary": secondary_results(),
        "per_network_delta_pp": {m: round(v, 3) for m, (v, _) in sorted(bo.items())},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(stats, OUT.open("w"), indent=2)

    print(f"wrote {OUT.relative_to(ROOT)}\n")
    print(f"{'quantity':<36}{'manuscript':>14}{'computed':>14}   ")
    for q, m, c, a in checks:
        print(f"{q:<36}{m:>14}{c:>14}   {'ok' if a else 'DIFFERS'}")
    if flips:
        print(f"\nODE sign flips (Boolean pp -> ODE pp): {flips}")


if __name__ == "__main__":
    main()
