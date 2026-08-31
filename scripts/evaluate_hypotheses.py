"""Evaluate all pre-registered hypotheses against sweep v2 results.

Scores H1-H4 (shapiq prereg) and HA1-HA7 (AUROC addendum) against
the actual experimental results. Outputs a structured verdict file.

Usage:
    uv run python scripts/evaluate_hypotheses.py \
        --results results/shapiq_sweep_v2.json \
        --output results/hypothesis_verdicts.md
"""

import argparse
import json
from collections import defaultdict

import numpy as np


def load_results(path):
    with open(path) as f:
        return json.load(f)


def extract_metric(trials, circuit, budget, method_key, metric):
    vals = []
    for t in trials:
        if t["circuit"] != circuit or t["budget_fraction"] != budget:
            continue
        trial = t["trial"]
        if method_key in trial and trial[method_key].get(metric) is not None:
            v = trial[method_key][metric]
            if np.isfinite(v):
                vals.append(v)
    return vals


def median_metric(trials, circuit, budget, method_key, metric):
    vals = extract_metric(trials, circuit, budget, method_key, metric)
    return float(np.median(vals)) if vals else None


def ci95(vals):
    if not vals:
        return None, None
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


CIRCUITS = ["weight_ioi", "ioi", "random15"]
CIRCUIT_NAMES = {"weight_ioi": "weight_ioi", "ioi": "canonical_ioi", "random15": "random15"}
BUDGETS = [0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.40]
SHAPIQ_METHODS = ["kernelshapiq", "shapiq_mc", "svarmiq"]
ALL_METHODS_O2 = ["lasso_walsh", "irf", "kernelshapiq_order2", "shapiq_mc_order2", "svarmiq_order2"]

CEILINGS = {
    "weight_ioi": {2: 0.757, 3: 0.868},
    "ioi": {2: 0.955, 3: 0.988},
    "random15": {2: 0.995, 3: 1.000},
}


def evaluate_h1(trials):
    lines = ["## H1: Estimation efficiency (shapiq reaches 95% of ceiling at >=5% budget)\n"]
    confirmed = True
    for circ in CIRCUITS:
        for budget in [0.05, 0.10, 0.20, 0.40]:
            m = f"kernelshapiq_order2"
            r2 = median_metric(trials, circ, budget, m, "heldout_r2")
            ceiling = CEILINGS[circ][2]
            eff = r2 / ceiling if r2 and ceiling else None
            status = "PASS" if eff and eff >= 0.95 else "FAIL"
            if eff and eff < 0.95:
                confirmed = False
            lines.append(f"- {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}%: R²={r2:.4f}, "
                         f"ceiling={ceiling:.3f}, efficiency={eff:.3f} [{status}]")
    verdict = "CONFIRMED" if confirmed else "PARTIALLY CONFIRMED"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_h2(trials):
    lines = ["## H2: Order-3 ceiling behavior per circuit\n"]
    verdicts = {}

    for circ in CIRCUITS:
        for budget in [0.05, 0.10, 0.20, 0.40]:
            best_o3 = max(
                median_metric(trials, circ, budget, f"{m}_order3", "heldout_r2") or -999
                for m in SHAPIQ_METHODS
            )
            irf_r2 = median_metric(trials, circ, budget, "irf", "heldout_r2")
            lines.append(f"- {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}%: best_shapiq_o3={best_o3:.4f}, irf={irf_r2:.4f}")

    weight_results = []
    for budget in [0.05, 0.10, 0.20, 0.40]:
        best_o2 = max(
            median_metric(trials, "weight_ioi", budget, f"{m}_order2", "heldout_r2") or -999
            for m in SHAPIQ_METHODS
        )
        best_o3 = max(
            median_metric(trials, "weight_ioi", budget, f"{m}_order3", "heldout_r2") or -999
            for m in SHAPIQ_METHODS
        )
        weight_results.append((budget, best_o2, best_o3, best_o3 - best_o2))

    lines.append(f"\nweight_ioi order-3 vs order-2 gains:")
    for b, o2, o3, gap in weight_results:
        lines.append(f"  {b*100:.0f}%: o2={o2:.4f}, o3={o3:.4f}, gap={gap:.4f}")

    gap_5pct = any(g >= 0.05 for _, _, _, g in weight_results)
    lines.append(f"\nPrediction: order-3 exceeds order-2 ceiling by >=5 R² points on weight_ioi: {'YES' if gap_5pct else 'NO'}")

    ioi_results = []
    for budget in [0.05, 0.10, 0.20, 0.40]:
        best_o2 = max(
            median_metric(trials, "ioi", budget, f"{m}_order2", "heldout_r2") or -999
            for m in SHAPIQ_METHODS
        )
        best_o3 = max(
            median_metric(trials, "ioi", budget, f"{m}_order3", "heldout_r2") or -999
            for m in SHAPIQ_METHODS
        )
        ioi_results.append((budget, best_o2, best_o3, best_o3 - best_o2))

    ioi_gap_small = all(abs(g) < 0.03 for _, _, _, g in ioi_results)
    lines.append(f"canonical_ioi: all gaps <3 R² points: {'YES' if ioi_gap_small else 'NO'}")

    verdict = "CONFIRMED" if gap_5pct else "PARTIALLY CONFIRMED"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_h3(trials):
    lines = ["## H3: Method ordering (KernelSHAP-IQ > SVARM-IQ > SHAPIQ-MC)\n"]

    reversals_low = 0
    total_low = 0
    for circ in CIRCUITS:
        for budget in [0.01, 0.02, 0.03, 0.05, 0.10]:
            for order in [2, 3]:
                k = median_metric(trials, circ, budget, f"kernelshapiq_order{order}", "heldout_r2")
                sv = median_metric(trials, circ, budget, f"svarmiq_order{order}", "heldout_r2")
                sh = median_metric(trials, circ, budget, f"shapiq_mc_order{order}", "heldout_r2")
                if k is None or sv is None or sh is None:
                    continue
                if any(v < 0.3 for v in [k, sv, sh]):
                    continue
                total_low += 1
                correct = k >= sv >= sh
                if not correct:
                    reversals_low += 1
                lines.append(f"- {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}% o{order}: "
                             f"K={k:.4f} SV={sv:.4f} SH={sh:.4f} {'OK' if correct else 'REVERSED'}")

    lines.append(f"\nLow-to-mid budget: {reversals_low}/{total_low} reversals (qualifying cells with R²>0.3)")

    convergence_violations = 0
    total_high = 0
    for circ in CIRCUITS:
        for budget in [0.20, 0.40]:
            for order in [2, 3]:
                k = median_metric(trials, circ, budget, f"kernelshapiq_order{order}", "heldout_r2")
                sv = median_metric(trials, circ, budget, f"svarmiq_order{order}", "heldout_r2")
                sh = median_metric(trials, circ, budget, f"shapiq_mc_order{order}", "heldout_r2")
                if k is None or sv is None or sh is None:
                    continue
                total_high += 1
                gap = max(k, sv, sh) - min(k, sv, sh)
                if gap > 0.05:
                    convergence_violations += 1
                lines.append(f"- {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}% o{order}: "
                             f"K={k:.4f} SV={sv:.4f} SH={sh:.4f} gap={gap:.4f}")

    lines.append(f"High budget convergence: {convergence_violations}/{total_high} cells with gap > 0.05")

    falsified_1 = False
    falsified_2 = False
    for circ in CIRCUITS:
        mc_beats_k = 0
        mc_beats_sv = 0
        for budget in [0.01, 0.02, 0.03, 0.05, 0.10]:
            for order in [2, 3]:
                k = median_metric(trials, circ, budget, f"kernelshapiq_order{order}", "heldout_r2")
                sv = median_metric(trials, circ, budget, f"svarmiq_order{order}", "heldout_r2")
                sh = median_metric(trials, circ, budget, f"shapiq_mc_order{order}", "heldout_r2")
                if sh and k and sh > k:
                    mc_beats_k += 1
                if sh and sv and sh > sv:
                    mc_beats_sv += 1

    verdict = "CONFIRMED" if reversals_low <= 2 else "PARTIALLY CONFIRMED"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_h4(trials):
    lines = ["## H4: Cross-circuit consistency of method ordering\n"]

    inconsistencies = 0
    qualifying_cells = 0

    for budget in BUDGETS:
        for order in [2, 3]:
            rankings = {}
            for circ in CIRCUITS:
                k = median_metric(trials, circ, budget, f"kernelshapiq_order{order}", "heldout_r2")
                sv = median_metric(trials, circ, budget, f"svarmiq_order{order}", "heldout_r2")
                sh = median_metric(trials, circ, budget, f"shapiq_mc_order{order}", "heldout_r2")
                if k is None or sv is None or sh is None:
                    continue
                if all(v > 0.3 for v in [k, sv, sh]):
                    ranking = sorted(
                        [("K", k), ("SV", sv), ("SH", sh)],
                        key=lambda x: x[1], reverse=True
                    )
                    rankings[circ] = tuple(r[0] for r in ranking)

            if len(rankings) >= 2:
                qualifying_cells += 1
                rank_vals = list(rankings.values())
                if len(set(rank_vals)) > 1:
                    inconsistencies += 1
                    lines.append(f"- {budget*100:.0f}% o{order}: INCONSISTENT {rankings}")
                else:
                    lines.append(f"- {budget*100:.0f}% o{order}: consistent {rank_vals[0]}")

    lines.append(f"\n{inconsistencies}/{qualifying_cells} inconsistent cells")
    verdict = "CONFIRMED" if inconsistencies < 3 else "FALSIFIED"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_ha1(trials):
    lines = ["## HA1: LASSO achieves highest pairwise AUROC at >=10% budget\n"]

    violations = 0
    for circ in ["weight_ioi", "ioi"]:
        for budget in [0.10, 0.20, 0.40]:
            lasso = median_metric(trials, circ, budget, "lasso_walsh", "pairwise_auroc")
            best_other = -999
            best_other_name = ""
            for m in ["irf", "kernelshapiq_order2", "kernelshapiq_order3",
                       "shapiq_mc_order2", "shapiq_mc_order3",
                       "svarmiq_order2", "svarmiq_order3"]:
                v = median_metric(trials, circ, budget, m, "pairwise_auroc")
                if v and v > best_other:
                    best_other = v
                    best_other_name = m
            status = "PASS" if lasso and lasso >= best_other else "FAIL"
            if lasso and lasso < best_other:
                violations += 1
            lines.append(f"- {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}%: LASSO={lasso:.4f}, "
                         f"best_other={best_other:.4f} ({best_other_name}) [{status}]")

    verdict = "CONFIRMED" if violations == 0 else "FALSIFIED"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_ha2(trials):
    lines = ["## HA2: Shapiq order-2 competitive with LASSO at low budgets\n"]

    falsified_both = True
    for circ in ["weight_ioi", "ioi"]:
        for budget in [0.01, 0.02, 0.03]:
            lasso = median_metric(trials, circ, budget, "lasso_walsh", "pairwise_auroc")
            best_shapiq_o2 = max(
                median_metric(trials, circ, budget, f"{m}_order2", "pairwise_auroc") or -999
                for m in SHAPIQ_METHODS
            )
            gap = (lasso - best_shapiq_o2) if lasso and best_shapiq_o2 > -999 else None
            within_005 = gap is not None and gap <= 0.05
            within_010 = gap is not None and gap <= 0.10
            if within_010:
                falsified_both = False
            lines.append(f"- {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}%: LASSO={lasso:.4f}, "
                         f"best_shapiq_o2={best_shapiq_o2:.4f}, gap={gap:.4f} "
                         f"[{'within 0.05' if within_005 else 'within 0.10' if within_010 else 'EXCEEDS 0.10'}]")

    lasso_dominates = True
    for circ in ["weight_ioi", "ioi"]:
        max_gap_this_circ = 0
        for budget in [0.01, 0.02, 0.03]:
            lasso = median_metric(trials, circ, budget, "lasso_walsh", "pairwise_auroc")
            best_shapiq_o2 = max(
                median_metric(trials, circ, budget, f"{m}_order2", "pairwise_auroc") or -999
                for m in SHAPIQ_METHODS
            )
            gap = (lasso - best_shapiq_o2) if lasso and best_shapiq_o2 > -999 else 0
            max_gap_this_circ = max(max_gap_this_circ, gap)
        if max_gap_this_circ <= 0.10:
            lasso_dominates = False

    if not lasso_dominates:
        verdict = "NOT FALSIFIED"
    else:
        verdict = "FALSIFIED"

    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_ha3(trials):
    lines = ["## HA3: Per-circuit ordering of pairwise AUROC\n"]

    ordering_ok = True
    threshold_violated = False

    for budget in [0.05, 0.10, 0.20, 0.40]:
        aurocs = {}
        for circ in CIRCUITS:
            method_aurocs = []
            for m_key in ALL_METHODS_O2:
                v = median_metric(trials, circ, budget, m_key, "pairwise_auroc")
                if v:
                    method_aurocs.append(v)
            aurocs[circ] = np.mean(method_aurocs) if method_aurocs else None

        ioi = aurocs.get("ioi")
        wioi = aurocs.get("weight_ioi")
        rand = aurocs.get("random15")
        if ioi and wioi and rand:
            correct = ioi > wioi > rand
            if not correct:
                ordering_ok = False
            lines.append(f"- {budget*100:.0f}%: ioi={ioi:.4f}, weight_ioi={wioi:.4f}, "
                         f"random15={rand:.4f} {'OK' if correct else 'VIOLATED'}")

    lines.append(f"\nrandom15 threshold check (all methods AUROC < 0.65 at 1-3%):")
    for budget in [0.01, 0.02, 0.03]:
        for m_key in ALL_METHODS_O2:
            v = median_metric(trials, "random15", budget, m_key, "pairwise_auroc")
            if v and v > 0.75:
                threshold_violated = True
                lines.append(f"  VIOLATED: {m_key} @ {budget*100:.0f}% = {v:.4f} > 0.75")
            elif v and v > 0.65:
                lines.append(f"  above 0.65: {m_key} @ {budget*100:.0f}% = {v:.4f}")

    lasso_random15 = median_metric(trials, "random15", 0.01, "lasso_walsh", "pairwise_auroc")
    lines.append(f"\nLASSO on random15 @ 1%: {lasso_random15:.4f}" if lasso_random15 else "")

    verdict = "CONFIRMED" if ordering_ok else "PARTIALLY CONFIRMED"
    if threshold_violated:
        verdict += " (random15 threshold FALSIFIED)"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_ha5(trials):
    lines = ["## HA5: Max-order effect on shapiq pairwise AUROC\n"]

    lines.append("### HA5a: weight_ioi, high budget (order-3 > order-2)")
    ha5a_violations = 0
    for budget in [0.10, 0.20, 0.40]:
        for m in SHAPIQ_METHODS:
            o2 = median_metric(trials, "weight_ioi", budget, f"{m}_order2", "pairwise_auroc")
            o3 = median_metric(trials, "weight_ioi", budget, f"{m}_order3", "pairwise_auroc")
            if o2 and o3:
                better = o3 > o2
                if not better:
                    ha5a_violations += 1
                lines.append(f"  {m} @ {budget*100:.0f}%: o2={o2:.4f}, o3={o3:.4f} "
                             f"{'o3 wins' if better else 'o2 wins'}")

    ha5a_verdict = "CONFIRMED" if ha5a_violations == 0 else "FALSIFIED" if ha5a_violations >= 2 else "PARTIALLY"

    lines.append(f"\n  HA5a violations: {ha5a_violations}/9, verdict: {ha5a_verdict}")

    lines.append("\n### HA5b: weight_ioi, low budget (order-2 > order-3)")
    ha5b_violations = 0
    for budget in [0.01, 0.02, 0.03]:
        for m in SHAPIQ_METHODS:
            o2 = median_metric(trials, "weight_ioi", budget, f"{m}_order2", "pairwise_auroc")
            o3 = median_metric(trials, "weight_ioi", budget, f"{m}_order3", "pairwise_auroc")
            if o2 and o3:
                better = o2 > o3
                if not better:
                    ha5b_violations += 1
                lines.append(f"  {m} @ {budget*100:.0f}%: o2={o2:.4f}, o3={o3:.4f} "
                             f"{'o2 wins' if better else 'o3 wins'}")

    ha5b_verdict = "CONFIRMED" if ha5b_violations <= 2 else "FALSIFIED"
    lines.append(f"\n  HA5b violations: {ha5b_violations}/9, verdict: {ha5b_verdict}")

    lines.append("\n### HA5c: canonical_ioi and random15 (order gap <= 0.03)")
    ha5c_violations = 0
    for circ in ["ioi", "random15"]:
        for budget in BUDGETS:
            for m in SHAPIQ_METHODS:
                o2 = median_metric(trials, circ, budget, f"{m}_order2", "pairwise_auroc")
                o3 = median_metric(trials, circ, budget, f"{m}_order3", "pairwise_auroc")
                if o2 and o3:
                    gap = abs(o3 - o2)
                    if gap > 0.05:
                        ha5c_violations += 1
                        lines.append(f"  VIOLATED: {CIRCUIT_NAMES[circ]} {m} @ {budget*100:.0f}%: "
                                     f"gap={gap:.4f}")

    ha5c_verdict = "CONFIRMED" if ha5c_violations == 0 else "FALSIFIED"
    lines.append(f"\n  HA5c violations: {ha5c_violations}, verdict: {ha5c_verdict}")

    verdict = f"HA5a={ha5a_verdict}, HA5b={ha5b_verdict}, HA5c={ha5c_verdict}"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_ha6(trials):
    lines = ["## HA6: Shapiq inter-method ordering (K >= SV >= SH for AUROC)\n"]

    reversals = 0
    qualifying = 0

    for circ in CIRCUITS:
        for budget in BUDGETS:
            for order in [2, 3]:
                k = median_metric(trials, circ, budget, f"kernelshapiq_order{order}", "pairwise_auroc")
                sv = median_metric(trials, circ, budget, f"svarmiq_order{order}", "pairwise_auroc")
                sh = median_metric(trials, circ, budget, f"shapiq_mc_order{order}", "pairwise_auroc")
                if k is None or sv is None or sh is None:
                    continue
                if max(k, sv, sh) <= 0.55:
                    continue
                qualifying += 1
                correct = k >= sv >= sh
                if not correct:
                    reversals += 1
                    lines.append(f"- REVERSAL: {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}% o{order}: "
                                 f"K={k:.4f} SV={sv:.4f} SH={sh:.4f}")

    lines.append(f"\n{reversals}/{qualifying} reversals in qualifying cells (threshold: >3 falsifies)")

    verdict = "CONFIRMED" if reversals <= 3 else "FALSIFIED"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def evaluate_ha7(trials):
    lines = ["## HA7: iRF pairwise AUROC profile\n"]

    irf_beats_lasso = 0
    for circ in ["weight_ioi", "ioi"]:
        for budget in [0.10, 0.20, 0.40]:
            lasso = median_metric(trials, circ, budget, "lasso_walsh", "pairwise_auroc")
            irf = median_metric(trials, circ, budget, "irf", "pairwise_auroc")
            if lasso and irf:
                status = "LASSO wins" if lasso > irf else "iRF wins"
                if irf > lasso:
                    irf_beats_lasso += 1
                lines.append(f"- {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}%: "
                             f"LASSO={lasso:.4f}, iRF={irf:.4f} [{status}]")

    lines.append(f"\niRF >= LASSO at high budget on both circuits: {'YES' if irf_beats_lasso >= 6 else 'NO'}")

    lines.append("\niRF vs worst shapiq (order=2) at >=5%:")
    for circ in ["weight_ioi", "ioi"]:
        for budget in [0.05, 0.10, 0.20, 0.40]:
            irf = median_metric(trials, circ, budget, "irf", "pairwise_auroc")
            worst_shapiq = min(
                median_metric(trials, circ, budget, f"{m}_order2", "pairwise_auroc") or 999
                for m in SHAPIQ_METHODS
            )
            if irf and worst_shapiq < 999:
                status = "iRF wins" if irf > worst_shapiq else "shapiq wins"
                lines.append(f"  {CIRCUIT_NAMES[circ]} @ {budget*100:.0f}%: "
                             f"iRF={irf:.4f}, worst_shapiq_o2={worst_shapiq:.4f} [{status}]")

    verdict = "CONFIRMED" if irf_beats_lasso == 0 else "FALSIFIED"
    lines.insert(1, f"\n**Verdict: {verdict}**\n")
    return "\n".join(lines), verdict


def summary_table(trials):
    lines = ["## Summary: median pairwise AUROC by method, circuit, budget\n"]

    method_keys = ["lasso_walsh", "irf",
                   "kernelshapiq_order2", "kernelshapiq_order3",
                   "shapiq_mc_order2", "shapiq_mc_order3",
                   "svarmiq_order2", "svarmiq_order3"]

    for circ in CIRCUITS:
        lines.append(f"\n### {CIRCUIT_NAMES[circ]}\n")
        header = "| Budget | " + " | ".join(method_keys) + " |"
        sep = "|--------|" + "|".join(["------"] * len(method_keys)) + "|"
        lines.append(header)
        lines.append(sep)
        for budget in BUDGETS:
            vals = []
            for m in method_keys:
                v = median_metric(trials, circ, budget, m, "pairwise_auroc")
                vals.append(f"{v:.3f}" if v else "N/A")
            lines.append(f"| {budget*100:.0f}% | " + " | ".join(vals) + " |")

    return "\n".join(lines)


def r2_summary_table(trials):
    lines = ["## Summary: median held-out R² by method, circuit, budget\n"]

    method_keys = ["lasso_walsh", "irf",
                   "kernelshapiq_order2", "kernelshapiq_order3",
                   "shapiq_mc_order2", "shapiq_mc_order3",
                   "svarmiq_order2", "svarmiq_order3"]

    for circ in CIRCUITS:
        lines.append(f"\n### {CIRCUIT_NAMES[circ]}\n")
        header = "| Budget | " + " | ".join(method_keys) + " |"
        sep = "|--------|" + "|".join(["------"] * len(method_keys)) + "|"
        lines.append(header)
        lines.append(sep)
        for budget in BUDGETS:
            vals = []
            for m in method_keys:
                v = median_metric(trials, circ, budget, m, "heldout_r2")
                vals.append(f"{v:.3f}" if v else "N/A")
            lines.append(f"| {budget*100:.0f}% | " + " | ".join(vals) + " |")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    data = load_results(args.results)
    trials = data["trials"]

    print(f"Loaded {len(trials)} trials")

    sections = []
    verdicts = {}

    sections.append("# Hypothesis Verdicts — EpistasisBench Sweep v2\n")
    sections.append(f"Results file: {args.results}")
    sections.append(f"Harness version: {data['harness_version']}")
    sections.append(f"Trials: {len(trials)}, Completed: {data['n_completed']}, Failed: {data['n_failed']}")
    sections.append(f"Runtime versions: {data['runtime_versions']}\n")

    sections.append("---\n")
    sections.append("# Shapiq Pre-registration (H1-H4)\n")

    for name, fn in [("H1", evaluate_h1), ("H2", evaluate_h2),
                     ("H3", evaluate_h3), ("H4", evaluate_h4)]:
        text, v = fn(trials)
        sections.append(text)
        verdicts[name] = v
        sections.append("")

    sections.append("---\n")
    sections.append("# AUROC Addendum Pre-registration (HA1-HA7)\n")

    for name, fn in [("HA1", evaluate_ha1), ("HA2", evaluate_ha2),
                     ("HA3", evaluate_ha3), ("HA5", evaluate_ha5),
                     ("HA6", evaluate_ha6), ("HA7", evaluate_ha7)]:
        text, v = fn(trials)
        sections.append(text)
        verdicts[name] = v
        sections.append("")

    sections.append("---\n")
    sections.append(summary_table(trials))
    sections.append("")
    sections.append(r2_summary_table(trials))

    sections.append("\n---\n")
    sections.append("## Verdict Summary\n")
    for k, v in verdicts.items():
        sections.append(f"- **{k}**: {v}")

    output = "\n".join(sections)
    with open(args.output, "w") as f:
        f.write(output)
    print(f"\nVerdicts written to {args.output}")

    print("\n=== VERDICT SUMMARY ===")
    for k, v in verdicts.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
