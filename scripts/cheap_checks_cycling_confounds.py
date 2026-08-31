"""Two cheap pre-finalization checks:
1. Is G4 Davidich's Δo3+ = -0.4pp meaningfully nonzero vs the energy floor?
2. Is cycling fraction confounded with n or mean_regs?

Uses only data from blind_experiment_summary.json + BUILTIN_MODELS rules.
"""
import json
import re
import sys
from scipy import stats
import numpy as np

sys.path.insert(0, "/Users/elliottower/Documents/GitHub/epistasis-bench")
from grn_coalition_sweep import BUILTIN_MODELS

# ── Load results ─────────────────────────────────────────────────────
with open("/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2/blind_experiment_summary.json") as f:
    results = json.load(f)

model_lookup = {r["model"]: r for r in results}

# ── Check 1: G4 delta vs energy floor ────────────────────────────────
print("=" * 60)
print("CHECK 1: G4 Davidich delta_o3+ vs energy floor")
print("=" * 60)
g4 = model_lookup["davidich_yeast"]
print(f"  global_o3+ = {g4['global_o3plus']:.6f} ({g4['global_o3plus']*100:.2f}%)")
print(f"  local_o3+  = {g4['local_o3plus']:.6f} ({g4['local_o3plus']*100:.2f}%)")
print(f"  delta_o3+  = {g4['delta_o3plus']:.6f} ({g4['delta_o3plus']*100:.2f} pp)")
print(f"  abs(delta) = {abs(g4['delta_o3plus']):.6f}")
print(f"  Quality gate threshold: abs_e3+ >= 1e-4")
print(f"  Global abs_e3+ = {g4['global_o3plus']:.6f} — clears gate: {g4['global_o3plus'] >= 1e-4}")
print(f"  Local abs_e3+  = {g4['local_o3plus']:.6f} — clears gate: {g4['local_o3plus'] >= 1e-4}")
print()

# Compare delta magnitude to other models
print("  All deltas for context:")
for r in sorted(results, key=lambda x: x["delta_o3plus"]):
    print(f"    {r['model']:25s}  delta = {r['delta_o3plus']*100:+7.2f} pp")
print()

# G4's delta is 10x smaller than the next-smallest nonzero delta
deltas = sorted([abs(r["delta_o3plus"]) for r in results])
print(f"  Sorted |delta|: {[f'{d*100:.2f}pp' for d in deltas]}")
print(f"  G4's |delta| = {abs(g4['delta_o3plus'])*100:.2f}pp is {deltas[1]/deltas[0]:.1f}x smaller than next")
print(f"  G4 Spearman p = {g4['spearman_p']:.3f} (not significant)")
print()
print("  VERDICT: G4 is borderline null — Δ = -0.4pp on a model where")
print("  pairwise rho is also non-significant (p=0.31). Could be 3/2/1.")

# ── Check 2: Cycling confounds ──────────────────────────────────────
print()
print("=" * 60)
print("CHECK 2: Cycling fraction vs n and mean regulators")
print("=" * 60)

# Compute mean regulators for each model
def count_regulators(rules):
    """Count unique gene names referenced in each rule's expression."""
    all_nodes = set(rules.keys())
    reg_counts = []
    for node, expr in rules.items():
        regs_found = set()
        for other_node in all_nodes:
            if re.search(r'\b' + re.escape(other_node) + r'\b', expr):
                if other_node != node:
                    regs_found.add(other_node)
        reg_counts.append(len(regs_found))
    return reg_counts


model_ids = ["faure_cellcycle", "tournier_apoptosis", "davidich_yeast",
             "drosophila_cellcycle", "fanconi_anemia", "arabidopsis_cellcycle"]

data = []
for mid in model_ids:
    r = model_lookup[mid]
    rules = BUILTIN_MODELS[mid]["rules"]
    reg_counts = count_regulators(rules)
    mean_regs = np.mean(reg_counts)
    n_edges = sum(reg_counts)
    density = n_edges / r["n"]
    data.append({
        "model": mid,
        "n": r["n"],
        "cycling": r["cycling_fraction"],
        "delta_o3plus": r["delta_o3plus"],
        "mean_regs": mean_regs,
        "n_edges": n_edges,
        "density": density,
    })

print(f"\n  {'Model':25s} {'n':>3s} {'cycling':>8s} {'mean_regs':>10s} {'n_edges':>8s} {'density':>8s} {'delta':>8s}")
print("  " + "-" * 72)
for d in sorted(data, key=lambda x: x["cycling"]):
    print(f"  {d['model']:25s} {d['n']:3d} {d['cycling']*100:7.1f}% {d['mean_regs']:10.2f} {d['n_edges']:8d} {d['density']:8.2f} {d['delta_o3plus']*100:+7.1f}pp")

ns = [d["n"] for d in data]
cycling = [d["cycling"] for d in data]
mean_regs = [d["mean_regs"] for d in data]
deltas = [d["delta_o3plus"] for d in data]
densities = [d["density"] for d in data]

print()
rho_cn, p_cn = stats.spearmanr(cycling, ns)
print(f"  Spearman(cycling, n)         = {rho_cn:+.3f}  (p = {p_cn:.3f})")

rho_cr, p_cr = stats.spearmanr(cycling, mean_regs)
print(f"  Spearman(cycling, mean_regs) = {rho_cr:+.3f}  (p = {p_cr:.3f})")

rho_cd, p_cd = stats.spearmanr(cycling, densities)
print(f"  Spearman(cycling, density)   = {rho_cd:+.3f}  (p = {p_cd:.3f})")

# Also check: do n or mean_regs themselves predict delta?
rho_dn, p_dn = stats.spearmanr(deltas, ns)
print(f"  Spearman(delta, n)           = {rho_dn:+.3f}  (p = {p_dn:.3f})")

rho_dr, p_dr = stats.spearmanr(deltas, mean_regs)
print(f"  Spearman(delta, mean_regs)   = {rho_dr:+.3f}  (p = {p_dr:.3f})")

rho_dd, p_dd = stats.spearmanr(deltas, densities)
print(f"  Spearman(delta, density)     = {rho_dd:+.3f}  (p = {p_dd:.3f})")

# cycling vs delta
rho_dcy, p_dcy = stats.spearmanr(deltas, cycling)
print(f"  Spearman(delta, cycling)     = {rho_dcy:+.3f}  (p = {p_dcy:.3f})")

print()
if abs(rho_cn) > 0.6:
    print("  WARNING: cycling is strongly correlated with n.")
    print("  Cannot attribute composition gap direction to cycling specifically.")
else:
    print(f"  Cycling-n correlation is weak (rho={rho_cn:+.3f}).")
    print("  Cycling is not a proxy for network size.")

if abs(rho_cr) > 0.6:
    print("  WARNING: cycling is strongly correlated with mean_regs.")
    print("  Confound with regulatory complexity.")
else:
    print(f"  Cycling-mean_regs correlation is moderate (rho={rho_cr:+.3f}).")

# Also: is the cycling threshold sensitive to counting self-regulators?
print()
print("  Input nodes per model:")
for mid in model_ids:
    rules = BUILTIN_MODELS[mid]["rules"]
    inputs = [node for node, expr in rules.items() if expr.strip() == node]
    r = model_lookup[mid]
    print(f"    {mid:25s}: {inputs} (n_dynamic = {r['n'] - len(inputs)})")

# Save results
output = {
    "check1_g4_delta": {
        "delta_o3plus_pp": g4["delta_o3plus"] * 100,
        "abs_delta": abs(g4["delta_o3plus"]),
        "clears_energy_gate": g4["global_o3plus"] >= 1e-4,
        "verdict": "borderline_null",
    },
    "check2_confounds": {
        "spearman_cycling_n": {"rho": rho_cn, "p": p_cn},
        "spearman_cycling_mean_regs": {"rho": rho_cr, "p": p_cr},
        "spearman_cycling_density": {"rho": rho_cd, "p": p_cd},
        "spearman_delta_n": {"rho": rho_dn, "p": p_dn},
        "spearman_delta_mean_regs": {"rho": rho_dr, "p": p_dr},
        "spearman_delta_density": {"rho": rho_dd, "p": p_dd},
        "spearman_delta_cycling": {"rho": rho_dcy, "p": p_dcy},
    },
    "per_model": data,
}

outpath = "/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2/cheap_checks_results.json"
with open(outpath, "w") as f:
    json.dump(output, f, indent=2, default=float)
print(f"\nSaved to {outpath}")
