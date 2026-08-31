"""Merge batch 1 (6 models) and batch 2 (21 models) and compute
cycling-vs-delta statistics on the full N=27 panel."""
import json
import numpy as np
from scipy import stats

batch1_path = "/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2/blind_experiment_summary.json"
batch2_path = "/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2/blind_batch2_summary.json"

with open(batch1_path) as f:
    batch1 = json.load(f)
with open(batch2_path) as f:
    batch2 = json.load(f)

batch2_results = batch2["results"]
all_models = batch1 + batch2_results

print(f"Batch 1: {len(batch1)} models")
print(f"Batch 2: {len(batch2_results)} models")
print(f"Total:   {len(all_models)} models")

# Extract arrays
names = [m["model"] for m in all_models]
ns = [m["n"] for m in all_models]
deltas = [m["delta_o3plus"] for m in all_models]
cycling = [m["cycling_fraction"] for m in all_models]
rhos = [m["spearman_rho"] for m in all_models]
ps = [m["spearman_p"] for m in all_models]
directions = [m["creation_or_destruction"] for m in all_models]

# Classify: creation / destruction / null
NULL_THRESHOLD = 0.005  # |delta| < 0.5pp = null
labels = []
for d in deltas:
    if abs(d) < NULL_THRESHOLD:
        labels.append("null")
    elif d > 0:
        labels.append("creation")
    else:
        labels.append("destruction")

print(f"\n{'='*80}")
print(f"ALL {len(all_models)} MODELS (sorted by cycling fraction)")
print(f"{'='*80}")
print(f"{'Model':30s} {'n':>3s} {'cycling':>8s} {'delta_pp':>9s} {'label':>12s} {'rho':>6s} {'p':>8s}")
print("-" * 80)
for i in np.argsort(cycling):
    m = all_models[i]
    d_pp = m["delta_o3plus"] * 100
    print(f"{m['model']:30s} {m['n']:3d} {m['cycling_fraction']*100:7.1f}% {d_pp:+8.2f}pp {labels[i]:>12s} {m['spearman_rho']:+6.3f} {m['spearman_p']:8.2e}")

# Count
n_creation = labels.count("creation")
n_destruction = labels.count("destruction")
n_null = labels.count("null")
print(f"\nCreation: {n_creation}, Destruction: {n_destruction}, Null: {n_null}")

# Spearman: delta vs cycling (all models)
rho_dc, p_dc = stats.spearmanr(deltas, cycling)
print(f"\n{'='*60}")
print(f"CYCLING HYPOTHESIS TEST (N={len(all_models)})")
print(f"{'='*60}")
print(f"Spearman(delta, cycling) = {rho_dc:+.3f}  p = {p_dc:.4f}")

# Exclude nulls and retest
non_null_idx = [i for i, l in enumerate(labels) if l != "null"]
nn_deltas = [deltas[i] for i in non_null_idx]
nn_cycling = [cycling[i] for i in non_null_idx]
nn_labels = [labels[i] for i in non_null_idx]
rho_nn, p_nn = stats.spearmanr(nn_deltas, nn_cycling)
print(f"Spearman(delta, cycling) excluding nulls (N={len(non_null_idx)}) = {rho_nn:+.3f}  p = {p_nn:.4f}")

# Sign separation test: do all creation precede all destruction in cycling order?
creation_cycling = [cycling[i] for i, l in enumerate(labels) if l == "creation"]
destruction_cycling = [cycling[i] for i, l in enumerate(labels) if l == "destruction"]
print(f"\nCreation cycling: min={min(creation_cycling)*100:.1f}%, max={max(creation_cycling)*100:.1f}%, median={np.median(creation_cycling)*100:.1f}%")
print(f"Destruction cycling: min={min(destruction_cycling)*100:.1f}%, max={max(destruction_cycling)*100:.1f}%, median={np.median(destruction_cycling)*100:.1f}%")

# Mann-Whitney U test: creation vs destruction cycling
u_stat, u_p = stats.mannwhitneyu(creation_cycling, destruction_cycling, alternative='two-sided')
print(f"Mann-Whitney U (creation vs destruction cycling): U={u_stat:.1f}, p={u_p:.4f}")

# Key counterexamples to the threshold
print(f"\n{'='*60}")
print(f"COUNTEREXAMPLES TO CYCLING THRESHOLD")
print(f"{'='*60}")
print("Low cycling + destruction:")
for i, m in enumerate(all_models):
    if labels[i] == "destruction" and m["cycling_fraction"] < 0.10:
        print(f"  {m['model']:30s} cycling={m['cycling_fraction']*100:.1f}% delta={m['delta_o3plus']*100:+.2f}pp")

print("High cycling + creation:")
for i, m in enumerate(all_models):
    if labels[i] == "creation" and m["cycling_fraction"] > 0.15:
        print(f"  {m['model']:30s} cycling={m['cycling_fraction']*100:.1f}% delta={m['delta_o3plus']*100:+.2f}pp")

# Universal positive rho check
print(f"\n{'='*60}")
print(f"UNIVERSAL POSITIVE RHO CHECK")
print(f"{'='*60}")
negative_rho = [(m["model"], m["spearman_rho"]) for m in all_models if m["spearman_rho"] < 0]
print(f"Models with negative rho: {len(negative_rho)}/{len(all_models)}")
for name, r in negative_rho:
    print(f"  {name}: rho = {r:.3f}")
print(f"Median rho: {np.median(rhos):.3f}")
print(f"Mean rho: {np.mean(rhos):.3f}")

# Significant rho count
sig_rho = sum(1 for p in ps if p < 0.05)
print(f"Significant (p<0.05): {sig_rho}/{len(all_models)}")

# Overall delta statistics
print(f"\n{'='*60}")
print(f"DELTA STATISTICS")
print(f"{'='*60}")
print(f"Median delta: {np.median(deltas)*100:+.2f}pp")
print(f"Mean delta: {np.mean(deltas)*100:+.2f}pp")
print(f"Creation bias: {n_creation}/{n_creation+n_destruction} non-null models show creation")

# Save merged
merged = {
    "n_models": len(all_models),
    "batch1_models": len(batch1),
    "batch2_models": len(batch2_results),
    "cycling_hypothesis": {
        "spearman_rho": float(rho_dc),
        "spearman_p": float(p_dc),
        "spearman_rho_no_nulls": float(rho_nn),
        "spearman_p_no_nulls": float(p_nn),
        "mann_whitney_U": float(u_stat),
        "mann_whitney_p": float(u_p),
        "verdict": "significant" if p_dc < 0.05 else "not_significant",
    },
    "rho_statistics": {
        "median": float(np.median(rhos)),
        "mean": float(np.mean(rhos)),
        "n_negative": len(negative_rho),
        "n_significant": sig_rho,
    },
    "delta_statistics": {
        "median_pp": float(np.median(deltas) * 100),
        "mean_pp": float(np.mean(deltas) * 100),
        "n_creation": n_creation,
        "n_destruction": n_destruction,
        "n_null": n_null,
    },
    "all_models": [
        {
            "model": m["model"],
            "n": m["n"],
            "delta_o3plus": m["delta_o3plus"],
            "cycling_fraction": m["cycling_fraction"],
            "spearman_rho": m["spearman_rho"],
            "spearman_p": m["spearman_p"],
            "label": labels[i],
            "batch": "batch1" if i < len(batch1) else "batch2",
        }
        for i, m in enumerate(all_models)
    ],
}

outpath = "/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2/merged_all_27_analysis.json"
with open(outpath, "w") as f:
    json.dump(merged, f, indent=2)
print(f"\nSaved to {outpath}")
