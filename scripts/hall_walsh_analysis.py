"""Walsh-Hadamard analysis of the Hall yeast biosynthetic knockout fitness landscape.

Computes the energy spectrum and composition gap for 6 biosynthetic gene
knockouts from Hall et al. (2010).

Local model: additive (the 6 genes lie in distinct biosynthetic pathways;
some share upstream precursors but the exact locus-to-gene mapping is unknown,
so the local model conservatively assumes no pairwise pathway interactions).
Local o3+ = 0.

Multiple fitness components are analyzed: haploid growth rate (primary),
diploid growth rate, mating efficiency, sporulation efficiency.

Usage:
    uv run python scripts/hall_walsh_analysis.py
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from data_utils import normalized_wht, energy_by_order


FITNESS_COMPONENTS = [
    "haploid_growth_rate",
    "diploid_growth_rate",
    "mating_efficiency",
    "sporulation_efficiency",
]


def main():
    data_path = Path(__file__).parent.parent / "data" / "hall_2010" / "fitness_data.json"

    with open(data_path) as f:
        data = json.load(f)

    n = data["metadata"]["n_loci"]
    locus_names = data["metadata"]["locus_names"]

    results_by_component = {}
    for component in FITNESS_COMPONENTS:
        values = np.zeros(2**n)
        valid = True
        for g in data["genotypes"]:
            idx = int(g["genotype"], 2)
            entry = g.get(component)
            if entry is None:
                valid = False
                break
            values[idx] = entry["mean"]

        if not valid:
            continue

        w = normalized_wht(values)
        energy = energy_by_order(w, n)
        total = energy.sum()
        spectrum = (energy / total).tolist() if total > 0 else energy.tolist()
        o3plus = sum(spectrum[3:])

        results_by_component[component] = {
            "spectrum": spectrum,
            "global_o3plus": o3plus,
            "value_range": [float(values.min()), float(values.max())],
            "value_std": float(values.std()),
        }

        print(f"  {component}:")
        print(f"    spectrum: {[f'{s:.4f}' for s in spectrum]}")
        print(f"    global o3+: {o3plus:.4f}")

    primary = "haploid_growth_rate"
    primary_gap = results_by_component[primary]["global_o3plus"]

    output = {
        "dataset": "hall_2010",
        "source": data["metadata"]["source"],
        "organism": "Saccharomyces cerevisiae (BY series)",
        "n_loci": n,
        "locus_names": locus_names,
        "candidate_genes": data["metadata"]["candidate_genes"],
        "locus_mapping_note": data["metadata"]["locus_name_note"],
        "n_genotypes": len(data["genotypes"]),
        "local_model": "additive (distinct biosynthetic pathways; locus-gene mapping unknown)",
        "local_o3plus": 0.0,
        "primary_component": primary,
        "results_by_component": results_by_component,
        "composition_gap": {
            "delta_o3plus": primary_gap,
            "sign": "creation" if primary_gap > 0.005 else ("destruction" if primary_gap < -0.005 else "null"),
            "interpretation": (
                "Higher-order epistasis beyond pairwise arises from whole-cell "
                "integration of independent biosynthetic pathways through shared "
                "precursor pools, growth-rate coupling, and import competition."
            ),
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = Path(__file__).parent.parent / "results" / "empirical" / "hall_2010_walsh.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nHall yeast Walsh analysis")
    print(f"  {n} loci, {len(data['genotypes'])} genotypes")
    print(f"  Primary ({primary}) global o3+: {primary_gap:.4f}")
    print(f"  Composition gap: {primary_gap:+.4f}")
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
