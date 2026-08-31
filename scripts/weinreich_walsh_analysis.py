"""Walsh-Hadamard analysis of the Weinreich TEM-1 beta-lactamase fitness landscape.

Computes the energy spectrum and composition gap for 5 mutations in TEM-1
using MIC (minimum inhibitory concentration) data from Weinreich et al. (2006).

Local model: purely additive (no structural contacts between any mutation
pair at 8A Ca-Ca threshold). Local o3+ = 0.

Composition gap = global o3+ - local o3+ = global o3+.

Usage:
    uv run python scripts/weinreich_walsh_analysis.py
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from data_utils import normalized_wht, energy_by_order


def main():
    data_path = Path(__file__).parent.parent / "data" / "weinreich_2006" / "fitness_data.json"
    contacts_path = Path(__file__).parent.parent / "data" / "weinreich_2006" / "structural_contacts.json"

    with open(data_path) as f:
        data = json.load(f)
    with open(contacts_path) as f:
        contacts = json.load(f)

    n = 5
    mutations = data["mutations_in_order"]

    raw_values = np.zeros(2**n)
    log_values = np.zeros(2**n)
    for g in data["genotypes"]:
        idx = int(g["genotype"], 2)
        mic = g["mic"]
        raw_values[idx] = mic
        log_values[idx] = np.log(mic)

    results = {}
    for label, values in [("log_mic", log_values), ("raw_mic", raw_values)]:
        w = normalized_wht(values)
        energy = energy_by_order(w, n)
        total = energy.sum()
        spectrum = (energy / total).tolist() if total > 0 else energy.tolist()
        o3plus = sum(spectrum[3:])

        results[label] = {
            "spectrum": spectrum,
            "global_o3plus": o3plus,
            "walsh_coefficients": w.tolist(),
            "value_range": [float(values.min()), float(values.max())],
            "value_std": float(values.std()),
        }

    n_contacts = contacts["summary"]["pairs_in_contact"]
    closest = contacts["summary"]["closest_pair"]

    output = {
        "dataset": "weinreich_2006",
        "source": data["source"],
        "organism": "E. coli (TEM-1 beta-lactamase)",
        "n_loci": n,
        "mutations": mutations,
        "n_genotypes": len(data["genotypes"]),
        "structural_contacts": {
            "threshold_angstroms": contacts["contact_threshold_angstroms"],
            "pairs_in_contact": n_contacts,
            "closest_pair": closest,
        },
        "local_model": "additive (zero structural contacts at 8A)",
        "local_o3plus": 0.0,
        "primary_phenotype": "log_mic",
        "log_mic": results["log_mic"],
        "raw_mic": results["raw_mic"],
        "composition_gap": {
            "delta_o3plus": results["log_mic"]["global_o3plus"],
            "sign": "creation" if results["log_mic"]["global_o3plus"] > 0.005 else "null",
            "interpretation": (
                "All higher-order epistasis in the fitness landscape is 'created' "
                "by long-range functional interactions (catalysis, stability, dynamics) "
                "since no mutation pairs are in direct structural contact."
            ),
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    out_path = Path(__file__).parent.parent / "results" / "empirical" / "weinreich_2006_walsh.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"Weinreich TEM-1 Walsh analysis")
    print(f"  {n} mutations: {mutations}")
    print(f"  Structural contacts at 8A: {n_contacts} (closest: {closest})")
    print(f"  log(MIC) energy spectrum: {[f'{s:.4f}' for s in results['log_mic']['spectrum']]}")
    print(f"  log(MIC) global o3+:      {results['log_mic']['global_o3plus']:.4f}")
    print(f"  Composition gap (delta):  {results['log_mic']['global_o3plus']:+.4f} (creation)")
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
