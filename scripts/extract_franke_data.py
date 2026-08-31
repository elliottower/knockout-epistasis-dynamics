"""Extract Franke et al. 2011 fitness landscape data from supplementary PDF.

Source: Franke et al. (2011) Evolutionary Accessibility of Mutational Pathways.
        PLOS Computational Biology 7(8): e1002134.
        Table S1: Mean fitness W (mycelium growth rate) of 186 segregants of
        A. niger relative to wildtype. Missing genotypes marked 'm'.

Usage:
    uv run python scripts/extract_franke_data.py
"""
import csv
import json
import sys
from pathlib import Path

import fitz

MARKERS = ["fwn", "arg", "pyr", "leu", "phe", "lys", "oli", "crn"]

PDF_PATH = Path(__file__).parent.parent / "data" / "franke2011_tableS1.pdf"
FALLBACK_PATH = Path(
    "/Users/elliottower/.claude/projects/"
    "-Users-elliottower-Documents-GitHub-factorization-circuits/"
    "1f8cce1a-7332-4537-9ab7-69983983aa12/tool-results/"
    "webfetch-1785298221668-5dbk44.pdf"
)


def extract_from_pdf(pdf_path):
    doc = fitz.open(str(pdf_path))
    all_rows = []

    for page in doc:
        tables = list(page.find_tables())
        if not tables:
            continue
        raw = tables[0].extract()
        for row in raw:
            cols = [col.split("\n") if col else [] for col in row]
            if not cols:
                continue
            n_entries = len(cols[0])
            for i in range(n_entries):
                entry = []
                for col in cols:
                    entry.append(col[i] if i < len(col) else "")
                all_rows.append(entry)

    return all_rows


def parse_rows(raw_rows):
    genotypes = []
    for row in raw_rows:
        if len(row) < 10:
            continue
        try:
            int(row[0].strip())
        except ValueError:
            continue

        mut_num = int(row[0].strip())
        markers = [int(row[i].strip()) for i in range(1, 9)]
        fitness_str = row[9].strip()

        if fitness_str == "m":
            fitness = None
        else:
            fitness = float(fitness_str)

        genotypes.append({
            "mutation_number": mut_num,
            "markers": dict(zip(MARKERS, markers)),
            "genotype_binary": markers,
            "fitness": fitness,
        })

    return genotypes


def to_fitness_array(genotypes):
    """Convert to a 256-element array indexed by binary genotype.

    Index = sum(marker_i * 2^i) where markers are in order fwn...crn.
    Missing genotypes get fitness 0.0 (non-viable).
    """
    import numpy as np
    values = np.zeros(256)
    for g in genotypes:
        idx = sum(b * (2 ** i) for i, b in enumerate(g["genotype_binary"]))
        values[idx] = g["fitness"] if g["fitness"] is not None else 0.0
    return values


def main():
    pdf_path = PDF_PATH if PDF_PATH.exists() else FALLBACK_PATH
    if not pdf_path.exists():
        print(f"PDF not found at {PDF_PATH} or {FALLBACK_PATH}")
        sys.exit(1)

    print(f"Reading PDF from {pdf_path}")
    raw_rows = extract_from_pdf(pdf_path)
    print(f"Extracted {len(raw_rows)} raw rows")

    genotypes = parse_rows(raw_rows)
    print(f"Parsed {len(genotypes)} genotypes")

    n_measured = sum(1 for g in genotypes if g["fitness"] is not None)
    n_missing = sum(1 for g in genotypes if g["fitness"] is None)
    print(f"  Measured: {n_measured}, Missing: {n_missing}")

    out_dir = Path(__file__).parent.parent / "data" / "empirical_landscapes"
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "franke2011_aspergillus.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(MARKERS + ["fitness", "missing"])
        for g in genotypes:
            row = g["genotype_binary"] + [
                g["fitness"] if g["fitness"] is not None else 0.0,
                1 if g["fitness"] is None else 0,
            ]
            writer.writerow(row)
    print(f"Saved CSV to {csv_path}")

    json_path = out_dir / "franke2011_aspergillus.json"
    with open(json_path, "w") as f:
        json.dump({
            "source": "Franke et al. (2011) PLOS Comp Bio 7(8):e1002134, Table S1",
            "organism": "Aspergillus niger",
            "phenotype": "mycelium growth rate relative to wildtype",
            "n_loci": 8,
            "n_genotypes": 256,
            "n_measured": n_measured,
            "n_missing": n_missing,
            "markers": MARKERS,
            "genotypes": [
                {
                    "binary": g["genotype_binary"],
                    "fitness": g["fitness"],
                    "mutation_number": g["mutation_number"],
                }
                for g in genotypes
            ],
        }, f, indent=2)
    print(f"Saved JSON to {json_path}")

    values = to_fitness_array(genotypes)
    print(f"\nFitness array shape: {values.shape}")
    print(f"Wildtype fitness: {values[0]:.3f}")
    print(f"Mean (non-zero): {values[values > 0].mean():.3f}")
    print(f"Min (non-zero): {values[values > 0].min():.3f}")
    print(f"Max: {values.max():.3f}")

    return genotypes, values


if __name__ == "__main__":
    main()
