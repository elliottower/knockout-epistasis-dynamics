"""The Boolean arm's Delta_3+ under each registered estimator, from its per-trajectory outputs.

The per-trajectory Boolean outputs are the ones the cycle-handling audit recomputed for all 28
networks with the production settings (512 initial states, seed 42). For the 27 networks with a
committed coalition table they reproduce it exactly. grieco_bladder (18 nodes) has no committed
table, and its recomputed plain value, 3.7655 pp, rounds to the published 3.77. They are stored
as chunk shards under results/audit/boolean_cycle_handling/.
The estimators are those of scripts/walsh_estimators.py; split-half uses initial states 0-255
against 256-511. Registered in experiments/2026-09-21_ode-primary-rerun/PREREG.md; not run
before that registration is frozen.

    uv run python -m scripts.boolean_estimators --out experiments/2026-09-21_ode-primary-rerun/results

Writes <out>/boolean_estimators.json.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from grn_coalition_sweep import extract_rule_fourier
from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.paths import AUDIT, RESULTS
from scripts.walsh_estimators import ENERGY, ESTIMATORS, EstimatorError, fractions, order3_fraction

SHARDS = AUDIT / "boolean_cycle_handling" / "shards"


class ShardAssemblyError(Exception):
    pass


def boolean_outputs(name: str, n: int, shards: Path = SHARDS) -> np.ndarray:
    """The (2^n, 512) per-trajectory outputs of one network, assembled from its chunk shards."""
    parts = []
    expected = 0
    for f in sorted((shards / name).glob("*.npz")):
        if f.name.endswith(".partial.npz"):
            continue
        with np.load(f) as z:
            if int(z["start"]) != expected:
                raise ShardAssemblyError(f"{name}: chunk at {int(z['start'])}, expected {expected}")
            parts.append(z["impl"])
            expected = int(z["end"])
    if expected != 2**n:
        raise ShardAssemblyError(f"{name}: chunks cover [0, {expected}) of [0, {2**n})")
    return np.concatenate(parts)


def delta_by_estimator(rules: dict, outputs: np.ndarray, n: int) -> dict:
    """Delta_3+ under each estimator, with the raw order energies it came from."""
    local = order3_fraction(np.asarray(extract_rule_fourier(rules)["local_energy_spectrum"], dtype=float))
    out = {"local_o3plus": local}
    for name in ESTIMATORS:
        energy = ENERGY[name](outputs, n)
        raw = {"energy_by_order": energy.tolist(), "total_energy": float(energy.sum())}
        try:
            out[name] = {"admissible": True, "delta_o3plus": order3_fraction(fractions(energy)) - local, **raw}
        except EstimatorError as err:
            out[name] = {"admissible": False, "reason": str(err), **raw}
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    table = json.loads((RESULTS / "paper_number_reconciliation.json").read_text())["per_network_table"]
    record = {"estimators": list(ESTIMATORS), "networks": {}}
    for name in sorted(table, key=lambda m: (table[m]["n"], m)):
        n = table[name]["n"]
        record["networks"][name] = {"n_nodes": n, **delta_by_estimator(ALL_MODELS[name]["rules"], boolean_outputs(name, n), n)}
    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / "boolean_estimators.json"
    tmp = path.with_name(path.stem + ".partial.json")
    tmp.write_text(json.dumps(record, indent=2))
    tmp.replace(path)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
