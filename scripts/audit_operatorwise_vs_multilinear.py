"""Where does the operator-wise ODE conversion depart from the multilinear extension?

The legacy ODE engine replaces AND, OR and NOT with product, probabilistic sum and
complement, recursively over each rule's formula. The multilinear extension of the
same rule's truth table (Wittmann et al. 2009, BooleCube) is unique and does not
depend on how the formula is written. The two agree on read-once formulas and can
differ when a regulator occurs more than once.

For every rule in the 28 published networks this measures the largest disagreement
on a fixed interior grid, and checks as a positive control that the two agree
exactly at every vertex of the unit cube. No dynamics are run.

Run from the repository root, in the project's pinned environment:
    uv run python -m scripts.audit_operatorwise_vs_multilinear

Writes results/audit/operatorwise_vs_multilinear.json.
"""

import itertools
import json
import re

import numpy as np

from grn_coalition_sweep import compile_network
from scripts.ode_coalition_sweep import ALL_MODELS, boolean_expr_to_continuous
from scripts.paths import AUDIT, PROJECT_ROOT, RESULTS

OUT = AUDIT / "operatorwise_vs_multilinear.json"
NETWORKS_FROM = RESULTS / "paper_number_reconciliation.json"

INTERIOR_LEVELS = (0.2, 0.5, 0.8)
SPARSE_LEVELS = (0.25, 0.75)
MAX_DENSE_POINTS = 20_000
TOLERANCE = 1e-9


class AuditError(Exception):
    pass


def interior_grid(k):
    """Deterministic interior points: a 3-level grid, or 2 levels plus the centre when large."""
    if len(INTERIOR_LEVELS) ** k <= MAX_DENSE_POINTS:
        pts = np.array(list(itertools.product(INTERIOR_LEVELS, repeat=k)), dtype=float)
        return pts, f"{INTERIOR_LEVELS}^{k}"
    pts = np.array(list(itertools.product(SPARSE_LEVELS, repeat=k)), dtype=float)
    pts = np.vstack([pts, np.full((1, k), 0.5)])
    return pts, f"{SPARSE_LEVELS}^{k} + centre"


def vertices(k):
    return np.array(list(itertools.product((0.0, 1.0), repeat=k)), dtype=float)


def multilinear(points, truth_table):
    """Unique multilinear extension of a truth table indexed with bit j = regulator j."""
    m, k = points.shape
    idx = np.arange(2**k)
    basis = np.ones((m, 2**k))
    for j in range(k):
        bit = ((idx >> j) & 1).astype(bool)
        xj = points[:, j:j + 1]
        basis *= np.where(bit, xj, 1.0 - xj)
    return basis @ truth_table.astype(float)


def operatorwise(points, fn, reg_names, all_names):
    x = {name: 0.0 for name in all_names}
    out = np.empty(points.shape[0])
    for i, p in enumerate(points):
        for name, v in zip(reg_names, p):
            x[name] = float(v)
        out[i] = fn(x)
    return out


def audit_network(model):
    if model not in ALL_MODELS:
        raise AuditError(f"network {model!r} is in the paper but not in ALL_MODELS")
    rules = ALL_MODELS[model]["rules"]
    compiled, node_names = compile_network(rules)
    rows = []
    for node, reg_indices, truth_table in compiled:
        expr = rules[node]
        regs = [node_names[i] for i in reg_indices]
        k = len(regs)
        occurrences = {r: len(re.findall(r"\b" + re.escape(r) + r"\b", expr)) for r in regs}
        repeated = any(c > 1 for c in occurrences.values())
        if k == 0:
            # A constant rule has one vertex and no interior; check the vertex.
            fn = boolean_expr_to_continuous(expr, node_names)
            diff = abs(fn({name: 0.0 for name in node_names}) - float(truth_table[0]))
            rows.append({"node": node, "k": 0, "repeated": False, "repeated_regulators": [],
                         "max_diff_interior": 0.0, "max_diff_vertices": diff, "grid": "constant"})
            continue
        fn = boolean_expr_to_continuous(expr, node_names)
        grid, grid_desc = interior_grid(k)
        verts = vertices(k)
        d_int = np.abs(operatorwise(grid, fn, regs, node_names) - multilinear(grid, truth_table))
        d_vert = np.abs(operatorwise(verts, fn, regs, node_names) - multilinear(verts, truth_table))
        rows.append({
            "node": node,
            "k": k,
            "repeated": repeated,
            "repeated_regulators": sorted(r for r, c in occurrences.items() if c > 1),
            "max_diff_interior": float(d_int.max()),
            "max_diff_vertices": float(d_vert.max()),
            "grid": grid_desc,
        })
    return rows


def main():
    networks = sorted(json.loads(NETWORKS_FROM.read_text())["per_network_table"])
    if len(networks) != 28:
        raise AuditError(f"expected 28 networks, found {len(networks)}")

    per_network = {m: audit_network(m) for m in networks}
    rules = [r for rows in per_network.values() for r in rows]
    differs = [r for r in rules if r["max_diff_interior"] > TOLERANCE]

    summary = {
        "n_networks": len(networks),
        "n_rules": len(rules),
        "n_rules_with_repeated_regulator": sum(r["repeated"] for r in rules),
        "n_rules_differing": len(differs),
        "n_rules_differing_without_repetition": sum(not r["repeated"] for r in differs),
        "n_rules_repeated_but_agreeing": sum(r["repeated"] and r["max_diff_interior"] <= TOLERANCE
                                             for r in rules),
        "max_diff_interior": max(r["max_diff_interior"] for r in rules),
        "vertex_control_max_diff": max(r["max_diff_vertices"] for r in rules),
        "n_networks_with_a_differing_rule": sum(
            any(r["max_diff_interior"] > TOLERANCE for r in rows) for rows in per_network.values()),
        "per_network": {
            m: {"n_rules": len(rows),
                "n_differing": sum(r["max_diff_interior"] > TOLERANCE for r in rows),
                "max_diff": max(r["max_diff_interior"] for r in rows)}
            for m, rows in per_network.items()},
    }
    result = {
        "question": "Where does the operator-wise conversion depart from the multilinear extension?",
        "inputs": {"networks_from": str(NETWORKS_FROM.relative_to(PROJECT_ROOT)),
                   "interior_levels": INTERIOR_LEVELS, "sparse_levels": SPARSE_LEVELS,
                   "max_dense_points": MAX_DENSE_POINTS, "tolerance": TOLERANCE},
        "summary": summary,
        "rules": per_network,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2))

    print(f"wrote {OUT.relative_to(PROJECT_ROOT)}")
    for key in ("n_networks", "n_rules", "n_rules_with_repeated_regulator", "n_rules_differing",
                "n_rules_differing_without_repetition", "n_rules_repeated_but_agreeing",
                "max_diff_interior", "vertex_control_max_diff", "n_networks_with_a_differing_rule"):
        print(f"  {key:<40} {summary[key]}")


if __name__ == "__main__":
    main()
