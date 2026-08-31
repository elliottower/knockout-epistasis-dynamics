"""Structural analysis of Boolean GRN rules for blind prediction.

Extracts: local energy spectra, feedback loops, AND-gate arity,
pathway depth, and competing pathways from Boolean rules alone.
NO dynamics. NO coalition sweeps. Pure rule-level analysis.
"""

import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from grn_coalition_sweep import (
    BUILTIN_MODELS,
    _find_regulators,
    compile_network,
    extract_interaction_graph,
    extract_rule_fourier,
    identify_input_nodes,
)


SURVIVING_MODELS = [
    "faure_cellcycle",
    "tournier_apoptosis",
    "davidich_yeast",
    "drosophila_cellcycle",
    "fanconi_anemia",
    "arabidopsis_cellcycle",
]


def find_feedback_loops(rules, max_length=4):
    """Find all feedback loops up to max_length by DFS on the regulatory graph.

    Returns list of (loop_nodes, sign) tuples.
    sign: +1 if even number of inhibitions, -1 if odd.
    """
    node_names = list(rules.keys())
    graph = extract_interaction_graph(rules)

    # Build adjacency: target -> [(regulator, sign)]
    # We need edge direction: regulator -> target with sign
    adj = defaultdict(list)
    for target, edges in graph.items():
        for edge in edges:
            reg = edge["regulator"]
            sign = edge["sign"]
            adj[reg].append((target, sign))

    loops = []
    visited_loops = set()

    for start in node_names:
        # DFS from start, looking for paths that return to start
        stack = [(start, [start], 1)]  # (current, path, cumulative_sign)
        while stack:
            current, path, cum_sign = stack.pop()
            for neighbor, edge_sign in adj[current]:
                new_sign = cum_sign * edge_sign if edge_sign != 0 else 0
                if neighbor == start and len(path) >= 2:
                    # Found a loop
                    loop_key = frozenset(path)
                    if loop_key not in visited_loops and len(path) <= max_length:
                        visited_loops.add(loop_key)
                        loops.append((list(path), new_sign))
                elif neighbor not in path and len(path) < max_length:
                    stack.append((neighbor, path + [neighbor], new_sign))

    return loops


def max_and_arity(rules):
    """Find the maximum number of literals in any product term across all rules."""
    max_arity = 0
    for node, expr in rules.items():
        # Split by OR (|) to get product terms
        terms = expr.split("|")
        for term in terms:
            term = term.strip()
            # Count literals (variables and negated variables)
            # A literal is a variable name possibly preceded by !
            parts = term.replace("&", " ").split()
            n_literals = 0
            for p in parts:
                p = p.strip().lstrip("(").rstrip(")")
                if p and p not in ("", "(", ")"):
                    n_literals += 1
            if n_literals > max_arity:
                max_arity = n_literals
    return max_arity


def compute_pathway_depth(rules, output_nodes):
    """Compute longest non-feedback path from any node to the output.

    Uses BFS backwards from output, tracking the longest simple path.
    """
    node_names = list(rules.keys())
    graph = extract_interaction_graph(rules)

    # Build reverse adjacency: target -> list of regulators
    reverse_adj = defaultdict(list)
    for target, edges in graph.items():
        for edge in edges:
            reg = edge["regulator"]
            reverse_adj[target].append(reg)

    # BFS/DFS backwards from output to find longest path
    max_depth = 0
    for output in output_nodes:
        # DFS tracking path to avoid cycles
        stack = [(output, 0, {output})]
        while stack:
            node, depth, visited = stack.pop()
            if depth > max_depth:
                max_depth = depth
            for reg in reverse_adj.get(node, []):
                if reg not in visited:
                    stack.append((reg, depth + 1, visited | {reg}))

    return max_depth


def find_redundant_nodes(rules):
    """Find pairs of nodes with identical update rules (modulo variable names)."""
    redundant = []
    nodes = list(rules.keys())
    for i, n1 in enumerate(nodes):
        for n2 in nodes[i+1:]:
            if rules[n1] == rules[n2]:
                redundant.append((n1, n2))
    return redundant


def count_competing_pathways(rules, output_nodes):
    """Count independent paths to the output node from different sources."""
    graph = extract_interaction_graph(rules)
    reverse_adj = defaultdict(list)
    for target, edges in graph.items():
        for edge in edges:
            reverse_adj[target].append(edge["regulator"])

    for output in output_nodes:
        direct_regs = reverse_adj.get(output, [])
        return len(direct_regs)
    return 0


def analyze_model(model_name):
    """Full structural analysis of one model."""
    model_info = BUILTIN_MODELS[model_name]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]
    node_names = list(rules.keys())
    n = len(node_names)

    # Input nodes
    inputs = identify_input_nodes(rules)

    # Interaction graph
    graph = extract_interaction_graph(rules)
    n_edges = sum(len(edges) for edges in graph.values())

    # Local Fourier
    rf = extract_rule_fourier(rules)

    # Feedback loops
    loops = find_feedback_loops(rules, max_length=4)
    pos_loops = [l for l in loops if l[1] > 0]
    neg_loops = [l for l in loops if l[1] < 0]
    mixed_loops = [l for l in loops if l[1] == 0]

    # AND-gate arity
    arity = max_and_arity(rules)

    # Pathway depth
    depth = compute_pathway_depth(rules, output_nodes)

    # Redundant nodes
    redundant = find_redundant_nodes(rules)

    # Output rule analysis
    output_rule = rules[output_nodes[0]]
    output_regs = _find_regulators(output_rule, node_names)

    # Per-node regulator count
    reg_counts = {}
    for node in node_names:
        regs = _find_regulators(rules[node], node_names)
        reg_counts[node] = len(regs)

    result = {
        "model": model_name,
        "description": model_info["description"],
        "n_nodes": n,
        "node_names": node_names,
        "output_nodes": output_nodes,
        "output_rule": output_rule,
        "output_regulators": output_regs,
        "n_output_regulators": len(output_regs),
        "input_nodes": inputs,
        "n_edges": n_edges,
        "n_coalitions": 2**n,
        "max_and_arity": arity,
        "pathway_depth": depth,
        "redundant_nodes": redundant,
        "n_feedback_loops_le4": len(loops),
        "n_positive_loops": len(pos_loops),
        "n_negative_loops": len(neg_loops),
        "n_mixed_loops": len(mixed_loops),
        "positive_loops": [(l[0], l[1]) for l in pos_loops],
        "negative_loops": [(l[0], l[1]) for l in neg_loops],
        "local_energy_spectrum": rf["local_energy_spectrum"],
        "n_local_pairwise": rf["n_pairwise"],
        "n_local_triples": rf["n_triples"],
        "per_node_regulator_count": reg_counts,
        "rules": rules,
    }

    return result


def main():
    results = {}
    for model_name in SURVIVING_MODELS:
        print(f"\n{'='*60}")
        print(f"Analyzing: {model_name}")
        print(f"{'='*60}")
        result = analyze_model(model_name)

        print(f"  Nodes: {result['n_nodes']}, Edges: {result['n_edges']}")
        print(f"  Output: {result['output_nodes']} = {result['output_rule']}")
        print(f"  Output regulators ({result['n_output_regulators']}): {result['output_regulators']}")
        print(f"  Input nodes: {result['input_nodes']}")
        print(f"  Max AND arity: {result['max_and_arity']}")
        print(f"  Pathway depth: {result['pathway_depth']}")
        print(f"  Redundant nodes: {result['redundant_nodes']}")
        print(f"  Feedback loops (len <= 4): {result['n_feedback_loops_le4']}")
        print(f"    Positive: {result['n_positive_loops']}")
        print(f"    Negative: {result['n_negative_loops']}")
        print(f"    Mixed: {result['n_mixed_loops']}")
        print(f"  Local energy spectrum: {[f'{x:.3f}' for x in result['local_energy_spectrum']]}")
        print(f"  Local pairwise: {result['n_local_pairwise']}, triples: {result['n_local_triples']}")
        print(f"  Coalitions: {result['n_coalitions']}")

        # Print feedback loops
        print(f"\n  Positive feedback loops:")
        for loop_nodes, sign in result['positive_loops'][:10]:
            print(f"    {' -> '.join(loop_nodes)} -> {loop_nodes[0]} (sign={sign:+d})")
        print(f"\n  Negative feedback loops:")
        for loop_nodes, sign in result['negative_loops'][:10]:
            print(f"    {' -> '.join(loop_nodes)} -> {loop_nodes[0]} (sign={sign:+d})")

        results[model_name] = result

    # Save structural analysis (excluding rules dict for JSON serialization)
    output_path = Path(__file__).parent.parent / "results" / "grn_v2" / "structural_analysis.json"
    serializable = {}
    for model_name, result in results.items():
        r = dict(result)
        r.pop("rules", None)
        serializable[model_name] = r
    with open(output_path, "w") as f:
        json.dump(serializable, f, indent=2)
    print(f"\n\nSaved structural analysis: {output_path}")


if __name__ == "__main__":
    main()
