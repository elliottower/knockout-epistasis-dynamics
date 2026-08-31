"""Null model generators for Boolean network composition gap analysis.

Three nested nulls, ordered by what they preserve:

| Null                    | Preserves                        | Isolates                        |
|-------------------------|----------------------------------|---------------------------------|
| Kauffman NK             | n, mean k                        | Whether anything about biology matters |
| Degree-preserving       | in-degree sequence               | Topology vs rules jointly       |
| Rule-preserving rewire  | exact truth tables + in-degrees  | Wiring alone (composition claim)|

The third is the sharpest test: if real networks show creation and
rule-preserving rewires don't, composition topology drives the gap
with rules held constant.

Input nodes (f(x) = x) are kept fixed in all three nulls — they
serve as external signals and aren't part of the composition story.
"""

import numpy as np

from data_utils import energy_by_order, normalized_wht, popcount_array


def identify_input_compiled(compiled):
    """Find input node indices: nodes whose truth table is identity (f(x) = x).

    An input node has exactly one regulator (itself) and truth table [0, 1].
    """
    input_indices = set()
    for col, (node_name, reg_indices, truth_table) in enumerate(compiled):
        if len(reg_indices) == 1 and reg_indices[0] == col:
            if len(truth_table) == 2 and truth_table[0] == 0 and truth_table[1] == 1:
                input_indices.add(col)
    return input_indices


def rule_preserving_rewire(compiled, rng):
    """Rewire a Boolean network keeping truth tables and in-degrees fixed.

    For each non-input node, randomly reassign which k nodes serve as
    regulators while keeping the same truth table. Tests whether the
    specific wiring pattern (which nodes regulate which) matters for
    the composition gap, with rules held constant.
    """
    n = len(compiled)
    input_indices = identify_input_compiled(compiled)
    new_compiled = []

    for col, (node_name, reg_indices, truth_table) in enumerate(compiled):
        if col in input_indices:
            new_compiled.append((node_name, list(reg_indices), truth_table.copy()))
            continue

        k = len(reg_indices)
        if k == 0:
            new_compiled.append((node_name, [], truth_table.copy()))
            continue

        new_regs = sorted(rng.choice(n, size=k, replace=False).tolist())
        new_compiled.append((node_name, new_regs, truth_table.copy()))

    return new_compiled


def degree_preserving_rewire(compiled, rng):
    """Rewire with random truth tables, keeping in-degree sequence fixed.

    For each non-input node, randomly reassign regulators AND generate a
    new random truth table. Tests whether the in-degree distribution alone
    explains the gap.
    """
    n = len(compiled)
    input_indices = identify_input_compiled(compiled)
    new_compiled = []

    for col, (node_name, reg_indices, truth_table) in enumerate(compiled):
        if col in input_indices:
            new_compiled.append((node_name, list(reg_indices), truth_table.copy()))
            continue

        k = len(reg_indices)
        if k == 0:
            new_compiled.append((node_name, [], truth_table.copy()))
            continue

        new_regs = sorted(rng.choice(n, size=k, replace=False).tolist())
        new_tt = rng.integers(0, 2, size=2**k).astype(np.int8)
        new_compiled.append((node_name, new_regs, new_tt))

    return new_compiled


def kauffman_nk(compiled, rng):
    """Generate a Kauffman NK network matched to the real network's size and mean connectivity.

    Preserves n and mean in-degree (rounded). Input nodes are kept fixed.
    All other nodes get random wiring and random truth tables.
    """
    n = len(compiled)
    input_indices = identify_input_compiled(compiled)
    non_input = [i for i in range(n) if i not in input_indices]

    if non_input:
        mean_k = np.mean([len(compiled[i][1]) for i in non_input])
        k = max(1, int(round(mean_k)))
    else:
        k = 1

    new_compiled = []
    for col, (node_name, reg_indices, truth_table) in enumerate(compiled):
        if col in input_indices:
            new_compiled.append((node_name, list(reg_indices), truth_table.copy()))
            continue

        new_regs = sorted(rng.choice(n, size=k, replace=False).tolist())
        new_tt = rng.integers(0, 2, size=2**k).astype(np.int8)
        new_compiled.append((node_name, new_regs, new_tt))

    return new_compiled


def sweep_compiled(compiled, output_indices, n_init=128, max_steps=200, seed=42,
                   clamp_value=0):
    """Coalition sweep on a pre-compiled network.

    Like grn_coalition_sweep.sweep_coalitions but takes compiled rules
    directly, avoiding the string-expression pipeline. Skips input nodes
    from the player set.
    """
    from grn_coalition_sweep import simulate_sync_output

    n = len(compiled)
    input_indices = identify_input_compiled(compiled)
    player_indices = [i for i in range(n) if i not in input_indices]
    n_players = len(player_indices)
    N = 2**n_players

    rng = np.random.default_rng(seed)
    init_states = rng.integers(0, 2, size=(n_init, n), dtype=np.int8)

    values = np.zeros((N, n_init), dtype=np.float64)
    total_cycling = 0
    total_fixed = 0

    for coalition in range(N):
        clamp_mask = np.zeros(n, dtype=bool)
        for bit_pos, node_idx in enumerate(player_indices):
            if not (coalition & (1 << bit_pos)):
                clamp_mask[node_idx] = True

        output, info = simulate_sync_output(
            init_states, compiled, clamp_mask, clamp_value,
            output_indices, max_steps=max_steps,
        )
        total_cycling += info["n_cycling"]
        total_fixed += info["n_fixed_point"]
        values[coalition] = output

    v_mean = values.mean(axis=1)
    return v_mean, n_players, total_cycling / (N * n_init)


def compute_delta_3plus(v_mean, n_players):
    """Compute global order-3+ energy fraction from the Walsh decomposition."""
    w = normalized_wht(v_mean)
    pc = popcount_array(n_players)
    total_energy = float(np.sum(w[pc > 0] ** 2))
    if total_energy < 1e-15:
        return None
    energy_3plus = float(np.sum(w[pc >= 3] ** 2))
    return 100.0 * energy_3plus / total_energy


def quality_gate(v_mean, min_unique=3, max_single_frac=0.95):
    """Check whether a value function is non-degenerate.

    Returns (passes, n_unique, max_frac).
    A value function with <3 unique values or >95% dominated by one value
    is degenerate — the Walsh spectrum is technically computable but
    reflects trivial dynamics (output always 0 or always 1).

    Using a lenient gate: the rejection rate itself is a key finding.
    """
    unique_vals = np.unique(v_mean)
    n_unique = len(unique_vals)
    counts = np.array([np.sum(v_mean == v) for v in unique_vals])
    max_frac = float(counts.max() / len(v_mean))
    passes = n_unique >= min_unique and max_frac <= max_single_frac
    return passes, n_unique, max_frac
