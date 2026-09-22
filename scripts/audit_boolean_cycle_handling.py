"""Does the implemented Boolean cycle handling give the outputs the registered procedure would?

The registrations (prereg_grn_v3.md, prereg_grn_v4.md) specify, for a trajectory that has
not converged by step 200: compare the final state with the states of the previous 50 steps;
if it recurs, average the output over one full period; if not, average the output over the
last 50 steps.

`grn_coalition_sweep.simulate_sync_output` does something else: from the step-200 state it
walks forward up to 200 further steps until that state recurs, and averages over the cycle;
if the state does not recur, it keeps the instantaneous output.

The two agree whenever the step-200 state lies on a cycle of period at most 50. They can
disagree for longer cycles and for trajectories still in their transient. This script
re-simulates every coalition of every network, computes both procedures from the same
trajectory, and counts where they disagree.

Positive control: the implemented branch is recomputed here, not called, so its per-initial-
state outputs must equal the committed coalition tables (results/grn_v2/*_coalition_blind.npz)
exactly. Grieco's table was never committed; for it the check is the published Delta_3+.

    uv run python -m scripts.audit_boolean_cycle_handling [--networks a b ...] [--workers 8]

Writes one JSON per network under results/audit/boolean_cycle_handling/, a coalition-value
file per network, and a summary. Work is checkpointed in append-only chunk shards and
resumes by skipping finished chunks.
"""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from data_utils import energy_spectrum, normalized_wht
from grn_coalition_sweep import compile_network, extract_rule_fourier, update_all
from scripts.ode_coalition_sweep import ALL_MODELS
from scripts.paths import AUDIT, PROJECT_ROOT, RESULTS

OUT = AUDIT / "boolean_cycle_handling"
PUBLISHED = RESULTS / "paper_number_reconciliation.json"
COMMITTED_TABLE = RESULTS / "grn_v2" / "{name}_coalition_blind.npz"

# The published sweep's settings (scripts/run_batch2_blind_sweep.py).
N_INIT = 512
SEED = 42
MAX_STEPS = 200
CLAMP_VALUE = 0
# The registered procedure's window.
WINDOW = 50
CHUNK = 2048
TOLERANCE = 1e-12


class AuditError(Exception):
    pass


def step_outputs(states: np.ndarray, output_indices: list[int]) -> np.ndarray:
    """Per-state output, with the production arithmetic: sum over output nodes, then divide."""
    out = np.zeros(states.shape[0], dtype=np.float64)
    for oi in output_indices:
        out += states[:, oi].astype(np.float64)
    return out / len(output_indices)


def both_procedures(init, compiled, clamp_mask, output_indices):
    """Return (implemented, registered, counters) for one coalition's N_INIT trajectories."""
    states = init.copy()
    states[:, clamp_mask] = CLAMP_VALUE
    converged = np.zeros(states.shape[0], dtype=bool)
    history = np.empty((WINDOW + 1,) + states.shape, dtype=states.dtype)  # steps 150..200
    for step in range(1, MAX_STEPS + 1):
        new = update_all(states, compiled)
        new[:, clamp_mask] = CLAMP_VALUE
        converged |= np.all(new == states, axis=1)
        states = new
        if step >= MAX_STEPS - WINDOW:
            history[step - (MAX_STEPS - WINDOW)] = states
    final = states

    implemented = step_outputs(final, output_indices)
    registered = implemented.copy()
    cyc = np.flatnonzero(~converged)
    counters = {"n_trajectories": int(final.shape[0]), "n_cycling": int(cyc.size),
                "n_unclosed_implemented": 0, "n_no_cycle_in_window_registered": 0,
                "n_period_over_window": 0}
    if cyc.size == 0:
        return implemented, registered, counters

    # Implemented: walk forward from the step-200 state until it recurs.
    anchor = final[cyc]
    walk = anchor.copy()
    vals = [step_outputs(walk, output_indices)]
    period = np.zeros(cyc.size, dtype=np.int64)
    open_ = np.ones(cyc.size, dtype=bool)
    clamp_rows = clamp_mask
    for t in range(1, MAX_STEPS + 1):
        walk = update_all(walk, compiled)
        walk[:, clamp_rows] = CLAMP_VALUE
        back = open_ & np.all(walk == anchor, axis=1)
        period[back] = t
        open_ &= ~back
        if not open_.any():
            break
        vals.append(step_outputs(walk, output_indices))
    vals = np.stack(vals)  # (steps, n_cyc), row 0 is the anchor
    closed = period > 0
    if len(output_indices) == 1:
        # Values are exactly 0 or 1, so any summation order gives production's mean exactly.
        for p in np.unique(period[closed]):
            cols = closed & (period == p)
            implemented[cyc[cols]] = vals[:p, cols].sum(axis=0) / p
    else:
        # Production takes np.mean of a Python list per trajectory; reproduce it exactly.
        for j in np.flatnonzero(closed):
            implemented[cyc[j]] = np.mean(vals[:period[j], j].tolist())
    counters["n_unclosed_implemented"] = int((~closed).sum())
    counters["n_period_over_window"] = int((period > WINDOW).sum())

    # Registered: look back WINDOW steps for the final state; else average the last WINDOW steps.
    hist_out = np.stack([step_outputs(h, output_indices) for h in history])[:, cyc]  # (51, n_cyc)
    eq = np.all(history[:WINDOW, cyc, :] == final[cyc][None, :, :], axis=2)         # steps 150..199
    found = eq.any(axis=0)
    latest = WINDOW - 1 - np.argmax(eq[::-1], axis=0)
    p_reg = WINDOW - latest
    reg = np.empty(cyc.size)
    reg[~found] = hist_out[1:, ~found].mean(axis=0)                                   # steps 151..200
    for p in np.unique(p_reg[found]):
        cols = found & (p_reg == p)
        reg[cols] = hist_out[WINDOW - p + 1:, cols].mean(axis=0)                     # one full period
    registered[cyc] = reg
    counters["n_no_cycle_in_window_registered"] = int((~found).sum())
    return implemented, registered, counters


def compute_chunk(name, start, end):
    """Both procedures for coalitions [start, end) of one network. Returns (impl, reg, totals)."""
    info = ALL_MODELS[name]
    compiled, node_names = compile_network(info["rules"])
    n = len(node_names)
    output_indices = [node_names.index(o) for o in info["output_nodes"]]
    init = np.random.default_rng(SEED).integers(0, 2, size=(N_INIT, n), dtype=np.int8)
    impl = np.empty((end - start, N_INIT))
    reg = np.empty((end - start, N_INIT))
    totals = {}
    for k, coalition in enumerate(range(start, end)):
        mask = np.array([not (coalition >> j) & 1 for j in range(n)], dtype=bool)
        a, b, c = both_procedures(init, compiled, mask, output_indices)
        impl[k], reg[k] = a, b
        for key, v in c.items():
            totals[key] = totals.get(key, 0) + v
    return impl, reg, totals


def run_chunk(args):
    name, start, end, shard = args
    if shard.exists():
        return str(shard)
    impl, reg, totals = compute_chunk(name, start, end)
    tmp = shard.with_name(shard.stem + ".partial.npz")
    np.savez_compressed(tmp, impl=impl, reg=reg, start=start, end=end,
                        counters=json.dumps(totals))
    tmp.replace(shard)
    return str(shard)


def delta_3plus(rules, values, n):
    spectrum = energy_spectrum(normalized_wht(values), n)
    local = np.asarray(extract_rule_fourier(rules)["local_energy_spectrum"], dtype=float)
    o3 = lambda s: float(s[3:].sum()) if s.shape[0] > 3 else 0.0
    return (o3(spectrum) - o3(local)) * 100


def audit_network(name, workers):
    info = ALL_MODELS[name]
    n = len(info["rules"])
    shard_dir = OUT / "shards" / name
    shard_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(name, s, min(s + CHUNK, 2**n), shard_dir / f"{s:07d}.npz") for s in range(0, 2**n, CHUNK)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        list(pool.map(run_chunk, jobs))

    impl, reg, totals = [], [], {}
    for _, s, e, shard in jobs:
        d = np.load(shard)
        impl.append(d["impl"]); reg.append(d["reg"])
        for key, v in json.loads(str(d["counters"])).items():
            totals[key] = totals.get(key, 0) + v
    impl = np.concatenate(impl); reg = np.concatenate(reg)

    control = {}
    table = Path(str(COMMITTED_TABLE).format(name=name))
    if table.exists():
        committed = np.load(table)["target_logits"]
        control = {"against": str(table.relative_to(PROJECT_ROOT)),
                   "exact_equal": bool(np.array_equal(impl, committed)),
                   "n_values_differing": int((impl != committed).sum()),
                   "max_abs_diff": float(np.max(np.abs(impl - committed)))}
    published = json.loads(PUBLISHED.read_text())["per_network_table"][name]["delta_3plus_pp"]
    v_impl, v_reg = impl.mean(axis=1), reg.mean(axis=1)
    d_impl, d_reg = delta_3plus(info["rules"], v_impl, n), delta_3plus(info["rules"], v_reg, n)
    differ = np.abs(impl - reg) > TOLERANCE
    record = {
        "network": name, "n_nodes": n, "n_coalitions": 2**n,
        "counters": totals,
        "n_trajectories_procedures_differ": int(differ.sum()),
        "n_coalitions_procedures_differ": int(differ.any(axis=1).sum()),
        "max_abs_trajectory_diff": float(np.max(np.abs(impl - reg))),
        "max_abs_coalition_value_diff": float(np.max(np.abs(v_impl - v_reg))),
        "delta_3plus_pp": {"published": published, "implemented": d_impl, "registered": d_reg},
        "positive_control_committed_table": control or "no committed table; compare delta_3plus_pp.implemented with published",
    }
    np.savez_compressed(OUT / f"{name}_coalition_values.npz", implemented=v_impl, registered=v_reg)
    (OUT / f"{name}.json").write_text(json.dumps(record, indent=2))
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--networks", nargs="*")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    table = json.loads(PUBLISHED.read_text())["per_network_table"]
    names = args.networks or sorted(table, key=lambda m: table[m]["n"])
    unknown = [m for m in names if m not in table]
    if unknown:
        raise AuditError(f"not paper networks: {unknown}")
    OUT.mkdir(parents=True, exist_ok=True)
    for name in names:
        r = audit_network(name, args.workers)
        pc = r["positive_control_committed_table"]
        pc_s = f"table exact={pc['exact_equal']}" if isinstance(pc, dict) else "no table"
        c = r["counters"]
        print(f"{name:<28} n={r['n_nodes']:<3} cycling {c['n_cycling']:>8}  unclosed {c['n_unclosed_implemented']:>5}  "
              f"period>50 {c['n_period_over_window']:>5}  procedures differ {r['n_trajectories_procedures_differ']:>6}  "
              f"d3+ pub {r['delta_3plus_pp']['published']:+.2f} impl {r['delta_3plus_pp']['implemented']:+.2f} "
              f"reg {r['delta_3plus_pp']['registered']:+.2f}  {pc_s}")

    records = [json.loads((OUT / f"{m}.json").read_text()) for m in sorted(table) if (OUT / f"{m}.json").exists()]
    summary = {
        "question": "Does the implemented Boolean cycle handling give the registered procedure's outputs?",
        "settings": {"n_init": N_INIT, "seed": SEED, "max_steps": MAX_STEPS, "window": WINDOW,
                     "clamp_value": CLAMP_VALUE, "tolerance": TOLERANCE},
        "n_networks_audited": len(records),
        "networks_where_procedures_differ": [r["network"] for r in records if r["n_trajectories_procedures_differ"]],
        "networks_failing_positive_control": [r["network"] for r in records
                                              if isinstance(r["positive_control_committed_table"], dict)
                                              and not r["positive_control_committed_table"]["exact_equal"]],
        "total_unclosed_implemented": sum(r["counters"]["n_unclosed_implemented"] for r in records),
        "total_period_over_window": sum(r["counters"]["n_period_over_window"] for r in records),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {(OUT / 'summary.json').relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
