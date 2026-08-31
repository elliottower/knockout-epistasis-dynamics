"""ODE coalition sweep: convert Boolean GRN rules to Hill-function ODEs.

For each Boolean rule f_i(x_1, ..., x_k) -> {0,1}, the ODE is:
    dx_i/dt = (1/tau) * (sigma(f_i_continuous(x)) - x_i)

where sigma is a sigmoid, f_i_continuous replaces AND/OR/NOT with
smooth approximations:
    NOT(x)   -> 1 - x
    AND(a,b) -> a * b
    OR(a,b)  -> a + b - a*b  (probabilistic OR)

The steady-state of this ODE approximates the Boolean fixed points
when the sigmoid is steep enough. Limit cycles in the Boolean system
may become stable oscillations in the ODE.

For the coalition sweep, knocked-out genes are clamped to 0 (or 1),
and the ODE is integrated from random initial conditions until
convergence or timeout. The output is time-averaged over the last
portion of the trajectory.

Usage:
    uv run python scripts/ode_coalition_sweep.py \
        --model faure_cellcycle --n-init 128 \
        --output results/grn_v2/ode_pilot/faure_cellcycle_ode.json
"""

import argparse
import json
import re
import sys
import time as _time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from tqdm import tqdm


class _SolveTimeout(Exception):
    pass

sys.path.insert(0, str(Path(__file__).parent.parent))

from data_utils import normalized_wht, energy_by_order
from grn_coalition_sweep import BUILTIN_MODELS, identify_input_nodes, extract_rule_fourier

# Import WEB_MODELS from batch scripts
try:
    from scripts.run_batch2b_extra_models import EXTRA_MODELS
except ImportError:
    EXTRA_MODELS = {}

try:
    from scripts.run_batch2_blind_sweep import WEB_MODELS
except ImportError:
    WEB_MODELS = {}

ALL_MODELS = {**BUILTIN_MODELS, **WEB_MODELS, **EXTRA_MODELS}


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def _find_regulators(expression, node_names):
    regs = []
    for name in node_names:
        if re.search(r"\b" + re.escape(name) + r"\b", expression):
            regs.append(name)
    return regs


def boolean_expr_to_continuous(expression, node_names):
    """Convert a Boolean expression string to a continuous function.

    Returns a callable f(x_dict) -> float in [0, 1], where x_dict maps
    node names to continuous values in [0, 1].

    AND -> product, OR -> probabilistic OR, NOT -> 1 - x.
    """
    def make_fn(expr_str):
        tokens = _tokenize(expr_str, node_names)
        tree = _parse_or(tokens, 0)
        return tree[0]

    return make_fn(expression)


def _tokenize(expr, node_names):
    tokens = []
    i = 0
    sorted_names = sorted(node_names, key=len, reverse=True)
    while i < len(expr):
        if expr[i].isspace():
            i += 1
            continue
        if expr[i] == '!':
            tokens.append(('NOT', None))
            i += 1
            continue
        if expr[i] == '&':
            tokens.append(('AND', None))
            i += 1
            continue
        if expr[i] == '|':
            tokens.append(('OR', None))
            i += 1
            continue
        if expr[i] == '(':
            tokens.append(('LPAREN', None))
            i += 1
            continue
        if expr[i] == ')':
            tokens.append(('RPAREN', None))
            i += 1
            continue
        matched = False
        for name in sorted_names:
            if expr[i:i+len(name)] == name:
                if i + len(name) < len(expr) and (expr[i+len(name)].isalnum() or expr[i+len(name)] == '_'):
                    continue
                tokens.append(('VAR', name))
                i += len(name)
                matched = True
                break
        if not matched:
            i += 1
    return tokens


def _parse_or(tokens, pos):
    left_fn, pos = _parse_and(tokens, pos)
    while pos < len(tokens) and tokens[pos][0] == 'OR':
        pos += 1
        right_fn, pos = _parse_and(tokens, pos)
        prev_left = left_fn
        prev_right = right_fn
        left_fn = lambda x, l=prev_left, r=prev_right: l(x) + r(x) - l(x) * r(x)
    return left_fn, pos


def _parse_and(tokens, pos):
    left_fn, pos = _parse_not(tokens, pos)
    while pos < len(tokens) and tokens[pos][0] == 'AND':
        pos += 1
        right_fn, pos = _parse_not(tokens, pos)
        prev_left = left_fn
        prev_right = right_fn
        left_fn = lambda x, l=prev_left, r=prev_right: l(x) * r(x)
    return left_fn, pos


def _parse_not(tokens, pos):
    if pos < len(tokens) and tokens[pos][0] == 'NOT':
        pos += 1
        inner_fn, pos = _parse_atom(tokens, pos)
        return lambda x, f=inner_fn: 1.0 - f(x), pos
    return _parse_atom(tokens, pos)


def _parse_atom(tokens, pos):
    if pos < len(tokens) and tokens[pos][0] == 'VAR':
        name = tokens[pos][1]
        return lambda x, n=name: x[n], pos + 1
    if pos < len(tokens) and tokens[pos][0] == 'LPAREN':
        pos += 1
        fn, pos = _parse_or(tokens, pos)
        if pos < len(tokens) and tokens[pos][0] == 'RPAREN':
            pos += 1
        return fn, pos
    return lambda x: 0.0, pos


def build_ode_system(rules, node_names, clamp_mask, clamp_value, tau=1.0,
                     hill_n=10.0, hill_k=0.5, deadline=None):
    """Build an ODE system from Boolean rules.

    Returns a function dydt = f(t, y) suitable for scipy.integrate.solve_ivp.
    deadline: list with one float element [wall_time_limit]. Checked every step;
              raises _SolveTimeout if exceeded. Use a mutable container so the
              caller can reset it between solve_ivp calls.
    """
    n = len(node_names)
    continuous_fns = []
    for node in node_names:
        expr = rules[node]
        fn = boolean_expr_to_continuous(expr, node_names)
        continuous_fns.append(fn)

    k_n = hill_k**hill_n
    unclamped = [(i, continuous_fns[i]) for i in range(n) if not clamp_mask[i]]
    clamped = np.array([i for i in range(n) if clamp_mask[i]], dtype=np.intp)
    x_dict = {name: 0.0 for name in node_names}

    def dydt(t, y):
        if deadline is not None and _time.monotonic() > deadline[0]:
            raise _SolveTimeout()
        for i in range(n):
            v = y[i]
            if v < 0.0:
                v = 0.0
            elif v > 1.0:
                v = 1.0
            x_dict[node_names[i]] = v
        dy = np.empty(n)
        for i, fn in unclamped:
            target = fn(x_dict)
            target = target**hill_n / (target**hill_n + k_n + 1e-30)
            dy[i] = (target - y[i]) / tau
        if len(clamped) > 0:
            dy[clamped] = (clamp_value - y[clamped]) * 100.0
        return dy

    return dydt


def simulate_ode_output(rules, node_names, clamp_mask, clamp_value,
                        output_indices, init_states, t_max=10.0,
                        t_tail=5.0, tau=1.0, hill_n=10.0, hill_k=0.5,
                        per_solve_timeout=0.1):
    """Simulate ODE dynamics for multiple initial conditions.

    Returns: (n_init,) array of time-averaged output values.
    """
    n = len(node_names)
    n_init = init_states.shape[0]
    deadline = [float('inf')]
    dydt = build_ode_system(rules, node_names, clamp_mask, clamp_value,
                            tau=tau, hill_n=hill_n, hill_k=hill_k,
                            deadline=deadline)

    outputs = np.zeros(n_init)
    for i in range(n_init):
        y0 = init_states[i].astype(np.float64)
        y0[clamp_mask] = clamp_value

        deadline[0] = _time.monotonic() + per_solve_timeout
        try:
            sol = solve_ivp(dydt, [0, t_max], y0, method='RK45',
                            max_step=2.0, rtol=1e-3, atol=1e-5,
                            dense_output=True)

            if sol.success:
                t_eval = np.linspace(max(0, t_max - t_tail), t_max, 50)
                y_tail = sol.sol(t_eval)

                out_val = 0.0
                for oi in output_indices:
                    out_val += np.mean(np.clip(y_tail[oi], 0, 1))
                out_val /= len(output_indices)
                outputs[i] = out_val
            else:
                outputs[i] = 0.0
        except _SolveTimeout:
            outputs[i] = 0.0
        except Exception:
            outputs[i] = 0.0

    return outputs


def sweep_coalitions_ode(rules, output_nodes, n_init=128, seed=42,
                         clamp_value=0, t_max=10.0, t_tail=5.0,
                         tau=1.0, hill_n=10.0, hill_k=0.5,
                         exclude_inputs=False,
                         checkpoint_path=None, checkpoint_every=1000,
                         on_checkpoint=None,
                         per_solve_timeout=0.1,
                         coalition_range=None):
    """Run exhaustive coalition sweep under ODE dynamics.

    Args:
        checkpoint_path: If set, save partial results as NPZ every
            checkpoint_every coalitions, and resume from this file if it exists.
        checkpoint_every: How often to save (in coalitions computed).
        on_checkpoint: Optional callback called after each checkpoint save,
            e.g. to commit a volume.
        coalition_range: Optional (start, end) tuple to sweep only a subset.
    """
    node_names = list(rules.keys())
    n_total = len(node_names)
    name_to_idx = {name: i for i, name in enumerate(node_names)}
    output_indices = [name_to_idx[o] for o in output_nodes]

    input_nodes = identify_input_nodes(rules) if exclude_inputs else []
    input_indices = [name_to_idx[inp] for inp in input_nodes]
    player_indices = [i for i in range(n_total) if i not in input_indices]
    player_names = [node_names[i] for i in player_indices]
    n_players = len(player_indices)
    N = 2**n_players

    rng = np.random.default_rng(seed)
    init_states = rng.random(size=(n_init, n_total))

    range_start = 0 if coalition_range is None else coalition_range[0]
    range_end = N if coalition_range is None else coalition_range[1]

    start_coalition = range_start
    values = np.zeros((range_end - range_start, n_init))

    if checkpoint_path is not None and Path(checkpoint_path).exists():
        ckpt = np.load(checkpoint_path)
        done = int(ckpt["completed"])
        values[:done] = ckpt["values"][:done]
        start_coalition = range_start + done
        print(f"  Resuming from checkpoint: {done}/{range_end - range_start} coalitions done")

    count = start_coalition - range_start
    for coalition in tqdm(range(start_coalition, range_end),
                          initial=count, total=range_end - range_start,
                          desc=f"ODE coalitions (n={n_players})"):
        clamp_mask = np.zeros(n_total, dtype=bool)
        for bit_pos, node_idx in enumerate(player_indices):
            if not (coalition & (1 << bit_pos)):
                clamp_mask[node_idx] = True

        outputs = simulate_ode_output(
            rules, node_names, clamp_mask, clamp_value,
            output_indices, init_states,
            t_max=t_max, t_tail=t_tail, tau=tau,
            hill_n=hill_n, hill_k=hill_k,
            per_solve_timeout=per_solve_timeout,
        )
        local_idx = coalition - range_start
        values[local_idx] = outputs
        count += 1

        if checkpoint_path is not None and count % checkpoint_every == 0:
            np.savez_compressed(
                checkpoint_path,
                values=values,
                completed=count,
                n_players=n_players,
                player_names=player_names,
                node_names=node_names,
                range_start=range_start,
                range_end=range_end,
            )
            if on_checkpoint is not None:
                on_checkpoint()

    return {
        "values": values,
        "node_names": node_names,
        "player_names": player_names,
        "n_players": n_players,
        "n_init": n_init,
        "range_start": range_start,
        "range_end": range_end,
    }


def score_ode_result(values, n_players, node_names):
    """Compute composition gap from ODE coalition sweep."""
    v_mean = values.mean(axis=1)
    w = normalized_wht(v_mean)
    energy = energy_by_order(w, n_players)
    total = energy.sum()
    if total > 0:
        spectrum = (energy / total).tolist()
    else:
        spectrum = energy.tolist()

    o3plus = sum(spectrum[3:]) if len(spectrum) > 3 else 0.0
    return {
        "global_spectrum": spectrum,
        "global_o3plus": o3plus,
        "n_unique_values": len(np.unique(np.round(v_mean, 6))),
        "v_mean_range": [float(v_mean.min()), float(v_mean.max())],
        "v_mean_std": float(v_mean.std()),
    }


def main():
    parser = argparse.ArgumentParser(description="ODE coalition sweep for Boolean GRNs")
    parser.add_argument("--model", required=True, help="Model name")
    parser.add_argument("--output", required=True, help="Output JSON path")
    parser.add_argument("--n-init", type=int, default=128)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--t-max", type=float, default=50.0)
    parser.add_argument("--t-tail", type=float, default=10.0)
    parser.add_argument("--tau", type=float, default=1.0)
    parser.add_argument("--hill-n", type=float, default=10.0)
    parser.add_argument("--hill-k", type=float, default=0.5)
    args = parser.parse_args()

    if args.model not in ALL_MODELS:
        print(f"Unknown model: {args.model}")
        print(f"Available: {sorted(ALL_MODELS.keys())}")
        sys.exit(1)

    model_info = ALL_MODELS[args.model]
    rules = model_info["rules"]
    output_nodes = model_info["output_nodes"]

    print(f"[{timestamp()}] ODE pilot: {args.model}")
    print(f"  {model_info['description']}")
    print(f"  n_init={args.n_init}, t_max={args.t_max}, hill_n={args.hill_n}")

    result = sweep_coalitions_ode(
        rules, output_nodes,
        n_init=args.n_init, seed=args.seed,
        t_max=args.t_max, t_tail=args.t_tail,
        tau=args.tau, hill_n=args.hill_n, hill_k=args.hill_k,
    )

    scores = score_ode_result(result["values"], result["n_players"],
                              result["player_names"])

    # Compute local rule energy (same for Boolean and ODE — rules don't change)
    rf = extract_rule_fourier(rules)
    local_spectrum_raw = rf.get("local_energy_spectrum", [])
    local_o3plus = sum(local_spectrum_raw[3:]) if len(local_spectrum_raw) > 3 else 0.0

    ode_delta = scores["global_o3plus"] - local_o3plus

    # Load Boolean results for comparison
    bool_path = Path(__file__).parent.parent / "results" / "grn_v2" / "merged_all_27_analysis.json"
    bool_delta = None
    bool_global_o3plus = None
    if bool_path.exists():
        with open(bool_path) as f:
            merged = json.load(f)
        for m in merged["all_models"]:
            if m["model"] == args.model:
                bool_delta = m.get("delta_o3plus")
                bool_global_o3plus = m.get("global_o3plus")
                break

    sign_preserved = None
    if bool_delta is not None:
        bool_sign = "creation" if bool_delta > 0.005 else ("destruction" if bool_delta < -0.005 else "null")
        ode_sign = "creation" if ode_delta > 0.005 else ("destruction" if ode_delta < -0.005 else "null")
        sign_preserved = (bool_sign == ode_sign)

    output = {
        "model": args.model,
        "method": "ode_hill",
        "description": model_info["description"],
        "n_players": result["n_players"],
        "n_init": args.n_init,
        "ode_params": {
            "t_max": args.t_max,
            "t_tail": args.t_tail,
            "tau": args.tau,
            "hill_n": args.hill_n,
            "hill_k": args.hill_k,
        },
        "local_o3plus": local_o3plus,
        "ode_global_o3plus": scores["global_o3plus"],
        "ode_delta_o3plus": ode_delta,
        "boolean_delta_o3plus": bool_delta,
        "boolean_global_o3plus": bool_global_o3plus,
        "sign_preserved": sign_preserved,
        "ode_global_spectrum": scores["global_spectrum"],
        "n_unique_values": scores["n_unique_values"],
        "v_mean_range": scores["v_mean_range"],
        "v_mean_std": scores["v_mean_std"],
        "timestamp": timestamp(),
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n[{timestamp()}] Saved: {out_path}")

    print(f"\n  Local o3+:       {local_o3plus:.4f}")
    print(f"  ODE global o3+:  {scores['global_o3plus']:.4f}")
    print(f"  ODE delta o3+:   {ode_delta:+.4f}")
    if bool_delta is not None:
        print(f"  Boolean delta:   {bool_delta:+.4f}")
        print(f"  Sign preserved:  {sign_preserved}")


if __name__ == "__main__":
    main()
