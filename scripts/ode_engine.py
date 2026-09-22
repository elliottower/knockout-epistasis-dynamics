"""Continuous dynamics for the composition-gap ODE arms.

This engine replaces the solver loop of `scripts/ode_coalition_sweep.py`. That loop had four
defects, found in review on 2026-09-21:

  1. A solver failure, a timeout and any exception were all recorded as output 0.0, which
     cannot be told apart from a genuine "off" steady state.
  2. The timeout was 0.1 s of wall-clock time, so whether a solve succeeded depended on the
     machine and its load.
  3. Nothing recorded how often either happened.
  4. `solve_ivp`'s success flag, which means only that integration reached t_max, was read
     as convergence.

Here, every trajectory carries a status and a failure is never turned into a value. The
wall-clock timeout is replaced by a deterministic cap on right-hand-side evaluations.
Solver attempts run under declared conditions: a larger budget only after the cap was hit,
a stiff method after any failure. When a classifier is supplied, a trajectory is accepted
only once its tail is classified as a fixed point or a sustained oscillation; otherwise the
horizon is extended deterministically, and a trajectory still unclassified at the last
horizon fails.

Three continuous extensions of each Boolean rule are available, all built from the truth
tables the Boolean arm uses (`grn_coalition_sweep.compile_network`):

  hillcube_normalized       Wittmann et al. 2009, eq. 5. A normalized Hill function on
                            each input, then the multilinear extension of the truth table.
                            Determined by the truth table, so it does not depend on how the
                            rule is written, and it equals the Boolean rule at every vertex.
                            Compiled to straight-line code by Shannon expansion.
  operatorwise_normalized   The legacy substitution (AND -> ab, OR -> a + b - ab,
                            NOT -> 1 - a) with one Hill sigmoid on the rule's output,
                            divided by its value at 1 so that 1 maps to 1.
  operatorwise_legacy       The legacy construction exactly, unnormalized, including its
                            1e-30 guard and its stiff restoring term for clamped nodes.
                            Kept only to audit the published arm.

For the two current constructions, clamped nodes are removed from the integrated state and
their value substituted directly, so a clamp is exact and adds no stiffness.

Measurement is kept apart from interpretation where it can be. Classification thresholds
and extension horizons come from a `ClassifierConfig` the registration freezes; the engine
records the diagnostics each classification rests on.

An accepted oscillation also carries two summaries that do not depend on the phase at which
its tail starts (`cycle_summary`): its period, and each node's mean over a whole number of
periods. `scripts/attractors.py` uses them to tell limit cycles apart.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Literal

import numpy as np
import scipy
from scipy.integrate import solve_ivp

import grn_coalition_sweep
import scripts.ode_coalition_sweep as legacy_sweep
from grn_coalition_sweep import compile_network
from scripts.ode_coalition_sweep import _tokenize

Construction = Literal["hillcube_normalized", "operatorwise_normalized", "operatorwise_legacy"]
CONSTRUCTIONS: tuple[Construction, ...] = (
    "hillcube_normalized", "operatorwise_normalized", "operatorwise_legacy")

Status = Literal["ok", "eval_cap_exceeded", "solver_failed", "did_not_reach_end", "nonfinite", "unclassified"]
STATUSES: tuple[Status, ...] = (
    "ok", "eval_cap_exceeded", "solver_failed", "did_not_reach_end", "nonfinite", "unclassified")

TrajectoryClass = Literal["none", "fixed", "oscillatory", "unclassified"]
# "none" means no classifier was supplied (the legacy audit), or the trajectory failed.
TRAJECTORY_CLASSES: tuple[TrajectoryClass, ...] = ("none", "fixed", "oscillatory", "unclassified")

When = Literal["first", "after_cap", "after_failure"]

# The legacy engine's Hill guard, reproduced exactly for the audit construction only.
LEGACY_HILL_GUARD = 1e-30


class ODEEngineError(Exception):
    """Base class for errors raised by this engine."""


class EvaluationCapExceeded(ODEEngineError):
    """The right-hand side was evaluated more times than the attempt allows."""


class NonFiniteDerivative(ODEEngineError):
    """The right-hand side produced NaN or inf. Raised at the evaluation itself, because the
    solvers do not agree on what to do with one: RK45 shrinks its step until it fails, and
    Radau raises ValueError from its LU factorization."""


class RuleCompileError(ODEEngineError):
    """A rule could not be compiled into a continuous function that matches its truth table."""


class IncompleteSweepError(ODEEngineError):
    """A sweep contains coalitions without an accepted value, so it cannot be scored."""


@dataclass(frozen=True)
class SolverAttempt:
    method: str
    max_rhs_evals: int
    when: When = "first"


@dataclass(frozen=True)
class SolverConfig:
    """How each trajectory is integrated. Every field is declared by the caller."""

    attempts: tuple[SolverAttempt, ...]
    rtol: float
    atol: float
    max_step: float
    t_max: float
    t_tail: float
    n_tail_samples: int

    def __post_init__(self) -> None:
        if not self.attempts:
            raise ODEEngineError("a SolverConfig needs at least one attempt")
        if self.attempts[0].when != "first":
            raise ODEEngineError("the first attempt must have when='first'")
        if any(a.when == "first" for a in self.attempts[1:]):
            raise ODEEngineError("only the first attempt may have when='first'")
        if any(a.max_rhs_evals < 1 for a in self.attempts):
            raise ODEEngineError("every attempt needs max_rhs_evals >= 1")
        if not 0 < self.t_tail <= self.t_max:
            raise ODEEngineError(f"t_tail {self.t_tail} must lie in (0, t_max={self.t_max}]")
        if self.n_tail_samples < 2:
            raise ODEEngineError("n_tail_samples must be at least 2")
        if self.rtol <= 0 or self.atol <= 0 or self.max_step <= 0:
            raise ODEEngineError("rtol, atol and max_step must be positive")


# The published arm's settings, as executed by modal_ode_sweep.py, minus the wall-clock
# timeout. One attempt, no fallback, and a cap high enough never to bind on these systems.
LEGACY_SOLVER = SolverConfig(
    attempts=(SolverAttempt("RK45", 10_000_000),),
    rtol=1e-3, atol=1e-5, max_step=2.0,
    t_max=30.0, t_tail=10.0, n_tail_samples=50,
)


@dataclass(frozen=True)
class ClassifierConfig:
    """Thresholds for accepting a trajectory's tail, and the horizons to extend through."""

    fixed_tail_range_max: float
    fixed_derivative_max: float
    oscillatory_tail_range_min: float
    halves_envelope_tol: float
    bounded_margin: float
    extension_horizons: tuple[float, ...]

    def __post_init__(self) -> None:
        if not self.fixed_tail_range_max < self.oscillatory_tail_range_min:
            raise ODEEngineError("the fixed-point range must lie below the oscillatory range")
        if list(self.extension_horizons) != sorted(set(self.extension_horizons)):
            raise ODEEngineError("extension horizons must be strictly increasing")


@dataclass(frozen=True)
class DynamicsConfig:
    construction: Construction
    hill_n: float
    hill_k: float = 0.5
    tau: float = 1.0
    clamp_stiffness: float = 100.0  # legacy construction only

    def __post_init__(self) -> None:
        if self.construction not in CONSTRUCTIONS:
            raise ODEEngineError(f"unknown construction {self.construction!r}")
        if self.hill_n <= 0 or not 0 < self.hill_k < 1 or self.tau <= 0:
            raise ODEEngineError(f"need hill_n > 0, 0 < hill_k < 1, tau > 0; got "
                                 f"{self.hill_n}, {self.hill_k}, {self.tau}")


# ---------------------------------------------------------------------------------------
# Rule compilation
# ---------------------------------------------------------------------------------------

def _exec_function(lines: list[str], result: str, label: str) -> Callable[[np.ndarray], float]:
    source = "def f(x):\n" + ("\n".join(lines) + "\n" if lines else "") + f"    return {result}\n"
    namespace: dict = {}
    exec(compile(source, label, "exec"), namespace)  # source built from parsed tokens and constants only
    return namespace["f"]


def compile_operatorwise(expression: str, node_names: list[str]) -> tuple[Callable[[np.ndarray], float], bool]:
    """Compile a rule to f(x) -> float with the legacy engine's exact arithmetic.

    The legacy parser (`_parse_or`, `_parse_and`, `_parse_not`, `_parse_atom` in
    ode_coalition_sweep.py) builds nested closures over a dict. This mirrors its recursive
    descent and associativity, so each operation happens in the same order on the same
    operands, but reads an array and evaluates each subexpression once. Tokenization is the
    legacy function itself.

    Returns the function and whether the legacy parser's silent fallback was reached: an
    unexpected token makes `_parse_atom` return the constant 0.0 without consuming it.
    """
    tokens = _tokenize(expression, node_names)
    index = {name: i for i, name in enumerate(node_names)}
    lines: list[str] = []
    counter = [0]
    fallback = [False]

    def temp(expr: str) -> str:
        counter[0] += 1
        name = f"t{counter[0]}"
        lines.append(f"    {name} = {expr}")
        return name

    def parse_or(pos: int) -> tuple[str, int]:
        left, pos = parse_and(pos)
        while pos < len(tokens) and tokens[pos][0] == "OR":
            right, pos = parse_and(pos + 1)
            left = temp(f"{left} + {right} - {left} * {right}")
        return left, pos

    def parse_and(pos: int) -> tuple[str, int]:
        left, pos = parse_not(pos)
        while pos < len(tokens) and tokens[pos][0] == "AND":
            right, pos = parse_not(pos + 1)
            left = temp(f"{left} * {right}")
        return left, pos

    def parse_not(pos: int) -> tuple[str, int]:
        if pos < len(tokens) and tokens[pos][0] == "NOT":
            inner, pos = parse_atom(pos + 1)
            return temp(f"1.0 - {inner}"), pos
        return parse_atom(pos)

    def parse_atom(pos: int) -> tuple[str, int]:
        if pos < len(tokens) and tokens[pos][0] == "VAR":
            return temp(f"x[{index[tokens[pos][1]]}]"), pos + 1
        if pos < len(tokens) and tokens[pos][0] == "LPAREN":
            inner, pos = parse_or(pos + 1)
            if pos < len(tokens) and tokens[pos][0] == "RPAREN":
                pos += 1
            return inner, pos
        fallback[0] = True
        return temp("0.0"), pos

    result, _ = parse_or(0)
    return _exec_function(lines, result, f"<operatorwise {expression!r}>"), fallback[0]


def compile_hillcube(reg_indices: np.ndarray, table: np.ndarray, hill_n: float,
                     hill_k: float) -> Callable[[np.ndarray], float]:
    """Compile a rule's normalized HillCube to straight-line code.

    Each deduplicated regulator's value passes through the normalized Hill function
    u = z^n (1 + K^n) / (z^n + K^n); the multilinear extension of the truth table is then
    evaluated by Shannon expansion on the highest remaining bit,
        ML(T) = (1 - u) * ML(T | bit = 0) + u * ML(T | bit = 1),
    folding any sub-table that is constant. Bit j of the table index is regulator j.
    """
    kn = hill_k ** hill_n
    c1 = 1.0 + kn
    lines: list[str] = []
    counter = [0]

    def temp(expr: str) -> str:
        counter[0] += 1
        name = f"t{counter[0]}"
        lines.append(f"    {name} = {expr}")
        return name

    us = []
    for r in reg_indices:
        zn = temp(f"x[{int(r)}] ** {hill_n!r}")
        us.append(temp(f"{zn} * {c1!r} / ({zn} + {kn!r})"))

    def expand(sub: np.ndarray, j: int) -> str:
        if np.all(sub == sub[0]):
            return repr(float(sub[0]))
        half = sub.shape[0] // 2
        low, high = expand(sub[:half], j - 1), expand(sub[half:], j - 1)
        u = us[j]
        if low == "0.0":
            return u if high == "1.0" else temp(f"{u} * {high}")
        if high == "0.0":
            return temp(f"1.0 - {u}") if low == "1.0" else temp(f"(1.0 - {u}) * {low}")
        return temp(f"(1.0 - {u}) * {low} + {u} * {high}")

    result = expand(np.asarray(table, dtype=float), len(reg_indices) - 1)
    return _exec_function(lines, result, "<hillcube>")


@dataclass(frozen=True)
class CompiledNetwork:
    node_names: tuple[str, ...]
    output_indices: tuple[int, ...]
    expressions: tuple[str, ...]
    reg_indices: tuple[np.ndarray, ...]
    truth_tables: tuple[np.ndarray, ...]
    bit_masks: tuple[np.ndarray, ...]
    operatorwise: tuple[Callable[[np.ndarray], float], ...]
    operatorwise_fallback: tuple[str, ...]
    operatorwise_vertex_mismatch: tuple[str, ...]


def compile_model(rules: dict[str, str], output_nodes: list[str]) -> CompiledNetwork:
    compiled, node_names = compile_network(rules)
    name_to_idx = {n: i for i, n in enumerate(node_names)}
    if not output_nodes:
        raise RuleCompileError("a network needs at least one output node")
    if len(set(output_nodes)) != len(output_nodes):
        raise RuleCompileError(f"duplicate output nodes: {output_nodes}")
    unknown = [o for o in output_nodes if o not in name_to_idx]
    if unknown:
        raise RuleCompileError(f"output nodes not in the network: {unknown}")

    reg_indices, tables, masks, opwise, fallback, mismatch = [], [], [], [], [], []
    for node, regs, table in compiled:
        k = len(regs)
        reg_indices.append(np.asarray(regs, dtype=np.intp))
        tables.append(table.astype(float))
        masks.append(np.array([[(b >> j) & 1 for b in range(2**k)] for j in range(k)], dtype=bool).reshape(k, 2**k))
        fn, used_fallback = compile_operatorwise(rules[node], node_names)
        opwise.append(fn)
        if used_fallback:
            fallback.append(node)
        # The operator-wise function must reproduce the truth table at every vertex.
        x = np.zeros(len(node_names))
        for b in range(2**k):
            for j, r in enumerate(regs):
                x[r] = float((b >> j) & 1)
            if abs(fn(x) - table[b]) > 1e-12:
                mismatch.append(node)
                break
            x[:] = 0.0
    return CompiledNetwork(
        node_names=tuple(node_names),
        output_indices=tuple(name_to_idx[o] for o in output_nodes),
        expressions=tuple(rules[node] for node in node_names),
        reg_indices=tuple(reg_indices),
        truth_tables=tuple(tables),
        bit_masks=tuple(masks),
        operatorwise=tuple(opwise),
        operatorwise_fallback=tuple(fallback),
        operatorwise_vertex_mismatch=tuple(mismatch),
    )


def model_digest(net: CompiledNetwork) -> str:
    """SHA-256 of everything that defines a network: node order, output nodes, and each rule's
    text, regulators and truth table."""
    canonical = {
        "node_names": list(net.node_names),
        "output_nodes": [net.node_names[i] for i in net.output_indices],
        "rules": [{"expression": e, "regulators": [int(r) for r in regs], "table": [float(v) for v in table]}
                  for e, regs, table in zip(net.expressions, net.reg_indices, net.truth_tables)],
    }
    return hashlib.sha256(json.dumps(canonical, sort_keys=True).encode()).hexdigest()


# The modules whose source determines a trajectory: this engine, the rule compiler and the
# legacy tokenizer. Their digest is part of every run's configuration.
TRAJECTORY_CODE = (Path(__file__), Path(grn_coalition_sweep.__file__), Path(legacy_sweep.__file__))


def code_digest() -> str:
    h = hashlib.sha256()
    for path in TRAJECTORY_CODE:
        h.update(path.name.encode())
        h.update(path.read_bytes())
    return h.hexdigest()


# ---------------------------------------------------------------------------------------
# Continuous extensions
# ---------------------------------------------------------------------------------------

def hill(z: np.ndarray | float, n: float, k: float) -> np.ndarray | float:
    zn = np.power(z, n)
    return zn / (zn + k**n)


def hill_normalized(z: np.ndarray | float, n: float, k: float) -> np.ndarray | float:
    """Hill function divided by its value at 1, so that 0 -> 0 and 1 -> 1 for every n.

    K is the parameter of the underlying Hill function, not the half-maximum of this one:
    h(K) = (1 + K^n) / 2, and the half-maximum lies at K (1 + 2 K^n)^(-1/n).
    """
    zn = np.power(z, n)
    return zn * (1.0 + k**n) / (zn + k**n)


def multilinear(u: np.ndarray, table: np.ndarray, masks: np.ndarray) -> float:
    """Multilinear extension of `table` at point u, with bit j of the index set by input j."""
    w = np.ones(table.shape[0])
    for j in range(u.shape[0]):
        w *= np.where(masks[j], u[j], 1.0 - u[j])
    return float(w @ table)


def rule_targets_reference(net: CompiledNetwork, dyn: DynamicsConfig, x: np.ndarray,
                           nodes: np.ndarray) -> np.ndarray:
    """Plain NumPy reference for each construction. The compiled path is tested against it."""
    out = np.empty(nodes.shape[0])
    n, k = dyn.hill_n, dyn.hill_k
    if dyn.construction == "hillcube_normalized":
        for m, i in enumerate(nodes):
            regs = net.reg_indices[i]
            if regs.shape[0] == 0:
                out[m] = net.truth_tables[i][0]
                continue
            out[m] = multilinear(hill_normalized(x[regs], n, k), net.truth_tables[i], net.bit_masks[i])
    elif dyn.construction == "operatorwise_normalized":
        for m, i in enumerate(nodes):
            out[m] = hill_normalized(net.operatorwise[i](x), n, k)
    else:
        kn = k**n
        for m, i in enumerate(nodes):
            t = net.operatorwise[i](x)
            out[m] = t**n / (t**n + kn + LEGACY_HILL_GUARD)
    return out


@dataclass(frozen=True)
class Prepared:
    """A network with one construction's per-node target functions, compiled once."""

    net: CompiledNetwork
    dyn: DynamicsConfig
    targets: tuple[Callable[[np.ndarray], float], ...]


def prepare(net: CompiledNetwork, dyn: DynamicsConfig) -> Prepared:
    n, k = dyn.hill_n, dyn.hill_k
    kn, c1 = k**n, 1.0 + k**n
    fns: list[Callable[[np.ndarray], float]] = []
    for i in range(len(net.node_names)):
        if dyn.construction == "hillcube_normalized":
            fns.append(compile_hillcube(net.reg_indices[i], net.truth_tables[i], n, k))
        elif dyn.construction == "operatorwise_normalized":
            def f(x, g=net.operatorwise[i]):
                tn = g(x) ** n
                return tn * c1 / (tn + kn)
            fns.append(f)
        else:
            def f(x, g=net.operatorwise[i]):
                t = g(x)
                return t**n / (t**n + kn + LEGACY_HILL_GUARD)
            fns.append(f)
    return Prepared(net=net, dyn=dyn, targets=tuple(fns))


def rule_targets(prep: Prepared, x: np.ndarray, nodes: np.ndarray) -> np.ndarray:
    """Target value of each listed node's rule at state x, which lies in [0, 1] when x does."""
    return np.array([prep.targets[i](x) for i in nodes], dtype=float)


# ---------------------------------------------------------------------------------------
# Right-hand sides
# ---------------------------------------------------------------------------------------

@dataclass
class _EvalCounter:
    cap: int
    n: int = 0


@dataclass(frozen=True)
class System:
    """The integrated system for one coalition: which components move, and how to rebuild the
    full state from them. The legacy construction integrates every node; the others integrate
    only the free nodes."""

    free: np.ndarray
    full_state_integrated: bool
    template: np.ndarray

    def expand(self, y: np.ndarray) -> np.ndarray:
        """Full state from the integrated components; y may carry a trailing sample axis."""
        if self.full_state_integrated:
            return y
        shape = (self.template.shape[0],) + y.shape[1:]
        full = np.empty(shape)
        full[...] = self.template.reshape((-1,) + (1,) * (y.ndim - 1))
        full[self.free] = y
        return full


def make_system(prep: Prepared, clamp_mask: np.ndarray, clamp_value: float) -> System:
    template = np.zeros(len(prep.net.node_names))
    template[clamp_mask] = clamp_value
    legacy = prep.dyn.construction == "operatorwise_legacy"
    return System(free=np.flatnonzero(~clamp_mask), full_state_integrated=legacy, template=template)


def make_rhs(prep: Prepared, system: System, clamp_mask: np.ndarray, clamp_value: float,
             counter: _EvalCounter | None) -> Callable[[float, np.ndarray], np.ndarray]:
    tau = prep.dyn.tau
    free = system.free
    targets = [prep.targets[i] for i in free]

    if system.full_state_integrated:
        clamped = np.flatnonzero(clamp_mask)
        stiffness = prep.dyn.clamp_stiffness

        def rhs(t: float, y: np.ndarray) -> np.ndarray:
            if counter is not None:
                counter.n += 1
                if counter.n > counter.cap:
                    raise EvaluationCapExceeded(f"{counter.cap} right-hand-side evaluations")
            x = np.clip(y, 0.0, 1.0)
            dy = np.empty_like(y)
            dy[free] = (np.array([f(x) for f in targets], dtype=float) - y[free]) / tau
            dy[clamped] = (clamp_value - y[clamped]) * stiffness
            if not np.all(np.isfinite(dy)):
                raise NonFiniteDerivative(f"non-finite derivative at t={t}")
            return dy
        return rhs

    x = system.template.copy()

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        if counter is not None:
            counter.n += 1
            if counter.n > counter.cap:
                raise EvaluationCapExceeded(f"{counter.cap} right-hand-side evaluations")
        x[free] = np.clip(y, 0.0, 1.0)
        dy = (np.array([f(x) for f in targets], dtype=float) - y) / tau
        if not np.all(np.isfinite(dy)):
            raise NonFiniteDerivative(f"non-finite derivative at t={t}")
        return dy
    return rhs


# ---------------------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------------------

@dataclass
class TailStats:
    tail_range: float
    final_derivative: float
    # Largest disagreement, over moving nodes, between the two halves of the tail in their
    # maxima and minima, each as a fraction of that node's range over the whole tail. A node
    # whose range is within the fixed-point bound counts as static and is not compared, since
    # its range is solver jitter. A sustained oscillation spanning at least one period per
    # half repeats its envelope. Half-means are not compared: each half cuts through a partial
    # period, which moves its mean by a fraction of the amplitude even on an exact limit cycle.
    halves_envelope_diff: float
    within_bounds: bool


def tail_stats(tail: np.ndarray, final_derivative: float, margin: float, static_below: float) -> TailStats:
    half = tail.shape[1] // 2
    first, second = tail[:, :half], tail[:, half:]
    node_range = tail.max(axis=1) - tail.min(axis=1)
    envelope = np.maximum(np.abs(first.max(axis=1) - second.max(axis=1)),
                          np.abs(first.min(axis=1) - second.min(axis=1)))
    moving = node_range > static_below
    return TailStats(
        tail_range=float(np.max(node_range)),
        final_derivative=final_derivative,
        halves_envelope_diff=float(np.max(envelope[moving] / node_range[moving])) if moving.any() else 0.0,
        within_bounds=bool(np.all(tail >= -margin) and np.all(tail <= 1.0 + margin)),
    )


def cycle_summary(t: np.ndarray, tail: np.ndarray) -> tuple[float, np.ndarray]:
    """Period of an oscillating tail, and each node's mean over a whole number of periods.

    tail is (n_nodes, n_samples) at times t, clipped to [0, 1]. Each node whose range is at
    least half the largest node range gives a period estimate: the mean interval between its
    upward crossings of its own midrange, interpolated linearly between samples. A node that
    rises once per cycle gives the period and one that rises twice gives half of it, so the
    period is the largest estimate. The mean of every node is taken over the samples between
    the first and last crossing of the node giving that estimate. Both are NaN if no such node
    crosses its midrange upward at least twice.
    """
    node_range = tail.max(axis=1) - tail.min(axis=1)
    best_period, window = float("nan"), None
    for i in np.flatnonzero(node_range >= 0.5 * node_range.max()):
        y = tail[i]
        mid = 0.5 * (y.max() + y.min())
        up = np.flatnonzero((y[:-1] < mid) & (y[1:] >= mid))
        if up.size < 2:
            continue
        times = t[up] + (mid - y[up]) / (y[up + 1] - y[up]) * (t[up + 1] - t[up])
        period = float((times[-1] - times[0]) / (up.size - 1))
        if window is None or period > best_period:
            best_period, window = period, (up[0] + 1, up[-1] + 1)
    if window is None:
        return float("nan"), np.full(tail.shape[0], np.nan)
    return best_period, tail[:, window[0]:window[1]].mean(axis=1)


def classify(stats: TailStats, cfg: ClassifierConfig) -> TrajectoryClass:
    if stats.tail_range <= cfg.fixed_tail_range_max and stats.final_derivative <= cfg.fixed_derivative_max:
        return "fixed"
    if (stats.tail_range > cfg.oscillatory_tail_range_min and stats.within_bounds
            and stats.halves_envelope_diff <= cfg.halves_envelope_tol):
        return "oscillatory"
    return "unclassified"


# ---------------------------------------------------------------------------------------
# Trajectories
# ---------------------------------------------------------------------------------------


@dataclass
class TrajectoryResult:
    status: Status
    output: float
    klass: TrajectoryClass
    attempt: int
    method: str
    nfev: int
    message: str
    horizon_index: int
    stats: TailStats | None  # kept for "unclassified" failures too, as the diagnostic
    # One entry per solver attempt run, across all horizons: horizon, attempt index, method,
    # condition, cap, right-hand-side evaluations, status, solver message.
    attempt_log: list[dict] = field(default_factory=list, repr=False)
    # For an "unclassified" failure only: its tail-window output at the last horizon. Never part
    # of a coalition value; kept so that a labeled sensitivity analysis needs no rerun.
    provisional_output: float = float("nan")
    final_state: np.ndarray | None = field(default=None, repr=False)
    # Per node over the accepting tail, clipped to [0, 1]: its mean, and its envelope. A fixed
    # point has all three equal; trajectories on one limit cycle share an envelope at any phase.
    tail_mean: np.ndarray | None = field(default=None, repr=False)
    tail_min: np.ndarray | None = field(default=None, repr=False)
    tail_max: np.ndarray | None = field(default=None, repr=False)
    tail: np.ndarray | None = field(default=None, repr=False)
    # For an accepted oscillation only (`cycle_summary`): its period, and each node's mean over a
    # whole number of periods. NaN, and None, otherwise.
    period: float = float("nan")
    cycle_mean: np.ndarray | None = field(default=None, repr=False)


def _should_run(when: When, previous: Status | None) -> bool:
    if when == "first":
        return previous is None
    if when == "after_cap":
        return previous == "eval_cap_exceeded"
    return previous is not None and previous != "ok"


def _segment(prep: Prepared, system: System, solver: SolverConfig, clamp_mask: np.ndarray,
             clamp_value: float, t0: float, t1: float, y0: np.ndarray):
    """Integrate [t0, t1] under the declared attempts.

    Returns (status, sol, attempt index, method, log), with one log entry per attempt run.
    Evaluations are counted by the engine's own counter, which sees every call: SciPy's
    `sol.nfev` omits the calls Radau and BDF make to estimate the Jacobian.
    """
    status: Status | None = None
    method, used = solver.attempts[0].method, 0
    log: list[dict] = []
    for a, attempt in enumerate(solver.attempts):
        if not _should_run(attempt.when, status):
            continue
        used, method = a, attempt.method
        counter = _EvalCounter(cap=attempt.max_rhs_evals)
        rhs = make_rhs(prep, system, clamp_mask, clamp_value, counter)
        sol = None
        try:
            sol = solve_ivp(rhs, (t0, t1), y0.copy(), method=attempt.method, max_step=solver.max_step,
                            rtol=solver.rtol, atol=solver.atol, dense_output=True)
        except EvaluationCapExceeded as err:
            status, message = "eval_cap_exceeded", str(err)
        except NonFiniteDerivative as err:
            status, message = "nonfinite", str(err)
        else:
            if not sol.success:
                status, message = "solver_failed", str(sol.message)
            elif abs(sol.t[-1] - t1) > 1e-9 * max(1.0, abs(t1)):
                status, message = "did_not_reach_end", f"stopped at t={sol.t[-1]}"
            elif not np.all(np.isfinite(sol.y)):
                status, message = "nonfinite", "non-finite state during integration"
            else:
                status, message = "ok", str(sol.message)
        log.append({"horizon": t1, "attempt": a, "method": attempt.method, "when": attempt.when,
                    "cap": attempt.max_rhs_evals, "nfev": counter.n, "status": status, "message": message})
        if status == "ok":
            return "ok", sol, used, method, log
    return status, None, used, method, log


def integrate(prep: Prepared, solver: SolverConfig, clamp_mask: np.ndarray, clamp_value: float,
              y0: np.ndarray, classifier: ClassifierConfig | None = None,
              keep_states: bool = False, keep_tail: bool = False) -> TrajectoryResult:
    """Integrate one trajectory. A failure of any kind returns output NaN and its status.

    y0 is a full state whose clamped components equal clamp_value. Without a classifier the
    trajectory is integrated to t_max once and accepted if the solver succeeded, as the
    legacy engine did. With one, it is accepted only when its tail classifies as a fixed point
    or an oscillation, and is otherwise continued to each extension horizon in turn.
    """
    net = prep.net
    system = make_system(prep, clamp_mask, clamp_value)
    horizons = (solver.t_max,) + (classifier.extension_horizons if classifier else ())
    if classifier and any(h <= solver.t_max for h in classifier.extension_horizons):
        raise ODEEngineError("extension horizons must exceed t_max")

    log: list[dict] = []

    def fail(status, a, m, msg, h, stats=None, tail=None, provisional=float("nan")):
        return TrajectoryResult(status=status, output=float("nan"), klass="none", attempt=a, method=m,
                                nfev=sum(e["nfev"] for e in log), message=msg, horizon_index=h, stats=stats,
                                tail=tail.astype(np.float32) if keep_tail and tail is not None else None,
                                attempt_log=log, provisional_output=provisional)

    if system.free.size == 0:  # every node clamped: the state is constant and exact
        template = system.template
        return TrajectoryResult(
            status="ok", output=float(np.mean([template[i] for i in net.output_indices])),
            klass="fixed" if classifier else "none", attempt=0, method="none", nfev=0,
            message="all nodes clamped", horizon_index=0, stats=TailStats(0.0, 0.0, 0.0, True),
            final_state=template.astype(np.float32) if keep_states else None,
            tail_mean=template.astype(np.float32) if keep_states else None,
            tail_min=template.astype(np.float32) if keep_states else None,
            tail_max=template.astype(np.float32) if keep_states else None)

    y = y0 if system.full_state_integrated else y0[system.free]
    t0 = 0.0
    stats, tail, output = None, None, float("nan")
    for h, horizon in enumerate(horizons):
        status, sol, a, method, segment_log = _segment(prep, system, solver, clamp_mask, clamp_value, t0, horizon, y)
        log.extend(segment_log)
        msg = segment_log[-1]["message"]
        if status != "ok":
            return fail(status, a, method, msg, h)
        # The tail is the same fraction of every horizon, sampled at the same density; at t_max
        # it is exactly [t_max - t_tail, t_max] with n_tail_samples points.
        scale = horizon / solver.t_max
        t_eval = np.linspace(horizon - solver.t_tail * scale, horizon, int(round(solver.n_tail_samples * scale)))
        tail = system.expand(sol.sol(t_eval))
        if not np.all(np.isfinite(tail)):
            return fail("nonfinite", a, method, "non-finite state in the tail", h)
        output = 0.0
        for i in net.output_indices:
            output += np.mean(np.clip(tail[i], 0.0, 1.0))
        output = float(output / len(net.output_indices))
        if not np.isfinite(output):
            return fail("nonfinite", a, method, "non-finite output", h)
        y_final = sol.y[:, -1]
        try:
            derivative = float(np.max(np.abs(make_rhs(prep, system, clamp_mask, clamp_value, None)(horizon, y_final))))
        except NonFiniteDerivative as err:
            return fail("nonfinite", a, method, str(err), h)
        stats = (tail_stats(tail, derivative, classifier.bounded_margin, classifier.fixed_tail_range_max)
                 if classifier else tail_stats(tail, derivative, 0.0, 0.0))
        klass: TrajectoryClass = classify(stats, classifier) if classifier else "none"
        if klass != "unclassified":
            clipped = np.clip(tail, 0.0, 1.0)
            period, cycle_mean = cycle_summary(t_eval, clipped) if klass == "oscillatory" else (float("nan"), None)
            return TrajectoryResult(
                status="ok", output=output, klass=klass, attempt=a, method=method,
                nfev=sum(e["nfev"] for e in log), message=msg, horizon_index=h, stats=stats, attempt_log=log,
                final_state=np.clip(system.expand(y_final), 0.0, 1.0).astype(np.float32) if keep_states else None,
                tail_mean=clipped.mean(axis=1).astype(np.float32) if keep_states else None,
                tail_min=clipped.min(axis=1).astype(np.float32) if keep_states else None,
                tail_max=clipped.max(axis=1).astype(np.float32) if keep_states else None,
                tail=tail.astype(np.float32) if keep_tail else None,
                period=period,
                cycle_mean=cycle_mean.astype(np.float32) if keep_states and cycle_mean is not None else None)
        t0, y = horizon, y_final
    return fail("unclassified", a, method, "not classifiable at the last horizon",
                len(horizons) - 1, stats, tail, output)


# ---------------------------------------------------------------------------------------
# Coalitions and sweeps
# ---------------------------------------------------------------------------------------

def initial_states(n_init: int, n_nodes: int, seed: int) -> np.ndarray:
    """The legacy engine's initial states: default_rng(seed).random((n_init, n_nodes))."""
    return np.random.default_rng(seed).random(size=(n_init, n_nodes))


def clamp_mask_for(coalition: int, n_nodes: int) -> np.ndarray:
    """Nodes outside the coalition are clamped. Bit j of the coalition is node j (legacy order)."""
    return np.array([not (coalition >> j) & 1 for j in range(n_nodes)], dtype=bool)


# Per-trajectory arrays a sweep stores, with their dtypes. Status and class are indices into
# STATUSES and TRAJECTORY_CLASSES.
TRAJECTORY_FIELDS = {
    "status": np.int8,
    "class": np.int8,
    "attempt": np.int8,
    "horizon_index": np.int8,
    "attempts_run": np.int8,
    "nfev": np.int32,
    "tail_range": np.float32,
    "final_derivative": np.float32,
    "unclassified_output": np.float64,
    "period": np.float32,
}


@dataclass
class CoalitionResult:
    value: float
    outputs: np.ndarray
    trajectories: dict[str, np.ndarray]
    final_states: np.ndarray | None
    tail_means: np.ndarray | None
    tail_mins: np.ndarray | None
    tail_maxs: np.ndarray | None
    cycle_means: np.ndarray | None
    results: list[TrajectoryResult] = field(default_factory=list, repr=False)

    def count(self, status: Status) -> int:
        return int((self.trajectories["status"] == STATUSES.index(status)).sum())


def simulate_coalition(prep: Prepared, solver: SolverConfig, coalition: int, init: np.ndarray,
                       clamp_value: float = 0.0, classifier: ClassifierConfig | None = None,
                       keep_states: bool = False, keep_tails: bool = False) -> CoalitionResult:
    """A coalition's value is the mean output over trajectories, and exists only if all succeed."""
    if init.shape[0] == 0:
        raise ODEEngineError("a coalition needs at least one initial state")
    n_nodes = len(prep.net.node_names)
    mask = clamp_mask_for(coalition, n_nodes)
    results = [integrate(prep, solver, mask, clamp_value, np.where(mask, clamp_value, y0),
                         classifier, keep_states, keep_tails) for y0 in init]
    outputs = np.array([r.output for r in results])
    value = float("nan")
    if all(r.status == "ok" for r in results):
        value = float(outputs.mean())
        if not np.isfinite(value):
            raise ODEEngineError(f"coalition {coalition}: every trajectory succeeded but the mean is {value}")
    nan = float("nan")
    trajectories = {
        "status": [STATUSES.index(r.status) for r in results],
        "class": [TRAJECTORY_CLASSES.index(r.klass) for r in results],
        "attempt": [r.attempt for r in results],
        "horizon_index": [r.horizon_index for r in results],
        "attempts_run": [len(r.attempt_log) for r in results],
        "nfev": [r.nfev for r in results],
        "tail_range": [r.stats.tail_range if r.stats else nan for r in results],
        "final_derivative": [r.stats.final_derivative if r.stats else nan for r in results],
        "unclassified_output": [r.provisional_output for r in results],
        "period": [r.period for r in results],
    }

    def stack(attr: str) -> np.ndarray | None:
        if not keep_states:
            return None
        blank = np.full(n_nodes, np.nan, dtype=np.float32)
        return np.stack([getattr(r, attr) if getattr(r, attr) is not None else blank for r in results])

    return CoalitionResult(
        value=value,
        outputs=outputs,
        trajectories={k: np.asarray(v, dtype=TRAJECTORY_FIELDS[k]) for k, v in trajectories.items()},
        final_states=stack("final_state"),
        tail_means=stack("tail_mean"),
        tail_mins=stack("tail_min"),
        tail_maxs=stack("tail_max"),
        cycle_means=stack("cycle_mean"),
        results=results,
    )


def provenance(prep: Prepared, network: str, solver: SolverConfig, classifier: ClassifierConfig | None,
               n_init: int, seed: int, clamp_value: float) -> dict:
    return {
        "network": network,
        "model_sha256": model_digest(prep.net),
        "code_sha256": code_digest(),
        "python": platform.python_version(),
        "machine": platform.machine(),
        "system": platform.system(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        # SciPy's solvers do their linear algebra through this library. On Apple Accelerate the
        # legacy engine differs from itself by ~1e-12 relative on identical inputs
        # (results/audit/solver_determinism.json).
        "blas": np.show_config(mode="dicts")["Build Dependencies"]["blas"].get("name"),
        "container_image": os.environ.get("MODAL_IMAGE_ID"),
        "dynamics": asdict(prep.dyn),
        "solver": asdict(solver),
        "classifier": asdict(classifier) if classifier else None,
        "n_init": n_init,
        "seed": seed,
        "clamp_value": clamp_value,
        "clamping": "restoring term" if prep.dyn.construction == "operatorwise_legacy" else "substituted",
        "coalition_acceptance": "all trajectories must have status ok",
        "statuses": list(STATUSES),
        "trajectory_classes": list(TRAJECTORY_CLASSES),
        "operatorwise_fallback_nodes": list(prep.net.operatorwise_fallback),
        "operatorwise_vertex_mismatch_nodes": list(prep.net.operatorwise_vertex_mismatch),
    }


def sweep(net: CompiledNetwork, dyn: DynamicsConfig, solver: SolverConfig, *, network: str, n_init: int,
          seed: int, clamp_value: float = 0.0, classifier: ClassifierConfig | None = None,
          coalition_range: tuple[int, int] | None = None, checkpoint: Path | None = None,
          checkpoint_every: int = 256, keep_states: bool = False,
          on_checkpoint: Callable[[], None] | None = None) -> dict:
    """Exhaustive coalition sweep, checkpointed within the unit of work and resumable.

    `on_checkpoint` runs after every checkpoint write, for example to commit a Modal volume.

    Refuses to run the normalized operator-wise construction on a network whose operator-wise
    rules do not reproduce their truth tables. The legacy construction exists to audit, so it
    runs regardless.
    """
    if dyn.construction == "operatorwise_normalized" and net.operatorwise_vertex_mismatch:
        raise RuleCompileError(
            f"operator-wise rules disagree with their truth tables at a vertex: "
            f"{list(net.operatorwise_vertex_mismatch)}")
    prep = prepare(net, dyn)
    n_nodes = len(net.node_names)
    start, end = coalition_range if coalition_range is not None else (0, 2**n_nodes)
    if not 0 <= start < end <= 2**n_nodes:
        raise ODEEngineError(f"coalition range [{start}, {end}) is not inside [0, {2**n_nodes})")
    size = end - start
    init = initial_states(n_init, n_nodes, seed)

    arrays = {
        "values": np.full(size, np.nan),
        "outputs": np.full((size, n_init), np.nan),
        **{f"trajectory_{k}": np.zeros((size, n_init), dtype=dt) for k, dt in TRAJECTORY_FIELDS.items()},
    }
    if keep_states:
        arrays["final_states"] = np.full((size, n_init, n_nodes), np.nan, dtype=np.float32)
        arrays["tail_mean_states"] = np.full((size, n_init, n_nodes), np.nan, dtype=np.float32)
        arrays["tail_min_states"] = np.full((size, n_init, n_nodes), np.nan, dtype=np.float32)
        arrays["tail_max_states"] = np.full((size, n_init, n_nodes), np.nan, dtype=np.float32)
        arrays["cycle_mean_states"] = np.full((size, n_init, n_nodes), np.nan, dtype=np.float32)

    prov = provenance(prep, network, solver, classifier, n_init, seed, clamp_value)
    fingerprint = json.dumps({**prov, "coalition_range": [start, end], "keep_states": keep_states},
                             sort_keys=True)
    done = 0
    if checkpoint is not None and checkpoint.exists():
        saved = np.load(checkpoint)
        if str(saved["fingerprint"]) != fingerprint:
            raise ODEEngineError(f"{checkpoint} was written under a different configuration or range")
        done = int(saved["completed"])
        for key in arrays:
            arrays[key][:done] = saved[key][:done]

    def save(completed: int) -> None:
        if checkpoint is None:
            return
        tmp = checkpoint.with_name(checkpoint.stem + ".partial.npz")
        np.savez_compressed(tmp, completed=completed, fingerprint=np.array(fingerprint), **arrays)
        tmp.replace(checkpoint)
        if on_checkpoint is not None:
            on_checkpoint()

    for offset in range(done, size):
        r = simulate_coalition(prep, solver, start + offset, init, clamp_value, classifier, keep_states)
        arrays["values"][offset] = r.value
        arrays["outputs"][offset] = r.outputs
        for k, v in r.trajectories.items():
            arrays[f"trajectory_{k}"][offset] = v
        if keep_states:
            arrays["final_states"][offset] = r.final_states
            arrays["tail_mean_states"][offset] = r.tail_means
            arrays["tail_min_states"][offset] = r.tail_mins
            arrays["tail_max_states"][offset] = r.tail_maxs
            arrays["cycle_mean_states"][offset] = r.cycle_means
        if (offset + 1) % checkpoint_every == 0:
            save(offset + 1)
    save(size)
    return {"arrays": arrays, "coalition_range": (start, end), "provenance": prov, "fingerprint": fingerprint}


def require_complete(values: np.ndarray) -> None:
    """Scoring a sweep with missing coalitions would silently change its Walsh spectrum."""
    missing = int(np.isnan(values).sum())
    if missing:
        raise IncompleteSweepError(f"{missing} of {values.size} coalitions have no accepted value")
