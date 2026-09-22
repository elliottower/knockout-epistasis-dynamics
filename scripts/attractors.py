"""The attractors a coalition's trajectories reach, and the entropy of its basins.

Two accepted trajectories of one coalition can reach the same attractor only if they have the
same class. Within a class they are grouped by complete linkage, cut at the tolerance: no group
holds two trajectories farther apart than the tolerance, so a chain of trajectories, each close
to the next, cannot join two distinct attractors into one.

  fixed points    distance is the largest difference between two tail means over nodes. The
                  classifier accepts a fixed point only if its tail range is at most 1e-5.
  oscillations    each trajectory is summarized by its period and, for every node, its mean
                  over a whole number of periods, its tail minimum and its tail maximum
                  (scripts/ode_engine.cycle_summary). None of these depends on the phase at
                  which the tail starts, while a tail mean does: on the 3-node repressor ring at
                  n_H = 4 and 10, 24 random starts on one limit cycle disagreed in tail mean by up
                  to 0.017 and in envelope by at most 0.0022. Two oscillations are within
                  tolerance if all three summaries agree within CYCLE_TOL in every node and their
                  periods agree within PERIOD_RTOL of the longer. An oscillation whose period is
                  undefined is compared on its tail mean and envelope, and is never grouped with
                  one whose period is defined.

A coalition's basin entropy is -sum p_a log2 p_a over its attractors, where p_a is the fraction
of its trajectories that reach attractor a. `network_basins` reports the statistics at the
primary tolerances and at every combination of FIXED_TOL_GRID and CYCLE_TOL_GRID, the registered
sensitivity analysis. Registered in experiments/2026-09-21_ode-primary-rerun/PREREG.md, and
validated before registration by scripts/validate_attractor_grouping.py.
"""

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from scripts.ode_engine import STATUSES, TRAJECTORY_CLASSES

FIXED_TOL = 1e-3
CYCLE_TOL = 0.05
PERIOD_RTOL = 0.05
# The grids lie inside the range the validation found to give one partition. Trajectories on one
# fixed point differed by at most 4e-7 and distinct fixed points by at least 0.49; trajectories on
# one cycle differed by up to 0.012 in their sampled summaries and distinct cycles by at least
# 0.34, so a cycle tolerance of 0.01 split single cycles
# (results/audit/attractor_grouping_validation_cycle_grid_0.01-0.05.json).
FIXED_TOL_GRID = (1e-4, 3e-4, 1e-3)
CYCLE_TOL_GRID = (0.025, 0.05, 0.1)
FIXED = TRAJECTORY_CLASSES.index("fixed")
OSCILLATORY = TRAJECTORY_CLASSES.index("oscillatory")
# Stands in for "never the same attractor" in the linkage, which needs finite distances.
FAR = 1e6


class AttractorError(Exception):
    pass


def chebyshev(states: np.ndarray) -> np.ndarray:
    """(m, n_nodes) -> (m, m) largest absolute difference over nodes."""
    return np.abs(states[:, None, :] - states[None, :, :]).max(axis=2)


def cycle_distance(tail_mean: np.ndarray, tail_min: np.ndarray, tail_max: np.ndarray, cycle_mean: np.ndarray,
                   period: np.ndarray, cycle_tol: float, period_rtol: float = PERIOD_RTOL) -> np.ndarray:
    """Distance between oscillations, scaled so that 1 is the tolerance."""
    undefined = np.isnan(period)
    level = np.where(undefined[:, None], tail_mean, cycle_mean)
    states = np.maximum.reduce([chebyshev(level), chebyshev(tail_min), chebyshev(tail_max)]) / cycle_tol
    with np.errstate(invalid="ignore"):
        periods = np.abs(period[:, None] - period[None, :]) / (period_rtol * np.fmax(period[:, None], period[None, :]))
    periods = np.where(undefined[:, None] & undefined[None, :], 0.0, periods)
    periods = np.where(undefined[:, None] ^ undefined[None, :], np.inf, periods)
    return np.maximum(states, periods)


def tree(distance: np.ndarray, method: str) -> np.ndarray | None:
    if distance.shape[0] < 2:
        return None
    d = np.minimum(distance, FAR)
    np.fill_diagonal(d, 0.0)
    return linkage(squareform(d, checks=False), method=method)


def cut(z: np.ndarray | None, m: int, threshold: float) -> np.ndarray:
    """Labels 0..k-1. Under complete linkage every group's diameter is at most `threshold`."""
    if z is None:
        return np.zeros(m, dtype=np.int64)
    return fcluster(z, t=threshold, criterion="distance").astype(np.int64) - 1


def labels_grid(klass: np.ndarray, tail_mean: np.ndarray, tail_min: np.ndarray, tail_max: np.ndarray,
                cycle_mean: np.ndarray, period: np.ndarray, fixed_tols=(FIXED_TOL,), cycle_tols=(CYCLE_TOL,),
                method: str = "complete") -> dict[tuple[float, float], np.ndarray]:
    """Attractor label of each of one coalition's trajectories, for every pair of tolerances.
    States are (n_init, n_nodes); `method="single"` exists only for the validation."""
    if not np.isin(klass, (FIXED, OSCILLATORY)).all():
        raise AttractorError("every trajectory must be classified fixed or oscillatory")
    fixed, osc = np.flatnonzero(klass == FIXED), np.flatnonzero(klass == OSCILLATORY)
    period = np.asarray(period, dtype=float)
    fixed_distance = chebyshev(tail_mean[fixed].astype(float))
    if fixed.size < 2 or fixed_distance.max() <= min(fixed_tols):
        fixed_labels = {ft: np.zeros(fixed.size, dtype=np.int64) for ft in fixed_tols}
    else:
        z = tree(fixed_distance, method)
        fixed_labels = {ft: cut(z, fixed.size, ft) for ft in fixed_tols}
    cycle_labels = {}
    for ct in cycle_tols:
        d = cycle_distance(*(a[osc].astype(float) for a in (tail_mean, tail_min, tail_max, cycle_mean)), period[osc], ct)
        cycle_labels[ct] = cut(tree(d, method), osc.size, 1.0)
    grid = {}
    for ft in fixed_tols:
        offset = int(fixed_labels[ft].max()) + 1 if fixed.size else 0
        for ct in cycle_tols:
            labels = np.empty(klass.size, dtype=np.int64)
            labels[fixed] = fixed_labels[ft]
            labels[osc] = cycle_labels[ct] + offset
            grid[(ft, ct)] = labels
    return grid


def attractor_labels(klass, tail_mean, tail_min, tail_max, cycle_mean, period, fixed_tol=FIXED_TOL,
                     cycle_tol=CYCLE_TOL, method="complete") -> np.ndarray:
    return labels_grid(klass, tail_mean, tail_min, tail_max, cycle_mean, period, (fixed_tol,), (cycle_tol,),
                       method)[(fixed_tol, cycle_tol)]


def basin_entropy_bits(labels: np.ndarray) -> float:
    p = np.bincount(labels) / labels.size
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def grid_key(fixed_tol: float, cycle_tol: float) -> str:
    return f"fixed_tol={fixed_tol:g},cycle_tol={cycle_tol:g}"


def network_basins(arrays: dict) -> dict:
    """Basin statistics over every coalition of a scored sweep that kept its states, at the
    primary tolerances and at every point of the tolerance grid."""
    if (arrays["trajectory_status"] != STATUSES.index("ok")).any():
        raise AttractorError("basin statistics need every trajectory accepted")
    keys = [(ft, ct) for ft in FIXED_TOL_GRID for ct in CYCLE_TOL_GRID]
    if (FIXED_TOL, CYCLE_TOL) not in keys:
        raise AttractorError("the primary tolerances must lie on the grid")
    n_coalitions = arrays["values"].shape[0]
    entropy = {k: np.empty(n_coalitions) for k in keys}
    count = {k: np.empty(n_coalitions, dtype=np.int64) for k in keys}
    klass, period = arrays["trajectory_class"], arrays["trajectory_period"]
    for c in range(n_coalitions):
        grid = labels_grid(klass[c], arrays["tail_mean_states"][c], arrays["tail_min_states"][c],
                           arrays["tail_max_states"][c], arrays["cycle_mean_states"][c], period[c],
                           FIXED_TOL_GRID, CYCLE_TOL_GRID)
        for k, labels in grid.items():
            entropy[k][c] = basin_entropy_bits(labels)
            count[k][c] = labels.max() + 1

    def stats(k: tuple[float, float]) -> dict:
        return {"mean_entropy_bits": float(entropy[k].mean()),
                "fraction_multistable": float((count[k] > 1).mean()),
                "max_attractors_in_a_coalition": int(count[k].max()),
                "attractor_count_histogram": {str(a): int(b) for a, b in zip(*np.unique(count[k], return_counts=True))}}

    return {
        "fixed_tol": FIXED_TOL, "cycle_tol": CYCLE_TOL, "period_rtol": PERIOD_RTOL, "linkage": "complete",
        **stats((FIXED_TOL, CYCLE_TOL)),
        "oscillatory_trajectories_without_a_period": int(((klass == OSCILLATORY) & np.isnan(period)).sum()),
        "sensitivity": {grid_key(*k): stats(k) for k in keys},
    }
