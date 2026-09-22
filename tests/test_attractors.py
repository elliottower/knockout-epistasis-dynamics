import numpy as np
import pytest

from scripts.attractors import (
    CYCLE_TOL_GRID,
    FIXED,
    FIXED_TOL_GRID,
    OSCILLATORY,
    AttractorError,
    attractor_labels,
    basin_entropy_bits,
    labels_grid,
    network_basins,
)
from scripts.ode_engine import STATUSES, TRAJECTORY_CLASSES, DynamicsConfig, compile_model, cycle_summary, integrate, prepare
from scripts.paths import PROJECT_ROOT
from scripts.run_ode_arm import load_classifier, load_solver

SOLVER = load_solver(PROJECT_ROOT / "experiments/ode_v2/solver_primary.json")
CLASSIFIER = load_classifier(PROJECT_ROOT / "experiments/ode_v2/classifier_primary.json")
TOGGLE = {"A": "!B", "B": "!A"}
RING3 = {"A": "!C", "B": "!A", "C": "!B"}
RING4 = {"A": "!D", "B": "!A", "C": "!B", "D": "!C"}
# A self-sustaining switch S selects a 3-node ring when on and a 5-node ring when off: two
# distinct limit cycles in one system.
SWITCHED_RINGS = {"S": "S", "A": "S & !C", "B": "S & !A", "C": "S & !B",
                  "D": "!S & !H", "E": "!S & !D", "F": "!S & !E", "G": "!S & !F", "H": "!S & !G"}


def trajectories(rules, n_init, hill_n=10.0):
    net = compile_model(rules, [next(iter(rules))])
    prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=hill_n))
    free = np.zeros(len(rules), dtype=bool)
    results = [integrate(prep, SOLVER, free, 0.0, y0, CLASSIFIER, keep_states=True)
               for y0 in np.random.default_rng().random((n_init, len(rules)))]
    assert all(r.status == "ok" for r in results)
    klass = np.array([FIXED if r.klass == "fixed" else OSCILLATORY for r in results])
    blank = np.full(len(rules), np.nan, dtype=np.float32)
    stack = lambda attr: np.stack([getattr(r, attr) if getattr(r, attr) is not None else blank for r in results])
    return (klass, stack("tail_mean"), stack("tail_min"), stack("tail_max"), stack("cycle_mean"),
            np.array([r.period for r in results]))


@pytest.mark.parametrize("rules,states", [
    (TOGGLE, [(1, 0), (0, 1)]),
    (RING4, [(1, 0, 1, 0), (0, 1, 0, 1)]),
])
def test_a_bistable_network_has_two_fixed_point_attractors_at_its_two_stable_states(rules, states):
    klass, mean, low, high, cycle, period = trajectories(rules, 32)
    assert (klass == FIXED).all()
    labels = attractor_labels(klass, mean, low, high, cycle, period)
    found = {tuple(np.round(mean[labels == a].mean(axis=0)).astype(int)) for a in np.unique(labels)}
    assert found <= set(states)
    # With 32 random starts, both basins are reached unless one holds almost none of the cube.
    assert len(found) == 2


def test_trajectories_on_one_limit_cycle_are_one_attractor_whatever_their_phase():
    klass, mean, low, high, cycle, period = trajectories(RING3, 16)
    assert (klass == OSCILLATORY).all()
    assert np.isfinite(period).all()
    for labels in labels_grid(klass, mean, low, high, cycle, period, FIXED_TOL_GRID, CYCLE_TOL_GRID).values():
        assert labels.max() == 0
    assert basin_entropy_bits(attractor_labels(klass, mean, low, high, cycle, period)) == 0.0


def test_two_limit_cycles_selected_by_a_switch_are_two_attractors():
    klass, mean, low, high, cycle, period = trajectories(SWITCHED_RINGS, 16)
    assert (klass == OSCILLATORY).all()
    on = mean[:, 0] > 0.5
    # 16 random starts all on one side of the switch has probability 2^-15.
    assert on.any() and (~on).any()
    labels = attractor_labels(klass, mean, low, high, cycle, period)
    assert labels.max() == 1
    assert len(set(labels[on])) == 1 and len(set(labels[~on])) == 1 and labels[on][0] != labels[~on][0]
    assert np.median(period[~on]) > np.median(period[on])


def fixed(means):
    means = np.asarray(means, dtype=float)
    return np.full(len(means), FIXED), means, means, means, np.full_like(means, np.nan), np.full(len(means), np.nan)


def test_a_chain_of_close_fixed_points_is_not_joined_into_one_attractor():
    # Each neighbour lies within the tolerance, the ends do not.
    args = fixed([[0.2], [0.2008], [0.2016]])
    assert attractor_labels(*args, fixed_tol=1e-3, method="single").max() == 0
    complete = attractor_labels(*args, fixed_tol=1e-3)
    assert complete.max() == 1
    assert complete[0] != complete[2]


def test_fixed_points_farther_apart_than_the_tolerance_are_distinct():
    labels = attractor_labels(*fixed([[0.2, 0.8], [0.2, 0.8 + 2e-3], [0.2, 0.8 + 5e-4]]))
    assert labels[0] == labels[2] != labels[1]


def cycles(periods, level=0.5):
    n = len(periods)
    state = np.full((n, 3), level)
    return np.full(n, OSCILLATORY), state, state - 0.4, state + 0.4, state, np.asarray(periods, dtype=float)


def test_oscillations_with_one_envelope_and_different_periods_are_distinct():
    assert attractor_labels(*cycles([5.0, 7.0])).max() == 1
    assert attractor_labels(*cycles([5.0, 5.1])).max() == 0


def test_an_oscillation_without_a_period_is_never_grouped_with_one_that_has_one():
    assert attractor_labels(*cycles([5.0, np.nan])).max() == 1
    assert attractor_labels(*cycles([np.nan, np.nan])).max() == 0


def test_a_fixed_point_and_an_oscillation_are_never_the_same_attractor():
    state = np.full((2, 3), 0.5)
    labels = attractor_labels(np.array([FIXED, OSCILLATORY]), state, state, state, state, np.array([np.nan, 5.0]))
    assert labels[0] != labels[1]


@pytest.mark.parametrize("labels,bits", [([0] * 32, 0.0), ([0] * 16 + [1] * 16, 1.0), ([0, 1, 2, 3] * 8, 2.0)])
def test_basin_entropy_in_bits(labels, bits):
    assert basin_entropy_bits(np.array(labels)) == pytest.approx(bits, abs=1e-12)


def test_an_unclassified_trajectory_has_no_attractor():
    state = np.zeros((2, 2))
    with pytest.raises(AttractorError):
        attractor_labels(np.array([FIXED, 0]), state, state, state, state, np.full(2, np.nan))


def test_network_basins_reports_every_point_of_the_tolerance_grid():
    # Two coalitions of 4 trajectories: one monostable, one with fixed points 5e-4 apart, which
    # are one attractor at the loosest fixed tolerance and two at the others.
    means = np.array([[[0.1, 0.1]] * 4, [[0.1, 0.1], [0.1, 0.1], [0.1, 0.1005], [0.1, 0.1005]]], dtype=np.float32)
    arrays = {"values": np.zeros(2), "trajectory_status": np.full((2, 4), STATUSES.index("ok")),
              "trajectory_class": np.full((2, 4), TRAJECTORY_CLASSES.index("fixed")),
              "trajectory_period": np.full((2, 4), np.nan, dtype=np.float32),
              "tail_mean_states": means, "tail_min_states": means, "tail_max_states": means,
              "cycle_mean_states": np.full_like(means, np.nan)}
    basins = network_basins(arrays)
    assert len(basins["sensitivity"]) == len(FIXED_TOL_GRID) * len(CYCLE_TOL_GRID)
    assert basins["mean_entropy_bits"] == pytest.approx(0.0)
    assert basins["sensitivity"]["fixed_tol=0.0001,cycle_tol=0.05"]["mean_entropy_bits"] == pytest.approx(0.5)
    assert basins["sensitivity"]["fixed_tol=0.0003,cycle_tol=0.1"]["fraction_multistable"] == pytest.approx(0.5)


def sampled(period, phase, span=20.0, rate=5.0):
    t = np.linspace(0.0, span, int(span * rate))
    return t, 0.5 + 0.4 * np.sin(2 * np.pi * t / period + phase)


def test_the_period_and_whole_period_mean_do_not_depend_on_phase():
    periods, means, window_means = [], [], []
    for phase in np.random.default_rng().uniform(0, 2 * np.pi, 200):
        t, y = sampled(6.0, phase)
        period, cycle_mean = cycle_summary(t, y[None, :])
        periods.append(period)
        means.append(cycle_mean[0])
        window_means.append(y.mean())
    assert np.max(np.abs(np.array(periods) / 6.0 - 1)) < 0.01
    assert np.max(np.abs(np.array(means) - 0.5)) < 0.01
    # The mean over the whole window, which ends part-way through a period, moves with phase.
    assert np.ptp(window_means) > 5 * np.ptp(means)


def test_a_node_that_rises_twice_per_cycle_does_not_halve_the_period():
    t, once = sampled(6.0, 0.3)
    _, twice = sampled(3.0, 1.1)
    period, _ = cycle_summary(t, np.stack([twice, once]))
    assert period == pytest.approx(6.0, rel=0.01)


def test_a_tail_with_fewer_than_two_upward_crossings_has_no_period():
    t, y = sampled(30.0, 0.0)
    period, cycle_mean = cycle_summary(t, y[None, :])
    assert np.isnan(period) and np.isnan(cycle_mean).all()
