"""Order-by-order Walsh energy of a coalition value function estimated from sampled trajectories.

A coalition's value is the mean output over m initial states, so it carries sampling error. The
plain estimator squares the Walsh coefficients of that mean, and the squared error adds to every
coefficient's energy. Most coefficients are high order, so the error inflates the order-3+
fraction. The same initial states are used for every coalition, so the error is correlated
across coalitions and no constant noise floor can be subtracted.

Let w_i(T) be the Walsh coefficients of the value function computed from initial state i alone.
The initial states are independent draws, so for i ≠ j, E[w_i(T) w_j(T)] = w(T)^2, the squared
coefficient of the value function that the full distribution of initial states defines.

  all-pairs    U(T) = ((Σ_i w_i(T))^2 − Σ_i w_i(T)^2) / (m (m − 1)), the average over every
               ordered pair i ≠ j: a U-statistic, unbiased for w(T)^2, with no choice of split.
  split-half   w_A(T) w_B(T), with A and B the means over the first and second halves of the
               states: unbiased for w(T)^2, with one fixed split. Its variance shares all-pairs'
               leading term and exceeds it in the term that is pure noise.
  plain        w(T)^2 of the mean over all states: biased upward by the sampling error.

Order energies are sums of these over |T| = k, so under all-pairs and split-half they are
unbiased for the unnormalized order energies. They can be negative at an order with little true
energy, and are reported as estimated. Fractions divide the order energies by the estimated total
(order 0 included, as the published spectra are). A fraction is a ratio of estimates and is not
unbiased: it reduces the plain estimator's upward bias, and its finite-sample behavior is checked
in tests/test_walsh_estimators.py and by the registered comparisons between estimators.

Admissibility, predeclared:

  MIN_TOTAL_ENERGY     A spectrum is admissible only if its estimated total energy is at least
                       1e-6. The total is the value function's mean square, so this requires a
                       root-mean-square value of at least 0.001; below that the output is zero
                       for practical purposes and no fraction is defined.
  MIN_RATIO_FRACTION   sigma_3+/sigma_2+ is admissible only if the order-3+ and order-2+
                       energies are each at least 1e-4 of the total, the published arm's gate.
                       Both are then positive, so the ratio cannot come from a nonpositive or
                       near-zero denominator.
"""

import numpy as np

from data_utils import normalized_wht, popcount_array

MIN_TOTAL_ENERGY = 1e-6
MIN_RATIO_FRACTION = 1e-4
ESTIMATORS = ("all_pairs", "split_half", "plain")


class EstimatorError(Exception):
    """The estimator is undefined or inadmissible for these inputs."""


def _check(outputs: np.ndarray, n: int, min_init: int) -> None:
    n_coalitions, n_init = outputs.shape
    if n_coalitions != 2**n:
        raise EstimatorError(f"{n_coalitions} coalitions for n = {n}")
    if n_init < min_init:
        raise EstimatorError(f"at least {min_init} initial states needed, got {n_init}")
    if not np.all(np.isfinite(outputs)):
        raise EstimatorError("outputs contain non-finite values")


def by_order(per_coefficient: np.ndarray, n: int) -> np.ndarray:
    pc = popcount_array(n)
    return np.array([per_coefficient[pc == k].sum() for k in range(n + 1)])


def plain_energy(outputs: np.ndarray, n: int) -> np.ndarray:
    _check(outputs, n, 1)
    return by_order(normalized_wht(outputs.mean(axis=1)) ** 2, n)


def split_half_energy(outputs: np.ndarray, n: int) -> np.ndarray:
    _check(outputs, n, 2)
    if outputs.shape[1] % 2:
        raise EstimatorError(f"split-half needs an even number of initial states, got {outputs.shape[1]}")
    half = outputs.shape[1] // 2
    return by_order(normalized_wht(outputs[:, :half].mean(axis=1)) * normalized_wht(outputs[:, half:].mean(axis=1)), n)


def all_pairs_energy(outputs: np.ndarray, n: int) -> np.ndarray:
    _check(outputs, n, 2)
    m = outputs.shape[1]
    total = np.zeros(2**n)
    squares = np.zeros(2**n)
    for i in range(m):
        w = normalized_wht(outputs[:, i])
        total += w
        squares += w * w
    return by_order((total * total - squares) / (m * (m - 1)), n)


ENERGY = {"all_pairs": all_pairs_energy, "split_half": split_half_energy, "plain": plain_energy}


def fractions(energy: np.ndarray) -> np.ndarray:
    total = float(energy.sum())
    if not total >= MIN_TOTAL_ENERGY:
        raise EstimatorError(f"total energy {total:.3g} is below the admissible {MIN_TOTAL_ENERGY:g}")
    return energy / total


def order3_fraction(spectrum: np.ndarray) -> float:
    return float(spectrum[3:].sum()) if spectrum.shape[0] > 3 else 0.0


def higher_order_ratio(energy: np.ndarray) -> float:
    """sigma_3+ / sigma_2+, the graded arm's normalized higher-order fraction, if admissible."""
    total = float(energy.sum())
    e3, e2 = float(energy[3:].sum()), float(energy[2:].sum())
    if not total >= MIN_TOTAL_ENERGY:
        raise EstimatorError(f"total energy {total:.3g} is below the admissible {MIN_TOTAL_ENERGY:g}")
    if not (e3 >= MIN_RATIO_FRACTION * total and e2 >= MIN_RATIO_FRACTION * total):
        raise EstimatorError(f"order-3+ ({e3 / total:.3g}) or order-2+ ({e2 / total:.3g}) fraction is below "
                             f"the admissible {MIN_RATIO_FRACTION:g}")
    return e3 / e2
