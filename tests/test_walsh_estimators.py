import numpy as np
import pytest

from data_utils import energy_by_order, energy_spectrum, normalized_wht
from scripts.boolean_estimators import ShardAssemblyError, boolean_outputs
from scripts.walsh_estimators import (
    MIN_RATIO_FRACTION,
    MIN_TOTAL_ENERGY,
    EstimatorError,
    all_pairs_energy,
    fractions,
    higher_order_ratio,
    order3_fraction,
    plain_energy,
    split_half_energy,
)

N = 6


def walsh_matrix(n):
    return np.array([[(-1) ** bin(s & t).count("1") for t in range(2**n)] for s in range(2**n)], dtype=float)


H = walsh_matrix(N)
# A value function with energy at orders 0 to 3, built from its Walsh coefficients.
W = np.zeros(2**N)
W[[0, 1, 3, 7]] = [0.5, 0.12, 0.08, 0.05]
V = H @ W
TRUE_ENERGY = energy_by_order(W, N)


def sampled_outputs(rng, n_init=32, scale=0.25):
    """Each initial state adds its own random function of the coalition, shared across all
    coalitions, as a fixed set of initial states does in the sweeps."""
    noise = rng.normal(size=(n_init, 2**N)) @ H.T * (scale / np.sqrt(2**N))
    return V[:, None] + noise.T


def draws(reps=3000):
    rng = np.random.default_rng()
    out = {name: {"e3": np.empty(reps), "f3": np.empty(reps)} for name in ("all_pairs", "split_half", "plain")}
    for r in range(reps):
        outputs = sampled_outputs(rng)
        for name, estimator in (("all_pairs", all_pairs_energy), ("split_half", split_half_energy), ("plain", plain_energy)):
            energy = estimator(outputs, N)
            out[name]["e3"][r] = energy[3:].sum()
            out[name]["f3"][r] = order3_fraction(fractions(energy))
    return out


@pytest.fixture(scope="module")
def sampled():
    return draws()


def z(samples, truth):
    return (samples.mean() - truth) / (samples.std() / np.sqrt(samples.size))


def test_the_value_function_has_the_coefficients_it_was_built_from():
    np.testing.assert_allclose(normalized_wht(V), W, atol=1e-15)


def test_without_sampling_error_every_estimator_gives_the_exact_spectrum():
    outputs = np.repeat(V[:, None], 32, axis=1)
    for estimator in (all_pairs_energy, split_half_energy, plain_energy):
        np.testing.assert_allclose(fractions(estimator(outputs, N)), energy_spectrum(W, N), atol=1e-12)


def test_with_two_initial_states_all_pairs_is_split_half():
    outputs = np.random.default_rng().random((2**N, 2))
    np.testing.assert_allclose(all_pairs_energy(outputs, N), split_half_energy(outputs, N), atol=1e-15)


def test_all_pairs_and_split_half_are_unbiased_for_the_unnormalized_order_3_energy_and_plain_is_not(sampled):
    truth = TRUE_ENERGY[3:].sum()
    assert abs(z(sampled["all_pairs"]["e3"], truth)) < 4
    assert abs(z(sampled["split_half"]["e3"], truth)) < 4
    assert z(sampled["plain"]["e3"], truth) > 20


def test_the_normalized_fractions_are_ratio_estimates_close_to_the_truth_here_and_plain_overshoots(sampled):
    # A fraction is not unbiased in general; in this example its bias is below the sampling noise.
    truth = order3_fraction(energy_spectrum(W, N))
    assert abs(z(sampled["all_pairs"]["f3"], truth)) < 4
    assert z(sampled["plain"]["f3"], truth) > 20


def test_all_pairs_varies_less_than_split_half(sampled):
    # The two share their leading variance term; all-pairs reduces only the term that is pure
    # noise, so the gain is modest when the signal dominates, as it does here.
    assert sampled["all_pairs"]["f3"].std() < sampled["split_half"]["f3"].std()


@pytest.mark.parametrize("total", [0.0, -1e-3, MIN_TOTAL_ENERGY / 2])
def test_a_spectrum_whose_total_is_nonpositive_or_near_zero_is_inadmissible(total):
    energy = np.zeros(N + 1)
    energy[0] = total
    with pytest.raises(EstimatorError, match="total energy"):
        fractions(energy)


@pytest.mark.parametrize("order2,order3", [
    (0.1, -1e-3),                    # order 3+ negative, order 2+ adequate
    (0.1, MIN_RATIO_FRACTION / 2),   # order 3+ low, order 2+ adequate
    (-0.2, 0.05),                    # order 2+ negative, order 3+ adequate
    (0.0, MIN_RATIO_FRACTION / 2),   # both low
])
def test_a_higher_order_ratio_with_a_nonpositive_or_near_zero_part_is_inadmissible(order2, order3):
    energy = np.zeros(N + 1)
    energy[0], energy[2], energy[3] = 1.0, order2, order3
    with pytest.raises(EstimatorError, match="fraction is below"):
        higher_order_ratio(energy)


def test_an_admissible_higher_order_ratio_is_order_3_over_order_2_and_up():
    energy = np.zeros(N + 1)
    energy[0], energy[2], energy[3], energy[4] = 1.0, 0.2, 0.05, 0.01
    assert higher_order_ratio(energy) == pytest.approx(0.06 / 0.26)


@pytest.mark.parametrize("estimator,outputs,match", [
    (split_half_energy, np.ones((2**N, 3)), "even"),
    (all_pairs_energy, np.ones((2**N, 1)), "at least 2"),
    (all_pairs_energy, np.ones((2**N - 1, 4)), "coalitions"),
    (all_pairs_energy, np.full((2**N, 4), np.nan), "non-finite"),
])
def test_estimators_refuse_inputs_they_are_undefined_for(estimator, outputs, match):
    with pytest.raises(EstimatorError, match=match):
        estimator(outputs, N)


def write_chunks(directory, name, outputs, chunk):
    (directory / name).mkdir(parents=True)
    for start in range(0, outputs.shape[0], chunk):
        np.savez(directory / name / f"{start:07d}.npz", impl=outputs[start:start + chunk],
                 start=start, end=min(start + chunk, outputs.shape[0]))


def test_boolean_outputs_are_assembled_from_their_chunks_in_order(tmp_path):
    outputs = np.random.default_rng().random((2**N, 8))
    write_chunks(tmp_path, "net", outputs, chunk=16)
    np.testing.assert_array_equal(boolean_outputs("net", N, tmp_path), outputs)


def test_a_missing_boolean_chunk_is_refused(tmp_path):
    outputs = np.random.default_rng().random((2**N, 8))
    write_chunks(tmp_path, "net", outputs, chunk=16)
    (tmp_path / "net" / "0000016.npz").unlink()
    with pytest.raises(ShardAssemblyError, match="expected 16"):
        boolean_outputs("net", N, tmp_path)
