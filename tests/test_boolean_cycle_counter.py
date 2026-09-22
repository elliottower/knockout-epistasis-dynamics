import numpy as np
import pytest

from grn_coalition_sweep import compile_network, simulate_sync_output

RING = {"A": "C", "B": "A", "C": "B"}  # synchronous period-3 rotation


def run_ring(max_steps):
    compiled, names = compile_network(RING)
    states = np.array([[1, 0, 0]], dtype=np.int8)
    no_clamp = np.zeros(3, dtype=bool)
    return simulate_sync_output(states, compiled, no_clamp, 0, [names.index("A")], max_steps=max_steps)


def test_a_cycle_longer_than_the_walk_is_counted_as_unclosed():
    output, info = run_ring(max_steps=2)
    assert info["n_cycling"] == 1
    assert info["n_cycle_unclosed"] == 1


def test_a_cycle_that_closes_is_averaged_and_not_counted():
    output, info = run_ring(max_steps=5)
    assert info["n_cycling"] == 1
    assert info["n_cycle_unclosed"] == 0
    assert output[0] == pytest.approx(1 / 3)
