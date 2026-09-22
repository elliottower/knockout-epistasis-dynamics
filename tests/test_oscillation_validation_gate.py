import pytest

from scripts.ode_engine import ClassifierConfig
from scripts.validate_oscillation_classifier import summarize

CLASSIFIER = ClassifierConfig(
    fixed_tail_range_max=1e-5, fixed_derivative_max=1e-6, oscillatory_tail_range_min=1e-4,
    halves_envelope_tol=0.025, bounded_margin=1e-3, extension_horizons=(120.0, 240.0, 480.0, 960.0))


def row(fate, status, klass):
    return {"fate": fate, "status": status, "class": klass,
            "per_horizon": [{"halves_envelope_diff": 0.001, "tail_range": 0.5}]}


@pytest.mark.parametrize("rows,passed", [
    ([row("sustained", "ok", "oscillatory"), row("settles", "ok", "fixed")], True),
    ([row("sustained", "unclassified", "none"), row("settles", "ok", "fixed")], False),
    ([row("sustained", "eval_cap_exceeded", "none"), row("sustained", "ok", "oscillatory")], False),
    ([row("sustained", "ok", "fixed")], False),
    ([row("sustained", "ok", "oscillatory"), row("settles", "ok", "oscillatory")], False),
    ([row("settles", "ok", "fixed"), row("undetermined", "unclassified", "none")], False),
])
def test_the_validation_passes_only_if_every_sustained_trajectory_is_accepted_as_oscillatory(rows, passed):
    assert summarize(rows, CLASSIFIER)["passed"] is passed
