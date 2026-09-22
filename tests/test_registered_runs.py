import pytest

from scripts.registered_runs import CLASSIFIER_PRIMARY, N_INIT, RUNS, SEED, expected_networks, launch_mismatches, paper_networks


def fields_for(run: str, **changes) -> dict:
    s = RUNS[run]
    fields = {"network": "lambda_phage", "construction": s["construction"], "hill_n": s["hill_n"], "hill_k": s["hill_k"],
              "solver_path": s["solver"], "classifier_path": s["classifier"], "n_init": N_INIT, "seed": SEED,
              "clamp_value": s["clamp_value"], "keep_states": s["keep_states"]}
    return {**fields, **changes}


def test_every_registered_launch_passes_and_is_refused_as_any_other_run():
    for run in RUNS:
        assert launch_mismatches(run, fields_for(run), expected_networks(run)) == []
        for other in RUNS:
            if other != run:
                assert launch_mismatches(other, fields_for(run), expected_networks(run)), (run, other)


@pytest.mark.parametrize("run,changes,differing", [
    ("sensitivity-c", {"hill_k": 0.5}, ["hill_k"]),                          # --hill-k 0.3 dropped
    ("graded-f0.5", {"clamp_value": 0.25}, ["clamp_value"]),
    ("sensitivity-a", {"keep_states": False}, ["keep_states"]),
    ("legacy-audit", {"classifier_path": CLASSIFIER_PRIMARY}, ["classifier"]),
    ("primary-hillcube-n10", {"n_init": 2}, ["n_init"]),
])
def test_a_launch_with_a_wrong_flag_is_refused_naming_the_flag(run, changes, differing):
    assert launch_mismatches(run, fields_for(run, **changes), expected_networks(run)) == differing


def test_an_unregistered_run_name_is_refused():
    assert launch_mismatches("graded-f1.0", fields_for("graded-f1"), expected_networks("graded-f1")) == ["run_name"]


def test_a_launch_with_a_network_dropped_or_added_is_refused():
    all28 = list(paper_networks())
    assert len(all28) == 28 and len(expected_networks("sensitivity-a")) == 17
    assert launch_mismatches("primary-hillcube-n10", fields_for("primary-hillcube-n10"), all28[1:]) == ["networks"]
    assert launch_mismatches("sensitivity-a", fields_for("sensitivity-a"), all28) == ["networks"]
    assert launch_mismatches("graded-f0.25", fields_for("graded-f0.25"), list(reversed(all28))) == []
