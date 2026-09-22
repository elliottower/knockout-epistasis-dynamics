"""The eleven registered runs, and the settings each is launched with.

scripts/modal_ode_arm.py refuses a registered launch whose run name or settings are not listed
here, and scripts/analyze_ode_rerun.py refuses a record whose settings are not its run's. The
registrations are experiments/2026-09-21_ode-primary-rerun/PREREG.md (primary-hillcube-n10,
legacy-audit), experiments/2026-09-21_graded-perturbation-rerun/PREREG.md (graded-f*) and
experiments/2026-09-21_ode-sensitivity/PREREG.md (sensitivity-*).
"""

import json

from scripts.paths import RESULTS

GRADED_LEVELS = (0.0, 0.25, 0.5, 0.75, 1.0)
# The sensitivity arm runs the networks with at most this many nodes; every other run all 28.
SENSITIVITY_MAX_NODES = 12
N_INIT, SEED = 32, 42
SOLVER_PRIMARY = "experiments/ode_v2/solver_primary.json"
SOLVER_LEGACY = "experiments/ode_v2/solver_legacy.json"
CLASSIFIER_PRIMARY = "experiments/ode_v2/classifier_primary.json"


def run_settings(construction: str, hill_n: float, hill_k: float = 0.5, clamp_value: float = 0.0,
                 solver: str = SOLVER_PRIMARY, classifier: str | None = CLASSIFIER_PRIMARY,
                 keep_states: bool = True) -> dict:
    return {"construction": construction, "hill_n": hill_n, "hill_k": hill_k, "clamp_value": clamp_value,
            "solver": solver, "classifier": classifier, "keep_states": keep_states}


RUNS = {
    "primary-hillcube-n10": run_settings("hillcube_normalized", 10.0),
    "legacy-audit": run_settings("operatorwise_legacy", 10.0, solver=SOLVER_LEGACY, classifier=None,
                                 keep_states=False),
    **{f"graded-f{f:g}": run_settings("hillcube_normalized", 10.0, clamp_value=f) for f in GRADED_LEVELS[1:]},
    "sensitivity-a": run_settings("hillcube_normalized", 2.0),
    "sensitivity-b": run_settings("hillcube_normalized", 4.0),
    "sensitivity-c": run_settings("hillcube_normalized", 10.0, hill_k=0.3),
    "sensitivity-d": run_settings("hillcube_normalized", 10.0, hill_k=0.7),
    "sensitivity-e": run_settings("operatorwise_normalized", 10.0),
}


def paper_networks() -> dict[str, int]:
    """The paper's 28 networks and their sizes, smallest first."""
    table = json.loads((RESULTS / "paper_number_reconciliation.json").read_text())["per_network_table"]
    return {name: table[name]["n"] for name in sorted(table, key=lambda m: (table[m]["n"], m))}


def expected_networks(run: str) -> list[str]:
    sizes = paper_networks()
    if run.startswith("sensitivity-"):
        return [n for n in sizes if sizes[n] <= SENSITIVITY_MAX_NODES]
    return list(sizes)


def launch_mismatches(run_name: str, fields: dict, networks: list[str]) -> list[str]:
    """The settings of a launch that differ from its registered run's, or ["run_name"] if the run
    is not registered. `fields` is one network's spec as scripts/modal_ode_arm.py main builds it,
    and `networks` the networks it launches."""
    if run_name not in RUNS:
        return ["run_name"]
    differing = [] if sorted(networks) == sorted(expected_networks(run_name)) else ["networks"]
    launched = {"construction": fields["construction"], "hill_n": fields["hill_n"], "hill_k": fields["hill_k"],
                "clamp_value": fields["clamp_value"], "solver": fields["solver_path"],
                "classifier": fields["classifier_path"], "keep_states": fields["keep_states"],
                "n_init": fields["n_init"], "seed": fields["seed"]}
    expected = {**RUNS[run_name], "n_init": N_INIT, "seed": SEED}
    return sorted(differing + [k for k in expected if launched[k] != expected[k]])
