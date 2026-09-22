"""Replay the oscillation validation panel at other envelope tolerances.

    uv run python -m scripts.replay_oscillation_tolerance

Reads the per-horizon diagnostics in results/audit/oscillation_classifier_validation.json, the
panel as run at a 5% envelope tolerance, and reclassifies every trajectory at each tolerance,
horizon by horizon, with the classifier's other thresholds unchanged. Each horizon's diagnostics
come from an independent integration from t = 0 to that horizon. The diagnostics carry no bounds
check; none is needed, because (F(x) - x)/tau with F in [0, 1] keeps every state in the unit cube.

Writes results/audit/oscillation_tolerance_replay.json.
"""

import json

from scripts.paths import AUDIT, PROJECT_ROOT

SOURCE = AUDIT / "oscillation_classifier_validation.json"
OUT = AUDIT / "oscillation_tolerance_replay.json"
TOLERANCES = (0.01, 0.02, 0.025, 0.03, 0.04, 0.05)


def replay(per_horizon: list[dict], classifier: dict, tol: float) -> tuple[str, int | None]:
    """The class and accepting horizon a trajectory gets at envelope tolerance `tol`."""
    for h, stats in enumerate(per_horizon):
        if (stats["tail_range"] <= classifier["fixed_tail_range_max"]
                and stats["final_derivative"] <= classifier["fixed_derivative_max"]):
            return "fixed", h
        if stats["tail_range"] > classifier["oscillatory_tail_range_min"] and stats["halves_envelope_diff"] <= tol:
            return "oscillatory", h
    return "unclassified", None


def main() -> None:
    panel = json.loads(SOURCE.read_text())
    by_tolerance = {}
    for tol in TOLERANCES:
        counts = {"sustained -> oscillatory": 0, "sustained accepted at horizon 0": 0,
                  "settles -> fixed": 0, "undetermined -> unclassified": 0}
        for row in panel["trajectories"]:
            klass, horizon = replay(row["per_horizon"], row["classifier"], tol)
            if row["fate"] == "sustained" and klass == "oscillatory":
                counts["sustained -> oscillatory"] += 1
                counts["sustained accepted at horizon 0"] += horizon == 0
            elif row["fate"] == "settles" and klass == "fixed":
                counts["settles -> fixed"] += 1
            elif row["fate"] == "undetermined" and klass == "unclassified":
                counts["undetermined -> unclassified"] += 1
        by_tolerance[f"{tol:g}"] = counts
    OUT.write_text(json.dumps({
        "source": f"{SOURCE.relative_to(PROJECT_ROOT).as_posix()} (per-horizon diagnostics, replayed)",
        "note": "each horizon's statistics come from an independent integration from t = 0 to that horizon",
        "by_tolerance": by_tolerance,
    }, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
