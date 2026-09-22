import json

import pytest

from scripts import analyze_ode_rerun as analysis
from scripts.attractors import CYCLE_TOL_GRID, FIXED_TOL_GRID, grid_key
from scripts.run_ode_arm import ArmSpec, configuration_of, expected_fingerprint

COMMIT = "f" * 40
LOCAL_O3PLUS = 0.3


def energy_for(delta_pp, total=0.1):
    """Order energies whose order-3+ fraction exceeds LOCAL_O3PLUS by delta_pp."""
    o3 = LOCAL_O3PLUS + delta_pp / 100
    return [total * (0.9 - o3), 0.0, total * 0.1, total * o3]


def record(name, delta_pp, scored=True, mismatches=0, o3plus=0.1, ratio=0.4, oscillatory=2, entropy=None,
           energy_3plus=0.01, admissible=True, ratio_admissible=True, total=0.1, commit=COMMIT, run=None):
    estimate = {"admissible": admissible, "delta_o3plus": delta_pp / 100, "global_o3plus": o3plus,
                "higher_order_ratio": ratio if ratio_admissible else None, "ratio_admissible": ratio_admissible,
                "ratio_reason": None if ratio_admissible else "order-3+ fraction is below the admissible 0.0001",
                "energy_3plus": energy_3plus, "energy_2plus": energy_3plus / ratio, "total_energy": total,
                "energy_by_order": energy_for(delta_pp, total)}
    r = {"network": name, "scored": scored, "launch_manifest": {"freeze": {"commit": commit}, "run_name": run},
         "configuration": analysis.expected_configuration(run) if run else None,
         "replicate": {"status_mismatches": mismatches, "class_mismatches": 0, "value_mismatches": 0},
         "summary": {"classes_of_accepted": {"none": 0, "fixed": 32 - oscillatory, "oscillatory": oscillatory, "unclassified": 0},
                     "coalitions_with_an_oscillating_trajectory": 1, "n_coalitions": 8}}
    if scored:
        r["scores"] = {"local_o3plus": LOCAL_O3PLUS, "estimators": {e: dict(estimate) for e in analysis.ESTIMATORS}}
    if entropy is not None:
        r["basins"] = {"mean_entropy_bits": entropy,
                       "sensitivity": {grid_key(ft, ct): {"mean_entropy_bits": entropy}
                                       for ft in FIXED_TOL_GRID for ct in CYCLE_TOL_GRID}}
    return r


@pytest.fixture
def boolean(tmp_path, monkeypatch):
    def write(values, plain=None, published=None, total=0.1, inadmissible=()):
        plain = plain or values
        rows = {}
        for n, v in values.items():
            row = {e: {"admissible": (n, e) not in inadmissible, "delta_o3plus": v / 100, "total_energy": total,
                       "energy_by_order": energy_for(v, total)} for e in ("all_pairs", "split_half")}
            row["plain"] = {"admissible": True, "delta_o3plus": plain[n] / 100, "total_energy": total,
                            "energy_by_order": energy_for(plain[n], total)}
            rows[n] = {"local_o3plus": LOCAL_O3PLUS, **row}
        (tmp_path / "boolean_estimators.json").write_text(json.dumps({"networks": rows}))
        registered = published if published is not None else plain
        (tmp_path / "identities.json").write_text(json.dumps(
            {"networks": {n: {"boolean_delta_3plus_pp_published": v} for n, v in registered.items()}}))
        monkeypatch.setattr(analysis, "BOOLEAN", tmp_path / "boolean_estimators.json")
        monkeypatch.setattr(analysis, "IDENTITIES", tmp_path / "identities.json")
    return write


def networks(k):
    return [f"net{i}" for i in range(k)]


def test_class_preservation_holds_only_if_every_non_null_network_keeps_its_class(boolean):
    names = networks(6)
    boolean({"net0": 5.0, "net1": -3.0, "net2": 0.2, "net3": 12.0, "net4": -1.0, "net5": 2.0})
    same = {n: record(n, d) for n, d in zip(names, [6.0, -2.5, 3.0, 11.0, -1.2, 2.5])}
    result = analysis.evaluate_primary(same, names)["split_half"]
    assert (result["n_non_null"], result["n_preserved"], result["hypotheses"]["H1"]) == (5, 5, True)
    flipped = {**same, "net4": record("net4", 0.3)}
    result = analysis.evaluate_primary(flipped, names)["split_half"]
    assert result["changed_class"] == ["net4"]
    assert result["hypotheses"]["H1"] is False


def test_the_primary_hypotheses_are_void_if_more_than_three_networks_cannot_be_scored(boolean):
    names = networks(8)
    boolean({n: 5.0 for n in names})
    records = {n: record(n, 5.0, scored=i >= 4) for i, n in enumerate(names)}
    assert analysis.evaluate_primary(records, names)["split_half"]["hypotheses"]["void"] is True
    records["net3"] = record("net3", 5.0)
    assert analysis.evaluate_primary(records, names)["split_half"]["hypotheses"]["void"] is False


def test_a_network_with_a_replicate_mismatch_is_withheld(boolean):
    names = networks(4)
    boolean({n: 5.0 for n in names})
    records = {n: record(n, 5.0) for n in names}
    records["net2"] = record("net2", -5.0, mismatches=1)
    result = analysis.evaluate_primary(records, names)["split_half"]
    assert result["excluded"] == {"net2": "withheld for a replicate mismatch"}
    assert result["hypotheses"]["H1"] is True


@pytest.mark.parametrize("n_larger,holds", [(18, True), (17, False)])
def test_the_magnitude_sign_test_needs_18_of_25(boolean, n_larger, holds):
    names = networks(25)
    boolean({n: 5.0 for n in names})
    records = {n: record(n, 6.0 if i < n_larger else 4.0) for i, n in enumerate(names)}
    assert analysis.evaluate_primary(records, names)["split_half"]["hypotheses"]["H4"] is holds


def test_the_total_energy_gate_sensitivity_rescores_from_raw_energies(boolean):
    names = networks(5)
    boolean({n: 5.0 for n in names})
    records = {n: record(n, 5.0) for n in names}
    records["net1"] = record("net1", 5.0, total=5e-6)
    rows = analysis.evaluate_primary(records, names)["all_pairs"]["admissibility_sensitivity"]
    table = {row["min_total_energy"]: row for row in rows}
    assert list(table) == [pytest.approx(1e-7), pytest.approx(1e-6), pytest.approx(1e-5)]
    loose, registered, strict = rows
    assert registered["excluded"] == {} and registered["all_preserved"] is True
    assert list(strict["excluded"]) == ["net1"]
    assert loose["excluded"] == {}


def test_graded_counts_gate_and_inversion(boolean):
    names = networks(5)
    by_level = {
        0.0: {n: record(n, d) for n, d in zip(names, [5, 6, 7, -2, 0.1])},
        0.25: {n: record(n, d) for n, d in zip(names, [3, 4, -1, -2, 0.2])},
        0.5: {n: record(n, d) for n, d in zip(names, [-3, -4, -1, 2, 0.3])},
        0.75: {n: record(n, d, ratio_admissible=n != "net0") for n, d in zip(names, [-3, -4, -1, 2, 0.3])},
        1.0: {n: record(n, d) for n, d in zip(names, [-3, -4, -1, 2, 0.3])},
    }
    result = analysis.evaluate_graded(by_level, names)["split_half"]
    assert result["levels"]["0.0"]["counts"] == {"creation": 3, "destruction": 1, "null": 1}
    assert list(result["levels"]["0.75"]["gated"]) == ["net0"]
    assert result["n_admissible_at_both"] == 4
    assert result["lowest_level_destruction_outnumbers_creation"] == 0.5
    assert result["first_class_change_by_network"]["net2"] == 0.25
    hypotheses = result["hypotheses"]
    assert (hypotheses["G1"], hypotheses["G2"], hypotheses["G3"]) == (True, True, True)
    # G4 needs at least 10 networks admissible at both levels; these are 5.
    assert hypotheses["G4"] is None and "only 4 networks" in hypotheses["G4_void_because"]


def test_sensitivity_is_void_at_a_setting_with_more_than_two_unscored_networks(boolean):
    names = networks(6)
    boolean({n: 5.0 for n in names})
    good = {n: record(n, 5.0) for n in names}
    three_missing = {n: record(n, 5.0, scored=i >= 3) for i, n in enumerate(names)}
    result = analysis.evaluate_sensitivity({"a": good, "b": good, "c": good, "d": three_missing, "e": good}, names)
    assert result["split_half"]["hypotheses"]["S1_by_setting"]["d"] is None
    assert result["split_half"]["hypotheses"]["S1"] is None
    assert result["split_half"]["hypotheses"]["S2"] is True


def test_the_boolean_robustness_check_holds_only_if_no_published_class_moves(boolean):
    boolean({"a": 5.0, "b": 0.4}, plain={"a": 6.0, "b": 0.6})
    check = analysis.boolean_estimator_check()
    assert (check["B1"], check["all_pairs_class_differs_from_published"]) == (False, ["b"])
    boolean({"a": 5.0, "b": 0.7}, plain={"a": 6.0, "b": 0.6})
    assert analysis.boolean_estimator_check()["B1"] is True


def test_the_boolean_robustness_check_is_void_when_a_registered_network_is_missing(boolean):
    boolean({"a": 5.0, "b": 3.0}, published={"a": 5.0, "b": 3.0, "c": -2.0})
    check = analysis.boolean_estimator_check()
    assert check["n_registered"] == 3
    assert (check["void"], check["B1"], check["missing_or_inadmissible"]) == (True, None, ["c"])


def test_the_boolean_robustness_check_is_void_unless_plain_reproduces_the_published_values(boolean):
    boolean({"a": 5.0, "b": 3.0}, plain={"a": 5.004, "b": 3.0}, published={"a": 5.0, "b": 3.0})
    assert analysis.boolean_estimator_check()["B1"] is True
    boolean({"a": 5.0, "b": 3.0}, plain={"a": 5.006, "b": 3.0}, published={"a": 5.0, "b": 3.0})
    check = analysis.boolean_estimator_check()
    assert (check["void"], check["plain_does_not_reproduce_publication"]) == (True, ["a"])


def test_an_inadmissible_split_half_estimate_does_not_void_the_boolean_robustness_check(boolean):
    boolean({"a": 5.0, "b": 3.0}, inadmissible={("b", "split_half")})
    check = analysis.boolean_estimator_check()
    assert (check["void"], check["B1"], check["split_half_inadmissible"]) == (False, True, ["b"])
    boolean({"a": 5.0, "b": 3.0}, inadmissible={("b", "all_pairs")})
    assert analysis.boolean_estimator_check()["missing_or_inadmissible"] == ["b"]


def test_oscillation_falling_with_the_gap_passes_the_one_sided_test_and_rising_does_not():
    names = networks(10)
    falling = {n: record(n, 10.0 - i, oscillatory=i) for i, n in enumerate(names)}
    rising = {n: record(n, float(i), oscillatory=i) for i, n in enumerate(names)}
    o1 = analysis.mechanism_tests(falling, names)["all_pairs"]["O1"]
    assert o1["holds"] is True
    assert o1["rho_interval_95"]["high"] <= 0
    assert analysis.mechanism_tests(rising, names)["all_pairs"]["O1"]["holds"] is False


def test_the_basin_test_is_void_when_fewer_than_five_networks_are_multistable():
    names = networks(10)
    few = {n: record(n, float(i), entropy=0.5 if i < 4 else 0.0) for i, n in enumerate(names)}
    enough = {n: record(n, float(i), entropy=0.1 * i) for i, n in enumerate(names)}
    assert analysis.mechanism_tests(few, names)["split_half"]["A1"]["void"] is True
    tests = analysis.mechanism_tests(enough, names)["split_half"]
    assert tests["A1"]["holds"] is True
    assert len(tests["A1_tolerance_sensitivity"]) == len(FIXED_TOL_GRID) * len(CYCLE_TOL_GRID)
    assert all(t["holds"] is True for t in tests["A1_tolerance_sensitivity"].values())


def test_the_tolerance_sensitivity_keeps_each_networks_value_at_each_grid_point():
    names = networks(10)
    records = {n: record(n, float(i), entropy=0.1 * i) for i, n in enumerate(names)}
    grid = analysis.mechanism_tests(records, names)["all_pairs"]["A1_tolerance_sensitivity"]
    assert all(entry["rows"]["net3"]["feature"] == pytest.approx(0.3) for entry in grid.values())
    table = analysis.basin_table(records, names)
    assert set(table) == set(names) and len(table["net0"]) == len(FIXED_TOL_GRID) * len(CYCLE_TOL_GRID)


def png(path):
    return path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_scatterplots_are_written_for_every_outcome_and_are_byte_identical_on_a_rerun(tmp_path, monkeypatch):
    monkeypatch.setattr(analysis, "PROJECT_ROOT", tmp_path)
    names = networks(10)
    # O1 holds here and A1 is void: only 4 networks are multistable.
    records = {n: record(n, 10.0 - i, oscillatory=i, entropy=0.5 if i < 4 else 0.0) for i, n in enumerate(names)}
    mechanisms = analysis.mechanism_tests(records, names)
    assert mechanisms["all_pairs"]["A1"]["void"] is True
    first = analysis.scatterplots(mechanisms, tmp_path / "figures")
    assert set(first) == set(analysis.ESTIMATORS)
    for estimator in analysis.ESTIMATORS:
        for test in ("O1", "A1"):
            assert first[estimator][test] == f"figures/{test}_{estimator}.png"
            assert png(tmp_path / first[estimator][test])
    before = (tmp_path / "figures/A1_all_pairs.png").read_bytes()
    analysis.scatterplots(mechanisms, tmp_path / "figures")
    assert (tmp_path / "figures/A1_all_pairs.png").read_bytes() == before


def test_no_scatterplot_is_written_when_no_network_is_usable(tmp_path):
    names = networks(4)
    records = {n: record(n, 1.0, scored=False) for n in names}
    assert analysis.scatterplot(analysis.mechanism_tests(records, names)["all_pairs"]["O1"], "O1", "all_pairs",
                                tmp_path / "O1.png") is None
    assert not (tmp_path / "O1.png").exists()


def test_the_mechanism_tests_are_void_if_more_than_three_networks_cannot_be_scored():
    names = networks(12)
    records = {n: record(n, float(i), scored=i >= 4, oscillatory=i, entropy=0.1 * i) for i, n in enumerate(names)}
    tests = analysis.mechanism_tests(records, names)["split_half"]
    assert tests["O1"]["void"] is True and tests["A1"]["void"] is True


def test_four_scored_networks_without_basins_void_the_basin_test():
    names = networks(12)
    records = {n: record(n, float(i), entropy=None if i < 4 else 0.1 * i) for i, n in enumerate(names)}
    a1 = analysis.mechanism_tests(records, names)["all_pairs"]["A1"]
    assert a1["excluded"] == {f"net{i}": "no mean_basin_entropy" for i in range(4)}
    assert a1["void"] is True
    records["net3"] = record("net3", 3.0, entropy=0.3)
    assert analysis.mechanism_tests(records, names)["all_pairs"]["A1"]["void"] is False


def test_four_scored_networks_without_an_oscillation_summary_void_the_oscillation_test():
    names = networks(12)
    records = {n: record(n, 12.0 - i, oscillatory=i + 1) for i, n in enumerate(names)}
    for i in range(4):
        del records[f"net{i}"]["summary"]
    o1 = analysis.mechanism_tests(records, names)["all_pairs"]["O1"]
    assert o1["excluded"] == {f"net{i}": "no oscillation_fraction" for i in range(4)}
    assert o1["void"] is True


def test_holm_adjusts_the_pair_and_leaves_a_void_test_void():
    assert analysis.holm({"O1": 0.01, "A1": 0.04}) == {"O1": pytest.approx(0.02), "A1": pytest.approx(0.04)}
    assert analysis.holm({"O1": 0.03, "A1": 0.02}) == {"O1": pytest.approx(0.04), "A1": pytest.approx(0.04)}
    assert analysis.holm({"O1": 0.01, "A1": None}) == {"O1": pytest.approx(0.01), "A1": None}


def test_an_inadmissible_spectrum_keeps_a_network_out_of_that_estimators_tests(boolean):
    names = networks(4)
    boolean({n: 5.0 for n in names})
    records = {n: record(n, 5.0) for n in names}
    records["net1"] = record("net1", 5.0, admissible=False)
    assert analysis.evaluate_primary(records, names)["all_pairs"]["excluded"] == {"net1": "inadmissible under all_pairs"}


def twelve_networks_at_five_levels(null_at_half, unscored_at_quarter=()):
    names = networks(12)
    base = {n: record(n, 5.0) for n in names}
    half = {n: record(n, 0.1 if i < null_at_half else 5.0) for i, n in enumerate(names)}
    quarter = {n: record(n, 5.0, scored=n not in unscored_at_quarter) for n in names}
    return names, {0.0: base, 0.25: quarter, 0.5: half, 0.75: base, 1.0: base}


@pytest.mark.parametrize("null_at_half,holds", [(3, True), (4, False)])
def test_the_gap_persists_at_every_partial_level_only_if_at_most_a_quarter_is_null(null_at_half, holds):
    names, by_level = twelve_networks_at_five_levels(null_at_half)
    assert analysis.evaluate_graded(by_level, names)["all_pairs"]["hypotheses"]["G5"] is holds


def test_g5_on_the_common_set_counts_only_networks_scored_at_every_level():
    names, by_level = twelve_networks_at_five_levels(3, unscored_at_quarter=("net0", "net1"))
    common = analysis.evaluate_graded(by_level, names)["all_pairs"]["G5_on_networks_scored_at_every_level"]
    assert common["n"] == 10 and "net0" not in common["networks"]
    # net2 is null at f = 0.5 and scored everywhere: 1 of 10.
    assert common["null_counts"]["0.5"] == 1 and common["holds"] is True


def test_g4_is_evaluated_once_ten_networks_are_admissible_at_both_levels():
    names, by_level = twelve_networks_at_five_levels(0)
    result = analysis.evaluate_graded(by_level, names)["all_pairs"]
    assert result["hypotheses"]["G4"] is True and result["hypotheses"]["G4_void_because"] is None
    assert set(result["ratio_gate_sensitivity"]) == {"1e-05", "0.0001", "0.001"}


def test_a_record_launched_from_another_commit_is_refused(tmp_path):
    run = tmp_path / "primary-hillcube-n10"
    run.mkdir()
    (run / "net0.json").write_text(json.dumps(record("net0", 5.0, run=run.name)))
    # A container's sidecar sits beside the records and is not one.
    (run / "net0.shard_0000000-0004096.task_ta-01.json").write_text(json.dumps({"task_id": "ta-01"}))
    assert set(analysis.load_records(run, COMMIT, ["net0"])) == {"net0"}
    (run / "net1.json").write_text(json.dumps(record("net1", 5.0, commit="e" * 40, run=run.name)))
    with pytest.raises(analysis.AnalysisError, match="not from the frozen commit"):
        analysis.load_records(run, COMMIT, ["net0", "net1"])


def test_a_run_missing_a_networks_record_or_holding_an_extra_one_is_refused(tmp_path):
    # A fetch made before every network finished would otherwise read as unscorable networks.
    run = tmp_path / "primary-hillcube-n10"
    run.mkdir()
    (run / "net0.json").write_text(json.dumps(record("net0", 5.0, run=run.name)))
    with pytest.raises(analysis.AnalysisError, match=r"lacks records for \['net1'\]"):
        analysis.load_records(run, COMMIT, ["net0", "net1"])
    with pytest.raises(analysis.AnalysisError, match=r"unregistered networks \['net0'\]"):
        analysis.load_records(run, COMMIT, ["net1"])


def test_a_record_whose_settings_are_not_its_runs_is_refused(tmp_path):
    # A dropped --hill-k 0.3 would give setting c the reference's records.
    run = tmp_path / "sensitivity-c"
    run.mkdir()
    (run / "net0.json").write_text(json.dumps(record("net0", 5.0, run="sensitivity-c")))
    assert set(analysis.load_records(run, COMMIT, ["net0"])) == {"net0"}
    dropped = record("net0", 5.0, run="sensitivity-c")
    dropped["configuration"] = analysis.expected_configuration("sensitivity-a")
    (run / "net0.json").write_text(json.dumps(dropped))
    with pytest.raises(analysis.AnalysisError, match="dynamics"):
        analysis.load_records(run, COMMIT, ["net0"])


def test_a_run_directory_that_is_not_registered_missing_or_empty_is_refused(tmp_path):
    (tmp_path / "graded-f1.0").mkdir()
    with pytest.raises(analysis.AnalysisError, match="not a registered run"):
        analysis.load_records(tmp_path / "graded-f1.0", COMMIT, ["net0"])
    # A run fetched under the wrong name leaves the registered directory missing.
    with pytest.raises(analysis.AnalysisError, match="does not exist"):
        analysis.load_records(tmp_path / "graded-f1", COMMIT, ["net0"])
    (tmp_path / "graded-f1").mkdir()
    with pytest.raises(analysis.AnalysisError, match="lacks records"):
        analysis.load_records(tmp_path / "graded-f1", COMMIT, ["net0"])


def test_a_record_launched_under_another_run_name_is_refused(tmp_path):
    run = tmp_path / "graded-f0.5"
    run.mkdir()
    r = record("net0", 5.0, run="graded-f0.5")
    r["launch_manifest"]["run_name"] = "graded-f0.50"
    (run / "net0.json").write_text(json.dumps(r))
    with pytest.raises(analysis.AnalysisError, match="launched as run graded-f0.50"):
        analysis.load_records(run, COMMIT, ["net0"])


def spec_for(run):
    s = analysis.RUNS[run]
    return ArmSpec(network="lambda_phage", construction=s["construction"], hill_n=s["hill_n"],
                   solver_path=analysis.REPO / s["solver"],
                   classifier_path=analysis.REPO / s["classifier"] if s["classifier"] else None,
                   n_init=32, seed=42, clamp_value=s["clamp_value"], keep_states=s["keep_states"], hill_k=s["hill_k"])


def test_the_configuration_the_pipeline_writes_passes_its_own_run_and_no_other():
    written = {run: {"configuration": json.loads(json.dumps(configuration_of(expected_fingerprint(spec_for(run)))))}
               for run in analysis.RUNS}
    for run, rec in written.items():
        assert analysis.configuration_mismatches(rec, run) == []
        for other in analysis.RUNS:
            if other != run:
                assert analysis.configuration_mismatches(rec, other), (run, other)


def test_every_registered_run_expects_its_own_settings():
    assert set(analysis.RUNS) == {"primary-hillcube-n10", "legacy-audit", "graded-f0.25", "graded-f0.5",
                                  "graded-f0.75", "graded-f1", *analysis.SETTINGS.values()}
    configs = {run: analysis.expected_configuration(run) for run in analysis.RUNS}
    assert len({json.dumps(c, sort_keys=True) for c in configs.values()}) == len(configs)
    assert configs["legacy-audit"]["classifier"] is None
    assert configs["graded-f0.75"]["clamp_value"] == 0.75


def test_an_unscored_network_that_is_not_like_for_like_stays_out_of_the_counts(monkeypatch):
    standard = {"t_max": 30.0, "t_tail": 10.0, "hill_n": 10.0, "hill_k": 0.5, "tau": 1.0}
    monkeypatch.setattr(analysis, "published_ode", lambda: {
        "net0": {"pp": 5.0, "settings": standard, "file": "a"},
        "net1": {"pp": 5.0, "settings": {**standard, "t_max": 10.0}, "file": "b"},
        "net2": {"pp": 5.0, "settings": standard, "file": "c"}})
    audit = analysis.legacy_audit({"net0": record("net0", 6.0), "net1": record("net1", 0.0, scored=False),
                                   "net2": record("net2", 5.2)}, ["net0", "net1", "net2"])
    assert audit["unscored"] == {"net1": {"reason": "cannot be scored", "like_for_like": False}}
    assert (audit["n_unscored_like_for_like"], audit["n_like_for_like"], audit["n_moved_over_0p5pp"]) == (0, 2, 1)


@pytest.mark.parametrize("settings,same", [
    ({"t_max": 30.0, "t_tail": 10.0, "hill_n": 10.0, "hill_k": 0.5, "tau": 1.0}, True),
    ({"t_max": 10.0, "t_tail": 5.0, "hill_n": 10.0, "hill_k": 0.5, "tau": 1.0, "n_init": 32}, False),
    (None, False),
])
def test_a_legacy_row_counts_only_if_the_published_run_used_the_reconstructed_settings(settings, same):
    assert analysis.legacy_like_for_like(settings)[0] is same


def test_each_arm_reads_its_registered_runs_and_writes_every_estimator(tmp_path, monkeypatch, boolean):
    sizes = {f"net{i}": 7 if i < 10 else 13 for i in range(12)}
    boolean({n: 5.0 for n in sizes})
    monkeypatch.setattr(analysis, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(analysis, "paper_networks", lambda: sizes)
    monkeypatch.setattr(analysis, "frozen_commit", lambda: COMMIT)
    for attr in ("PRIMARY", "GRADED", "SENSITIVITY"):
        monkeypatch.setattr(analysis, attr, tmp_path / attr.lower())
    runs = {analysis.PRIMARY: ["primary-hillcube-n10", "legacy-audit"],
            analysis.GRADED: [f"graded-f{f:g}" for f in analysis.GRADED_LEVELS[1:]],
            analysis.SENSITIVITY: list(analysis.SETTINGS.values())}
    for arm_dir, names in runs.items():
        for run in names:
            (arm_dir / "results" / run).mkdir(parents=True)
            # The sensitivity arm runs only the networks with at most 12 nodes.
            for n in (x for x in sizes if arm_dir != analysis.SENSITIVITY or sizes[x] <= 12):
                (arm_dir / "results" / run / f"{n}.json").write_text(json.dumps(record(n, 6.0, run=run)))
    for arm in ("primary", "graded", "sensitivity"):
        analysis.main([arm])
    primary = json.loads((analysis.PRIMARY / "results" / "analysis.json").read_text())
    graded = json.loads((analysis.GRADED / "results" / "analysis.json").read_text())
    sensitivity = json.loads((analysis.SENSITIVITY / "results" / "analysis.json").read_text())
    assert primary["frozen_commit"] == graded["frozen_commit"] == sensitivity["frozen_commit"] == COMMIT
    assert set(primary["tests"]) == set(graded["tests"]) == set(analysis.ESTIMATORS)
    assert primary["tests"]["all_pairs"]["n_preserved"] == 12
    assert primary["boolean_estimator_check"]["B1"] is True
    # These records carry no basin statistics, so A1 has no usable network and no plot.
    assert all(png(tmp_path / tests["O1"]) and tests["A1"] is None for tests in primary["scatterplots"].values())
    assert "basins_by_network_and_tolerance" in primary
    assert set(primary["legacy_audit"]["rows"]) == set(sizes)
    assert [graded["tests"]["all_pairs"]["levels"][str(f)]["n_scored"] for f in analysis.GRADED_LEVELS] == [12] * 5
    assert sensitivity["networks"] == [f"net{i}" for i in range(10)]
    assert set(sensitivity["reference"]) == set(sensitivity["tests"]) == set(analysis.ESTIMATORS)
    assert sensitivity["reference"]["all_pairs"]["n_scored"] == 10
    assert sensitivity["tests"]["all_pairs"]["hypotheses"]["S1"] is True
