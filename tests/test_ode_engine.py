import dataclasses
import inspect
import json
import shutil

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from scripts import ode_engine as engine
from scripts import run_ode_arm
from scripts.ode_coalition_sweep import ALL_MODELS, boolean_expr_to_continuous, simulate_ode_output
from scripts.ode_engine import (
    LEGACY_SOLVER,
    STATUSES,
    TRAJECTORY_CLASSES,
    ClassifierConfig,
    DynamicsConfig,
    IncompleteSweepError,
    ODEEngineError,
    SolverAttempt,
    SolverConfig,
    clamp_mask_for,
    compile_model,
    hill,
    hill_normalized,
    initial_states,
    integrate,
    model_digest,
    prepare,
    require_complete,
    rule_targets,
    rule_targets_reference,
    simulate_coalition,
    sweep,
)
from scripts.paths import PROJECT_ROOT

PAPER_NETWORKS = sorted([
    "albert_segment_polarity", "arabidopsis_cellcycle", "arellano_rootstem",
    "asymmetric_cell_division", "blood_stem_cell", "calzone_cell_fate",
    "calzone_cellfate_reduced", "cell_cycle_transcription", "davidich_yeast", "drosophila_cellcycle",
    "emt_switch", "fanconi_anemia", "faure_cellcycle", "fumia_cellcycle", "grieco_bladder",
    "hematopoiesis_aging", "irons_cardiac", "lac_operon", "lambda_phage", "li_budding_yeast",
    "mendoza_thelper", "morphogenetic_checkpoint", "myeloid_progenitors", "pair_rule_module",
    "remy_p53_mdm2", "saadatpour_guardcell", "tournier_apoptosis", "zanudo_tlgl",
])
HILL_NS = (1.0, 2.0, 4.0, 10.0)
# The legacy engine differs from itself by up to 1.14e-12 relative on identical inputs
# (scripts/audit_solver_determinism.py). Its right-hand side is scalar Python, so the variation
# arises inside SciPy's solver, most likely its BLAS-backed stage arithmetic; the call has not
# been isolated. Comparisons allow |x - y| <= 1e-12 + 1e-9 * max(|x|, |y|), and nothing more.
REPRODUCTION_ATOL = 1e-12
REPRODUCTION_RTOL = 1e-9

PRIMARY = SolverConfig(
    attempts=(SolverAttempt("RK45", 1_000_000),), rtol=1e-6, atol=1e-8,
    max_step=2.0, t_max=60.0, t_tail=20.0, n_tail_samples=100)
CLASSIFIER = ClassifierConfig(
    fixed_tail_range_max=1e-5, fixed_derivative_max=1e-6, oscillatory_tail_range_min=1e-4,
    halves_envelope_tol=0.05, bounded_margin=1e-3, extension_horizons=(120.0, 240.0))
RING = {"A": "!C", "B": "!A", "C": "!B"}
# X holds its own value and Y copies X, so with X clamped at 0, Y decays as exp(-t / tau).
RELAY = {"X": "X", "Y": "X"}


def assert_reproduces(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    bound = REPRODUCTION_ATOL + REPRODUCTION_RTOL * np.maximum(np.abs(x), np.abs(y))
    assert np.all(np.abs(x - y) <= bound), float(np.max(np.abs(x - y) - bound))


def model(name):
    info = ALL_MODELS[name]
    return info["rules"], info["output_nodes"]


def vertex_states(net, node):
    regs = net.reg_indices[node]
    for b in range(2 ** len(regs)):
        x = np.zeros(len(net.node_names))
        for j, r in enumerate(regs):
            x[r] = float((b >> j) & 1)
        yield b, x


def test_every_paper_network_is_available():
    assert set(PAPER_NETWORKS) <= set(ALL_MODELS)
    assert len(PAPER_NETWORKS) == 28


@pytest.mark.parametrize("name", PAPER_NETWORKS)
def test_compiled_operatorwise_rule_equals_legacy_closure_on_real_rules(name):
    rules, outputs = model(name)
    net = compile_model(rules, outputs)
    names = list(net.node_names)
    points = np.random.default_rng().random((200, len(names)))
    for i, node in enumerate(names):
        legacy = boolean_expr_to_continuous(rules[node], names)
        for x in points:
            assert net.operatorwise[i](x) == legacy(dict(zip(names, x)))


@pytest.mark.parametrize("name", PAPER_NETWORKS)
def test_operatorwise_rules_reproduce_their_truth_tables_at_every_vertex(name):
    net = compile_model(*model(name))
    assert net.operatorwise_vertex_mismatch == ()
    assert net.operatorwise_fallback == ()


@pytest.mark.parametrize("name", PAPER_NETWORKS)
@pytest.mark.parametrize("n", HILL_NS)
def test_hillcube_equals_the_boolean_rule_at_every_vertex(name, n):
    net = compile_model(*model(name))
    prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=n))
    for i in range(len(net.node_names)):
        for b, x in vertex_states(net, i):
            assert rule_targets(prep, x, np.array([i]))[0] == pytest.approx(net.truth_tables[i][b], abs=1e-12)


@pytest.mark.parametrize("name", PAPER_NETWORKS)
@pytest.mark.parametrize("construction", ["hillcube_normalized", "operatorwise_normalized", "operatorwise_legacy"])
def test_compiled_targets_equal_the_numpy_reference_inside_the_cube(name, construction):
    net = compile_model(*model(name))
    nodes = np.arange(len(net.node_names))
    for n in HILL_NS:
        dyn = DynamicsConfig(construction, hill_n=n)
        prep = prepare(net, dyn)
        for x in np.random.default_rng().random((50, len(nodes))):
            assert_reproduces(rule_targets(prep, x, nodes), rule_targets_reference(net, dyn, x, nodes))


def test_hillcube_does_not_depend_on_how_the_rule_is_written_and_operatorwise_does():
    forms = ["A & (B | C)", "(A & B) | (A & C)", "!(!A | (!B & !C))"]
    nets = [compile_model({"A": "A", "B": "B", "C": "C", "Y": f}, ["Y"]) for f in forms]
    y = nets[0].node_names.index("Y")
    points = np.random.default_rng().random((500, 4))

    def targets(construction):
        preps = [prepare(net, DynamicsConfig(construction, hill_n=2.0)) for net in nets]
        return [[rule_targets(p, x, np.array([y]))[0] for x in points] for p in preps]

    hillcube, operatorwise = targets("hillcube_normalized"), targets("operatorwise_normalized")
    assert hillcube[1] == pytest.approx(hillcube[0], abs=1e-12)
    assert hillcube[2] == pytest.approx(hillcube[0], abs=1e-12)
    assert max(abs(a - b) for a, b in zip(operatorwise[0], operatorwise[1])) > 1e-3


@pytest.mark.parametrize("n", HILL_NS)
def test_normalized_hill_maps_0_to_0_and_1_to_1_and_legacy_hill_does_not_reach_1(n):
    assert hill_normalized(0.0, n, 0.5) == 0.0
    assert hill_normalized(1.0, n, 0.5) == pytest.approx(1.0, abs=1e-15)
    assert hill(1.0, n, 0.5) == pytest.approx(1.0 / (1.0 + 0.5**n))


@pytest.mark.parametrize("n", HILL_NS)
@pytest.mark.parametrize("k", [0.3, 0.5, 0.7])
def test_normalized_hill_is_one_half_at_the_stated_half_maximum_and_not_at_k(n, k):
    assert hill_normalized(k * (1 + 2 * k**n) ** (-1 / n), n, k) == pytest.approx(0.5, abs=1e-12)
    assert hill_normalized(k, n, k) == pytest.approx((1 + k**n) / 2, abs=1e-12)


def legacy_outputs(rules, outputs, coalition, init):
    names = list(rules)
    mask = clamp_mask_for(coalition, len(names))
    return simulate_ode_output(
        rules, names, mask, 0, [names.index(o) for o in outputs], init,
        t_max=LEGACY_SOLVER.t_max, t_tail=LEGACY_SOLVER.t_tail, tau=1.0,
        hill_n=10.0, hill_k=0.5, per_solve_timeout=1e9)


@pytest.mark.parametrize("name", ["lambda_phage", "arellano_rootstem", "davidich_yeast"])
def test_legacy_construction_reproduces_the_legacy_engine_without_its_timeout(name):
    rules, outputs = model(name)
    prep = prepare(compile_model(rules, outputs), DynamicsConfig("operatorwise_legacy", hill_n=10.0))
    init = initial_states(4, len(rules), seed=42)
    for coalition in (0, 2 ** len(rules) - 1, 2 ** (len(rules) - 1) + 3):
        new = simulate_coalition(prep, LEGACY_SOLVER, coalition, init)
        assert new.count("ok") == len(init)
        assert_reproduces(new.outputs, legacy_outputs(rules, outputs, coalition, init))


def test_a_failed_trajectory_is_nan_and_the_legacy_engine_records_it_as_zero():
    rules, outputs = model("lambda_phage")
    prep = prepare(compile_model(rules, outputs), DynamicsConfig("operatorwise_legacy", hill_n=10.0))
    init = initial_states(3, len(rules), seed=42)
    starved = dataclasses.replace(LEGACY_SOLVER, attempts=(SolverAttempt("RK45", 5),))
    coalition = 2 ** len(rules) - 1
    new = simulate_coalition(prep, starved, coalition, init)
    assert new.count("eval_cap_exceeded") == len(init)
    assert np.isnan(new.outputs).all()
    assert np.isnan(new.value)

    names = list(rules)
    old = simulate_ode_output(rules, names, clamp_mask_for(coalition, len(names)), 0,
                              [names.index(o) for o in outputs], init, t_max=30.0, t_tail=10.0,
                              per_solve_timeout=1e-9)
    assert (old == 0.0).all()


@pytest.mark.parametrize("method", ["RK45", "LSODA", "Radau", "BDF"])
def test_the_evaluation_cap_binds_inside_every_solver(method):
    rules, outputs = model("lambda_phage")
    prep = prepare(compile_model(rules, outputs), DynamicsConfig("hillcube_normalized", hill_n=4.0))
    mask = clamp_mask_for(2 ** len(rules) - 1, len(rules))
    y0 = initial_states(1, len(rules), seed=7)[0]
    result = integrate(prep, dataclasses.replace(PRIMARY, attempts=(SolverAttempt(method, 5),)), mask, 0.0, y0)
    assert result.status == "eval_cap_exceeded"
    assert result.nfev == 6
    assert np.isnan(result.output)


def lambda_trajectory(attempts):
    rules, outputs = model("lambda_phage")
    prep = prepare(compile_model(rules, outputs), DynamicsConfig("hillcube_normalized", hill_n=4.0))
    mask = clamp_mask_for(2 ** len(rules) - 1, len(rules))
    y0 = initial_states(1, len(rules), seed=7)[0]
    return integrate(prep, dataclasses.replace(PRIMARY, attempts=attempts), mask, 0.0, y0)


def test_a_capped_first_attempt_runs_the_larger_budget_and_not_the_stiff_solver():
    r = lambda_trajectory((SolverAttempt("RK45", 5), SolverAttempt("RK45", 1_000_000, "after_cap"),
                           SolverAttempt("Radau", 1_000_000, "after_failure")))
    assert (r.status, r.attempt, r.method) == ("ok", 1, "RK45")


def test_a_stiff_attempt_after_failure_also_runs_after_a_cap():
    r = lambda_trajectory((SolverAttempt("RK45", 5), SolverAttempt("Radau", 1_000_000, "after_failure")))
    assert (r.status, r.attempt, r.method) == ("ok", 1, "Radau")


def test_a_later_attempt_never_runs_when_the_first_succeeds():
    r = lambda_trajectory((SolverAttempt("RK45", 1_000_000), SolverAttempt("Radau", 1_000_000, "after_failure")))
    assert (r.status, r.attempt, r.method) == ("ok", 0, "RK45")


@pytest.mark.parametrize("previous,when,runs", [
    (None, "first", True), ("ok", "after_cap", False), ("ok", "after_failure", False),
    ("eval_cap_exceeded", "after_cap", True), ("eval_cap_exceeded", "after_failure", True),
    ("solver_failed", "after_cap", False), ("solver_failed", "after_failure", True),
    ("nonfinite", "after_cap", False), ("did_not_reach_end", "after_failure", True),
])
def test_attempt_conditions(previous, when, runs):
    assert engine._should_run(when, previous) is runs


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
@pytest.mark.parametrize("threshold", [-1.0, 0.6])
def test_a_non_finite_derivative_fails_the_trajectory_under_every_solver(bad, threshold):
    net = compile_model(RELAY, ["Y"])
    prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=4.0))
    # From the start (threshold -1) or once Y passes 0.6 on its way up from 0.5.
    broken = dataclasses.replace(prep, targets=(
        prep.targets[0], lambda x, good=prep.targets[1]: bad if x[1] > threshold else good(x)))
    solver = dataclasses.replace(PRIMARY, attempts=(
        SolverAttempt("RK45", 100_000), SolverAttempt("Radau", 100_000, "after_failure")))
    r = integrate(broken, solver, np.zeros(2, dtype=bool), 0.0, np.array([0.9, 0.5]))
    assert r.status == "nonfinite"
    assert r.method == "Radau"
    assert np.isnan(r.output)


def test_the_engine_never_reads_the_clock():
    source = inspect.getsource(engine)
    for clock in ("import time", "from time", "monotonic", "perf_counter", "datetime"):
        assert clock not in source


def test_repeated_runs_agree_to_platform_precision():
    rules, outputs = model("davidich_yeast")
    prep = prepare(compile_model(rules, outputs), DynamicsConfig("hillcube_normalized", hill_n=2.0))
    init = initial_states(4, len(rules), seed=3)
    a = simulate_coalition(prep, PRIMARY, 2 ** len(rules) - 1, init)
    b = simulate_coalition(prep, PRIMARY, 2 ** len(rules) - 1, init)
    assert_reproduces(a.outputs, b.outputs)


@pytest.mark.parametrize("name", ["lambda_phage", "davidich_yeast"])
def test_substituting_clamped_nodes_equals_integrating_them_with_a_restoring_term(name):
    rules, outputs = model(name)
    net = compile_model(rules, outputs)
    prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=4.0))
    # A tail spanning the whole run, so the transient is compared as well as where it ends.
    tight = dataclasses.replace(PRIMARY, rtol=1e-10, atol=1e-12, t_max=5.0, t_tail=5.0, n_tail_samples=50)
    n = len(net.node_names)
    rng = np.random.default_rng()
    # Two or three nodes clamped, so most free nodes read other free nodes.
    mask = np.zeros(n, dtype=bool)
    while not 2 <= mask.sum() <= 3:
        mask = clamp_mask_for(int(rng.integers(0, 2**n)), n)
    y0 = np.where(mask, 0.0, rng.random(n))

    def full_state(t, y):
        x = np.clip(y, 0.0, 1.0)
        dy = np.array([f(x) for f in prep.targets]) - y
        dy[mask] = 100.0 * (0.0 - y[mask])
        return dy

    sol = solve_ivp(full_state, (0.0, tight.t_max), y0, method="RK45", max_step=tight.max_step,
                    rtol=tight.rtol, atol=tight.atol, dense_output=True)
    expected = sol.sol(np.linspace(0.0, tight.t_max, tight.n_tail_samples))
    reduced = integrate(prep, tight, mask, 0.0, y0, keep_tail=True).tail  # stored as float32
    assert reduced == pytest.approx(expected, abs=1e-6)


def test_clamped_nodes_stay_exactly_at_the_clamp_value():
    rules, outputs = model("lambda_phage")
    prep = prepare(compile_model(rules, outputs), DynamicsConfig("hillcube_normalized", hill_n=10.0))
    coalition = 0b1010101
    mask = clamp_mask_for(coalition, len(rules))
    r = simulate_coalition(prep, PRIMARY, coalition, initial_states(3, len(rules), seed=5),
                           classifier=CLASSIFIER, keep_states=True)
    assert (r.final_states[:, mask] == 0.0).all()
    assert (r.tail_means[:, mask] == 0.0).all()


def test_a_coalition_with_every_node_clamped_needs_no_integration():
    rules, outputs = model("lambda_phage")
    prep = prepare(compile_model(rules, outputs), DynamicsConfig("hillcube_normalized", hill_n=10.0))
    r = simulate_coalition(prep, PRIMARY, 0, initial_states(3, len(rules), seed=5), classifier=CLASSIFIER)
    assert r.value == 0.0
    assert (r.trajectories["nfev"] == 0).all()
    assert (r.trajectories["class"] == TRAJECTORY_CLASSES.index("fixed")).all()


def test_a_three_node_repressor_ring_oscillates_at_high_hill_coefficient_and_settles_at_low():
    net = compile_model(RING, ["A"])
    free = np.zeros(3, dtype=bool)
    y0 = np.array([0.9, 0.2, 0.4])
    high = integrate(prepare(net, DynamicsConfig("hillcube_normalized", hill_n=10.0)), PRIMARY, free, 0.0, y0, CLASSIFIER)
    low = integrate(prepare(net, DynamicsConfig("hillcube_normalized", hill_n=2.0)), PRIMARY, free, 0.0, y0, CLASSIFIER)
    assert (high.status, high.klass) == ("ok", "oscillatory")
    assert high.stats.tail_range > 0.5
    assert (low.status, low.klass) == ("ok", "fixed")


def test_a_slow_decay_is_extended_until_it_classifies_as_fixed():
    # Y = 0.5 exp(-t / 8): tail ranges ~3e-3 at t_max 60 and ~2e-5 at 120, both above the fixed-point bound; ~1e-9 at 240.
    net = compile_model(RELAY, ["Y"])
    prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=4.0, tau=8.0))
    r = integrate(prep, PRIMARY, np.array([True, False]), 0.0, np.array([0.0, 0.5]), CLASSIFIER)
    assert (r.status, r.klass, r.horizon_index) == ("ok", "fixed", 2)


def test_a_decay_too_slow_for_the_last_horizon_fails_as_unclassified_and_keeps_its_diagnostics():
    net = compile_model(RELAY, ["Y"])
    prep = prepare(net, DynamicsConfig("hillcube_normalized", hill_n=4.0, tau=40.0))
    r = integrate(prep, PRIMARY, np.array([True, False]), 0.0, np.array([0.0, 0.5]), CLASSIFIER)
    assert (r.status, r.klass, r.horizon_index) == ("unclassified", "none", 2)
    assert np.isnan(r.output)
    assert r.stats.tail_range == pytest.approx(0.5 * (np.exp(-160 / 40) - np.exp(-240 / 40)), rel=1e-3)
    # Its tail-window output survives only as a provisional value, for a labeled sensitivity analysis.
    samples = np.linspace(160.0, 240.0, round(PRIMARY.n_tail_samples * 240 / 60))
    assert r.provisional_output == pytest.approx(np.mean(0.5 * np.exp(-samples / 40)), rel=1e-4)


class Interrupted(Exception):
    pass


def test_a_sweep_interrupted_midway_resumes_to_the_uninterrupted_result(tmp_path, monkeypatch):
    rules, outputs = model("lambda_phage")
    net = compile_model(rules, outputs)
    dyn = DynamicsConfig("hillcube_normalized", hill_n=4.0)
    kwargs = dict(network="lambda_phage", n_init=2, seed=11, classifier=CLASSIFIER, coalition_range=(0, 10),
                  checkpoint_every=2)
    whole = sweep(net, dyn, PRIMARY, **kwargs)["arrays"]

    real = engine.simulate_coalition
    calls = []

    def dies_on_the_fifth_call(*args, **kw):
        calls.append(1)
        if len(calls) == 5:
            raise Interrupted
        return real(*args, **kw)

    ckpt = tmp_path / "sweep.npz"
    monkeypatch.setattr(engine, "simulate_coalition", dies_on_the_fifth_call)
    with pytest.raises(Interrupted):
        sweep(net, dyn, PRIMARY, checkpoint=ckpt, **kwargs)
    assert int(np.load(ckpt)["completed"]) == 4

    monkeypatch.setattr(engine, "simulate_coalition", real)
    resumed = sweep(net, dyn, PRIMARY, checkpoint=ckpt, **kwargs)["arrays"]
    assert set(resumed) == set(whole)
    for key in whole:
        np.testing.assert_array_equal(resumed[key], whole[key])


def test_a_checkpoint_from_a_different_configuration_is_refused(tmp_path):
    net = compile_model(*model("lambda_phage"))
    ckpt = tmp_path / "sweep.npz"
    sweep(net, DynamicsConfig("hillcube_normalized", hill_n=4.0), PRIMARY, network="lambda_phage",
          n_init=2, seed=11, classifier=CLASSIFIER, coalition_range=(0, 4), checkpoint=ckpt)
    with pytest.raises(ODEEngineError):
        sweep(net, DynamicsConfig("hillcube_normalized", hill_n=2.0), PRIMARY, network="lambda_phage",
              n_init=2, seed=11, classifier=CLASSIFIER, coalition_range=(0, 4), checkpoint=ckpt)


@pytest.mark.parametrize("coalition_range", [(0, 129), (5, 5), (-1, 3), (10, 4)])
def test_a_coalition_range_outside_the_network_is_refused(coalition_range):
    net = compile_model(*model("lambda_phage"))
    with pytest.raises(ODEEngineError, match="not inside"):
        sweep(net, DynamicsConfig("hillcube_normalized", hill_n=4.0), PRIMARY, network="lambda_phage",
              n_init=1, seed=11, classifier=CLASSIFIER, coalition_range=coalition_range)


def test_evaluation_counts_include_every_call_across_attempts_and_jacobian_estimates():
    rules, outputs = model("lambda_phage")
    prep = prepare(compile_model(rules, outputs), DynamicsConfig("hillcube_normalized", hill_n=4.0))
    calls = []

    def counted_target(x, f=prep.targets[0]):
        calls.append(1)
        return f(x)

    counted = dataclasses.replace(prep, targets=(counted_target,) + prep.targets[1:])
    solver = dataclasses.replace(PRIMARY, attempts=(SolverAttempt("RK45", 5),
                                                    SolverAttempt("Radau", 1_000_000, "after_failure")))
    r = integrate(counted, solver, np.zeros(len(rules), dtype=bool), 0.0, initial_states(1, len(rules), seed=7)[0])
    assert r.status == "ok"
    assert [(e["method"], e["when"], e["status"]) for e in r.attempt_log] == [
        ("RK45", "first", "eval_cap_exceeded"), ("Radau", "after_failure", "ok")]
    assert r.attempt_log[0]["nfev"] == 6
    assert r.nfev == sum(e["nfev"] for e in r.attempt_log)
    # Node 0 is free, so its target runs once per right-hand-side call, except the call that
    # exceeded the cap, which is counted and raises before evaluating; and once more for the final
    # derivative, a diagnostic that is not counted.
    capped, radau = r.attempt_log[0]["nfev"], r.attempt_log[1]["nfev"]
    assert len(calls) == (capped - 1) + radau + 1


def test_a_sweep_with_a_missing_coalition_cannot_be_scored():
    with pytest.raises(IncompleteSweepError):
        require_complete(np.array([0.1, np.nan, 0.3]))


SOLVER_FILE = PROJECT_ROOT / "experiments/ode_v2/solver_primary.json"
CLASSIFIER_FILE = PROJECT_ROOT / "experiments/ode_v2/classifier_primary.json"


def lambda_spec(**changes):
    fields = dict(network="lambda_phage", construction="hillcube_normalized", hill_n=10.0,
                  solver_path=SOLVER_FILE, classifier_path=CLASSIFIER_FILE, n_init=2, seed=42)
    return run_ode_arm.ArmSpec(**{**fields, **changes})


def test_the_runner_writes_a_scored_record_whose_totals_and_replicates_add_up(tmp_path):
    run_ode_arm.main(["--network", "lambda_phage", "--construction", "hillcube_normalized", "--hill-n", "10",
                      "--solver", str(SOLVER_FILE), "--classifier", str(CLASSIFIER_FILE), "--n-init", "2",
                      "--seed", "42", "--out", str(tmp_path), "--coalitions-per-shard", "50", "--keep-states"])
    stem = "lambda_phage__hillcube_normalized__n10__k0.5__c0"
    record = json.loads((tmp_path / f"{stem}.json").read_text())
    shards = sorted(tmp_path.glob(f"{stem}.shard_*.npz"))
    assert [f.name.split(".shard_")[1] for f in shards] == ["0000000-0000050.npz", "0000050-0000100.npz",
                                                            "0000100-0000128.npz"]
    status = np.concatenate([np.load(f)["trajectory_status"] for f in shards])
    s, rep = record["summary"], record["replicate"]
    assert s["n_coalitions"] == 2**7 == status.shape[0]
    assert sum(s["trajectories"].values()) == 2**7 * 2
    assert sum(s["classes_of_accepted"].values()) == s["trajectories"]["ok"]
    assert record["scored"] == (s["n_coalitions_without_value"] == 0)
    assert (rep["status_mismatches"], rep["class_mismatches"], rep["value_mismatches"]) == (0, 0, 0)
    assert rep["n_coalitions"] >= 7 + 1
    if record["scored"]:
        assert set(record["scores"]["estimators"]) == {"all_pairs", "split_half", "plain"}
        assert record["basins"]["attractor_count_histogram"]
        assert sum(record["basins"]["attractor_count_histogram"].values()) == 2**7
    assert record["configuration"]["network"] == "lambda_phage"
    assert record["configuration"]["model_sha256"] == model_digest(compile_model(*model("lambda_phage")))
    replicate = np.load(tmp_path / rep["files"][0])
    assert json.loads(str(replicate["fingerprint"]))["model_sha256"] == record["configuration"]["model_sha256"]
    np.testing.assert_array_equal(replicate["status"], status[replicate["coalitions"]])
    last = json.loads(CLASSIFIER_FILE.read_text())["extension_horizons"][-1]
    n_tails = replicate["tail_coalitions"].size
    assert 4 <= n_tails <= 4 + 12
    assert replicate["tails"].shape == (n_tails, 2, 7, round(100 * last / 60))


def test_shards_merge_to_the_unsharded_sweep(tmp_path):
    spec = lambda_spec()
    for start, end in run_ode_arm.plan_shards(7, 128):
        run_ode_arm.run_shard(spec, start, end, tmp_path / "whole")
    for start, end in run_ode_arm.plan_shards(7, 40):
        run_ode_arm.run_shard(spec, start, end, tmp_path / "split")
    whole, _, _ = run_ode_arm.merge_shards(spec, tmp_path / "whole")
    split, _, _ = run_ode_arm.merge_shards(spec, tmp_path / "split")
    assert set(whole) == set(split)
    for key in whole:
        np.testing.assert_array_equal(whole[key], split[key])


def shards_on_disk(tmp_path, ranges, spec=None):
    spec = spec or lambda_spec()
    for start, end in ranges:
        run_ode_arm.run_shard(spec, start, end, tmp_path)
    return spec


def test_a_missing_range_of_coalitions_is_refused(tmp_path):
    spec = shards_on_disk(tmp_path, [(0, 64), (96, 128)])
    with pytest.raises(run_ode_arm.ShardError, match="missing"):
        run_ode_arm.merge_shards(spec, tmp_path)


def test_overlapping_shards_are_refused(tmp_path):
    spec = shards_on_disk(tmp_path, [(0, 64), (60, 128)])
    with pytest.raises(run_ode_arm.ShardError, match="overlap"):
        run_ode_arm.merge_shards(spec, tmp_path)


def test_an_incomplete_shard_is_refused(tmp_path, monkeypatch):
    spec = shards_on_disk(tmp_path, [(64, 128)])
    real = engine.simulate_coalition
    calls = []

    def dies_on_the_fiftieth_call(*args, **kw):
        calls.append(1)
        if len(calls) == 50:
            raise Interrupted
        return real(*args, **kw)

    monkeypatch.setattr(engine, "simulate_coalition", dies_on_the_fiftieth_call)
    with pytest.raises(Interrupted):
        run_ode_arm.run_shard(spec, 0, 64, tmp_path, checkpoint_every=16)
    with pytest.raises(run_ode_arm.ShardError, match="48 of 64"):
        run_ode_arm.merge_shards(spec, tmp_path)


@pytest.mark.parametrize("changes", [{"seed": 7}, {"n_init": 4}, {"hill_n": 4.0}])
def test_a_shard_from_another_configuration_is_refused(tmp_path, changes):
    spec = shards_on_disk(tmp_path, [(0, 64)])
    other = lambda_spec(**changes)
    run_ode_arm.run_shard(other, 64, 128, tmp_path)
    # Rename the foreign shard into this spec's namespace, as a mislabeled upload would.
    foreign = run_ode_arm.shard_file(tmp_path, other, 64, 128)
    foreign.rename(run_ode_arm.shard_file(tmp_path, spec, 64, 128))
    with pytest.raises(run_ode_arm.ShardError, match="differs from the requested run"):
        run_ode_arm.merge_shards(spec, tmp_path)


def test_shards_run_under_another_solver_file_are_refused(tmp_path):
    spec = shards_on_disk(tmp_path, [(0, 128)])
    relabeled = dataclasses.replace(spec, solver_path=PROJECT_ROOT / "experiments/ode_v2/solver_legacy.json")
    with pytest.raises(run_ode_arm.ShardError, match="solver"):
        run_ode_arm.merge_shards(relabeled, tmp_path)


def test_a_shard_from_another_network_of_the_same_size_is_refused(tmp_path):
    spec = lambda_spec(network="arellano_rootstem")
    other = lambda_spec(network="asymmetric_cell_division")
    assert len(spec.net().node_names) == len(other.net().node_names) == 9
    run_ode_arm.run_shard(spec, 0, 256, tmp_path)
    run_ode_arm.run_shard(other, 256, 512, tmp_path)
    run_ode_arm.shard_file(tmp_path, other, 256, 512).rename(run_ode_arm.shard_file(tmp_path, spec, 256, 512))
    with pytest.raises(run_ode_arm.ShardError, match="model_sha256.*network"):
        run_ode_arm.merge_shards(spec, tmp_path)


def negated_rule(name):
    """The model with its first rule negated, which changes that rule's truth table at every vertex."""
    info = ALL_MODELS[name]
    rules = dict(info["rules"])
    node = next(iter(rules))
    rules[node] = f"!({rules[node]})"
    return {**info, "rules": rules}


def test_a_checkpoint_is_not_resumed_after_a_rule_changes(tmp_path, monkeypatch):
    spec = lambda_spec()
    real = engine.simulate_coalition
    calls = []

    def dies_on_the_twentieth_call(*args, **kw):
        calls.append(1)
        if len(calls) == 20:
            raise Interrupted
        return real(*args, **kw)

    monkeypatch.setattr(engine, "simulate_coalition", dies_on_the_twentieth_call)
    with pytest.raises(Interrupted):
        run_ode_arm.run_shard(spec, 0, 40, tmp_path, checkpoint_every=8)
    monkeypatch.setattr(engine, "simulate_coalition", real)
    assert int(np.load(run_ode_arm.shard_file(tmp_path, spec, 0, 40))["completed"]) == 16
    monkeypatch.setitem(ALL_MODELS, "lambda_phage", negated_rule("lambda_phage"))
    with pytest.raises(ODEEngineError, match="different configuration"):
        run_ode_arm.run_shard(spec, 0, 40, tmp_path, checkpoint_every=8)


@pytest.fixture(scope="module")
def finished_lambda_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("finished_lambda_run")
    spec = lambda_spec()
    run_ode_arm.run_shard(spec, 0, 128, out)
    for i, chunk in enumerate(run_ode_arm.plan_replicates(spec, out, 1)):
        run_ode_arm.run_replicate(spec, i, chunk, run_ode_arm.replicate_file(out, spec, i))
    return spec, out


def with_changed_json(tmp_path, source, key, value):
    data = json.loads(source.read_text())
    data[key] = value
    changed = tmp_path / source.name
    changed.write_text(json.dumps(data))
    return changed


@pytest.mark.parametrize("change", ["none", "solver tolerance", "classifier threshold", "seed", "output nodes",
                                    "truth table"])
def test_a_replicate_from_another_configuration_is_refused(finished_lambda_run, tmp_path, monkeypatch, change):
    spec, finished = finished_lambda_run
    out = tmp_path / "run"
    shutil.copytree(finished, out)
    path = run_ode_arm.replicate_file(out, spec, 0)
    coalitions = np.load(path)["coalitions"]
    variant = spec
    if change == "solver tolerance":
        variant = dataclasses.replace(spec, solver_path=with_changed_json(tmp_path, SOLVER_FILE, "rtol", 1e-7))
    elif change == "classifier threshold":
        variant = dataclasses.replace(spec, classifier_path=with_changed_json(tmp_path, CLASSIFIER_FILE,
                                                                              "halves_envelope_tol", 0.06))
    elif change == "seed":
        variant = dataclasses.replace(spec, seed=spec.seed + 1)
    elif change == "output nodes":
        info = ALL_MODELS["lambda_phage"]
        other_node = next(n for n in info["rules"] if n not in info["output_nodes"])
        monkeypatch.setitem(ALL_MODELS, "lambda_phage", {**info, "output_nodes": [other_node]})
    elif change == "truth table":
        monkeypatch.setitem(ALL_MODELS, "lambda_phage", negated_rule("lambda_phage"))
    run_ode_arm.run_replicate(variant, 0, coalitions, path)
    monkeypatch.undo()
    if change == "none":
        assert run_ode_arm.finalize(spec, out)["replicate"]["value_mismatches"] == 0
    else:
        with pytest.raises(run_ode_arm.ShardError, match="replicate_000.npz differs from the requested run"):
            run_ode_arm.finalize(spec, out)


def test_a_replicate_is_reused_only_under_the_same_fingerprint_chunk_and_coalitions(finished_lambda_run):
    spec, out = finished_lambda_run
    path = run_ode_arm.replicate_file(out, spec, 0)
    coalitions = np.load(path)["coalitions"]
    assert run_ode_arm.replicate_is_current(spec, 0, coalitions, path)
    assert not run_ode_arm.replicate_is_current(dataclasses.replace(spec, seed=spec.seed + 1), 0, coalitions, path)
    assert not run_ode_arm.replicate_is_current(spec, 1, coalitions, path)
    assert not run_ode_arm.replicate_is_current(spec, 0, coalitions[:-1], path)


def synthetic_run(n_nodes, n_init=2):
    shape = (2**n_nodes, n_init)
    return {"values": np.full(2**n_nodes, 0.5), "trajectory_status": np.zeros(shape, dtype=np.int8),
            "trajectory_class": np.full(shape, TRAJECTORY_CLASSES.index("fixed"), dtype=np.int8),
            "trajectory_horizon_index": np.zeros(shape, dtype=np.int8),
            "trajectory_attempts_run": np.ones(shape, dtype=np.int8),
            "trajectory_nfev": np.full(shape, 100, dtype=np.int32)}


def test_the_replicate_sample_takes_the_sentinels_and_every_flag_under_the_cap():
    arrays = synthetic_run(5)
    arrays["trajectory_horizon_index"][9, 0] = 1
    arrays["trajectory_attempts_run"][9, 0] = 2      # extended once, no retry
    arrays["trajectory_attempts_run"][12, 1] = 2     # retried within its first horizon
    arrays["trajectory_class"][[3, 7], 0] = TRAJECTORY_CLASSES.index("oscillatory")
    arrays["trajectory_nfev"][11, 0] = 10**6
    sample, composition = run_ode_arm.replicate_coalitions(5, {"any": "configuration"}, arrays)
    assert {0, 31, 11, 9, 12} <= set(sample.tolist())
    assert {3, 7} & set(sample.tolist())
    assert composition["flag_cap"] == 2
    assert (composition["retried_total"], composition["retried_sampled"]) == (1, 1)
    assert (composition["extended_total"], composition["extended_sampled"]) == (1, 1)
    assert (composition["oscillating_total"], composition["oscillating_sampled"]) == (2, 1)


def test_flags_beyond_the_cap_are_subsampled_the_same_way_every_time():
    arrays = synthetic_run(6)
    retried = [2, 5, 8, 13, 21, 34, 40, 41, 50, 60]
    arrays["trajectory_attempts_run"][retried, 0] = 2
    first, composition = run_ode_arm.replicate_coalitions(6, {"any": "configuration"}, arrays)
    again, _ = run_ode_arm.replicate_coalitions(6, {"any": "configuration"}, arrays)
    np.testing.assert_array_equal(first, again)
    assert composition["flag_cap"] == 3
    assert (composition["retried_total"], composition["retried_sampled"]) == (10, 3)


def test_a_run_with_a_coalition_lacking_a_value_is_not_replicated():
    arrays = synthetic_run(4)
    arrays["values"][6] = np.nan
    arrays["trajectory_status"][6, 1] = STATUSES.index("unclassified")
    sample, composition = run_ode_arm.replicate_coalitions(4, {"any": "configuration"}, arrays)
    assert sample.size == 0
    assert composition == {"skipped": "1 of 16 coalitions have no value; the run cannot be scored"}


def test_an_unscoreable_run_is_recorded_unscored_and_unreplicated(tmp_path):
    starved = with_changed_json(tmp_path, SOLVER_FILE, "attempts", [{"method": "RK45", "max_rhs_evals": 5, "when": "first"}])
    run_ode_arm.main(["--network", "lambda_phage", "--construction", "hillcube_normalized", "--hill-n", "10",
                      "--solver", str(starved), "--classifier", str(CLASSIFIER_FILE), "--n-init", "2",
                      "--seed", "42", "--out", str(tmp_path / "run")])
    record = json.loads((tmp_path / "run" / "lambda_phage__hillcube_normalized__n10__k0.5__c0.json").read_text())
    assert record["scored"] is False
    assert "scores" not in record
    assert record["summary"]["trajectories"]["eval_cap_exceeded"] == 2 * (2**7 - 1)   # all but the all-clamped coalition
    assert record["replicate"]["sample_composition"]["skipped"].startswith(f"{2**7 - 1} of {2**7} coalitions")
    assert not list((tmp_path / "run").glob("*.replicate_*.npz"))


@pytest.mark.parametrize("n_init", [1, 3])
def test_a_run_spec_needs_an_even_number_of_initial_states_for_the_split_half_estimator(n_init):
    with pytest.raises(ODEEngineError, match="even"):
        lambda_spec(n_init=n_init)


def test_the_threshold_reaches_the_dynamics_and_names_the_run(tmp_path):
    spec = lambda_spec(hill_k=0.3)
    assert spec.dynamics().hill_k == 0.3
    assert "__k0.3__" in spec.stem
    run_ode_arm.run_shard(spec, 0, 8, tmp_path)
    fingerprint = json.loads(str(np.load(run_ode_arm.shard_file(tmp_path, spec, 0, 8))["fingerprint"]))
    assert fingerprint["dynamics"]["hill_k"] == 0.3
