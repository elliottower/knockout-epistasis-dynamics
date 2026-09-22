# Shared context for the ODE rerun registrations

Three registrations govern the rerun of the paper's continuous-dynamics arms on the corrected
engine, one per arm:

| registration | arm |
|---|---|
| `experiments/2026-09-21_ode-primary-rerun/PREREG.md` | HillCube conversion of all 28 networks at knockout; the legacy-engine audit |
| `experiments/2026-09-21_graded-perturbation-rerun/PREREG.md` | clamp levels 0.25, 0.5, 0.75 and 1 |
| `experiments/2026-09-21_ode-sensitivity/PREREG.md` | Hill coefficient, threshold and conversion syntax, networks with at most 12 nodes |

This file holds what they share. Each registration freezes with this file in the same commit,
and each names it; a change to this file after the freeze is a change to all three.

## Engine and frozen configuration

The trajectories come from `scripts/ode_engine.py`, run by `scripts/run_ode_arm.py` and, on
Modal, `scripts/modal_ode_arm.py`. Their design was reviewed in six rounds before this
registration (`~/Documents/backup/knockout-ode-review-{3,4,5,6}-2026-09-21/`, and two earlier
rounds recorded in `HOUSE_STYLE/docs/oversight.log`).

`experiments/ode_v2/registered_identities.json`, written by
`scripts/record_registered_identities.py` at the freeze, records:

- `code_sha256`, the digest of the three modules whose source determines a trajectory
  (`scripts/ode_engine.py`, `grn_coalition_sweep.py`, `scripts/ode_coalition_sweep.py`);
- the SHA-256 of `solver_primary.json` and `classifier_primary.json`;
- each network's `model_sha256`: ordered node names, output nodes, and each rule's text,
  regulators and truth table.

Every shard, replicate and record carries the same digests, and a run whose digests differ from
this file is not a run under these registrations. This trajectory-code digest is distinct from
the launch manifest's digest, which covers every file shipped to Modal, the wrappers included,
and is recorded per run.

Solver (`solver_primary.json`): RK45 with a cap of 200,000 right-hand-side evaluations; RK45
with 2,000,000 only if the cap was hit; Radau with 2,000,000 after a failure of any kind.
Radau and BDF both succeeded on all 24 forced trajectories of the benchmark. BDF was about twice
as fast. Radau was chosen because it is an implicit Runge-Kutta method of order 5 that is A- and
L-stable, while BDF above order 2 is not A-stable; the benchmark did not decide between them.
rtol 1e-8, atol 1e-10, max step 2, t_max 60, tail the final 20 time units sampled at 100 points.

Classifier (`classifier_primary.json`): a trajectory is fixed if its tail range is at most 1e-5
and its final derivative at most 1e-6; oscillatory if its tail range exceeds 1e-4, it stays in
[0, 1] to within 1e-3, and each moving node's tail halves repeat their maximum and minimum to
within 2.5% of that node's range. Otherwise it is continued to t = 120, 240, 480 and 960, the
tail kept at the final third and at constant sampling density. A trajectory still unclassified
at 960 fails. An accepted oscillation also records its period and each node's mean over a whole
number of periods (`scripts/ode_engine.cycle_summary`), which do not depend on the phase at which
the tail starts.

## Networks, dynamics and initial states

The 28 paper networks, as defined by `scripts/ode_coalition_sweep.ALL_MODELS` and pinned by the
model digests. Node j is bit j of the coalition index. Genes outside a coalition are clamped:
they are removed from the integrated state and their clamp value substituted.

Each free node follows dx_i/dt = (F_i(x) - x_i)/τ with τ = 1, where F_i is the **normalized
HillCube** of its Boolean rule (Wittmann et al. 2009, eq. 5): each regulator's value passes
through h(z) = z^n (1 + K^n)/(z^n + K^n), and F_i is the multilinear interpolation of the rule's
truth table at those values. F_i equals the Boolean rule at every vertex and does not depend on
how the rule is written. K is the parameter of the Hill function, not its half-maximum; the
half-maximum of h is K (1 + 2K^n)^(-1/n):

| n_H | K = 0.3 | K = 0.5 | K = 0.7 |
|---|---|---|---|
| 2 | 0.276 | 0.408 | 0.498 |
| 4 | 0.299 | 0.486 | 0.635 |
| 10 | 0.300 | 0.500 | 0.696 |

Initial states: `numpy.random.default_rng(42).random((32, n))`, the published ODE arm's, shared
by every coalition of a network. Clamped components are set to the clamp value.

## Value function

A trajectory's output is the mean, over the tail window of the horizon at which it was
accepted, of its output nodes clipped to [0, 1], averaged over output nodes. A coalition's value
is the mean over its 32 trajectories and exists only if all 32 are accepted. A network is scored
only if every coalition has a value.

## Estimators

Walsh coefficients are w = WHT(v)/2^n. Energy at order k is summed over coalitions T with
popcount k, and fractions are normalized by the total including order 0, as in the published
spectra. Δ3+ = E3+(global) − E3+(local). The local spectrum comes from the rules' truth tables
and has no sampling error.

The plain estimator squares the coefficients of the mean over the initial states. Their sampling
error then adds to every coefficient's energy, mostly at high orders, and inflates the order-3+
fraction. Because the same initial states are used for every coalition, the error is correlated
across coalitions and no constant floor can be subtracted. With w_i the coefficients of the
value function computed from initial state i alone, the states are independent draws, so
E[w_i(T) w_j(T)] = w(T)^2 for i ≠ j.

- **all-pairs** (primary). U(T) = ((Σ_i w_i(T))^2 − Σ_i w_i(T)^2)/(m(m − 1)), the average of
  w_i w_j over every ordered pair of distinct states: a U-statistic.
- **split-half** (secondary). w_A(T) w_B(T), with A and B the means over initial states 0-15
  and 16-31 (0-255 and 256-511 for the Boolean arm).
- **plain** (secondary), the published estimator.

Summed over an order, all-pairs and split-half estimate the **unnormalized** order energies
without bias. A normalized fraction is a ratio of two such estimates and is not unbiased; it is
used to reduce the plain estimator's upward sampling bias. Its finite-sample behavior is checked
by the registered comparisons between the three estimators, and in
`tests/test_walsh_estimators.py` on a value function of known spectrum with state-shared noise.
Both unbiased estimates matched the true unnormalized order-3+ energy there, while the plain
estimate overshot it by more than 100 standard errors. All-pairs varied somewhat less than
split-half: the two share their leading variance term, and all-pairs reduces the term that is
pure noise. Estimated order energies can be negative at an order with little true energy, and are
reported as estimated. Implemented in `scripts/walsh_estimators.py`.

**Admissibility**, predeclared in that module. A spectrum is admissible only if its estimated
total energy, the value function's mean square, is at least 1e-6, a root-mean-square value of
0.001. A network whose spectrum is inadmissible under an estimator is left out of that estimator's
tests, and counts toward the void rules as a network that cannot be scored. The ratio σ3+/σ2+ is
admissible only if the order-3+ and order-2+ energies are each at least 1e-4 of the total, so it
never divides by a nonpositive or near-zero quantity. Admissibility decides whether a normalized
spectrum or ratio is numerically interpretable. It excludes no network on biological grounds.
Both thresholds are judgment, not calibration. Every record reports each estimator's raw order
energies, total, and ratio numerator and denominator, and every exclusion is listed by network,
estimator and level. As a sensitivity analysis, each arm's class-preservation or G4 result is
recomputed from the raw energies with each gate one decade looser and one decade stricter
(`scripts/analyze_ode_rerun.py`, `TOTAL_ENERGY_GRID` and `RATIO_FRACTION_GRID`).

For the Boolean reference, the same three estimators are applied to the 512 per-trajectory
Boolean outputs (`scripts/boolean_estimators.py`). These come from the cycle-handling audit,
which recomputed all 28 networks. For the 27 with a committed coalition table it reproduces the
table exactly. `grieco_bladder` (18 nodes) has no committed table; its recomputed plain value,
3.7655 pp, rounds to the published 3.77. The plain estimator reproduces all 28 published values to
within 0.005 pp, half their last published digit (largest gap 0.00495, `asymmetric_cell_division`).
The published values are registered in `registered_identities.json`.

**Classes.** Creation if Δ3+ > +0.5 pp, destruction if Δ3+ < −0.5 pp, null otherwise: the
paper's null band. A claim that a class is preserved always states its denominator and names
every excluded network. It is never phrased as holding across the paper's networks while
exclusions remain.

## Attractors

A coalition's trajectories are grouped into attractors by `scripts/attractors.py`, within class
and by complete linkage, so that no two trajectories in one attractor are farther apart than the
tolerance. Fixed points are compared on their tail means, at 1e-3. Oscillations are compared on
their period, within 5% of the longer, and on each node's whole-period mean, tail minimum and
tail maximum, at 0.05. Basin entropy is −Σ p_a log2 p_a over a coalition's attractors, in bits.
Every record reports the basin statistics at these tolerances and at each of the nine
combinations of fixed-point tolerance 1e-4, 3e-4, 1e-3 and cycle tolerance 0.025, 0.05, 0.1.

The grouping was validated before registration (`scripts/validate_attractor_grouping.py`,
`results/audit/attractor_grouping_validation.json`). The panel had four controls with known
attractors (toggle switch, 4-node ring, 3-node ring, and a switch selecting a 3- or 5-node ring)
and 892 coalitions, about 32 from each of the 28 networks. Its gate was set before it ran, and
the grouping passed it:

| measured | value |
|---|---|
| controls with their known number of attractors, at every grid point | 4 of 4 |
| coalitions with more than one fixed-point attractor | 547 of 892 |
| coalitions with two or more oscillating trajectories | 19 |
| partition identical at all nine grid points | 892 of 892 |
| complete and single linkage identical | 892 of 892 |
| largest distance within one fixed-point attractor | 4.2e-7 |
| smallest distance between two fixed-point attractors | 0.49 |
| largest spread of a cycle summary within one attractor | 0.010 |
| smallest distance between two cycle attractors | 0.34 |
| largest relative period spread within one cycle attractor | 0.0029 |

A first run used cycle tolerances 0.01, 0.025 and 0.05 and failed its first criterion: at 0.01
the switched-rings control split each of its two cycles in two
(`results/audit/attractor_grouping_validation_cycle_grid_0.01-0.05.json`). The sampled cycle
summaries differ by up to 0.012 on one cycle, so the grid was moved to 0.025–0.1, inside the
range where the partition does not change.

## Runs and acceptance

Runs execute on Modal. Each network is split into shards of 4,096 coalitions, checkpointed every
64 coalitions and resumable. A shard, replicate or record is accepted only if its
configuration, including the network, model and code digests, equals the one the requested
spec produces from the files as frozen. A deterministic replicate sample is recomputed in
separate containers: the three sentinel coalitions, 5% of each coalition size, retried and
extended coalitions up to a cap of 5% of the network, and 5% of coalitions with an
oscillating trajectory. It must agree with the sweep in status, class and value to within
1e-12 + 1e-9·max(|x|, |y|). A network that cannot be scored is not replicated.

## Foreknowledge common to all three

- **The published continuous-dynamics results, all read.** They come from the legacy engine:
  an operator-wise conversion, rtol 1e-3, t_max 30, and a 0.1 s wall-clock timeout under which
  failures and exceptions were recorded as output 0. They are the paper's Supplement ODE table,
  class preservation 25 of 25 non-null networks, r = 0.988, median |ΔODE − ΔBool| 1.4 pp,
  magnitudes larger under ODE in 17 of 25, and the graded table with its sentences.
- **The published Boolean results for all 28 networks, read**, and the 28-network audit of
  their cycle handling (`results/audit/boolean_cycle_handling/`), which found the implemented
  and registered procedures identical everywhere.
- **The static conversion audit** (`results/audit/operatorwise_vs_multilinear.json`): the
  operator-wise conversion differs from the multilinear one in 97 of 339 rules across 24
  networks, so the published and the rerun arms measure different continuous extensions.
- **Engine development integrated trajectories without computing any Walsh quantity.** This
  covers the calibration of the classifier, the oscillation validation panel, the solver
  benchmarks and the timing of every network (`results/audit/`), each sampling at most 200
  coalitions per network, most 8 to 24. These runs showed trajectory classes: some coalitions
  oscillate, and some approach fixed points slowly.
- **The attractor validation panel integrated 892 coalitions and wrote only pooled
  statistics**, listed under Attractors. No value, attractor count or basin entropy was written
  or read per network.
- **One Δ3+ has been computed and not read.** The Modal smoke test finalized `lambda_phage`
  at 2 initial states, and the local test suite finalizes it in temporary directories. The
  smoke report omits scores; nobody has opened the records that hold them.

The predictions in the three registrations are therefore made with the published
legacy-engine numbers in view, and cannot be blind to them. What the registrations fix before
any rerun number exists is the estimator, the decision rules and the handling of networks that
cannot be scored.

## Freeze and launch

The freeze is one sequence, checked by `scripts/verify_freeze.py`:

1. Commit the code, this file, the configuration files and `registered_identities.json`,
   regenerated by `scripts/record_registered_identities.py` from that tree.
2. Commit the three registrations alone. This is the plan commit.
3. Run `prereg freeze` in each registration's directory. Each header names the plan commit.
4. Commit the three freeze headers alone.
5. Tag that commit `ode-rerun-freeze` with a signed annotated tag, and push the branch and the
   tag.

Every launch runs from a checkout of the tagged commit. `scripts/verify_freeze.py` must pass
there first, and `scripts/modal_ode_arm.py` runs the same check itself and refuses to launch
otherwise. The check requires:

- the three registrations frozen at one plan commit, each passing `prereg check`;
- HEAD equal to the tagged commit, the tag annotated with a signature that verifies, and origin
  holding the same tag object;
- nothing but the three registrations changed between the plan commit and HEAD;
- every file the registrations rely on tracked at the plan commit and unchanged since, including
  `pyproject.toml`, `uv.lock`, the Modal wrappers, which define the container image, the analysis
  code, the verifier itself, and `results/audit/boolean_cycle_handling/`, whose per-trajectory
  shards are the Boolean arm's input;
- no modified or untracked file under a relevant path;
- `registered_identities.json` equal to the identities computed from the tree.

The launch records the tag and both commits in every record's launch manifest. The analysis
refuses a record whose manifest names another commit.

## Compute

Measured single-core time per trajectory, 16 coalitions × 4 initial states per network
(`results/audit/engine_timing_all_networks.json`), gives about 270 core-hours for one full sweep
of the 28 networks. Grieco accounts for 128 of them and Calzone (full) for 48. The graded arm is
four such sweeps and the sensitivity arm five sweeps of the 17 networks with at most 12 nodes
(about 12.5 core-hours each). The legacy-engine audit is about a fifth of a full sweep, timed on
four networks (`results/audit/engine_timing_legacy_audit.json`). The Modal price read on
2026-09-21 was $0.0000131 per core-second and $0.00000222 per GiB-second.
