# Does the composition gap keep its class when each network is converted to normalized HillCube ODEs?

**Status:** DRAFT — not frozen.

Sections use the [OSF Preregistration](https://osf.io/prereg/) question titles verbatim, so
this maps onto a registration without being rewritten. A question that does not apply is
answered **N/A** with the reason, never deleted.

## Research questions or hypotheses

The paper states, in its abstract, Discussion and Conclusion, that the composition gap's sign
survives conversion of each Boolean network to Hill-function ODEs in all 25 non-null networks.
That result came from the legacy engine (`experiments/ode_v2/CONTEXT.md`). This registration
asks whether it holds under the corrected engine and the normalized HillCube conversion.

A network is **non-null** if its Boolean Δ3+ under the all-pairs estimator lies outside the
±0.5 pp null band. All Δ3+ values below are all-pairs, the primary estimator, unless stated
otherwise (`experiments/ode_v2/CONTEXT.md`, Estimators).

**H1.** Every scored non-null network has the same class, creation or destruction, under
HillCube ODE dynamics as under Boolean dynamics.

**H2.** Across the scored networks, Boolean and ODE Δ3+ have Pearson r ≥ 0.90, a predeclared
descriptive success criterion.

**H3.** Across the scored non-null networks, the median of |Δ3+(ODE) − Δ3+(Boolean)| is at most
2.0 pp, a predeclared descriptive success criterion.

**H4.** |Δ3+(ODE)| exceeds |Δ3+(Boolean)| in more scored non-null networks than not, by an
exact one-sided binomial sign test at p < 0.05.

**B1.** Every published Boolean class of the 28 networks (19 creation, 6 destruction, 3 null,
from the published values in `experiments/ode_v2/registered_identities.json`) is unchanged when
the Boolean Δ3+ is re-estimated from the same 512 trajectories with the all-pairs estimator.

**O1.** Across scored networks, the fraction of ODE trajectories classified oscillatory is
negatively associated with ODE Δ3+, by a one-sided Spearman test at p < 0.05. A secondary
mechanistic test.

**A1.** Across scored networks, the mean basin entropy of a coalition's ODE trajectories is
positively associated with ODE Δ3+, by a one-sided Spearman test at p < 0.05. A secondary
mechanistic test.

H1 carries the design: it is the claim the abstract makes. H1 holding would strengthen that
claim beyond what the published arm supports, because the HillCube conversion does not depend
on how a rule is written. The published conversion differed from the rules' multilinear
extension in 97 of 339 rules.

H4 tests the paper's reading of the published 17 of 25 as smooth dynamics amplifying
basin-boundary effects. With 25 non-null networks the test needs at least 18; the published 17
gives p = 0.054.

B1 checks the published Boolean counts against the sampling bias of the plain estimator. B1
holding would strengthen the paper's headline counts.

O1 and A1 test the two sources of the gap in the paper's Hypothesis 1: cycle averaging, which
compresses higher-order structure, and basin-boundary rerouting, which creates it. O1 is also a
replication in continuous dynamics of the Boolean association between cycling fraction and
Δ3+. The two test different mechanisms, and each is interpreted on its own one-sided p-value at
0.05; no claim rests on both holding. The Holm-adjusted p-values for the pair are reported
alongside.

Attractors are grouped as in `experiments/ode_v2/CONTEXT.md` (Attractors): within class, by
complete linkage, fixed points at 1e-3 on their tail means, oscillations at 0.05 on their
whole-period means and envelopes and within 5% in period. The grouping was validated before
this registration. A coalition's basin entropy is −Σ p_a log2 p_a over its attractors, with p_a
the fraction of its 32 trajectories reaching attractor a. A network's value is the mean over
its coalitions.

## Foreknowledge of data or evidence

Common to all three registrations: `experiments/ode_v2/CONTEXT.md`. Specific to this arm, all
read before writing it: the published ODE value for each of the 28 networks (Supplement,
"Full ODE comparison table"); class preservation in 25 of 25 non-null networks; r = 0.988;
median absolute difference 1.4 pp; larger ODE magnitude in 17 of 25. All come from the legacy
engine.

H2's threshold of 0.90 and H3's of 2.0 pp were chosen with the published 0.988 and 1.4 pp in
view. They are set looser than the published values, so that H2 and H3 test whether the
corrected engine leaves the published relation substantially intact, not whether it reproduces
it.

For O1: the published Boolean association between cycling fraction and Δ3+ (ρ = −0.83 on the 6
discovery networks, −0.37 on the 21 held out, −0.51 pooled) was read, and O1's direction is
taken from it. No ODE oscillation fraction exists for any network. Calibration only showed that
some coalitions oscillate. A1's direction comes from the paper's Hypothesis 1. The attractor
validation panel grouped the trajectories of about 32 coalitions per network and wrote only
pooled statistics: 547 of its 892 coalitions have more than one fixed-point attractor, and 19
have two or more oscillating trajectories. No basin entropy has been computed for any network.

## Explanation of foreknowledge and managing unintended influences

The rerun changes the conversion (97 rules), the solver, the clamp, the failure accounting and
the estimator, so it cannot reproduce the published values by construction. It could still be
tuned toward them. What prevents that: the engine, configuration and network definitions are
frozen by digest in `experiments/ode_v2/registered_identities.json`, and the analysis script
`scripts/analyze_ode_rerun.py` is frozen with this file. The runner prints no score. The only
ODE Δ3+ values that exist are those of the smoke test and the test suite (`lambda_phage`, 2
initial states), and nobody has read them.

## Study type

Computational experiment: deterministic simulation of fixed, published network models.

## Intention for causal interpretation

Within the models only: each coalition value is the outcome of simulated knockouts. No claim
about the organisms is intended.

## Blinding of experimental treatments

N/A — no participants or assigned treatments; every run is a deterministic computation.

## Additional blinding during research or analysis

No ODE Δ3+ is inspected until `scripts/analyze_ode_rerun.py primary` has run on all 28 records.
The runner and the Modal wrapper print only statuses and replicate agreement.

## Study design

Within-network comparison. For each of the 28 networks, the Boolean Δ3+ (512 initial states) is
compared with the ODE Δ3+ (normalized HillCube, n_H = 10, K = 0.5, clamp value 0, 32 initial
states, seed 42, frozen solver and classifier). Every one of the 2^n coalitions is computed.
The unit of analysis is the network.

## Randomization

N/A — nothing is assigned. The initial states are the fixed pseudo-random draw
`default_rng(42)`, and the replicate sample is drawn from a generator seeded by the run's
configuration digest.

## Data collection procedures

The freeze and every launch follow `experiments/ode_v2/CONTEXT.md` (Freeze and launch): the
signed tag `ode-rerun-freeze`, and a launch only from a checkout of the tagged commit where
`uv run python -m scripts.verify_freeze` passes. One Modal launch per arm:

```
uv run --with modal==1.4.3 modal run --detach -m scripts.modal_ode_arm::main \
    --run-name primary-hillcube-n10 --networks <the 28 networks> \
    --construction hillcube_normalized --hill-n 10 \
    --solver experiments/ode_v2/solver_primary.json \
    --classifier experiments/ode_v2/classifier_primary.json \
    --n-init 32 --seed 42 --keep-states
```

The legacy audit (Other planned analysis, item 2) uses `--run-name legacy-audit
--construction operatorwise_legacy --hill-n 10 --solver experiments/ode_v2/solver_legacy.json`
with no classifier. Each launch writes a local launch record (`results/ode_v2/launches/`).
Records are copied from the Modal volume into `results/` next to this file. The Boolean
estimates come from `scripts/boolean_estimators.py`.

## Data collection procedures - File upload

N/A — the procedure is the code named above, frozen by digest.

## Sample size

28 networks with 7 to 18 nodes: 128 to 262,144 coalitions each, 32 trajectories per coalition.
Under the published Boolean classes, 25 are non-null; the all-pairs Boolean classes fix the
count used here.

H1 and B1 are censuses of these networks and say nothing about networks outside the set. H4's
sign test at one-sided α = 0.05 needs at least 18 of 25, 17 of 24, or 16 of 22 or 23. H2 and H3
are predeclared descriptive success criteria, not tests. With 28 scored networks, the one-sided
Spearman tests of O1 and A1 have power 0.62 at |ρ| = 0.37, the held-out Boolean effect, and
0.86 at |ρ| = 0.5. A failure of O1 or A1 is therefore weak evidence against the mechanism, and
the paper says so.

## Sample size rationale

The networks are the paper's, all included. Each value function is computed on every coalition,
because the Walsh spectrum requires the complete function. 32 initial states match the
published ODE arm. The all-pairs estimator removes the bias their sampling error adds to the
unnormalized energies, and reduces it in the normalized fractions.

## Starting and stopping rules

Each network is run once to completion. An interrupted shard resumes from its checkpoint. A
finalized record is not rerun unless the acceptance checks refuse it, and any such rerun is
logged below the line with its cause.

## Manipulated variables

The dynamics: synchronous Boolean update against the normalized HillCube ODE. Within each, the
set of genes left unclamped: all 2^n coalitions.

## Measured variables

Per trajectory: status, class, solver attempts, accepting horizon, right-hand-side
evaluations, tail range, final derivative, output, provisional output if unclassified, and the
final state, each node's tail mean, minimum and maximum, and for an oscillation its period and
whole-period mean. Per network: the value function, the global spectra and Δ3+ under each
estimator, replicate agreement, the oscillatory fraction and the basin statistics at every point
of the tolerance grid.

## Measured variables - File upload

N/A — the record format is defined by `scripts/run_ode_arm.py`.

## Indices

Δ3+ and its class, as defined in `experiments/ode_v2/CONTEXT.md`; |Δ3+(ODE) − Δ3+(Boolean)|; a
network's oscillatory fraction (accepted trajectories classified oscillatory, over all
coalitions); and its mean basin entropy, as defined under the hypotheses.

## Indices - File upload

N/A — defined in `scripts/walsh_estimators.py` and `scripts/attractors.py`.

## Statistical models

H2: Pearson correlation of Boolean and ODE Δ3+ across scored networks, with Spearman reported
alongside. H4: exact one-sided binomial test of the number of scored non-null networks with
|Δ3+(ODE)| > |Δ3+(Boolean)| against p = 0.5. A network whose two magnitudes agree to 1e-12 counts
as a tie and is dropped from the test. O1 and A1: Spearman correlation across scored networks,
one-sided (negative for O1, positive for A1). Each is reported with its ρ, its exact number of
networks, the number of tied values on each axis, a 95% percentile interval from 10,000
network-level bootstrap resamples (seed 0), and a scatterplot (`results/figures/`), whatever its
outcome. H1, H3 and B1 are counts and a median.

## Statistical models - File upload

N/A — `scripts/analyze_ode_rerun.py` (`primary`), frozen with this file.

## Transformations

N/A — Δ3+ is used in percentage points as computed.

## Inference criteria

| hypothesis | holds when |
|---|---|
| H1 | every scored non-null network has the same class under ODE and Boolean dynamics |
| H2 | Pearson r ≥ 0.90 across scored networks |
| H3 | median absolute difference ≤ 2.0 pp across scored non-null networks |
| H4 | exact one-sided binomial p < 0.05 |
| B1 | no network's Boolean class differs between the all-pairs and plain estimators |
| O1 | one-sided Spearman p < 0.05, negative |
| A1 | one-sided Spearman p < 0.05, positive |

**H1–H4, O1 and A1 are void if more than 3 of the 28 networks cannot be scored.**

**O1 is void if fewer than 5 scored networks have any oscillating trajectory, and A1 if fewer
than 5 have any coalition with more than one attractor.**

**A scored network whose record lacks the oscillation summary or the basin statistics counts,
for O1 or A1, as a network that cannot be scored, and is named.**

**A network with any replicate mismatch is withheld from H1–H4 until its cause is found and
logged below the line.**

What the paper says under each outcome:

- **H1 holds:** the abstract's sentence stands for the corrected engine, with its denominator
  and every excluded network named ("every scored non-null network, x of y"). It is never
  stated as holding across the paper's networks while exclusions remain.
- **H1 fails:** the abstract, Discussion and Conclusion report the count and name each network
  that changes class. The claim that the effect is not an artifact of Boolean discretization is
  restricted to the networks whose class survives.
- **H2 or H3 fails:** the figure caption and the sentence reporting the median difference carry
  the new values, and the text no longer calls the two arms close.
- **H4 fails:** the sentence reading larger ODE magnitudes as amplification is removed.
- **B1 fails:** the Boolean counts are reported under the published and the all-pairs
  estimators, with each network whose class moves named. The all-pairs counts are then the
  Boolean result, corrected for sampling bias.
- **B1 is void** unless all 28 registered networks have admissible plain and all-pairs Boolean
  estimates and the plain estimate reproduces each published value to within 0.005 pp, half its
  last published digit, and in class. A void B1 is reported with the networks that caused it.
- **O1:** if it holds, the cycling association replicates under continuous dynamics; if it fails
  or is void, it does not, and the power is stated.
- **A1:** read the same way, for basin-boundary rerouting.

## Data inclusion and exclusion

All 28 networks are run. A network is excluded from H1–H4 only if it cannot be scored or is
withheld for a replicate mismatch, or if its spectrum is inadmissible under the estimator in
use (`experiments/ode_v2/CONTEXT.md`, Admissibility). Each exclusion is reported with its cause,
and an inadmissible network counts toward the void rules as one that cannot be scored. Networks in
the Boolean null band are excluded from H1, H3 and H4 and included in H2, O1 and A1. B1 covers
all 28 networks.

## Missing data

A coalition lacks a value when any of its trajectories fails, and a network with such a coalition
is not scored under H1–H4. If every failure in a network is an unclassified trajectory, a
sensitivity analysis scores it with each unclassified trajectory's provisional tail-window
output. It is labeled as such and never enters H1–H4.

## Other planned analysis

1. **Secondary estimators.** H1–H4, B1, O1 and A1 recomputed with the split-half and plain
   estimators on both sides.
2. **Legacy-engine audit.** The published conversion (`operatorwise_legacy`) with its published
   solver settings, run with the plain estimator and compared network by network with the
   published ODE table. It reports how many networks can no longer be scored because the legacy
   solver fails where the published arm recorded 0, and each scored network's change in Δ3+ and
   class. Its claim is narrow: it reconstructs the legacy engine without its timeout, and is not
   a validation of the corrected arm.
   Held fixed, and tested to reproduce the legacy engine to within 1e-12 + 1e-9 relative:
   - the operator-wise parser;
   - the Hill function with its 1e-30 guard;
   - the restoring-term clamp;
   - the initial states;
   - RK45 at rtol 1e-3, atol 1e-5, max step 2 and t_max 30;
   - the output as the clipped mean over [20, 30] at 50 samples.
   The departures, all deliberate:
   - there is no 0.1 s wall-clock timeout;
   - a deterministic cap of 10^7 right-hand-side evaluations takes its place;
   - a solver failure, a non-finite value or an exception fails the trajectory, and so leaves
     its network unscored, where the legacy engine recorded 0;
   - the right-hand side checks its own finiteness, and the solution must reach t_max;
   - it runs under Python 3.13, where the published arm ran under 3.12, with the same NumPy
     2.2.6 and SciPy 1.15.3; the BLAS may differ, which moves values by about 1e-12 relative
     (`results/audit/solver_determinism.json`).
3. **Oscillation and basins, descriptive.** Per network, the fraction of coalitions containing
   an oscillating trajectory, the fraction of multistable coalitions, and the distribution of
   attractor counts.
4. **Attractor tolerance.** A1 recomputed at each of the nine combinations of fixed-point
   tolerance (1e-4, 3e-4, 1e-3) and cycle tolerance (0.025, 0.05, 0.1), with each network's mean
   entropy and attractor-count distribution at each.
5. **Admissibility gate.** H1 recomputed from the raw energies with the total-energy gate at
   1e-7 and at 1e-5, with the networks each excludes.
6. **Provisional outputs.** The sensitivity analysis under Missing data.

Anything else is exploratory and labeled so.

## Context and additional information

**Maximum claim under this registration.** Confirmatory: whether, under the corrected engine
and the normalized HillCube conversion at n_H = 10 and K = 0.5, every non-null network among the
paper's 28 keeps its composition-gap class (H1); whether the Boolean and ODE gaps stay strongly
correlated (H2) and close in magnitude (H3); and whether ODE magnitudes are larger more often
than chance allows (H4); whether the published Boolean classes survive the correction for
sampling bias (B1); and, as secondary mechanistic tests, whether across networks the gap falls
with the fraction of oscillating trajectories (O1) and rises with basin entropy (A1). A claim of
robustness across Hill coefficients, thresholds or conversion syntax belongs to
`experiments/2026-09-21_ode-sensitivity/PREREG.md`, and a claim about perturbation strength to
`experiments/2026-09-21_graded-perturbation-rerun/PREREG.md`. A causal claim that oscillation or
basin geometry produces the gap is not confirmatory: O1 and A1 are associations across 28
networks.

This registration supersedes the published ODE arm for the rerun. The published values remain
the record of what the legacy engine produced, and the legacy-engine audit measures how they
move without its timeout.

---

## Log

Append only. Never edit above the line.

The last column is what distinguishes an amendment from a deviation, so you do not have to
decide which word to use: `nothing run`, `no results seen`, `results not opened`, `results seen`.

```
2026-09-21  created                              nothing run
```
