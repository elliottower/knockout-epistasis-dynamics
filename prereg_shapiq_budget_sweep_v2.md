# Pre-registration: shapiq approximators in EpistasisBench budget sweep

## Status: FROZEN — see prereg_shapiq_freeze_record.md for hash

## Contamination statement

This pre-registration is **sequential, not blind**. H1 and H2 are
informed by prior LASSO-Walsh and iRF results (results/full_sweep_v1.json)
and are explicitly labeled as such. The reasoning is comparative:
shapiq performance is predicted relative to established LASSO/iRF
baselines, so full blindness to those baselines is neither achievable
nor appropriate.

H3, H4, the max_order design decision, and the falsification criteria
were derived blind by an independent agent with access only to term-count
math and method descriptions — no LASSO/iRF numbers, no diagnostic
performance, no circuit-specific energy spectra.

This is a different standard from a fully blind pre-registration. We
state it explicitly so readers do not assume independence that does
not hold.

## Question

Do budget-limited Shapley interaction approximators (KernelSHAP-IQ,
SHAPIQ MC, SVARM-IQ) recover interaction structure from exhaustive
coalition tables more or less accurately than sparse Walsh recovery
(LASSO) and nonparametric tree methods (iRF)?

## max_order design decision (derived blind)

At n=15, the number of interaction terms by max_order:

| max_order | Terms | Budget 1% (327) | 2% (655) | 3% (983) | 5% (1638) | 10% (3276) | 20% (6553) | 40% (13107) |
|-----------|-------|------------------|----------|----------|-----------|------------|------------|-------------|
| 2         | 120   | 2.7              | 5.5      | 8.2      | 13.7      | 27.3       | 54.6       | 109.2       |
| 3         | 575   | 0.57             | 1.14     | 1.71     | 2.85      | 5.70       | 11.40      | 22.79       |
| 4         | 1940  | 0.17             | 0.34     | 0.51     | 0.84      | 1.69       | 3.38       | 6.76        |

Ratios are budget / n_terms. Below 1.0 = underdetermined. Stable
regression (WLS) typically requires ratio > 5. MC methods benefit
from evaluation sharing across terms, so the ratio is a pessimistic
bound for them, but ratios well below 1 still produce noise-dominated
estimates.

**Decision.** The core sweep uses max_order in {2, 3}:

- **max_order=2 (120 terms)**: Overdetermined at all budget levels
  (ratio >= 2.7). The regime where all methods should produce
  meaningful estimates.
- **max_order=3 (575 terms)**: Underdetermined at 1% (0.57), barely
  determined at 2-3% (1.1-1.7), comfortable at 5%+ (2.85+). The
  bias-variance crossover — where method differences should be
  most informative.
- **max_order=4 (1940 terms)**: Underdetermined at every budget
  through 10% (ratio 1.69 at 10%). First reaches ratio > 3 at 20%.
  Included at 20% and 40% budgets ONLY, as a boundary test for
  the feasibility ceiling. Excluded at lower budgets because
  running an experiment you expect to fail by construction is
  manufacturing a foregone result, not testing a hypothesis.
  **Confound note:** because max_order=4 is only observed at high
  budgets, any order-4 effect is fully confounded with budget.
  Results cannot distinguish "order 4 works well" from "order 4
  works well when you have 6553+ samples."

This decision was derived from first-principles sample complexity
(the ratio table above), not from any diagnostic performance data.

## Scoring

Primary metric: held-out R² (same as LASSO/iRF). For shapiq methods,
reconstruction uses the plug-in predictor:

  v_hat(S) = baseline + sum_{T ⊆ S, |T| > 0} phi_T

This is exact at max_order=n and approximate at lower orders.

Secondary: interaction-index MSE against exact k-SII values computed
via shapiq.ExactComputer. Scores approximators in their native space
rather than through reconstruction.

**Pre-designated held-out set.** To ensure all methods are scored
on identical test coalitions, a fixed held-out set (10% of 2^15 =
3277 coalitions, drawn once per seed before any method runs) is
withheld from all methods. LASSO and iRF draw their training sample
from the remaining 90%. Shapiq approximators may query any coalition
in the remaining 90% via smart sampling. All methods are scored on
the same pre-designated held-out coalitions. This removes the
confound that smart samplers choose informative coalitions, leaving
the held-out set as a systematically less-informative remainder.

## Predictions

### Analytically known ceilings

The k-SII reconstruction ceiling (exact k-SII values, plug-in
reconstruction, measured R² on all coalitions) is computable via
shapiq.ExactComputer and is NOT a prediction:

| Circuit | k-SII order≤2 | k-SII order≤3 |
|---------|---------------|---------------|
| weight_ioi | 0.757 | 0.868 |
| canonical_ioi | 0.955 | 0.988 |
| random15 | 0.995 | 1.000 |

For reference, iRF achieves R² ~0.93 at 5%+ budget (from prior
LASSO/iRF sweep). The ceiling exceeds iRF on canonical_ioi and
random15 at order 2, but falls well below iRF on weight_ioi.

### H1 (informed by prior LASSO/iRF results, not blind)

Shapiq approximators will reach within 5% of their circuit-specific
k-SII ceiling at sufficient budget (>= 5%), testing estimation
efficiency rather than the ceiling itself. Specifically:

- On weight_ioi (ceiling 0.757 at order 2): shapiq achieves
  R² ~0.72-0.76, well below iRF's ~0.93. The gap is dominated
  by the truncation ceiling, not estimation error.
- On canonical_ioi (ceiling 0.955 at order 2): shapiq achieves
  R² ~0.91-0.96, comparable to or exceeding iRF. Near-additive
  circuits have high ceilings at low truncation orders.
- On random15 (ceiling 0.995 at order 2): shapiq achieves
  R² ~0.95-0.99, exceeding iRF. Null circuits are trivially
  reconstructed at any order.

The prediction is that the estimation efficiency (achieved R² /
ceiling R²) exceeds 0.95 at 5%+ budget for KernelSHAP-IQ at
max_order=2.

### H2 (informed by prior results, not blind)

At max_order=3 with sufficient budget (>= 5%), the best shapiq
approximator will approach but not exceed iRF on weight_ioi
(ceiling 0.868 vs iRF ~0.93), and will match or exceed iRF on
canonical_ioi and random15 (ceilings 0.988 and 1.000).

On weight_ioi, truncation at order 3 costs 6.2 R² points
relative to iRF (ceiling 0.868 vs iRF ~0.93). The 13.2% of
reconstructable variance beyond order 3 (1 - 0.868) is
inaccessible to any order-3-truncated method.

### H3 (derived blind, with falsification criteria)

At low-to-mid budgets (1-10%), KernelSHAP-IQ achieves strictly
higher held-out R² than SVARM-IQ, which in turn achieves strictly
higher R² than SHAPIQ MC. At high budgets (20-40%), the three
methods converge (R² differences below 0.05).

Rationale: KernelSHAP-IQ solves a structured WLS problem that
imposes global consistency — estimated interactions must jointly
explain observed coalition values — functioning as implicit
regularization. SVARM-IQ applies stratification and variance
reduction over naive MC. At high budgets, all converge at standard
MC rates.

Qualifier: at max_order values where the system is underdetermined,
KernelSHAP-IQ may lose its advantage or underperform MC if it
overfits the sampled coalitions. This prediction applies to cells
with overdetermination ratio >= 1.

**Falsification criteria:**

1. H3 is falsified if SHAPIQ MC achieves higher held-out R² than
   KernelSHAP-IQ at budget <= 10% and max_order <= 3 on at least
   2 of 3 circuits.
2. H3 is falsified if SHAPIQ MC achieves higher held-out R² than
   SVARM-IQ at budget <= 10% and max_order <= 3 on at least 2 of
   3 circuits.
3. H3 is falsified if the R² gap between the best and worst method
   exceeds 0.05 at budget = 40% and max_order <= 3 on at least 2
   of 3 circuits.

### H4 (derived blind, with falsification criteria)

The ordinal ranking KernelSHAP-IQ > SVARM-IQ > SHAPIQ MC is
consistent across circuits at any fixed (budget, max_order) pair
where the system is determined. The magnitude of gaps varies —
circuits with more high-order energy show larger absolute R²
differences because the estimation problem is harder.

Rationale: estimator properties (bias, variance, convergence rate)
depend on budget, term count, and estimator design, not on specific
game values. Ordering should be invariant to the circuit.

**Qualifying cells.** The sweep contains 16 (budget, max_order)
cells total: 7 at order 2, 7 at order 3, 2 at order 4. Of these,
11 have overdetermination ratio > 3 (6 at order 2 for budgets 2-40%,
3 at order 3 for budgets 10-40%, 2 at order 4 for budgets 20-40%).
The R² > 0.3 filter is expected to pass most of these (the hardest
circuit, weight_ioi, has a k-SII ceiling of 0.757 at order 2),
yielding at least 9 qualifying cells. The threshold of 3 reversals
is achievable.

**Falsification criteria:**

1. H4 is falsified if the ordinal ranking of the three methods
   (by median held-out R²) differs across circuits at the same
   (budget, max_order) for at least 3 of the (budget, max_order)
   cells where all methods achieve R² > 0.3.
2. The ranking reversal must occur at a cell with overdetermination
   ratio > 3 (ruling out the underdetermined qualifier from H3).

## Methods

Three shapiq approximators, all using k-SII index:
- KernelSHAP-IQ (regression-based, Fumagalli et al. 2024)
- SHAPIQ MC (Monte Carlo, Fumagalli et al. 2023)
- SVARM-IQ (stratified Monte Carlo, Kolpaczki et al. 2024)

## Protocol

1. Finalize this prereg (incorporate review feedback)
2. SHA-freeze and timestamp
3. Run full sweep:
   - 7 budgets x 20 seeds x 3 circuits x 3 approximators
   - max_order in {2, 3} at all budgets
   - max_order = 4 at 20% and 40% budgets only
4. Score all methods on held-out R² and interaction-index MSE
5. Evaluate H1-H4 against falsification criteria

## Circuits

Same three as the LASSO/iRF sweep:
- weight_ioi (interactive)
- canonical IOI (near-additive)
- random15 (null)

## Risks and limitations

**Noise amplification in plug-in reconstruction.** The reconstruction
sums C(|S|,1) + ... + C(|S|, max_order) interaction estimates per
coalition. At max_order=3 and |S|=10, this is 175 terms. Independent
estimation errors with variance sigma² produce reconstruction variance
175 * sigma². Even unbiased estimates can yield poor R² through
aggregation. The optimal max_order is a non-trivial function of
budget that we cannot fully predict in advance.

**Universal failure at low budgets.** At 1% budget and max_order >= 3,
all methods may produce R² near or below zero. Method comparisons in
this regime are uninformative. Conclusions about method ordering are
drawn only from cells where at least one method achieves R² > 0.1.

**Metric mismatch.** The methods estimate k-SII interaction indices,
not coalition values. Held-out R² evaluates the latter through the
plug-in reconstruction pipeline. A method could estimate interactions
accurately but reconstruct poorly (if errors align constructively
when summed), or vice versa. Results evaluate the methods as
components of an estimate-then-reconstruct pipeline, not as
interaction estimators per se.

**Truncation bias confounds method comparison.** At max_order < n,
all methods hit the same R² ceiling set by truncation. For circuits
with significant energy above the truncation order, method differences
are compressed against this shared ceiling. Running max_order = n
(32767 terms) would remove this confound but is infeasible at any
budget level.

## What this does NOT test

- Whether shapiq methods are useful for their intended purpose (ranking
  interactions, feature attribution). They may be excellent at ranking
  while poor at value reconstruction.
- Mean-ablation ground truth (deferred).
- Other interaction indices (STII, FSII) — only k-SII tested.
