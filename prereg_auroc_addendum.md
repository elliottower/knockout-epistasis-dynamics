# Pre-registration addendum: Pairwise Interaction AUROC

## Status: FROZEN
SHA-256: 93ebe466c644a4d0fd7a5bd29044a4a6531787ddbfa42f1e4052a439ea4281a9
(Hash computed before this status line was added; verify on the DRAFT version.)

## Contamination statement

This addendum was written by an automated agent (Claude Opus 4.6)
that has NOT seen any pairwise AUROC results, any reconstruction R²
results, any shapiq sweep output, or any per-trial metrics from any
method on any circuit. The agent has seen:

- The metric implementation (`data_utils.py`, `compute_recovery_metrics`)
- The energy spectra and k-SII reconstruction ceilings (from the
  frozen shapiq pre-registration)
- The method source code (`lasso_walsh_oracle.py`,
  `irf_interaction_discovery.py`, `shapiq_approximators.py`,
  `budget_sweep_harness.py`, `shapiq_sweep_harness.py`)
- The existing pre-registrations (`prereg_shapiq_budget_sweep_v2.md`,
  `prereg_grn_v4.md`)

Predictions are derived from structural properties of the methods
(estimation target, parameter counts, basis transforms) and the
circuits (energy spectra, k-SII ceilings). No pilot data, no
diagnostic runs, no empirical tuning of thresholds.


## Metric definition

For each of the C(15,2) = 105 head pairs, extract the true and
recovered pairwise Walsh coefficient magnitudes:

    true_score[i,j]  = |w_true[(1 << i) | (1 << j)]|
    pred_score[i,j]  = |w_recovered[(1 << i) | (1 << j)]|

Binary labels: pairs with true score >= median(true scores) are
positive (truly interacting). This gives a ~53/52 split.

AUROC = area under the ROC curve of pred_score predicting binary
labels. AUROC = 1.0 means perfect separation; AUROC = 0.5 means
no detection ability.

This is a ranking metric. It is invariant to monotone
transformations of pred_score and insensitive to calibration
errors. It tests whether a method can distinguish strong pairwise
interactions from weak ones, not whether it recovers accurate
magnitudes.


## Structural reasoning

### What each method estimates

| Method | Native output | Path to pairwise Walsh |
|--------|--------------|----------------------|
| LASSO-Walsh | Walsh coefficients (order <= 4) | Direct: w_recovered[pair_idx] IS the LASSO coefficient |
| iRF | v_hat(S) for all 2^15 coalitions | WHT of forest predictions on full table |
| KernelSHAP-IQ | k-SII values (order <= max_order) | Reconstruct v_hat from k-SII, then WHT(v_hat) |
| SHAPIQ-MC | k-SII values (order <= max_order) | Same as KernelSHAP-IQ |
| SVARM-IQ | k-SII values (order <= max_order) | Same as KernelSHAP-IQ |

LASSO is the only method whose native output is the quantity being
scored. All other methods require one or two transforms.

### Parameter counts and sample ratios

At n=15 with budget fraction f, the available samples are
m = f * 32768 (minus 10% held out for LASSO/iRF but not shapiq).

| Method | Parameters | 1% (m=295*) | 5% (m=1474*) | 10% (m=2949*) | 40% (m=11796*) |
|--------|-----------|-------------|--------------|----------------|-----------------|
| LASSO (order<=4) | 1941 | 0.15 | 0.76 | 1.52 | 6.08 |
| shapiq (order=2) | 121 | 2.71** | 13.5** | 27.1** | 108.2** |
| shapiq (order=3) | 576 | 0.57** | 2.85** | 5.69** | 22.8** |

*LASSO/iRF receive 90% of the budget (held-out set); shapiq
receives 100% via smart sampling. Values shown are for LASSO/iRF
(0.9 * f * 32768).

**Shapiq ratios use full budget (f * 32768) because shapiq methods
choose their own coalitions.

LASSO is underdetermined (ratio < 1) at 1-3% budget. L1
regularization handles this by selecting a sparse subset of the
1941 features, but pairwise coefficients compete with 1365
order-3 and 455 order-4 features for selection. At low budgets,
some true pairwise signals may be zeroed out in favor of spurious
higher-order terms.

Shapiq at order=2 is overdetermined at all budgets (ratio >= 2.7).
The pairwise k-SII estimates are stable even at 1% budget.

### Basis transform distortion

Shapiq methods estimate k-SII values, not Walsh coefficients. The
pairwise Walsh coefficients used for AUROC scoring are derived
by: (1) reconstructing v_hat(S) from k-SII via plug-in prediction,
(2) computing WHT(v_hat). This two-step transform introduces a
specific kind of error: the order-truncated k-SII reconstruction
misses all k-SII terms above max_order. When the reconstruction is
transformed back to Walsh space, the truncation error projects
onto ALL Walsh orders, including pairwise.

The severity depends on the circuit's interaction structure:

| Circuit | k-SII ceiling order<=2 | k-SII ceiling order<=3 | Pairwise energy frac |
|---------|----------------------|----------------------|---------------------|
| weight_ioi | 0.757 | 0.868 | 0.7% |
| canonical_ioi | 0.955 | 0.988 | 1.5% |
| random15 | 0.995 | 0.999 | 0.2% |

For weight_ioi, the order-2 reconstruction misses 24.3% of
variance. This missing variance projects into pairwise Walsh
coefficients, potentially distorting the ranking. For canonical_ioi
and random15, the reconstruction is nearly complete at order 2.

### Signal-to-noise ratio per circuit

Pairwise energy as a fraction of total (non-intercept) energy:

- canonical_ioi: 1.5% / 5.5% = 27% of non-intercept energy is pairwise
- weight_ioi: 0.7% / 2.5% = 28% of non-intercept energy is pairwise
- random15: 0.2% / 0.7% = 29% of non-intercept energy is pairwise

These ratios are similar, but absolute magnitudes differ by an
order of magnitude. For AUROC (ranking), what matters is the
ratio of inter-class separation (gap between the 53rd and 54th
ranked pair) to estimation noise. Circuits with larger absolute
pairwise magnitudes tolerate more estimation noise before the
ranking is scrambled.


## Predictions

### HA1: LASSO achieves the highest pairwise AUROC at medium-high budgets

At budgets >= 10%, LASSO-Walsh achieves the highest median
pairwise AUROC among all five methods on canonical_ioi and
weight_ioi, for any max_order setting used by the shapiq methods.

**Rationale.** LASSO directly estimates Walsh coefficients. The
pairwise magnitudes used for AUROC scoring are literally the
LASSO's own output — no basis transform, no reconstruction, no
WHT inversion. At 10% budget, LASSO has ~2949 samples for 1941
features (ratio 1.52), enough for L1 regularization to reliably
identify the largest pairwise terms. The shapiq methods' two-step
transform (k-SII -> reconstruction -> WHT) introduces truncation
error that LASSO avoids entirely. iRF's forest predictions must
generalize to unseen coalitions before WHT extraction, adding
extrapolation noise.

**Qualifier.** This prediction excludes random15, where all methods
may cluster near AUROC = 0.5 (see HA3), making method ordering
uninformative.

**Falsification.** HA1 is falsified if any shapiq method or iRF
achieves higher median pairwise AUROC than LASSO at budget >= 10%
on both canonical_ioi and weight_ioi.


### HA2: Shapiq order-2 is competitive with LASSO at low budgets

At 1-3% budget, the best shapiq method at max_order=2 achieves
pairwise AUROC within 0.05 of LASSO on canonical_ioi and
weight_ioi.

**Rationale.** At 1% budget, LASSO has ~295 samples for 1941
features (ratio 0.15) — severely underdetermined. L1
regularization selects a sparse subset, but the 105 pairwise
features compete with 1365 order-3 features for selection. Some
true pairwise signals may be zeroed out or attenuated. Shapiq at
order=2 has ~328 samples for 121 parameters (ratio 2.7) — well
overdetermined. The pairwise k-SII estimates should be more
stable, and the truncation-induced bias in pairwise Walsh
estimates is small for canonical_ioi (ceiling 0.955) and moderate
for weight_ioi (ceiling 0.757).

The favorable sample complexity of shapiq order-2 compensates for
the basis transform indirection, at least partially.

**Falsification.** HA2 is falsified if LASSO exceeds the best
shapiq order-2 method by more than 0.10 AUROC at budget <= 3% on
both circuits. This would mean LASSO's direct Walsh estimation
dominates even in the regime where it is severely underdetermined.


### HA3: Per-circuit ordering of pairwise AUROC

Across all methods and budgets >= 5%, the median pairwise AUROC
follows:

    canonical_ioi > weight_ioi > random15

**Rationale.** canonical_ioi has the largest absolute pairwise
energy (1.5% of total) and the strongest non-intercept structure
(5.5% total). The pairwise magnitudes are large enough that
estimation noise at >= 5% budget should not scramble the ranking.

weight_ioi has moderate pairwise energy (0.7%), less than half of
canonical_ioi. The pairwise ranking is detectable but noisier.

random15 has very weak pairwise energy (0.2%). With 99.3% of
energy in the intercept, the pairwise coefficients are small in
absolute terms. The median split creates labels from what are
effectively noise-floor magnitudes. Most methods should achieve
AUROC near 0.5 on random15, especially at low budgets.

**Specific threshold for random15.** At 1-3% budget, all five
methods achieve median AUROC < 0.65 on random15. The weak signal
is indistinguishable from estimation noise at low sample counts.

**Falsification criteria.**

1. HA3 ordering is falsified if random15 achieves higher median
   AUROC than weight_ioi (averaged across methods) at any budget
   >= 5%.
2. The random15 threshold is falsified if any method achieves
   median AUROC > 0.75 on random15 at budget <= 3%. This would
   mean the pairwise signal in the null circuit is stronger than
   the energy spectrum implies.


### HA4: AUROC saturates faster than R-squared with budget

For each method, define the AUROC saturation budget as the
smallest budget at which median AUROC reaches 95% of the
method's AUROC at 40% budget. For R², define the analogous
saturation budget.

On canonical_ioi, the AUROC saturation budget is <= the R²
saturation budget for all five methods.

**Rationale.** AUROC is a ranking metric: once a method places the
top 53 pairs above the bottom 52, adding budget cannot improve
AUROC further. R² requires accurate magnitudes for all 2^15
coalitions, which demands progressively more samples. Ranking
the 105 pairwise coefficients is a lower bar than reconstructing
32768 values.

This prediction is strongest for LASSO, whose L1 regularization
produces a sparse solution — the largest pairwise coefficients
emerge first, and the ranking stabilizes before the magnitudes
converge. For shapiq methods, the basis transform complicates
saturation: reconstruction R² and pairwise AUROC may improve at
different rates because they depend on different aspects of the
k-SII estimates.

**Falsification.** HA4 is falsified if, for at least 3 of 5
methods on canonical_ioi, the AUROC saturation budget exceeds
the R² saturation budget. This would mean pairwise ranking is
harder to recover than function reconstruction.


### HA5: Max-order effect on shapiq pairwise AUROC

**HA5a (weight_ioi, high budget).** At budgets >= 10%, shapiq
methods at max_order=3 achieve higher pairwise AUROC than
max_order=2 on weight_ioi, for all three shapiq methods.

**Rationale.** The k-SII ceiling gap between order 2 and order 3
on weight_ioi is large (0.757 vs 0.868). The order-2
reconstruction misses 24.3% of variance; order-3 misses only
13.2%. This missing variance projects into the pairwise Walsh
coefficients through the WHT of the reconstruction. Reducing
the truncation error from 24.3% to 13.2% should reduce the
distortion of pairwise rankings.

At >= 10% budget, the order-3 system is well-determined (ratio
5.7 at 10%), so estimation variance does not overwhelm the bias
reduction.

**HA5b (weight_ioi, low budget).** At 1-3% budget, shapiq methods
at max_order=2 achieve higher pairwise AUROC than max_order=3 on
weight_ioi.

**Rationale.** At 1% budget, order-3 has 576 parameters for 328
samples (ratio 0.57) — underdetermined. The estimation noise from
the underdetermined order-3 k-SII system propagates through the
reconstruction and WHT into pairwise Walsh estimates, scrambling
the ranking. Order-2 with 121 parameters (ratio 2.7) produces
cleaner estimates despite higher truncation bias.

**HA5c (canonical_ioi and random15).** The AUROC difference
between max_order=2 and max_order=3 is <= 0.03 (in absolute
AUROC) at all budgets on canonical_ioi and random15.

**Rationale.** The ceiling gaps are small: canonical_ioi (0.955
vs 0.988 = 3.3 percentage points), random15 (0.995 vs 0.999 =
0.4 percentage points). The order-2 reconstruction already
captures enough of the function for the pairwise Walsh ranking
to be accurate. Adding order-3 terms provides minimal additional
information about pairwise structure.

**Falsification criteria.**

1. HA5a is falsified if max_order=2 achieves higher AUROC than
   max_order=3 for at least 2 of 3 shapiq methods on weight_ioi
   at budget >= 10%.
2. HA5b is falsified if max_order=3 achieves higher AUROC than
   max_order=2 for at least 2 of 3 shapiq methods on weight_ioi
   at budget <= 3%.
3. HA5c is falsified if the absolute AUROC difference between
   orders exceeds 0.05 for any shapiq method at any budget on
   canonical_ioi or random15.


### HA6: Shapiq inter-method ordering for pairwise AUROC

At fixed (budget, max_order, circuit), the ordinal ranking of the
three shapiq methods by median pairwise AUROC matches the
reconstruction R² ranking from the shapiq pre-registration (H3):

    KernelSHAP-IQ >= SVARM-IQ >= SHAPIQ-MC

**Rationale.** Pairwise AUROC is a monotone function of the quality
of the pairwise Walsh estimates, which are derived from the k-SII
estimates through reconstruction and WHT. Better k-SII estimation
(higher R²) should produce better pairwise Walsh estimates and
therefore higher AUROC. The ordering among shapiq methods reflects
their estimation efficiency (KernelSHAP-IQ's WLS > SVARM-IQ's
stratified MC > SHAPIQ-MC's naive MC), which should affect AUROC
and R² in the same direction.

**Qualifier.** AUROC is less sensitive than R² to systematic biases
(all pairs shifted by the same amount). If two shapiq methods
differ only in systematic bias (not in ranking quality), their R²
may differ but their AUROC may be equal. The prediction is that
the ordering holds, not that the AUROC gaps match the R² gaps.

**Falsification.** HA6 is falsified if the ordinal ranking differs
(e.g., SHAPIQ-MC > KernelSHAP-IQ) at more than 3 of the
(budget, max_order) cells where all three methods achieve AUROC
> 0.55 on at least one circuit.


### HA7: iRF pairwise AUROC profile

iRF achieves lower pairwise AUROC than LASSO at budgets >= 10% on
canonical_ioi and weight_ioi, and higher AUROC than the worst
shapiq method (at max_order=2) at budgets >= 5%.

**Rationale.** iRF's path to pairwise Walsh coefficients runs
through: (1) train a forest on sampled coalitions, (2) predict on
all 2^15 coalitions (extrapolation for unseen ones), (3) compute
WHT of predictions. The forest's predictions on unseen coalitions
are extrapolations that accumulate error. When the WHT is applied
to these noisy predictions, the pairwise coefficients pick up the
extrapolation noise.

At high budgets, the forest has seen enough coalitions that its
predictions are accurate on most of the space, and the WHT-derived
pairwise estimates are reliable. But LASSO's direct estimation
avoids the extrapolation entirely.

Relative to shapiq methods, iRF has the advantage of not being
truncated to any max_order. But it lacks the sample efficiency of
shapiq order-2 at low budgets (iRF has many more effective
parameters than 121). At medium-high budgets, iRF's flexible
nonparametric model should outperform the worst shapiq method.

**Falsification.** HA7 is falsified if iRF achieves higher AUROC
than LASSO at budgets >= 10% on both canonical_ioi and weight_ioi.
This would mean the forest's nonparametric flexibility outweighs
LASSO's direct estimation advantage for ranking.


## Summary of predicted ordering by budget regime

### Low budget (1-3%)

    canonical_ioi: LASSO ~ shapiq_o2 >> shapiq_o3 >> iRF
    weight_ioi:    shapiq_o2 ~ LASSO >> shapiq_o3 >> iRF
    random15:      all ~ 0.5 (no detectable signal)

Shapiq order-2 benefits from overdetermined estimation (121
parameters). LASSO is underdetermined (1941 parameters) but L1
selection rescues the strongest signals. iRF forests are poorly
trained with < 1000 samples. Shapiq order-3 is underdetermined
(576 parameters) and noise-dominated.

### Medium budget (5-10%)

    canonical_ioi: LASSO > iRF ~ shapiq_o2 ~ shapiq_o3
    weight_ioi:    LASSO > shapiq_o3 ~ iRF > shapiq_o2
    random15:      LASSO > rest (small absolute differences)

LASSO becomes well-determined and its direct estimation advantage
dominates. On weight_ioi, order-3 shapiq starts to outperform
order-2 as the bias reduction exceeds the variance cost. iRF
becomes competitive as the forest has sufficient training data.

### High budget (20-40%)

    canonical_ioi: LASSO > iRF ~ shapiq_o3 > shapiq_o2
    weight_ioi:    LASSO > iRF ~ shapiq_o3 >> shapiq_o2
    random15:      all methods may exceed 0.5 modestly

All methods converge, but LASSO retains its direct-estimation
advantage. Shapiq order-3 matches or exceeds order-2 as estimation
variance vanishes. The gap between methods shrinks because the
estimation problem is well-conditioned for everyone.


## Risks and limitations

**Median-split sensitivity.** The median split creates a ~53/52
partition of 105 pairs. If the true pairwise magnitude distribution
is heavy-tailed (a few very large pairs, many near-zero), the
AUROC task is easy: just detect the few large ones. If the
distribution is more uniform (gradual decline from the largest to
smallest), the task is harder because pairs near the median are
difficult to classify. The shape of this distribution has not been
examined.

**AUROC floor on random15.** If the true pairwise magnitudes on
random15 are all effectively zero (indistinguishable from
numerical noise in the WHT), then the median split creates
arbitrary labels and AUROC is meaningless for all methods. This is
a feature, not a bug — the null circuit SHOULD produce AUROC near
0.5 — but it makes method comparisons on random15 uninformative.

**Basis transform linearity.** The k-SII -> Walsh pipeline is
linear (reconstruction is a sum of k-SII terms, WHT is linear).
This means errors in k-SII estimates propagate linearly into Walsh
coefficient estimates. However, the AUROC scoring function is
nonlinear (it depends on the rank order of magnitudes, not
magnitudes themselves). Small linear errors in Walsh estimates can
produce discontinuous changes in ranking, making AUROC predictions
harder to calibrate than R² predictions.

**LASSO shrinkage bias.** L1 regularization shrinks all
coefficients toward zero. The pairwise magnitudes used for AUROC
scoring are the shrunken values. If shrinkage is approximately
uniform across pairs (all shrunk by roughly the same amount), the
ranking is preserved and AUROC is unaffected. If shrinkage is
differential (smaller coefficients shrunk proportionally more,
which is typical of L1), the ranking may be distorted. This
distortion tends to IMPROVE AUROC (large pairs are even more
separated from small ones), creating a mild bias in favor of LASSO
on this metric.

**Held-out set asymmetry.** LASSO and iRF receive a fixed random
sample (minus 10% held-out), while shapiq methods choose their own
coalitions via smart sampling. The smart sampling is optimized for
k-SII estimation, not for pairwise Walsh AUROC. This design is
correct for a fair budget comparison, but the shapiq methods are
solving a different optimization than pairwise Walsh detection.
