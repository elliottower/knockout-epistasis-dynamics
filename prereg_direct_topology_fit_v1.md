# Pre-registration: does topology predict the composition gap, without any predictor?

**Frozen before the regression is run. No fit has been computed.**

## Why this replaces a weaker test

The paper currently establishes "the gap is not predictable from standard
topological features" from the **failure of pre-registered predictions** --
10/21 = 48% direction accuracy, indistinguishable from chance.

That inference has a hole, stated in `PREREG_METHODOLOGY.md`: a
prediction-based null is ambiguous between *the features carry no information*
and *the predictor was poor*. Blinding cannot separate them.

A direct fit removes the predictor from the argument entirely.

## Features — taken from the paper, not chosen by me

The paper names its predictors: **feedback loop counts, rule complexity,
in-degree distribution**. Those are the feature families fitted here, so there
is no feature-selection freedom to exploit:

| family | concrete features |
|---|---|
| in-degree distribution | mean, max, variance of in-degree; edge density |
| feedback loops | count of positive and negative feedback cycles up to length 4; total cycle count |
| rule complexity | mean and max number of pairwise and triple Fourier terms per gene; fraction of genes with any nonlinear term |
| size (nuisance) | number of nodes |

## Target

$\Delta_{3+}$ = (global order-3+ spectral energy) − (local-rule order-3+
spectral energy), per network, from `energy_spectrum` in the committed
`results/grn_v2/*_composition_blind.json`.

## Predictions

**H1 (primary).** A multiple regression of $\Delta_{3+}$ on the full feature set
achieves **leave-one-out cross-validated $R^2 < 0.2$**. The features do not
carry usable information about the gap.

**H2 (primary).** No single feature reaches $|{\rm Spearman}\ \rho| \geq 0.5$
with $\Delta_{3+}$ after Holm correction across the feature set.

**H3, the falsifier that matters.** If cross-validated $R^2 \geq 0.4$, **the
published null is an artifact of poor predictions, not a property of the
systems.** The paper's claim that the gap is unpredictable from topology is
then wrong and must be retracted and replaced with the fitted relationship.
This outcome is stated first so it cannot be reinterpreted afterwards.

**H4 (secondary).** Node count alone does not predict the gap
($|\rho| < 0.3$), consistent with the matched random-network control already
reported.

## Decision rule, fixed in advance

| outcome | consequence |
|---|---|
| H1 and H2 hold | The null is established by data rather than by prediction failure. The prediction exercise is demoted to corroboration and the paper's claim is strengthened. |
| **H3 fires** | The published claim is wrong. Report the fit, retract the unpredictability claim, and state that the pre-registered prediction failure reflected the predictor rather than the systems. |
| H1 holds, H2 fails | One feature carries signal the joint fit cannot use. Report both; the unpredictability claim is narrowed to "no usable joint predictor". |

## Not varied

The 27 networks with committed composition and wiring JSONs, the $\Delta_{3+}$
definition already used in the paper, and the feature families as named in the
manuscript. Leave-one-out cross-validation because $n = 27$ is small.

## Script

`scripts/direct_topology_fit.py`. To be written; nothing has been run.
