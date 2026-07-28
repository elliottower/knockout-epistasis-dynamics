# Pre-registration: shapiq approximators in EpistasisBench budget sweep

## Status: DRAFT — awaiting Perplexity review before execution

## Question

Do budget-limited Shapley interaction approximators (KernelSHAP-IQ,
SHAPIQ Monte Carlo, SVARM-IQ) recover interaction structure from
exhaustive coalition tables more or less accurately than sparse Walsh
recovery (LASSO) and nonparametric tree methods (iRF)?

## Background

The existing LASSO vs iRF sweep (results/full_sweep_v1.json) established:
- LASSO-Walsh has near-zero absolute error that keeps dropping with budget
- iRF has an invariant ~0.05-0.11 error floor due to axis-aligned splits
- Both methods receive random coalition samples and predict held-out values

Shapiq approximators are a fundamentally different class: they choose
which coalitions to evaluate (smart sampling) and estimate Shapley
interaction indices rather than Walsh coefficients. The comparison tests
whether smart sampling + interaction-index estimation outperforms random
sampling + basis-aligned sparse recovery.

## Key design decision: max_order

At n=15, the number of interaction terms by max_order:
- order 2: 120 terms (15 + 105)
- order 3: 575 terms (+ 455)
- order 4: 1940 terms (+ 1365)

Budget at 5% = 1638 evaluations. At max_order=4, the system is
underdetermined (1940 terms > 1638 samples). Preliminary diagnostic
(not yet complete) is testing whether max_order=2 or 3 gives better
reconstruction than max_order=4.

LASSO-Walsh uses max_order=4 with L1 regularization, which handles
underdetermination via sparsity. Shapiq approximators have no explicit
sparsity prior, so they may need max_order matched to budget.

## Scoring

Primary metric: held-out R² (same as LASSO/iRF).

For shapiq methods, reconstruction uses the plug-in predictor:
  v_hat(S) = baseline + sum_{T subset S, |T|>0} phi_T

This is exact at max_order=n and approximate at lower orders.
Reconstruction quality depends on both the approximation accuracy of
individual interaction values AND the error amplification from summing
~100+ noisy terms per coalition.

Secondary: interaction-index MSE against exact k-SII values computed
via shapiq.ExactComputer. This scores the approximators in their native
space rather than through reconstruction.

## Predictions

H1: At max_order=2, shapiq approximators will score BETWEEN LASSO and
iRF on held-out R² at matched budgets. Rationale: smart sampling should
outperform random sampling when the budget is sufficient relative to
the number of terms, but the truncation to order 2 caps their ceiling
below LASSO's max_order=4.

H2: At max_order=4, shapiq approximators will score WORSE than iRF
(negative R²). Rationale: underdetermined system + error amplification
in reconstruction. The approximators estimate all 1940 terms with
1638 samples and no sparsity prior, then sum ~100 noisy terms per
coalition.

H3: KernelSHAP-IQ will outperform SHAPIQ MC and SVARM-IQ at low
budgets. Rationale: regression-based methods are more sample-efficient
than Monte Carlo methods when the number of terms is small relative
to budget.

H4: The ordering will be consistent across circuits. The gap between
methods reflects algorithmic capacity, not circuit-specific structure.

## Methods

Three shapiq approximators:
- KernelSHAP-IQ (regression-based, Fumagalli et al. 2024)
- SHAPIQ MC (Monte Carlo, Fumagalli et al. 2023)
- SVARM-IQ (stratified Monte Carlo, Kolpaczki et al. 2024)

Index: k-SII (matching existing shapiq analyses in the weight paper)
Max orders tested: 2, 3, 4

## Protocol

1. Determine feasible max_order from diagnostic (this prereg)
2. Get Perplexity review of the prereg
3. Run full sweep: 7 budgets x 20 seeds x 3 circuits x 3 approximators
4. Score all methods on the same held-out R² metric

## Circuits

Same three as the LASSO/iRF sweep:
- weight_ioi (interactive, k99=241)
- canonical IOI (near-additive, k99=142)
- random15 (null, k99=34)

## What this does NOT test

- Approximation accuracy of the interaction indices themselves (secondary)
- Whether shapiq methods are useful for their intended purpose (ranking
  interactions, feature attribution). They may be excellent at ranking
  interactions while being poor at value reconstruction.
- Mean-ablation ground truth (deferred, per weight paper)
