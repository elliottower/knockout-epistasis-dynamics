# Pre-registration: does 1% recovery of epistatic coefficients generalise?

**Frozen before any recovery sweep is run on the 28 systems.**

## The claim being tested, which is not ours

Faure, Lehner, Miró Pina, Serrano Colome and Weghorn (*PLOS Comput Biol*
20(5):e1012132, 2024) subsampled a combinatorially complete landscape **to 1% of
variants**, fitted LASSO on a Walsh basis, and recovered epistatic coefficients at
**$r = 0.99$** against known ground truth.

Their landscape was **simulated**, and it was one landscape. This project holds 28
systems with exact $2^n$ coalition tables, $n = 7$--18, spanning three orders of
magnitude in combinatorial complexity, derived from published models of real
pathways.

**The question is therefore not whether 1% recovery is possible -- that is
settled -- but whether it generalises, and where it fails.**

## Design

Replicate their protocol without modification: uniform random subsample to budget
$b \in \{1, 5, 10, 25, 50\}\%$ of the $2^n$ coalitions, LASSO on the Walsh basis,
correlate recovered coefficients against the exact spectrum. 50 repeats per
(network, budget) for intervals.

Their protocol is used **as published**, not improved. This is a generalisation
test, and changing the method would confound it.

## Predictions

**H1 (primary).** Median recovery at $b = 1\%$ across the 28 systems is **below
0.99**, and the spread across systems is wide -- interquartile range exceeding
0.10. A single-landscape result is unlikely to hold uniformly across three orders
of magnitude of complexity.

**H2 (primary).** Recovery **degrades with $n$**: negative rank correlation
between $n$ and recovery at fixed budget. Larger systems have exponentially more
coefficients at fixed sampling fraction.

**H3 (primary), the one that connects the two projects.** Recovery degrades with
the **order-3+ energy fraction**: systems with a large composition gap are harder
to recover at low budget. If this holds, the two projects are linked by a
mechanism neither predicted -- the systems where dynamics create the most
higher-order structure are exactly the ones where sparse recovery fails.

**H4 (secondary).** The additive baseline is not competitive at any budget. If
additivity matches LASSO-Walsh recovery, the recovery result is uninformative --
mandatory since the MULTI-evolve dispute.

**H5, the null outcome.** If recovery exceeds 0.95 on every system at 1%, Faure's
result generalises completely, there is nothing to add, and this arm is reported
as a successful replication and dropped.

## Decision rule, fixed in advance

| outcome | consequence |
|---|---|
| H1, H2, H3 hold | A characterisation of **when** sparse recovery works, with the failure mode tied to the composition gap. This is the version worth publishing, and it cites Faure as the result being extended rather than competing with it. |
| H1 and H2 hold, H3 fails | Recovery degrades with size but not with interaction structure. Useful and practical; no link to the composition gap, and the two projects stay separate. |
| **H5 fires** | Complete generalisation. Report as replication, drop the arm, and do not claim a budget contribution anywhere. |
| H4 fails | The comparison is uninformative at any budget and the whole arm is withdrawn. |

## Explicit non-claim

**No version of this arm may claim priority for recovery from 1% of
combinations.** That number is Faure's. The v12 title asserted it; that assertion
is being removed. Any wording implying we established low-budget recoverability is
a factual error about the literature.

## Not varied

The 28 coalition tables, the Walsh basis, Faure's published protocol, and the
budget grid.

## Script

To be written. Hash recorded in `docs/FROZEN_SHAS.md` before running.

---

# Amendment 1 — before anything ran

**No script exists, no hash was recorded, and no estimator has been run against
these coalition tables.** Nothing below is a revision in light of data; it is a
revision in light of a review, made while the arm is still unexecuted.

## A1.1 H1 was unfalsifiable and is replaced

As written, H1 predicted "median recovery below 0.99, IQR > 0.10". **0.99 is a
ceiling**, so any imperfection whatsoever confirms it. A prediction that cannot
lose is not a prediction.

Replaced with:

> **H1 (primary, revised).** Median Pearson $r$ at $b = 1\%$, over systems with
> $n \geq 13$, falls **below 0.85**.
>
> **Faure generalises cleanly** if median $r \geq 0.95$ with interquartile range
> below 0.10. That outcome is stated here so it cannot be reinterpreted later:
> it means the arm is a successful replication, contributes nothing beyond it,
> and is dropped.

## A1.2 The 1% budget is degenerate at small $n$ — a design flaw in both arms

A percentage budget over $2^n$ collapses to almost nothing on the small systems:

| $n$ | $2^n$ | coalitions at $b = 1\%$ |
|---|---|---|
| 7 | 128 | **1** |
| 10 | 1,024 | 10 |
| 13 | 8,192 | 81 |
| 18 | 262,144 | 2,621 |

At $n = 7$ a 1% budget is a **single coalition**. Any recovery statistic there is
noise, and pooling it into a median across 28 systems would drag the headline
number down for a reason that has nothing to do with recovery difficulty.

Corrections, applying to this arm **and** to
`prereg_estimator_headtohead_v1.md`, which inherits the same budget grid:

1. **Report absolute sample count alongside every relative budget.** The budget
   axis is $(b, 2^n \cdot b)$, never $b$ alone.
2. **The 1% condition is read only over $n \geq 13$** (81+ coalitions). Smaller
   systems are reported at that budget but excluded from the pooled 1% statistic,
   and the exclusion is stated in the figure.
3. **Add an absolute-budget axis** -- $\{64, 256, 1024, 4096\}$ coalitions -- run
   alongside the relative one. It is the axis a practitioner actually has ("I can
   afford 1,000 measurements"), and it makes systems of different $n$ comparable
   in a way percentages do not.

## A1.3 H3 is deliberately duplicated, and that is the point

H3 here (recovery degrades with the order-3+ energy fraction) is the **same
hypothesis** as H4 in `prereg_estimator_headtohead_v1.md`. This is designed
triangulation, not accidental duplication: H3 tests it with one estimator
(LASSO), H4 with a five-estimator panel.

If it holds under both, the link between the composition gap and recovery
difficulty is **estimator-independent** rather than an artifact of L1 penalisation
— which is the difference between a Discussion paragraph and a footnote. If it
holds under one and not the other, the estimator is doing the work and the
cross-project claim is withdrawn.
