# How does the composition gap change with perturbation strength under the corrected ODE engine?

**Status:** FROZEN at `78071efa918c`
**Plan sha256:** `37cf619ef64b277a52f9abaf348ec73a6412a02e628651561b13a355ea623701`
**Frozen:** 2026-09-22
**Log:** 3 entries, head `ec381430`

Sections use the [OSF Preregistration](https://osf.io/prereg/) question titles verbatim, so
this maps onto a registration without being rewritten. A question that does not apply is
answered **N/A** with the reason, never deleted.

## Research questions or hypotheses

The paper states that creation dominance is specific to complete knockout and inverts by
half-strength knockdown, that the gap persists at partial perturbation, and that perturbation
strength does not drive the normalized higher-order fraction. Those results came from the legacy
engine under `prereg_graded_perturbation_v1.md` and its two amendments. This registration asks
whether they hold under the corrected engine.

The clamp level f is the value perturbed genes are clamped to: f = 0 is complete knockout and
f = 1 full activation. Classes, the null band, the estimators and their admissibility are as in
`experiments/ode_v2/CONTEXT.md`. All Δ3+ values are all-pairs, the primary estimator, unless
stated otherwise.

**G1.** At f = 0, scored networks with creation outnumber those with destruction.

**G2.** At f = 0.5, scored networks with destruction outnumber those with creation.

**G3.** At f = 0.25, at most 2 scored networks lie in the null band.

**G4.** The median σ3+/σ2+, over networks whose ratio is admissible at both f = 0 and f = 0.75,
differs by at most 0.05 between the two levels.

**G5.** At every level other than f = 0, at most a quarter of the scored networks lie in the null
band.

G2 carries the design: it is the abstract's claim that creation dominance inverts by
half-strength knockdown. G1 and G2 holding together would strengthen that claim, now on a
conversion that does not depend on how rules are written and with 28 networks instead of 27.

G3 and G5 together test the paper's sentence that the gap persists at partial perturbation: G3
at quarter strength, where the published arm had no network in the null band, and G5 at every
level. G5's threshold was set with the published counts in view, so G5 holding shows that the
corrected engine keeps the published pattern. It does not independently confirm the published
persistence claim.

G4 is the claim that perturbation strength does not drive the higher-order fraction. Its
alternative is the v1 registration's H2, which predicted that σ3+/σ2+ rises by a factor of at
least 1.5 from f = 0.75 to f = 0. That factor is reported whatever G4's outcome.

## Foreknowledge of data or evidence

Common to all three registrations: `experiments/ode_v2/CONTEXT.md`. Specific to this arm, all
read before writing it, from the legacy engine (`results/graded_perturbation_scored.json`, the
paper's graded table and text):

- creating/destroying counts of 17/7, 15/12, 8/16, 7/15 and 6/15 at f = 0, 0.25, 0.5, 0.75 and 1;
- ten of 27 networks changing sign, eight of them by f = 0.25;
- 0, 3, 5 and 6 of 27 networks in the null band at f = 0.25, 0.5, 0.75 and 1;
- median σ3+/σ2+ of 0.384 at f = 0 and 0.387 at f = 0.75, a factor of 0.99;
- the v1 verdicts: H1, H2 and H3 failed, and the H5 falsifier did not fire.

G1–G5 restate the published claims. Their thresholds were set with these values in view. G3's
limit of 2, G4's of 0.05 and G5's quarter (7 of 28) are looser than the published 0, 0.003 and
at most 6 of 27. They are predeclared descriptive success criteria, not inferential tests.

## Explanation of foreknowledge and managing unintended influences

As for the primary arm: the engine, configuration and networks are frozen by digest;
`scripts/analyze_ode_rerun.py` is frozen with this file; the runner prints no score; and no
graded value exists under the corrected engine. The f = 0 sweep is the primary arm's run, so
this arm cannot choose its own reference level.

## Study type

Computational experiment: deterministic simulation of fixed, published network models.

## Intention for causal interpretation

Within the models only: the perturbation level is set by intervention in the simulation. No
claim about the organisms is intended.

## Blinding of experimental treatments

N/A — no participants or assigned treatments; every run is a deterministic computation.

## Additional blinding during research or analysis

No graded Δ3+ is inspected until `scripts/analyze_ode_rerun.py graded` has run on every record
at every level.

## Study design

Within-network comparison across five clamp levels, f ∈ {0, 0.25, 0.5, 0.75, 1}, with every
other setting as in the primary arm: normalized HillCube, n_H = 10, K = 0.5, 32 initial states,
seed 42, and the frozen solver and classifier. The f = 0 level is the primary arm's sweep, reused
and not recomputed. All 28 networks, every coalition. The unit of analysis is the network at a
level.

## Randomization

N/A — nothing is assigned; the initial states and the replicate sample are fixed as in the
primary arm.

## Data collection procedures

The freeze and every launch follow `experiments/ode_v2/CONTEXT.md` (Freeze and launch). A
launch runs only from a checkout of the tagged commit, after
`uv run python -m scripts.verify_freeze` passes there. Four Modal launches, one per level, with run
names `graded-f0.25`, `graded-f0.5`, `graded-f0.75` and `graded-f1`:

```
uv run --with modal==1.4.3 modal run --detach -m scripts.modal_ode_arm::main \
    --run-name graded-f<level> --networks <the 28 networks> \
    --construction hillcube_normalized --hill-n 10 --clamp-value <level> \
    --solver experiments/ode_v2/solver_primary.json \
    --classifier experiments/ode_v2/classifier_primary.json \
    --n-init 32 --seed 42 --keep-states
```

Each run's records, and nothing else, are copied from the Modal volume into
`results/<run name>/` next to this file, for example:

```
uv run --with modal==1.4.3 modal run -m scripts.modal_ode_arm::fetch --run-name graded-f0.25 \
    --out experiments/2026-09-21_graded-perturbation-rerun/results/graded-f0.25
```

## Data collection procedures - File upload

N/A — the procedure is the code named above, frozen by digest.

## Sample size

28 networks at five levels. The published arm had 27; `arabidopsis_cellcycle` is added because
the corrected engine covers it. G1 and G2 compare two counts at one level each, G3 and G5 are
counts, and G4 compares medians over at least 10 networks. None generalizes beyond these
networks.

## Sample size rationale

Every paper network and every coalition, as in the primary arm. The levels are the published
arm's, so the rerun answers the same questions.

## Starting and stopping rules

Each network at each level is run once to completion, with the primary arm's resume and rerun
rules.

## Manipulated variables

The clamp level f of perturbed genes: 0, 0.25, 0.5, 0.75, 1. Within each level, the coalition
of genes left unclamped: all 2^n.

## Measured variables

As in the primary arm, at each level.

## Measured variables - File upload

N/A — the record format is defined by `scripts/run_ode_arm.py`.

## Indices

Δ3+ and its class at each level; σ3+/σ2+, the order-3+ energy divided by the order-2+ energy of
the all-pairs spectrum, reported with its numerator and denominator; the null fraction at each
level; and the lowest level at which destroying networks outnumber creating ones.

## Indices - File upload

N/A — defined in `scripts/walsh_estimators.py`.

## Statistical models

N/A — G1–G3 and G5 are counts, and G4 compares two medians across networks; no test is applied.

## Statistical models - File upload

N/A — `scripts/analyze_ode_rerun.py` (`graded`), frozen with this file.

## Transformations

N/A — Δ3+ in percentage points and σ3+/σ2+ as computed.

## Inference criteria

| hypothesis | holds when |
|---|---|
| G1 | at f = 0, creating > destroying among scored networks |
| G2 | at f = 0.5, destroying > creating among scored networks |
| G3 | at f = 0.25, at most 2 scored networks lie in the null band |
| G4 | \|median σ3+/σ2+ at f = 0 − median at f = 0.75\| ≤ 0.05, over networks admissible at both levels |
| G5 | at each of f = 0.25, 0.5, 0.75 and 1, null networks ≤ a quarter of scored networks |

**G1–G5 are void at a level where more than 3 of the 28 networks cannot be scored.**

**G4 is void if fewer than 10 networks have an admissible ratio at both f = 0 and f = 0.75.**

**A network enters G4 at a level only if its ratio is admissible there:** its order-3+ and
order-2+ energies must each be at least 1e-4 of the total (`experiments/ode_v2/CONTEXT.md`,
Admissibility). Otherwise it is gated at that level, excluded from G4 there and nowhere else.
The order-3+ condition is the published arm's gate as its scoring code applied it
(`scripts/modal_graded_perturbation.py`). Amendment 1 of the v1 registration describes an
absolute-energy gate instead, which is inconsistent with that code. The absolute-energy version
is run as a sensitivity analysis.

If G1 and G2 hold, the abstract's sentence stands for the corrected engine. If G2 fails, the
abstract no longer says creation inverts by half-strength knockdown: it reports the lowest level
at which destruction outnumbers creation, or that it never does. If G3 fails, the sentence that
the gap persists at quarter strength is restated with the null-band count. If G5 fails, the
paper no longer says the gap persists throughout, and reports the null fraction at each level.
If G4 fails, the paper reports the factor, and if the factor is at least 1.5, that v1's H2 holds
under the corrected engine.

## Data inclusion and exclusion

All 28 networks at all five levels. A network is excluded at a level only if it cannot be scored
or is withheld for a replicate mismatch there, and from G4 also if it is gated. Every exclusion
is reported with its cause.

## Missing data

As in the primary arm: a network with a coalition lacking a value is not scored at that level.
If every failure is unclassified, a labeled sensitivity analysis scores it with provisional
outputs.

## Other planned analysis

1. **Secondary estimators.** G1–G5 recomputed with the split-half and plain estimators.
2. **Absolute-energy gate.** G4 recomputed over networks whose ratio is admissible and whose
   absolute order-3+ energy is at least 1e-4 at both levels, the gate as amendment 1 of the v1
   registration describes it.
3. **Ratio gate.** G4 recomputed with the ratio gate at 1e-5 and at 1e-3 of the total, with the
   networks each admits or excludes at each level.
4. **G5 on a common set.** G5 recomputed over the networks scored at all five levels, so that the
   denominator is the same at every level. The null count and the number scored are reported at
   every level, for G5 and for this set.
5. **The published table.** Creating, destroying and null counts and the median Δ3+ at each
   level, in the published table's format.
6. **Sign changes.** The number of networks whose class at each level differs from their class
   at f = 0, and for each network the lowest level at which its class changes.
7. **The v1 factor.** The median σ3+/σ2+ ratio between f = 0 and f = 0.75.

Anything else is exploratory and labeled so.

## Context and additional information

**Maximum claim under this registration.** Confirmatory: whether, under the corrected engine,
creation outnumbers destruction at complete knockout (G1) and destruction outnumbers creation at
half-strength knockdown (G2); whether the gap stays outside the null band at f = 0.25 in all but
at most two networks (G3) and in at least three quarters of the networks at every level (G5);
and whether the normalized higher-order fraction is flat between f = 0 and f = 0.75 (G4). The
level at which the balance inverts, the sign-change counts and the v1 factor are descriptive. A
mechanism for the inversion is exploratory.

This registration supersedes `prereg_graded_perturbation_v1.md` for the rerun. The published
values remain the record of what the legacy engine produced.

---

## Log

Append only. Never edit above the line.

The last column is what distinguishes an amendment from a deviation, so you do not have to
decide which word to use: `nothing run`, `no results seen`, `results not opened`, `results seen`.

```
2026-09-21  created                              nothing run
2026-09-22  frozen at 851b7a96154f                nothing run  ·09adf740
2026-09-22  frozen at 78071efa918c                nothing run  ·ec381430
```
