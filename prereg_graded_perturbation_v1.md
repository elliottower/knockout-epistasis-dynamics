# Pre-registration: is the composition gap a property of composition, or of complete removal?

**Frozen before any graded sweep is run. No graded result exists.**

## The question

Every result in the composition gap paper uses a **complete** knockout -- a binary
clamp to zero. That is the most violent perturbation available, and this project's
own measurement is that **zero-ablation inflates pairwise epistasis 10--50x
relative to mean ablation**.

So "dynamics create higher-order epistasis" may really be "removing a node
entirely creates higher-order epistasis". Those are different claims and only one
is currently supported.

## Design

The perturbation set $S$ stays binary -- a gene is perturbed or it is not -- so the
**Walsh--Hadamard decomposition still applies unchanged**. What varies is the
*level* $f$ to which perturbed genes are set, as a fraction of wild-type activity,
applied uniformly within a sweep.

Hill-function ODE conversions, already built for all 28 networks. Sweep
$f \in \{0, 0.25, 0.5, 0.75\}$ for loss of function and
$f \in \{1.25, 1.5, 2.0\}$ for gain of function. $f = 0$ reproduces the existing
binary result and is the consistency check; $f = 1$ is no perturbation and is
excluded, since the landscape is then constant by construction.

**The headline quantity is the normalised higher-order fraction**
$\sigma_{3+}/\sigma_{2+}$, not the raw gap. Raw magnitude must fall as
$f \to 1$ for trivial reasons; the fraction need not.

## Predictions

**H1 (primary).** The gap **sign** is preserved at every $f \in (0,1)$. Creation
remains creation under partial knockdown.

**H2 (primary).** The normalised higher-order fraction is **not constant in $f$**:
it increases as $f \to 0$, by at least a factor of 1.5 between $f = 0.75$ and
$f = 0$. This follows from the zero-versus-mean ablation measurement and, if it
fails, that measurement does not generalise beyond pairwise.

**H3 (secondary).** The dose-response is **monotone** in $f$. A threshold or
non-monotone shape would indicate a regime boundary that the current account does
not describe, and would be reported as such rather than smoothed.

**H4 (secondary).** Sweeping $f$ through 1 into over-expression, the
creation/destruction ratio **flips at some $f^\ast$**, and we report where. If
$f^\ast = 1$ exactly, the reversal is a statement about the reference point rather
than about biology.

**H5, the falsifier.** If the gap vanishes -- sign flip, or $|\Delta_{3+}| < 0.5$
pp -- for any $f \geq 0.25$, the published result is **specific to complete
knockout** and the paper must say so in the abstract.

## Decision rule, fixed in advance

| outcome | consequence |
|---|---|
| H1 and H2 hold | The gap survives graded perturbation and intensifies with perturbation strength. The claim generalises beyond complete removal, and the dose-response becomes a figure. |
| H1 holds, H2 fails | The higher-order fraction is scale-free in perturbation strength -- more surprising than the prediction, and it means the 10--50x ablation effect is confined to order 2. Report as the finding. |
| **H5 fires** | The result is about complete removal, not about composition. Title and abstract are restated. This is the outcome that would most change the paper. |
| H4 gives $f^\ast = 1$ | The loss/gain reversal is a reference-point artifact. Its prominence in the paper is reduced accordingly. |

## Not varied

The 28 networks, their wiring, the ODE conversion already validated (27/28
preserve gap sign), the attractor readout, the knockout coalition enumeration,
and the Walsh basis. Only the perturbation level changes.

## Script

To be written. Hash recorded in `docs/FROZEN_SHAS.md` before running.

---

# Amendment 1 — before anything ran

**No graded sweep has executed.** Revision follows a code audit of the existing
sweep machinery, not results.

## A1.1 The Boolean path cannot express a graded perturbation, and would fail silently

`simulate_sync_output` applies `current[:, clamp_mask] = clamp_value`, and the
update indexes a per-gene truth table via
`idx += states[:, ri].astype(np.int64) << j`.

A fractional clamp value therefore **casts to an integer**: `clamp_value = 0.5`
becomes `0`. The sweep would run without error and return results **identical to
full knockout**, so a graded Boolean sweep would silently reproduce the
published `blind` result at every level and read as "the gap is invariant to
perturbation strength".

**The Boolean path is excluded from this arm.** Graded perturbation is defined
only on the Hill-function ODE conversion, where node states are continuous.

## A1.2 $f$ is the clamp value, not a fraction of wild-type

The original text swept $f \in \{1.25, 1.5, 2.0\}$ for gain of function. That is
unphysical here: the Hill sigmoid
$\sigma(z) = z^{n_H}/(z^{n_H} + K^{n_H})$ bounds node states to $[0,1]$, so a
clamp above 1 lies outside the model's range.

Redefined: **$f$ is the value perturbed genes are clamped to**, sweeping
$f \in \{0, 0.25, 0.5, 0.75, 1\}$. The endpoints are the two regimes the paper
already reports -- $f = 0$ is loss of function (`blind`), $f = 1$ is gain of
function (`clamp1`) -- and the intermediate levels are the new measurement.

This makes H4 well-posed rather than unphysical: the loss/gain reversal is a
sweep between two published endpoints, and $f^\ast$ is where the
creation/destruction ratio crosses.

## A1.3 A reproduction gate, which the original omitted

Because $f = 0$ and $f = 1$ correspond to sweeps already computed for 27
networks, they are a free check on the implementation.

**Gate, evaluated before any intermediate $f$ is interpreted:** the graded
harness at $f = 0$ must reproduce each network's committed
`*_composition_blind.json` order-3+ energies, and at $f = 1$ its
`*_composition_clamp1.json`, to within Monte-Carlo tolerance over the 32 initial
conditions. **If the endpoints do not reproduce, the harness is wrong and no
intermediate result is read.**

This is stated now because a graded sweep that disagrees with the published
endpoints could otherwise be reported as a finding about perturbation strength.

## A1.4 Variance floor on the headline quantity

The headline is $\sigma_{3+}/\sigma_{2+}$. As the perturbation weakens the value
function flattens, so both terms approach zero and the ratio becomes a quotient
of small numbers.

The existing `spectrum_gate(v_mean, n, min_abs_energy_3plus=1e-4)` in
`composition_scorer.py` is applied unchanged at every $f$. Networks failing the
gate at a given $f$ are **reported as gated, not as zero**, and are excluded
from the pooled statistic at that $f$ with the exclusion stated on the figure.

## A1.5 Compute

ODE coalition tables were never saved locally -- `results/ode_full/` holds
summaries for 13 of 28 networks and the NPZs written by
`scripts/modal_ode_sweep_save_coalitions.py` live on the Modal volume. Each $f$
therefore requires re-integration: $2^n$ coalitions $\times$ 32 initial
conditions $\times$ RK45.

Runs on Modal. Three intermediate levels ($f = 0.25, 0.5, 0.75$) plus the two
reproduction endpoints.

## Predictions unchanged

H1, H2, H3 and H5 stand as written. H4 is now well-posed rather than
unphysical, and its prediction -- that the creation/destruction ratio flips at
some $f^\ast$, with the location reported -- is unchanged.
