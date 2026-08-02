# Pre-registration: does the composition gap require dynamics?

**Frozen before any static-system spectrum is computed.**

## Why this is not already answered

The paper's proposed mechanism is that composition **through time** generates
higher-order structure absent from the wiring. The evidence offered is a
comparison with proteins: Brookes, Aghazadeh and Listgarten (*PNAS* 2022) find
higher-order epistasis in protein landscapes is **localised on the pairwise
contact graph**, so there the wiring does predict it.

That is a cross-literature comparison, and proteins differ from gene regulatory
networks in many ways besides having static rather than dynamic phenotypes. The
mechanism is currently **asserted, not tested**.

## Design

For each of the 28 networks, hold the wiring fixed and define **two phenotypes on
the same graph**:

- **Dynamic:** the existing attractor-level readout.
- **Static:** a pairwise energy over the same edges,
  $E(x) = \sum_{(i,j) \in \mathcal{E}} J_{ij} x_i x_j$, with $J_{ij}$ signed by
  the activating or inhibiting character of the edge, passed through a saturating
  nonlinearity $\phi$ so that the phenotype is $\phi(E)$.

The nonlinearity is **required**, not incidental. Without it the static phenotype
is quadratic in the perturbation indicators and has exactly zero higher-order
terms by construction, which would make the test vacuous. With it, higher-order
terms appear -- and the question becomes whether they are **predictable from the
wiring**, exactly as Brookes et al. found for proteins.

Same nodes, same edges, same knockout coalitions, same Walsh basis.

## The measured quantity

Not the presence of higher-order energy -- both systems will have some -- but its
**predictability from local structure**: the correlation between order-3+ Walsh
coefficients and those predicted from the wiring alone.

## Predictions

**H1 (primary).** In the static systems, order-3+ coefficients are **predictable
from the wiring**: correlation $\geq 0.5$ between predicted and actual, in at
least 20 of 28 networks. This is the protein-like regime.

**H2 (primary).** In the dynamic systems, they are **not**: correlation $< 0.3$,
consistent with the paper's existing finding that topological predictors score at
chance.

**H3 (primary).** The paired difference across the 28 networks is significant, and
each network serves as its own control since the graph is identical.

**H4, the falsifier.** If static systems also show unpredictable higher-order
structure, **the dynamics mechanism is wrong**, the Brookes framing in the
Discussion is withdrawn, and the explanation for the composition gap becomes open
again.

## Decision rule, fixed in advance

| outcome | consequence |
|---|---|
| H1, H2, H3 hold | The mechanism is established **within one system** rather than inferred across two literatures. This replaces the Brookes citation with a result and is the strongest available version of the paper's central explanation. |
| **H4 fires** | Dynamics is not the mechanism. A major Discussion claim is deleted, and the paper reverts to reporting the gap without explaining it. |
| H1 fails but H2 holds | Neither system's higher-order structure is predictable from wiring, so the contrast with proteins is not about dynamics but about something else -- possibly that Brookes' contact graph carries more information than a regulatory edge list. Report honestly; do not claim the mechanism. |

## What this cannot show

The static phenotype is a construction, not a measured protein landscape. It
isolates the dynamic/static distinction while holding topology fixed, which is
what a cross-literature comparison cannot do -- but it is not evidence about real
proteins, and the paper must not present it as such.

The choice of $\phi$ matters. Report results for at least two saturating
nonlinearities and show the conclusion is invariant to it.

## Not varied

The 28 wiring graphs, node sets, knockout coalitions, and Walsh basis.

## Script

To be written. Hash recorded in `docs/FROZEN_SHAS.md` before running.
