# Pre-registration: a GNK positive control, and whether dynamics creates off-graph structure

**Frozen before the GNK control is run. Supersedes the static arm of
`prereg_static_vs_dynamic_v1.md`, whose design does not work.**

## Why the previous static arm was invalid

`prereg_static_vs_dynamic_v1.md` built its static phenotype as
$\phi(\sum_{(i,j)\in\mathcal{E}} J_{ij} x_i x_j)$ -- a saturating nonlinearity
over a pairwise energy.

A pairwise energy has **exactly zero** Walsh coefficients above order 2. So
every order-3+ term in that phenotype is produced by the nonlinearity. It is
100% scale artifact in the sense of Sailer and Harms (*PLoS Comput Biol*
13(5):e1005541, 2017), who state that "the most important confounding effect is
the scale of the map" and who **remove** the nonlinear scale by power transform
*before* decomposing. The previous design applied a nonlinearity and decomposed
without removing it, which is backwards, and linearising it would leave a
phenotype with no higher-order structure at all.

The static arm therefore tested nothing. Its H1/H3/H4 verdicts are withdrawn.
**H2 -- that dynamic systems are not predictable from wiring -- does not depend
on the static construction and stands.**

## The replacement, taken from the source rather than invented

Brookes, Aghazadeh and Listgarten (*PNAS* 119(1):e2109649118, 2022) define the
GNK model, read here from the full text:

> "Assign a 'subsequence fitness,' $f_j(s[j])$, to every possible subsequence,
> $s[j]$, by drawing a value from the normal distribution with mean equal to
> zero and variance equal to $1/L$... For every $s$, the subsequence fitness
> values are summed to produce the total fitness values
> $f(s) = \sum_{j=1}^{L} f_j(s[j])$."

with neighborhoods defined so that "the neighborhood of a position $j$ contains
all positions that are in structural contact with it."

**Construction here.** For each network, $V_j = \{j\} \cup \{$nodes sharing a
regulatory edge with $j\}$. Draw $f_j \sim \mathcal{N}(0, 1/L)$ independently for
each of the $2^{|V_j|}$ states of $V_j$. The phenotype of a knockout coalition is
$\sum_j f_j(x_{V_j})$, with $x_i = 0$ for knocked-out nodes and $1$ otherwise.

This produces **genuine** higher-order structure localised on the graph, rather
than manufactured by a scale transform.

## What this buys that the previous design did not

The star construction is a **theorem** about GNK, not an empirical claim about
it: given Structural neighborhoods, an $r$th-order coefficient is nonzero when
$r-1$ positions contact a central position. So GNK is a **positive control on
the predictor itself**, which this project has never had.

## Predictions

**H0 (positive control, must pass before anything else is read).** On GNK
phenotypes built from each network's regulatory graph, the star predictor
recovers order-3+ coefficients at median Spearman $\rho \geq 0.5$ across the 27
networks. **If H0 fails the predictor implementation is wrong**, every
conclusion drawn with it -- including the withdrawn static arm and any statement
about dynamic systems -- is unsupported, and nothing further is interpreted.

**H1 (primary).** Dynamic (attractor) phenotypes score materially lower than
GNK on the same predictor and the same graph: median $\rho$ at least $0.25$
below the GNK median.

**H2 (primary).** The dynamic median stays below $0.3$ in absolute terms,
reproducing the already-measured result under the corrected predictor.

**H3, the falsifier.** If dynamic and GNK score within $0.1$ of each other, the
attractor phenotype localises on the wiring just as a GNK model does, dynamics
generates no off-graph structure, and the mechanism claim is dead.

## Decision rule, fixed in advance

| outcome | consequence |
|---|---|
| H0 passes, H1 and H2 hold | Dynamics creates higher-order structure that is **not** localised on the wiring, where a graph-local model on the same graph is. This is the mechanism the paper currently asserts from a cross-literature comparison, measured inside one system. |
| **H0 fails** | The predictor is broken. Everything measured with it is withdrawn, including the existing statement that dynamic systems are unpredictable from wiring. |
| **H3 fires** | Dynamics is not the mechanism. The Discussion claim is deleted and the gap is reported without an explanation. |
| H1 fails, H2 holds | Neither localises well. Consistent with a regulatory edge list being a poorer object than a structural contact map, and no mechanism is claimed. |

## Controls

Ten independent GNK draws per network, since $f_j$ is random; report the median
across draws and the spread. Neighborhood sizes are recorded per network,
because a network whose $|V_j|$ approaches $n$ has a GNK model that is nearly
unconstrained and cannot support the star property.

## Not varied

The 27 networks, their regulatory graphs, the committed coalition tables, the
Walsh basis, and the star predictor as verified against the source.

## Script

`scripts/gnk_control.py`. To be written; nothing has been run.

---

# Post-run analysis note

Recorded after the run, prompted by an external review. No prediction is
changed; this documents a metric artifact and identifies which statistic should
carry the claim.

## The absolute correlation is confounded by class balance

GNK's Spearman $\rho$ correlates with the **fraction of coefficients the
predictor marks active** at $\rho = 0.975$ ($p < 0.001$), and with mean
neighborhood size at $0.755$. Networks scoring below the H0 threshold have a
mean active fraction of $0.04$; those at or above have $0.34$.

Correlating a magnitude against a binary indicator is limited by that
indicator's balance. So the absolute medians -- GNK $+0.583$, dynamic $+0.256$ --
partly measure how sparse each network's predictions are, not how well the
predictor works.

**H0's threshold was therefore the wrong gate.** It is balance-confounded.

## The gate that does validate the predictor

The direct theorem check is balance-independent: for GNK, **100.0% of order-3+
spectral energy sits on predicted-active coefficients, on every network tested**,
against 78.9% under the earlier incorrect predictor. That is what the theorem
requires and it is the check the positive control should have been built on.

## The statistic that should carry the claim

The **paired** within-network comparison, which holds graph, predictor and
active fraction fixed and varies only the phenotype:

- dynamic scores below GNK in **24 of 25 networks**
- dynamic reaches a median of **0.58** of its own network's GNK value
  (IQR $0.34$--$0.76$)

Report the paired result. The absolute medians may be quoted alongside with the
balance confound stated, and never on their own.

## Still open

Whether a regulatory edge list is a fair analogue of a structural contact map
remains unresolved. The GNK control establishes that the predictor is correct
**when the phenotype is built from the graph it is given**; it does not
establish that this graph is the right object for the attractor phenotype. Any
claim from this arm must be scoped to "structure created by dynamics is not
localised on the regulatory graph", not to contact maps in general.
