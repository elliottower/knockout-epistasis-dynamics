# Dynamical composition creates and destroys epistasis in gene regulatory networks

## Thesis (one sentence)

Feedback loops in Boolean gene regulatory networks generate emergent
higher-order epistasis absent from any single update rule, while linear
cascades destroy it — so measured phenotypic epistasis systematically
diverges from molecular interaction structure.

## Contribution

Biologists routinely interpret epistasis from double-knockout screens as
evidence of molecular interaction. This assumes a correspondence between
functional epistasis (non-additive knockout effects on phenotype) and
structural interactions (which genes appear together in regulatory rules).

We test this assumption directly using curated Boolean network models
where both objects are known exactly:

- **Local rule Fourier structure**: the Walsh-Hadamard transform of each
  gene's truth table gives exact, representation-independent interaction
  coefficients at every order.
- **Attractor-level epistasis**: exhaustive knockout sweeps over all 2^n
  gene combinations yield the global value function, whose Walsh
  decomposition gives attractor-level interaction coefficients.

Comparing the two reveals whether dynamics preserve, destroy, or create
interaction structure.

## Method

### Local rule extraction

Each gene's Boolean update rule defines a function f_k: {0,1}^{r_k} -> {0,1}
where r_k is the number of regulators. Since r_k <= 6 for all models,
the truth table has at most 64 entries. We take its Walsh-Hadamard transform
and normalize to get exact Fourier coefficients. A nonzero coefficient at
order >= 2 indicates a genuine multi-way interaction among regulators,
independent of how the rule is written (DNF, CNF, or any other form).

This replaces DNF clause-splitting, which is representation-dependent:
`a & b | a & !b` simplifies to `a` (no interaction with b), but clause
co-occurrence reports a spurious pair.

### Coalition sweep

For each of the 2^n subsets S of genes (coalitions), genes not in S are
clamped to 0 (knockout). The network runs synchronous Boolean dynamics
from N random initial states until attractor convergence. Limit cycles
are detected and the output is averaged over one full period. The value
function v(S) = mean output over initial states.

Asynchronous update (random-permutation sequential) is run for models
published under that scheme.

### Composition scoring

For all C(n,2) gene pairs:
- **Local magnitude**: sum of |local Walsh pairwise coefficient| across
  all genes where the pair co-occurs (continuous, not binary)
- **Global magnitude**: |global Walsh pairwise coefficient| from the
  attractor value function

Primary metric: Spearman rank correlation between local and global
magnitudes. No threshold, no class-imbalance problem.

Secondary metrics:
- **Creation rate**: fraction of top-k global pairs with zero local source
- **Destruction rate**: fraction of top-k local pairs with near-zero
  global magnitude
- **Order spectrum comparison**: energy fraction at each interaction order,
  local vs global

## Preliminary results (n_init=8 code test — NOT the registered experiment)

| Model | n | Spearman rho | p-value | Global o3+ | Local o3+ | Delta | Direction |
|-------|---|-------------|---------|------------|-----------|-------|-----------|
| faure_cellcycle (G1) | 10 | 0.33 | 0.028 | 0.9% | 3.8% | -2.9pp | Destroys |
| tournier_apoptosis (G2) | 12 | 0.28 | 0.025 | 17.5% | 4.2% | +13.3pp | **Creates** |
| davidich_yeast (G4) | 10 | 0.08 | 0.59 | 3.9% | 2.6% | +1.3pp | Creates (weak) |

### What the numbers mean

**G2 is the headline.** Three positive feedback loops (C3a-C8a executioner
loop, C3a-IAP double negative, C3a-CARP-C8a reinforcement) create a
bistable landscape with a sharp apoptosis threshold. This threshold
generates 17.5% order-3+ energy at the attractor level from rules that
contain only 4.2% locally. The extra 13.3 percentage points are emergent
— they exist in the phenotype but not in any single molecular rule.

**G1 is the control.** One positive and one negative feedback loop. The
negative feedback (CycB-Cdc20-cdh1-CycB) damps higher-order structure.
Composition destroys 2.9pp of order-3+ energy.

**G4 is the surprise.** Despite similar biology to G1 (CDK-driven cell
cycle), Spearman rho ~ 0. Local rule magnitudes have no predictive power
for global attractor magnitudes. 75% of top global pairs have no local
source. Two positive feedback loops reshuffle the interaction landscape
completely.

### Structural predictor: feedback loop count

| Model | Positive feedback loops | Negative feedback loops | Composition effect |
|-------|------------------------|------------------------|-------------------|
| G1 faure | 1 | 1 | Destroys o3+ |
| G2 tournier | 3 | 1 | Creates o3+ (13pp) |
| G4 davidich | 2 | 1 | Reshuffles (rho ~ 0) |

Positive feedback loop count predicts whether dynamics create higher-order
structure. This is testable: adding more networks with varying feedback
topology should show a monotonic relationship.

## What's needed for a full paper

### Data
- [ ] Full experiment at n_init=512 for G1, G2, G4 (with bootstrap CIs)
- [ ] 5-7 more networks selected on logic depth (survey in progress)
- [ ] Input nodes excluded from player set (CycD, TNF, Start fixed per-condition)
- [ ] Both clamp-to-0 and clamp-to-1 ablation
- [ ] Sync + async for all models
- [ ] Variance floor gate with absolute order-3+ energy threshold

### Analysis
- [ ] Spearman with bootstrap 95% CIs over initial states
- [ ] Regression: feedback loop count vs composition delta
- [ ] Creation/destruction rates across all models
- [ ] Per-order spectrum comparison (not just 3+)
- [ ] Triple-level composition (Spearman on order-3 coefficients)

### Framing
- [ ] Relate to experimental epistasis literature (double-knockout screens)
- [ ] Discuss implications for interpreting CRISPR combinatorial screens
- [ ] Connect to mech-interp: if dynamics reshuffle interactions in 10-node
      Boolean networks, what happens in 100M-parameter transformers?

### Venue
TBD. This could be:
- A standalone genetics methods paper (PLoS Comp Bio, Bioinformatics)
- A section of the EpistasisBench paper (but it's a different thesis)
- A short communication / letter (the result is crisp enough)

## Key sentence for the abstract

Across N curated Boolean gene regulatory networks, dynamical composition
creates up to X percentage points of higher-order epistasis absent from
any molecular rule, with positive feedback loop count as the structural
predictor — demonstrating that phenotypic epistasis from knockout screens
does not, in general, correspond to molecular interaction structure.
