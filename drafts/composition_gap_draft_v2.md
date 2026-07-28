# Order-by-order quantification of the gap between molecular interactions and phenotypic epistasis in gene regulatory networks

## Thesis (one sentence)

The divergence between molecular interaction structure and phenotypic
epistasis — established qualitatively by Gjuvsland et al. (2007) and
Phillips (2008) — is quantifiable order-by-order via paired Walsh
decompositions of local update rules and attractor-level knockout
landscapes. The composition gap is signed (dynamics create or destroy
structure at each order) and varies across network architectures in
ways consistent with feedback topology, though N=3 networks do not
permit a formal regression.

## Prior art and positioning

**What is known.** Statistical epistasis is a generic feature of gene
regulatory networks (Gjuvsland et al. 2007, Genetics). Nonlinear
dynamics generate non-additive phenotypic effects from architectures
with and without epistatic molecular interactions. Phillips (2008,
Nature Reviews Genetics) formalized the statistical-versus-functional
epistasis distinction. Baier et al. (2023, Science Advances) measured
pairwise and triplet epistasis in a synthetic 3-node GRN, finding
environment-dependent sign switching. The global epistasis literature
(Otwinowski et al. 2018; Reddy & Desai 2021) further establishes that
nonlinear genotype-phenotype maps reshape interaction spectra.

**What is not known.** Prior work characterizes the *existence* of the
gap between molecular and phenotypic epistasis. No study has:
1. Computed the complete Walsh-Fourier spectrum to order n for both
   the local update rules and the attractor-level knockout landscape
2. Compared the two spectra order-by-order on the same scale
3. Quantified whether composition creates or destroys structure at
   each interaction order, with magnitudes
4. Examined whether feedback topology covaries with the direction
   (requires N >> 3 for formal testing)

**Our contribution.** For each gene's Boolean update rule (a function
of at most 6 regulators), we take the Walsh-Hadamard transform of the
truth table — at most 64 entries — to obtain exact, representation-
independent interaction coefficients at every order. For the full
network, we run exhaustive 2^n knockout sweeps and take the Walsh
transform of the attractor-level value function. Comparing the two
spectra gives a signed, quantified composition gap at each interaction
order.

## Lead result: G4 (davidich_yeast)

Spearman correlation between local rule pairwise magnitudes and
attractor-level pairwise magnitudes: rho = 0.15 (p = 0.31, 95% CI
crosses zero). Knowing which gene pairs interact in the molecular
rules tells you nothing about which pairs show knockout epistasis.
75% of the strongest global pairwise interactions have no local
source — created entirely by dynamical composition.

This is the most practically relevant finding for anyone designing
double-knockout screens: the molecular wiring diagram is not a useful
prior for predicting which gene pairs will show epistasis.

## Second result: G2 (tournier_apoptosis)

Dynamics create higher-order epistasis. Local rules contain 4.2%
order-3+ energy; the attractor landscape contains 12.6% (+8.5pp,
absolute 5.3e-3). The pairwise Spearman correlation (rho=0.25, CI
[0.03, 0.45]) provides weak evidence for a positive association —
the CI barely excludes zero. Three positive feedback loops create a
bistable apoptosis threshold that generates multi-way knockout
interactions absent from any single molecular rule. 62% of top
global pairs have no local source.

## Third result: G1 (faure_cellcycle)

Dynamics DESTROY higher-order structure. The negative feedback loop
(CycB-Cdc20-cdh1-CycB) damps order-3+ energy from 3.8% to 0.8%
(-3.0pp). Moderate Spearman correlation (rho = 0.41, p < 0.01) —
local pairwise structure partially survives composition. This is the
only network where the molecular wiring diagram is a statistically
significant predictor of knockout epistasis.

## Method

### Local rule Fourier extraction

Each gene's Boolean update rule defines f_k: {0,1}^{r_k} -> {0,1}.
Since r_k <= 6 for all models, the truth table has at most 64 entries.
Its Walsh-Hadamard transform gives exact Fourier coefficients. A
nonzero coefficient at order >= 2 indicates a genuine multi-way
interaction among regulators, independent of how the rule is written.

This is representation-independent: `a & b | a & !b` simplifies to
`a`, and the WHT correctly reports zero pairwise interaction. DNF
clause-splitting would report a spurious pair.

### Coalition sweep

For each of 2^n coalitions S, genes not in S are clamped to 0
(knockout). Synchronous Boolean dynamics from N=512 random initial
states until attractor convergence. Limit cycles detected and output
averaged over one full period.

### Composition scoring

For all C(n,2) gene pairs:
- **Local magnitude**: sum of |local Walsh pairwise coefficient| across
  all genes where the pair co-occurs
- **Global magnitude**: |global Walsh pairwise coefficient| from the
  attractor value function

Primary metric: **Spearman rank correlation** between local and global
magnitudes. Bootstrap 95% CIs over initial states.

Secondary: creation rate (fraction of top-k global pairs with zero
local source), destruction rate, per-order energy spectrum comparison.

**Aggregation choice.** Local pairwise coefficients are summed (|w|)
across all genes in which the pair appears. Max-aggregation is reported
as a robustness check. The two give [TBD: similar / different] Spearman
values.

## Results (n_init=512, sync update, clamp-to-0)

| Model | n | Spearman rho [95% CI] | p | Global o3+ | Local o3+ | Delta | Direction |
|-------|---|----------------------|---|------------|-----------|-------|-----------|
| faure_cellcycle (G1) | 10 | 0.409 [0.129, 0.637] | 5.3e-3 | 0.8% | 3.8% | -3.0pp | Destroys |
| tournier_apoptosis (G2) | 12 | 0.247 [0.026, 0.450] | 4.6e-2 | 12.6% | 4.2% | +8.5pp | Creates |
| davidich_yeast (G4) | 10 | 0.154 [-0.148, 0.423] | 3.1e-1 | 2.2% | 2.6% | -0.4pp | Neutral |

### Statistical notes

**Correlation strength.** Each network is an independent biological
system; the three Spearman tests answer separate questions ("does
local predict global in this network?"), so no family-wise correction
is applied. Instead, bootstrap 95% CIs on rho convey the precision
of each estimate directly. G1's CI [0.13, 0.64] excludes zero — a
real positive association. G2's CI [0.03, 0.45] barely excludes
zero — weak evidence for a positive association. G4's CI [-0.15,
0.42] spans zero — no detectable association. All three CIs overlap
substantially, so the apparent ordering G1 > G2 > G4 cannot be
distinguished from a shared rho near 0.25.

**G1 energy floor.** G1's global order-3+ energy (0.8% fractional)
corresponds to an absolute squared-coefficient sum of 1.7e-3. This
clears the 1e-4 noise floor by 17x, so the "destruction" label
reflects real signal compression. For comparison, G2 and G4 have
absolute order-3+ energies of 5.3e-3 and 4.9e-3.

### Per-network case studies

Each network is treated as an independent case study. With N=3
networks, no cross-network trend or regression is claimed.

**G1 (faure_cellcycle, n=10, 45 pairs)**: The only network where
local pairwise structure significantly predicts attractor-level
epistasis after correction (rho=0.41, p=0.005). Dynamics compress
higher-order structure: order-3+ energy drops from 3.8% to 0.8%
(-3.0pp, absolute 1.7e-3). The negative feedback loop
(CycB-Cdc20-cdh1-CycB) damps multi-way interactions. Triple-level
correlation is moderate (rho=0.39, CI [0.23, 0.51]). 19.6% of
initial conditions reach limit cycles.

**G2 (tournier_apoptosis, n=12, 66 pairs)**: Pairwise correlation
is weak (rho=0.25, CI [0.03, 0.45] — barely excludes zero).
Dynamics create higher-order epistasis: order-3+ energy
increases from 4.2% to 12.6% (+8.5pp, absolute 5.3e-3). 62% of
top global pairs have no local source. Triple-level correlation is
near zero (rho=0.05, CI spans zero), indicating that the created
higher-order structure has no relationship to local triple
interactions. The bistable apoptosis threshold — maintained by
positive feedback — generates multi-way knockout interactions absent
from any single molecular rule. 7.4% cycling.

**G4 (davidich_yeast, n=10, 45 pairs)**: No significant pairwise
correlation (rho=0.15, CI [-0.15, 0.42], p=0.31). The molecular
wiring diagram provides no useful prior for predicting which gene
pairs show knockout epistasis. 75% of top global pairs have no local
source. Order-3+ energy is similar locally and globally (2.6% vs
2.2%, absolute 4.9e-3), so dynamics reshuffle pairwise structure
without net creation or destruction at higher orders. Triple-level
correlation is weak but positive (rho=0.21, CI [0.09, 0.32]).
17.7% cycling.

## What's still needed

### Data
- [x] n_init=512 results
- [ ] 5-7 more networks selected on logic depth (survey in progress)
- [ ] Input nodes excluded from player set, run per-condition
- [ ] Max-aggregation robustness check
- [ ] Both clamp-to-0 and clamp-to-1

### Analysis
- [ ] With N >= 8 networks: regress composition delta on feedback count
- [x] With N = 3: per-network case studies only (no regression)
- [x] Triple-level Spearman (order-3 coefficients) — G1: 0.39, G2: 0.05, G4: 0.21
- [x] Absolute-energy floor check on G1's order-3+ — 1.7e-3, clears 1e-4 by 17x
- [x] CI overlap analysis — all three rho CIs overlap; no correction needed (independent systems)

### Related work to cite
- Gjuvsland et al. 2007, Genetics (statistical epistasis in GRNs)
- Phillips 2008, Nat Rev Genet (statistical vs functional epistasis)
- Baier et al. 2023, Science Advances (triplet epistasis in synthetic GRN)
- Otwinowski et al. 2018 (global epistasis)
- Reddy & Desai 2021 (global epistasis on fitness landscapes)

### Venue
Genetics-first. PLoS Comp Bio or Genetics (same journal as Gjuvsland).
Short communication if N stays at 3.
