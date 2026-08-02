# Pre-Registered Predictions: Composition Gap in Epistatic Landscapes

**Status:** These predictions were made WITHOUT access to any computed results.
Only the simulation code (`grn_coalition_sweep.py`, `composition_scorer.py`,
`data_utils.py`) and theoretical framework (`theorem_scope.md`) were read.
No files in `results/`, no JSON output files, and no computed summary
statistics were consulted.

**Date:** 2026-07-28

**Basis for predictions:** The Walsh-Hadamard decomposition framework, the
algebraic/dynamical decomposition of the composition gap (Propositions 1--3
in `theorem_scope.md`), and general knowledge of Boolean network dynamics,
protein biophysics, and yeast metabolic genetics.

---

## Prediction 1: Cycling fraction vs. composition gap direction

**Question:** What is the Spearman correlation between `cycling_fraction`
(fraction of coalition-by-initial-state pairs entering limit cycles) and
`delta_o3plus` (global minus local order-3+ energy fraction) across the 27
Boolean gene regulatory networks?

### (a) Predicted sign: POSITIVE

### (b) Confidence: Medium

### (c) Reasoning

The composition gap decomposes as Delta_3+ = Delta_alg + Delta_dyn
(Proposition 1/3). DAG networks have zero cycling and Delta_dyn = 0, so their
composition gap is purely algebraic. Feedback networks have cycling > 0 and
can contribute a nonzero Delta_dyn through coalition-dependent attractor basin
geometry: different knockout sets route initial states into different basins
whose output statistics encode multi-way gene dependencies that exceed what
any single local rule contains. Because cycling requires feedback loops, and
feedback both lengthens compositional paths (increasing Delta_alg via
Proposition 2's chain-composition mechanism) and enables the dynamical
component Delta_dyn, networks with higher cycling fractions should have larger
composition gaps on average.

A competing effect exists: cycle averaging acts as a low-pass filter on the
output (a cycle alternating 0/1 yields 0.5, destroying coalition-dependent
structure for that trajectory). However, the dominant effect is attractor
multiplicity: feedback creates multiple basins with distinct outputs, and the
coalition-dependent partitioning of initial states across these basins
introduces sharp boundaries in v(S) that encode higher-order interactions. The
smoothing within basins is secondary to the structure created between basins.

---

## Prediction 2: ODE pilot (Hill function continuous dynamics)

**Question:** If the 4 smallest Boolean networks (n = 7--10) are converted to
continuous ODE dynamics using Hill functions, how does the composition gap
change?

### (a) Will the SIGN of Delta_3+ be preserved in at least 3/4 networks? YES

**Reasoning:** Hill function ODEs are monotone approximations of Boolean
dynamics. For sufficiently steep Hill coefficients, the Thomas--Snoussi
correspondence guarantees that the qualitative attractor structure (number and
stability of fixed points, presence of oscillations) is preserved. The sign
of Delta_3+ reflects whether the network's feedback topology and path
structure create or destroy higher-order epistasis relative to local rules.
This is a qualitative property of the network architecture, not a quantitative
artifact of the Boolean discretization. Converting to ODE dynamics changes
magnitudes but should preserve the direction of the composition gap in most
networks.

The 1/4 exception allowance accounts for networks near Delta_3+ = 0 where
small quantitative shifts under Hill dynamics could flip the sign, and for
networks where Boolean-specific attractor artifacts (e.g., synchronous update
creating spurious 2-cycles that vanish under ODE dynamics) qualitatively
change the attractor landscape.

### (b) Will the MAGNITUDE of Delta_3+ increase or decrease? DECREASE

**Reasoning:** Hill functions replace sharp Boolean switches with graded
sigmoidal responses. The coalition sweep still uses binary inputs (gene
present or knocked out), so the value function remains a pseudo-Boolean
function on {0,1}^n. However, the internal dynamics are now continuous.

Three mechanisms reduce |Delta_3+| under ODE dynamics:

1. **Smoother basin boundaries.** Boolean dynamics create sharp, digital
   transitions between attractor basins as a function of the coalition. ODE
   dynamics produce smoother transitions (gradual shifts in basin size and
   attractor position), reducing the higher-order Walsh energy that encodes
   these boundaries.

2. **Loss of Boolean amplification.** Boolean update rules amplify small
   differences into 0/1 outputs at each step, creating sharp nonlinearities
   that compound along paths. Hill functions with finite cooperativity
   (n = 2--4) partially preserve intermediate values, weakening the
   compositional amplification that generates higher-order terms along
   regulatory paths.

3. **Attractor value continuity.** ODE attractors produce continuous output
   values rather than {0,1}, reducing the effective dynamic range over which
   coalition-dependent variations can create higher-order spectral structure.

Both local and global higher-order energy decrease under Hill dynamics, but
the global decrease is larger because the compositional amplification along
multi-step paths (the dominant source of Delta_alg) is attenuated at each
step.

**Confidence:** Medium. The sign-preservation prediction is stronger than the
magnitude prediction. The magnitude could increase if Hill function dynamics
create richer attractor landscapes (e.g., stable limit cycles that were
absent under synchronous Boolean update), but I consider this unlikely for
the smallest networks.

---

## Prediction 3: Weinreich TEM-1 beta-lactamase fitness landscape

**Question:** The Weinreich et al. (2006) dataset measures antibiotic
resistance (MIC) for all 32 genotypes formed by 5 binary mutations in TEM-1
beta-lactamase. Structural contacts from the crystal structure define which
mutation pairs interact physically. What is the sign of the composition gap
(local structural interactions vs. global fitness landscape)?

### (a) Predicted sign: POSITIVE (creation of higher-order epistasis)

### (b) Reasoning

Structural contacts from the crystal structure capture direct pairwise
physical interactions between mutation sites: steric clashes, hydrogen bonds,
hydrophobic packing between neighboring side chains. These interactions are
inherently pairwise in the Walsh sense (one residue contacts one other
residue), so the local Walsh structure from structural contacts has order at
most 2 for most site pairs.

Antibiotic resistance, however, depends on collective enzyme properties that
couple all 5 sites simultaneously:

1. **Catalytic efficiency** requires precise geometric coordination of
   active-site residues, substrate, and the catalytic water molecule. A
   mutation at one site can be compensated only if mutations at multiple other
   sites jointly reposition the catalytic machinery. This creates irreducible
   3-way and higher interactions.

2. **Protein stability** depends on the global free energy landscape, which
   has many-body terms from cooperative folding, desolvation entropy, and
   cavity packing. Stability effects of mutations are modulated by the entire
   sequence context, not just structural neighbors.

3. **The Weinreich study itself demonstrated extensive sign epistasis**,
   where the fitness effect of a mutation reverses depending on the background
   genotype. Sign epistasis distributed across 5 sites generically produces
   higher-order Walsh terms, because the sign-switching pattern cannot be
   captured by sums of pairwise interactions.

The composition gap is positive because structural contacts provide only the
pairwise skeleton, while the fitness landscape integrates stability, catalysis,
and dynamics into a highly nonlinear phenotype that encodes multi-way
dependencies. This is the protein analog of Proposition 2: composing pairwise
contacts through the protein's 3D structure creates higher-order epistasis
in the fitness readout.

**Confidence:** High.

---

## Prediction 4: Hall et al. yeast biosynthetic gene knockouts

**Question:** Hall et al. (2010) measured fitness for all 64 combinations of
6 biosynthetic gene knockouts in yeast. The metabolic pathway defines which
genes interact locally. What is the sign of the composition gap?

### (a) Predicted sign: POSITIVE (creation), but small in magnitude

### (b) Reasoning

Biosynthetic gene knockouts remove the cell's ability to synthesize specific
metabolites (amino acids, nucleotides, or vitamins). The metabolic pathway
defines local interactions through substrate-product relationships: gene A
produces the substrate for gene B, creating a direct pairwise dependency.
In linear biosynthetic pathways, most enzymes have 1--2 substrates, so local
interactions are predominantly pairwise (Walsh order 2 or below). Higher-order
local Walsh structure requires a single enzyme to depend on 3+ substrates
from different upstream genes, which is uncommon in linear pathways.

The global fitness landscape adds higher-order structure through:

1. **Metabolic flux balance.** Cellular growth rate is a nonlinear function
   of all available metabolite pools. The growth-rate function couples all
   biosynthetic pathways through shared precursor pools, ATP/NADH budgets,
   and ribosome allocation. These whole-cell constraints create 3-way and
   higher dependencies that no single local enzyme rule contains.

2. **Import competition.** Multiple auxotrophies force the cell to import
   several metabolites simultaneously. Import systems share membrane
   transporters, proton motive force, and regulatory pathways, creating
   interactions between knockouts in pathways that are locally independent.

3. **Diminishing-returns epistasis.** Systematic studies of biosynthetic
   knockouts in yeast show that fitness costs of sequential knockouts are
   subadditive (buffering epistasis). This diminishing-returns pattern is
   primarily pairwise and contributes to order-2 structure in both local
   and global spectra. It does not strongly contribute to order-3+.

The composition gap is positive because the local metabolic pathway
captures only direct enzymatic connections (pairwise), while the global
fitness landscape adds whole-cell flux constraints that couple all 6
pathways simultaneously. However, the magnitude should be small because
biosynthetic knockouts are dominated by main effects (order-1) and pairwise
buffering (order-2), with limited higher-order interaction energy in either
the local or global spectrum.

**Confidence:** Medium-low. The positive direction is predicted by the generic
argument that pathway composition creates higher-order terms (Proposition 2
analog). However, if the 6 genes lie in completely independent pathways with
no shared precursors, the global landscape could be nearly additive and
the composition gap could be indistinguishable from zero. The prediction is
sensitive to the specific pathway architecture, which I have not inspected.

---

## Summary table

| Prediction | Quantity | Predicted value | Confidence |
|------------|----------|-----------------|------------|
| 1 | Spearman(cycling_frac, delta_o3plus) | Positive | Medium |
| 2a | ODE sign preservation (>= 3/4) | Yes | Medium-high |
| 2b | ODE magnitude change | Decrease | Medium |
| 3 | Weinreich TEM-1 composition gap | Positive (creation) | High |
| 4 | Hall yeast composition gap | Positive (creation), small | Medium-low |
