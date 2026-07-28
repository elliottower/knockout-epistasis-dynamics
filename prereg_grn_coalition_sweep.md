# Pre-registration: GRN Coalition Sweep for EpistasisBench

## Motivation

Boolean gene regulatory network models provide two independent ground
truths that no other EpistasisBench regime offers:
1. Exhaustive interaction structure from 2^n knockout combinations
   (same Walsh-Hadamard decomposition as the transformer circuits)
2. The published, curated wiring diagram: which genes regulate which,
   with what sign (activation or repression)

The second ground truth enables the faithfulness question: does the
recovered interaction structure correspond to the real regulatory edges?
This is the exact question mech-interp circuit discovery asks but cannot
answer, because transformers have no curated ground-truth wiring to
validate against. GRNs do.

Knockout = ablation. Fixing a gene's update rule to 0 is the same
operation as zero-ablating an attention head. The attractor is the logit.
Coalition S = set of genes whose update rules are intact (not knocked out).

## Models

Four curated Boolean networks from the PyBoolNet repository, selected for:
- Biological relevance (cell cycle, apoptosis, signaling)
- Size in the tractable exhaustive range (10-13 nodes)
- Known wiring with published citations
- Diversity of expected interaction structure

| ID | Model | Nodes | Coalitions | Output node | Citation |
|----|-------|-------|------------|-------------|----------|
| G1 | faure_cellcycle | 10 | 1,024 | CycB | Faure et al., Bioinformatics 2006 |
| G2 | tournier_apoptosis | 12 | 4,096 | C3a | Tournier & Chaves, J Theor Biol 2009 |
| G3 | saadatpour_guardcell | 13 | 8,192 | KEV | Saadatpour et al., J Theor Biol 2010 |
| G4 | davidich_yeast | 10 | 1,024 | Cdc2_Cdc13 | Davidich & Bornholdt, PLoS ONE 2008 |

### Output node selection rationale

Each output node is the biologically canonical readout for that system:
- **CycB** (G1): Cyclin B marks mitotic entry; its activity is the
  cell-cycle "decision point"
- **C3a** (G2): Active Caspase-3 is the executioner caspase; its
  activation commits the cell to apoptosis
- **KEV** (G3): Potassium efflux through the vacuolar channel; the
  final guard-cell-closure readout
- **Cdc2_Cdc13** (G4): CDK-cyclin complex; drives the yeast cell
  into mitosis (same role as CycB in G1, different organism)

### Value function

For each coalition S and each of N_init = 512 random initial states:
1. Fix knocked-out genes (not in S) to 0 at initialization
2. Run synchronous Boolean update with knocked-out genes clamped to 0
   at every step
3. Run until fixed point or max 200 steps
4. Record the output node's value (0 or 1)

The value function v(S, init) is binary per initial state. The
pooled-mean v(S) = mean over initial states gives a continuous
quantity in [0, 1] representing the probability of the target
phenotype under knockout pattern S.

This is directly analogous to the transformer setting: initial
states are prompts, the binary output is whether the phenotype
fires, and the mean over initial states is the logit-diff averaged
over prompts.

## Predictions

Written AFTER running code tests (n_init=16) on all four models.
The code tests confirmed the script produces non-degenerate output
and revealed the approximate energy spectra. All four models' spectra
below are therefore POST-HOC at the 16-init-state resolution. The
registered predictions are about what holds under the FULL 512-init
run, and the faithfulness metrics (P3, P4) which have NOT been computed.

Code-test spectra (n_init=16, exploratory, not the registered experiment):

| Model | order-0 | order-1 | order-2 | order-3+ | unique vals |
|-------|---------|---------|---------|----------|-------------|
| G1 faure_cellcycle | 34.5% | 49.3% | 15.1% | 1.1% | 15 |
| G2 tournier_apoptosis | 14.8% | 37.7% | 32.3% | 15.2% | 12 |
| G3 saadatpour_guardcell | 3.2% | 15.7% | 31.1% | 50.0% | 6 |
| G4 davidich_yeast | 41.5% | 46.3% | 7.6% | 4.5% | 16 |

G3's 50% order-3+ was a surprise — see prediction revision below.

### P1: Energy spectrum diversity

The four models span at least two distinct interaction regimes:
- At least one model with >55% order-1 energy (near-additive)
- At least one model with <45% order-1 energy (interactive)

**Falsification:** All four models have order-1 energy within 10
percentage points of each other. This would mean Boolean network
topology does not produce diverse interaction structure, and the
GRN regime adds no benchmark value beyond the transformer circuits.

### P2: Order-1 energy predictions per model

| Model | Predicted order-1 | Predicted order-3+ | Rationale |
|-------|-------------------|---------------------|-----------|
| G1 faure_cellcycle | 49.3% (OBSERVED) | 0.97% (OBSERVED) | Post-hoc. CycB depends on cdh1 and Cdc20, which depend on multiple cyclins; moderate interaction depth. |
| G2 tournier_apoptosis | 35-55% | 2-10% | Apoptosis has strong positive feedback (C3a-C8a loop) and multi-step inhibition (TNF-NFkB-IAP-C3a). Feedback loops create higher-order interactions. Wider range because feedback dynamics are harder to predict. |
| G3 saadatpour_guardcell | 55-75% | <2% | Predominantly a linear signaling cascade (Ca2 -> PLC -> InsP3 -> CIS -> Ca2 feedback). Fewer parallel regulatory paths than cell cycle or apoptosis. |
| G4 davidich_yeast | 40-55% | <2% | Same biological function as G1 (cell cycle), different organism. CDK regulation is structurally similar. |

**Falsification criteria:**
- G2 order-3+ < 1% would mean feedback loops do not generate
  detectable higher-order epistasis in Boolean models
- G3 order-1 < 40% would mean the signaling cascade is more
  interactive than predicted, invalidating the "linear cascade =
  near-additive" hypothesis
- G4 and G1 disagree by >20 percentage points on order-1, despite
  modeling the same biological process

### P3: Wiring diagram faithfulness (the crown result)

For each model, compute the AUROC for recovering true regulatory
edges from Walsh pairwise coefficients:
- For each unordered pair (i, j), score = |w_{i,j}| (absolute
  Walsh pairwise coefficient from the pooled-mean value function)
- Label = 1 if the wiring diagram contains edge i->j OR j->i, else 0
- Compute AUROC over all C(n, 2) pairs

**Prediction:** AUROC > 0.60 for at least 3 of 4 models. Walsh
pairwise coefficients should recover true edges better than chance
(0.50) because the Walsh coefficient for a pair measures their
joint contribution to the output, and pairs connected by a
regulatory edge contribute jointly by construction.

**Falsification:** AUROC < 0.55 for all 4 models. This would mean
Walsh pairwise coefficients do not correspond to regulatory edges,
which is a genuine negative finding (structural wiring does not
predict functional epistasis).

**Nuance:** Edge density affects AUROC interpretation. G1 has
density 31/90 = 0.34, meaning one-third of all pairs are connected.
High density compresses the AUROC toward 0.50. Report edge density
alongside AUROC, and also report precision@k for k = n_edges
(the biologically meaningful operating point).

### P4: Motif enrichment in higher-order coefficients

For each model with detectable order-3+ energy (>0.5%):
- Take the top-10 three-way Walsh coefficients by magnitude
- For each triple (i, j, k), check whether the wiring diagram
  contains a connected subgraph (path or motif) among the three nodes
- Compare the fraction of connected triples to a permutation null
  (1000 random draws of 10 triples from the same order-3 set)

**Prediction:** Connected triples are enriched in the top-10 vs
null (one-sided p < 0.05 by permutation test) for models with
order-3+ > 1%.

**Falsification:** No enrichment in any model. This would mean
higher-order Walsh structure does not correspond to regulatory
motifs.

### P5: Cross-regime rank stability (deferred)

This prediction covers the full benchmark, not just GRN. Deferred
until all regimes (transformer, protein, GRN) have coalition tables
and at least 4 methods are scored on all of them.

Prediction: method rankings (by aggregate detection score across
metrics) correlate at Spearman rho > 0.5 between GRN regimes and
transformer circuit regimes. Methods that detect real structure
should do so regardless of domain.

## Experimental protocol

### Step 1: Generate coalition tables

Run `grn_coalition_sweep.py` for G2, G3, G4 (G1 already done as
code test). Parameters:
- `--n-init 512` (512 random initial states)
- `--seed 42`
- `--max-steps 200`

Store results in `results/` within the epistasis-bench repo.

### Step 2: Walsh decomposition + wiring extraction

For each model:
1. Compute pooled-mean value function v(S) = mean over initial states
2. Compute normalized WHT coefficients w_T = (1/2^n) * WHT(v)[T]
3. Compute energy spectrum (fraction of variance at each order)
4. Extract k_99 sparsity (number of coefficients capturing 99% of energy)
5. Extract interaction graph from the Boolean rules (already done
   during sweep — stored in `*_wiring.json`)

### Step 3: Faithfulness scoring

1. Compute pairwise Walsh coefficients |w_{i,j}| for all C(n,2) pairs
2. Load wiring diagram edge set
3. Compute AUROC, precision@k, recall@k
4. Compute three-way motif enrichment (P4)

### Step 4: Report

Table: G1/G2/G3/G4 x order-1/order-2/order-3+/AUROC/precision@n_edges

## Confound controls

1. **Self-regulation.** Some nodes regulate themselves (CycD, TNF, Start).
   Self-edges are excluded from the wiring-diagram ground truth because
   the Walsh coefficient w_{i} (order-1) measures self-importance, not
   pairwise interaction.

2. **Input nodes.** CycD (G1), TNF (G2), Start (G4) are external inputs
   (their update rule is self-referential: f(x) = x). They persist
   indefinitely once set. Report results both with and without input
   nodes in the player set, since input nodes that persist create a
   different kind of "knockout" than nodes whose state depends on
   the network.

3. **Attractor multiplicity.** Synchronous update from different initial
   states may converge to different attractors, including limit cycles.
   We use a 200-step ceiling; states that have not converged are
   recorded as-is. Report the fraction of coalition-x-init-state
   combinations that converge to a fixed point vs oscillate.

4. **Knockout vs knockdown.** We fix knocked-out nodes to 0, which is
   a complete loss-of-function. An alternative is stochastic knockdown
   (node fires with probability p < 1). This is analogous to mean
   ablation vs zero ablation. We use zero-knockout in v1 for simplicity
   and note this as a limitation.

## What this pre-registration does NOT cover

- The budget-vs-accuracy experiment (requires method scoring framework)
- Specific detection methods scored against the GRN ground truth
- The Poelwijk 2019 protein landscape
- The Weinreich 2018 multi-landscape collection
- Transformer circuit regimes (separate pre-registration exists)
- Any cross-domain comparison

Those are separate pre-registrations that depend on this one producing
non-degenerate coalition tables with diverse interaction structure.

## Code freeze checklist

Before computing the SHA-256:

1. [ ] grn_coalition_sweep.py runs for all 4 models without error
2. [ ] Output npz files pass shape assertion (2^n coalitions)
3. [ ] Wiring JSON contains expected edge counts per model
4. [ ] data_utils.py Walsh functions verified against mobius_wht.py
   on a shared test vector
5. [ ] G1 results from code test match embedded predictions above
6. [ ] All four models have non-degenerate value functions (more than
   2 unique values in pooled mean)

Seeds: initial-state seed 42; no other randomization in the sweep.

## References

- Faure et al., Bioinformatics 22.14 (2006): e124-e131
- Tournier & Chaves, J Theor Biol 260.2 (2009): 196-209
- Saadatpour et al., J Theor Biol 266.4 (2010): 641-656
- Davidich & Bornholdt, PLoS ONE 3.2 (2008): e1672
- PyBoolNet: Klarner et al., Bioinformatics 33.5 (2017): 770-772
