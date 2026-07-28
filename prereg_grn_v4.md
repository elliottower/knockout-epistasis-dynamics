# Pre-registration: GRN Coalition Sweep v4

**Changes from v3:** Added quantitative gate criteria (three-part:
abs_e3+ >= 1e-4, unique >= 10, no single-value modal fraction > 90%).
Output-node rule locked: use original publication readout, no
substitution. Added cross-network hypothesis P9 (positive feedback
fraction predicts order-3+ creation). Separated exploratory set
(G1, G2, G4 — already observed) from confirmatory set (gate
survivors among G5-G10 — not yet run). Per-network framing: case
studies with mechanistic explanations; cross-network patterns
reported as exploratory/suggestive regardless of N.

## Motivation

Boolean gene regulatory network models provide two independent ground
truths that no other EpistasisBench regime offers:
1. **Exhaustive interaction structure** from 2^n knockout combinations
   (same Walsh-Hadamard decomposition as the transformer circuits)
2. **Boolean update rules** that define the exact functional
   dependencies between genes: if gene k's rule is `f_k = i AND j`,
   then i and j interact by construction

The second ground truth is the faithfulness question: does the
recovered interaction structure correspond to the real functional
dependencies? This is the exact question mech-interp circuit
discovery asks but cannot answer, because transformers have no wiring
diagram validated by independent experimental evidence. GRNs do.

Knockout = ablation. Fixing a gene's update rule to 0 is the same
operation as zero-ablating an attention head. Fixing a gene to 1
(constitutive activation) is the gain-of-function complement. The
attractor is the logit. Coalition S = set of genes whose update
rules are intact (not clamped).


## Models

Nine curated Boolean networks from PyBoolNet and biodivine repositories,
selected for:
- Biological relevance (cell cycle, apoptosis, signaling, cell fate)
- Size in the tractable exhaustive range (10-15 nodes)
- Known wiring with published citations
- Mean regulators per node >= 2.0, at least 4 nodes with 3+ regulators
- Diversity of expected interaction structure

| ID | Model | Dyn | Inp | Coalitions | Output | Citation |
|----|-------|-----|-----|------------|--------|----------|
| G1 | faure_cellcycle | 10 | 1 (CycD) | 1,024 | CycB | Faure et al., Bioinformatics 2006 |
| G2 | tournier_apoptosis | 12 | 1 (TNF) | 4,096 | C3a | Tournier & Chaves, J Theor Biol 2009 |
| G4 | davidich_yeast | 10 | 1 (Start) | 1,024 | Cdc2_Cdc13 | Davidich & Bornholdt, PLoS ONE 2008 |
| G5 | drosophila_cellcycle | 11 | 3 (Ago, CycD, Notch) | 16,384 | CycB | Faure et al., J Theor Biol 2008 |
| G6 | myeloid_progenitors | 11 | 0 | 2,048 | GATA1 | Krumsiek et al., PLoS ONE 2011 |
| G7 | blood_stem_cell | 11 | 0 | 2,048 | GATA1 | Bonzanni et al., Bioinformatics 2013 |
| G8 | emt_switch | 12 | 0 | 4,096 | Ecadherin_mRNA | Steinway et al., PLoS Comput Biol 2014 |
| G9 | fanconi_anemia | 15 | 0 | 32,768 | CHKREC | Rodriguez et al., BMC Syst Biol 2015 |
| G10 | arabidopsis_cellcycle | 14 | 0 | 16,384 | CYCB1_1 | Ortiz-Gutierrez et al., PLoS Comput Biol 2015 |

G3 (saadatpour_guardcell) was dropped: its 13 nodes are mostly
single-regulator copy chains (mean regs/node = 1.5), producing
degenerate interaction structure dominated by chain-AND artifacts
rather than genuine multi-input logic.

### Input node handling

Nodes with self-referential rules (f(x) = x) are external inputs
that persist indefinitely once set. These are identified automatically
by `identify_input_nodes()`. When `--exclude-inputs` is passed:
1. Input nodes are removed from the player set
2. Input nodes are clamped to a fixed configuration (0 or 1) for
   the entire sweep
3. The sweep runs separately for each input configuration of interest
4. Coalition bitmasks map to dynamic-node positions only

Models with inputs (G1, G2, G4, G5) are swept twice: once with all
inputs OFF (clamp=0) and once with all inputs ON (clamp=1). Models
without inputs (G6-G10) are swept once per ablation type.

### Output node selection rationale

Each output node is the biologically canonical readout for that system:
- **CycB** (G1, G5): Cyclin B marks mitotic entry
- **C3a** (G2): Active Caspase-3 commits the cell to apoptosis
- **Cdc2_Cdc13** (G4): CDK-cyclin complex drives yeast mitosis
- **GATA1** (G6, G7): Erythroid commitment transcription factor
- **Ecadherin_mRNA** (G8): Epithelial marker (high = epithelial state)
- **CHKREC** (G9): Checkpoint recovery (successful DNA repair)
- **CYCB1_1** (G10): B-type cyclin, mitotic commitment marker


## Value function

For each coalition S and each of N_init = 512 random initial states:
1. Fix clamped genes at initialization and at every subsequent step
2. Run Boolean dynamics until convergence or max steps
3. Record the output node's value (0 or 1)

The value function v(S, init) is binary per initial state. The
pooled-mean v(S) = mean over initial states gives a continuous
quantity in [0, 1] representing the probability of the target
phenotype under ablation pattern S.

This is directly analogous to the transformer setting: initial states
are prompts, the binary output is whether the phenotype fires, and the
mean over initial states is the logit-diff averaged over prompts.

### Ablation types

Two ablation modes are registered:

1. **Clamp-to-0 (knockout):** Genes not in the coalition have their
   update rule forced to 0. This is loss-of-function for activators
   but gain-of-function for downstream targets of repressors (removing
   the brake). Analogous to zero-ablation.

2. **Clamp-to-1 (constitutive activation):** Genes not in the
   coalition have their update rule forced to 1. This is
   gain-of-function for the clamped gene and loss-of-function for
   targets that depend on `!gene`. Analogous to the mean-ablation axis
   in the transformer setting, where the ablation value matters.

Both ablation types produce separate coalition tables, value functions,
Walsh decompositions, and energy spectra. Differences between the two
spectra reveal where the knockout/activation asymmetry is structurally
important.

### Update scheme

Two update schemes are registered:

1. **Synchronous:** All nodes update simultaneously at each step.
   This is the default for G1 (faure_cellcycle) and G4
   (davidich_yeast), which were published as synchronous models.

2. **Asynchronous (random-order sequential):** Nodes update one at a
   time in a random permutation, one full permutation per "step."
   This is the published update scheme for G2 (tournier_apoptosis).
   Under synchronous update, some models produce spurious limit
   cycles that do not appear in asynchronous dynamics.

Both update schemes are run for all models. The PRIMARY result uses
each model's published scheme. Disagreements between sync and async
results are flagged and reported. For G2, synchronous results are
labeled SECONDARY.


## Structural analysis of Boolean update rules

The predictions below are derived entirely from reading the Boolean
rules in `grn_coalition_sweep.py`. No data has been generated and no
code has been run. The analysis traces feedback loops, pathway depth,
AND-gate arity, and redundancy structure to predict the interaction
order spectrum.


### G1: faure_cellcycle (10 nodes)

**Output rule:** `CycB = !cdh1 & !Cdc20`

A 2-input AND of negations. CycB is active when both its inhibitors
(cdh1 and Cdc20) are off.

**Regulatory chain to output:**
- Cdc20 = CycB (direct copy, 1-step feedback)
- cdh1 = p27&!CycB | !CycB&!CycA | Cdc20 (depends on 4 variables)
- p27 = p27&!CycE&!CycD&!CycB | ... (4-input AND terms, self-loop)
- CycA has a complex 4-clause rule involving cdh1, Rb, E2F, UbcH10,
  Cdc20, and CycA itself
- E2F = p27&!Rb&!CycB | !Rb&!CycB&!CycA
- Rb = !CycE&!CycD&!CycB&!CycA | p27&!CycD&!CycB

**Feedback loops (3):**
1. CycB -> Cdc20 -> cdh1 (via Cdc20 term in cdh1) -> CycB. Sign:
   CycB=1 sets Cdc20=1, Cdc20=1 raises cdh1, cdh1=1 kills CycB.
   Negative feedback.
2. CycB -> cdh1 (CycB appears negated in cdh1 rule) -> CycB. Sign:
   CycB=1 suppresses cdh1, less cdh1 sustains CycB. Positive feedback.
3. CycB -> E2F (!CycB in E2F) -> Rb -> CycE -> ... -> CycB. Long
   chain through the Rb/E2F/cyclin module. Mixed sign.

**Nodes upstream of output:** All 10 nodes are in a single connected
component. CycD is an external input (self-referential rule: CycD =
CycD) that feeds through Rb -> E2F -> CycA -> cdh1 -> CycB.

**AND-gate analysis:** The output is a 2-input AND. The deepest
upstream ANDs are 4-input (p27's rule, Rb's first clause). Because
these are 2-4 steps upstream of the output, their higher-order
structure is attenuated by intervening OR gates and feedback before
reaching CycB.

**Depth from input to output:** CycD -> Rb -> E2F -> CycA -> cdh1 ->
CycB gives a longest non-feedback path of ~5 steps.

**Prediction (clamp-to-0):**
- order-0 (intercept): 25-40%
- order-1 (main effects): 40-55%
- order-2 (pairwise): 15-30%
- order-3+: 2-10%

**Rationale:** The 2-input AND at the output bottleneck limits direct
higher-order interactions at the output level to pairwise between
cdh1 and Cdc20. Upstream multi-input AND gates propagate some
higher-order structure, but the intervening OR clauses and the short
negative feedback loop (Cdc20 = CycB) dampen it. One external input
(CycD) contributes mainly order-1 energy through its long chain to
the output. The single negative and single positive feedback loop
create moderate bistability, enough for a non-trivial attractor
landscape but not the rich switching seen in systems with more
feedback.


### G2: tournier_apoptosis (12 nodes)

**Output rule:** `C3a = !IAP & C8a`

A 2-input AND. C3a fires when IAP is absent and C8a is present.

**Regulatory chain to output:**
- C8a = T2&!CARP | !CARP&C3a (feedback from C3a in second clause)
- IAP = !TNF&NFkBnuc | !TNF&!C3a | NFkBnuc&!C3a (feedback from C3a)
- T2 = TNF&!FLIP (2-input AND)
- FLIP = NFkBnuc (copy)
- CARP = !TNF&NFkB | !TNF&!C3a | NFkBnuc&!C3a
- IKKa = TNF&!C3a&!A20a (3-input AND)
- A20a = TNF&NFkBnuc (2-input AND)
- IkB = !TNF&NFkBnuc | !TNF&!IKKa | NFkBnuc&!IKKa
- NFkB = !IkB
- NFkBnuc = NFkB&!IkB (equivalent to !IkB & !IkB = !IkB = NFkB;
  both require IkB off)

**Feedback loops (4+):**
1. C3a -> C8a -> C3a (via `!CARP&C3a` clause). Positive: C3a
   directly sustains C8a, C8a sustains C3a. Self-amplifying
   executioner loop.
2. C3a -> IAP -> C3a. `!C3a` in IAP rule means C3a=1 suppresses IAP,
   less IAP frees C3a. Positive (double negative).
3. C3a -> CARP -> C8a -> C3a. `!C3a` in CARP rule means C3a removes
   CARP, helping C8a, helping C3a. Positive.
4. C3a -> IKKa (`!C3a` term) -> IkB -> NFkB -> NFkBnuc -> IAP ->
   C3a. Long chain: C3a suppresses IKKa, less IKKa raises IkB, IkB
   kills NFkB/NFkBnuc, less NFkBnuc can reduce IAP (but IAP also has
   `!C3a` terms). Mixed-sign, depth 6.

Three of four feedback loops are positive. This is a classic
bistability signature: once apoptosis fires (C3a=1), three
independent reinforcement mechanisms lock it in. This creates a
sharp threshold in the knockout landscape.

**Two competing pathways from TNF to C3a:**
- Activation: TNF -> T2 (via !FLIP) -> C8a -> C3a (depth 3)
- Survival: TNF -> IKKa -> IkB -> NFkB -> NFkBnuc -> FLIP/IAP/CARP
  (depth 5-6), which inhibits C3a through IAP and C8a through
  FLIP/CARP

These competing pathways create pathway-level redundancy and
interference. Knocking out nodes in one pathway vs the other has
opposite effects on C3a.

**AND-gate analysis:** IKKa = TNF&!C3a&!A20a is the only 3-input AND
in the network. It creates a genuine 3-way interaction among TNF,
C3a, and A20a. The output rule is 2-input AND, but the upstream
3-input AND and the multiple feedback loops propagate higher-order
structure to the output.

**Nodes upstream of output:** 11 of 12 nodes regulate C3a directly
or indirectly; only TNF is an external input (TNF = TNF).

**Depth from input to output:** 3 steps (activation path) to 6 steps
(survival path).

**Prediction (clamp-to-0):**
- order-0: 10-25%
- order-1: 25-40%
- order-2: 25-40%
- order-3+: 10-25%

**Rationale:** Three positive feedback loops create a bistable
landscape with a sharp threshold between survival and apoptosis
attractors. Different knockout patterns push the system to opposite
sides of this threshold, generating diverse attractor-dependent
output values and rich multi-way interactions. The two competing
pathways (activation vs survival) create strong pairwise interactions
between nodes in different pathways. The 3-input AND in IKKa
propagates genuine order-3 structure. The deep chain through the NFkB
module (depth 6) creates multi-step dependencies that generate
moderate higher-order interactions.

**Key difference from G1:** G2 has three positive feedback loops
(vs one in G1), two competing signal pathways (vs one dominant
chain), and a 3-input AND gate. All three features increase
higher-order interaction energy.


### G4: davidich_yeast (10 nodes)

**Output rule:** `Cdc2_Cdc13 = !Ste9&!Rum1&!Slp1 | !Ste9&!Rum1&Cdc2_Cdc13_A`

Equivalent to: `!Ste9 & !Rum1 & (!Slp1 | Cdc2_Cdc13_A)`.
A 4-variable function with a 3-input AND core (!Ste9 & !Rum1 &
!Slp1) OR'd with a positive-feedback clause.

**Regulatory chain to output:**
- Ste9 and Rum1 have IDENTICAL update rules:
  `!SK&!Cdc2_Cdc13&!Cdc2_Cdc13_A | !SK&!Cdc2_Cdc13&PP | !Cdc2_Cdc13&!Cdc2_Cdc13_A&PP`
- SK = Start (copy of external input)
- Start = Start (external input)
- Cdc2_Cdc13_A = Cdc2_Cdc13 & Cdc25 & !Wee1_Mik1 (3-input AND)
- Cdc25 = Cdc2_Cdc13 | Cdc2_Cdc13_A (OR)
- Wee1_Mik1 = !Cdc2_Cdc13 & !Cdc2_Cdc13_A (2-input AND of negations)
- Slp1 = Cdc2_Cdc13_A (copy)
- PP = Slp1 (copy of copy: PP = Cdc2_Cdc13_A with 2-step delay)

**Feedback loops (3):**
1. Cdc2_Cdc13 -> Ste9/Rum1 (inhibits via !Cdc2_Cdc13 terms) ->
   Cdc2_Cdc13. Positive: CDK suppresses its own inhibitors.
2. Cdc2_Cdc13 -> Cdc2_Cdc13_A (via 3-input AND with Cdc25,
   !Wee1_Mik1) -> Cdc2_Cdc13. Positive: CDK activates its activated
   form, which feeds back to sustain CDK.
3. Cdc2_Cdc13 -> Cdc2_Cdc13_A -> Slp1 -> PP -> Ste9/Rum1 ->
   Cdc2_Cdc13. Negative: activated CDK triggers Slp1/PP, which
   reactivates inhibitors, suppressing CDK. Depth 5.

Two positive feedback loops and one negative feedback loop. The
balance between positive (locking CDK on) and negative (delaying
then shutting CDK off) drives the cell cycle oscillation.

**AND-gate analysis:** The output's first clause is a 3-input AND
(!Ste9 & !Rum1 & !Slp1). However, Ste9 and Rum1 have identical
rules, so they are functionally redundant: their contributions to
the Walsh decomposition are symmetric. The effective interaction
between them is simple (knocking out either one partially relieves
inhibition, knocking out both fully relieves it). The 3-input AND
for Cdc2_Cdc13_A (Cdc2_Cdc13 & Cdc25 & !Wee1_Mik1) creates
additional 3-way structure upstream.

**Nodes upstream of output:** All 10 nodes participate in a single
connected component. Start is the sole external input.

**Depth from input to output:** Start -> SK -> (inhibits Ste9/Rum1)
-> Cdc2_Cdc13 gives depth 3. The negative feedback loop extends
depth to 5.

**Key structural feature: Ste9/Rum1 redundancy.** Because their
rules are identical, the Walsh spectrum treats them symmetrically.
The pairwise interaction between Ste9 and Rum1 reflects the
redundancy of having two copies of the same inhibitor: knocking out
one partially relieves inhibition, knocking out both fully relieves
it. This is an AND-like interaction (both must be removed for full
effect).

**Prediction (clamp-to-0):**
- order-0: 30-45%
- order-1: 40-55%
- order-2: 15-30%
- order-3+: 3-12%

**Rationale:** The output's 3-input AND clause creates up to order-3
interactions, but the Ste9/Rum1 symmetry means two of the three
variables in that AND are interchangeable, reducing the effective
interaction complexity. The positive feedback loops create moderate
bistability (CDK on vs off attractors) but less extreme than G2
(only two positive loops vs three, and no competing activation/
survival pathways). Start as an external input contributes mainly
order-1 energy. The 3-input AND in Cdc2_Cdc13_A adds some order-3
structure but it is attenuated by the downstream OR in the output
rule. Overall spectrum should resemble G1, consistent with both
models implementing CDK-based cell cycle control in different
organisms.


### Prediction summary (clamp-to-0, synchronous update)

Structural analysis predictions for the original models (G1, G2, G4).
New models (G5-G10) have not been structurally analyzed at this level
of detail; their spectra are exploratory.

| Model | order-0 | order-1 | order-2 | order-3+ | Primary driver |
|-------|---------|---------|---------|----------|----------------|
| G1 faure_cellcycle | 25-40% | 40-55% | 15-30% | 2-10% | 2-input AND output, moderate feedback |
| G2 tournier_apoptosis | 10-25% | 25-40% | 25-40% | 10-25% | 3 positive feedback loops, competing pathways |
| G4 davidich_yeast | 30-45% | 40-55% | 15-30% | 3-12% | 3-input AND with redundant inhibitors |


### Clamp-to-1 predictions

Constitutive activation (clamp-to-1) should produce a different
spectrum from knockout (clamp-to-0) because clamping a repressor to
1 permanently blocks its target, while clamping it to 0 relieves the
target.

**G1:** CycB = !cdh1 & !Cdc20. Clamping cdh1 or Cdc20 to 1
permanently kills CycB. Clamping them to 0 frees CycB. The asymmetry
is large because the output rule uses only negated inputs. Prediction:
clamp-to-1 shifts energy toward order-1 (clamping any single
inhibitor is sufficient to kill output) and reduces order-2 (the
AND interaction between cdh1 and Cdc20 becomes redundant when either
one being ON kills output).

**G2:** The most complex asymmetry. The survival pathway (IAP, CARP,
FLIP) contains repressors of C3a. Clamping IAP to 1 permanently
blocks C3a regardless of C8a. Clamping IAP to 0 frees C3a. The
activation pathway (T2, C8a) contains activators. Clamping T2 to 1
forces activation regardless of FLIP. The mix of activators and
repressors on different pathways makes clamp-to-1 vs clamp-to-0 most
divergent for G2. Prediction: order-1 increases (clamping any single
repressor ON is sufficient to block the output), order-3+ may
increase because constitutive activation of upstream nodes breaks
the normal feedback-loop dynamics.

**G4:** Similar asymmetry to G1. The output depends on !Ste9, !Rum1,
!Slp1 (all negated). Clamping any of these to 1 blocks the output.
Prediction: clamp-to-1 shifts energy toward order-1, similar to G1.


## Ground truths

### PRIMARY: Boolean update rules (CONFIRMATORY)

The Boolean rules define exact functional dependencies between genes.
Two genes i and j functionally interact in gene k's regulation if
they co-occur in a product term (AND clause) of gene k's rule. This
ground truth is:
- **Exact**: derived analytically from the rules, not from noisy
  measurements
- **Multi-order**: AND clauses with 3+ inputs define genuine
  higher-order interactions
- **Asymmetric**: interactions are per-target-gene, not symmetric
  pairwise

**Interaction extraction protocol (v3: exact local Fourier).** For
each gene k, compute the Walsh-Hadamard transform of k's truth table
(a 2^m vector where m = number of regulators of k). Each Walsh
coefficient w_T gives the exact interaction strength among the
regulators indexed by T. This is representation-independent: it does
not depend on how the Boolean expression is written (DNF, CNF, or
arbitrary). The coefficient magnitudes give continuous interaction
strengths, replacing the binary AND-clause labeling from v2.

Example: IKKa = TNF & !C3a & !A20a. The WHT of IKKa's truth table
(8 entries for 3 regulators) produces coefficients at all 2^3
indices. The order-3 coefficient (TNF, C3a, A20a) = 1/8, and the
three pairwise coefficients each = 1/8. This gives continuous
strengths that sum across all genes to produce global pairwise and
triple scores.

For pairwise faithfulness (P3), the score for pair (i,j) is the sum
of |local w_{i,j}| across all genes where both i and j are
regulators. For triple faithfulness (P4), the score for triple
(i,j,k) is analogous.

### SECONDARY: Wiring diagram (edge list)

The edge list records which genes regulate which, extracted from the
Boolean rules via `extract_interaction_graph()`. This is a structural
measure: an edge from i to j means i appears in j's Boolean rule.

The wiring diagram is a SUPERSET of the functional interaction
structure. Two genes can share an edge (gene i appears in gene j's
rule) without functionally interacting (they appear in different OR
clauses, so they contribute independently). The wiring diagram is
kept as a secondary result to answer: "how well does functional
epistasis track structural wiring?"


## Predictions

### Status labels

- **CONFIRMATORY:** The faithfulness test (P3, P4) using rules-based
  ground truth. This analysis has never been computed.
- **EXPLORATORY:** The energy spectrum predictions (P1, P2) and all
  clamp-to-1 results. Code tests were run on earlier versions of this
  protocol (results deleted, but the analysis has been seen). Any
  spectrum prediction, even if it differs from prior observations, is
  exploratory because the analyst has been exposed to pilot data.


### P1: Energy spectrum diversity (EXPLORATORY)

The nine models span at least two distinct interaction regimes:
- At least one model with >50% order-1 energy (near-additive)
- At least one model with <30% order-1 energy (interactive)

**Falsification:** All nine models have order-1 energy within 10
percentage points of each other. This would mean Boolean network
topology does not produce diverse interaction structure, and the
GRN regime adds no benchmark value beyond the transformer circuits.


### P2: Order spectrum predictions per model (EXPLORATORY)

Predictions from the structural analysis above (clamp-to-0,
synchronous update):

| Model | order-1 | order-3+ | Structural basis |
|-------|---------|----------|------------------|
| G1 faure_cellcycle | 40-55% | 2-10% | 2-input AND output, moderate feedback |
| G2 tournier_apoptosis | 25-40% | 10-25% | 3 positive feedback loops, competing pathways, 3-input AND |
| G4 davidich_yeast | 40-55% | 3-12% | 3-input AND with redundant inhibitors |

**Falsification criteria:**
- G2 order-3+ < 5% would mean multiple positive feedback loops do
  not generate detectable higher-order epistasis in Boolean models,
  falsifying the "feedback drives higher-order" hypothesis.
- G4 and G1 differ by > 25 percentage points on order-1, despite
  modeling the same biological process (CDK-driven cell cycle entry)
  in different organisms.


### P3: Faithfulness to Boolean rules (CONFIRMATORY)

For each model, compute the Spearman rank correlation between
local-rule pairwise interaction strengths and global (attractor-level)
pairwise Walsh coefficient magnitudes:

- **Local pairwise strengths:** For each unordered pair (i, j),
  compute the exact Walsh-Hadamard transform of each gene's truth
  table. Sum |w_{i,j}| across all genes where both i and j appear
  as regulators. This gives continuous local interaction strength,
  replacing the binary AND-clause labeling from v2.
- **Global pairwise magnitudes:** |w_{i,j}| from the normalized WHT
  of the pooled-mean value function.
- Compute Spearman rho over all C(n, 2) pairs, with 2000-bootstrap
  95% confidence intervals.

**Prediction:** Spearman rho > 0.3 for at least 6 of 9 models.
Local Fourier pairwise strengths should rank-correlate with
attractor-level pairwise structure because genes with strong local
multiplicative interactions propagate that structure through dynamics.
The correlation need not be perfect: composition can destroy local
interactions (via feedback cancellation) or create new ones (via
multi-gene pathway convergence).

**Falsification:** Spearman rho < 0.1 and CI includes 0 for all 9
models. This would mean local rule interaction structure has no
detectable relationship to knockout-derived epistasis.

**Advantages over v2 AUROC:** Continuous strengths avoid the base-rate
problem (models with many AND clauses had compressed AUROC). Spearman
is rank-based and insensitive to monotone rescaling. Bootstrap CIs
quantify estimation uncertainty directly.


### P4: Motif enrichment in higher-order coefficients (CONFIRMATORY)

For each model with detectable order-3+ energy passing the spectrum
gate (absolute order-3+ energy > 1e-4):
- Take the top-10 three-way Walsh coefficients by magnitude
- For each triple (i, j, k), check whether the triple co-occurs in
  any single AND clause of any gene's Boolean rule (rule-derived
  three-way interaction), computed via WHT of each gene's truth table
- Compare the fraction of rule-matched triples to a permutation null
  (1000 random draws of 10 triples from the same order-3 set)

**Spectrum gate:** Absolute order-3+ energy (sum of squared
coefficients at order >= 3) must exceed 1e-4 to qualify. This
prevents testing on models where order-3+ structure is numerical
noise.

**Prediction:** Rule-matched triples are enriched in the top-10 vs
null (one-sided p < 0.05 after Benjamini-Yekutieli correction) for
models passing the spectrum gate.

**Falsification:** No enrichment in any qualifying model.


### P5: Sync vs async disagreement (EXPLORATORY)

For all models, compute the energy spectrum under both synchronous
and asynchronous update. Report the L1 distance between the two
spectra.

**Prediction:** For G2 (tournier_apoptosis, published async), sync
and async spectra differ by < 5 percentage points. The positive
feedback loops in G2 drive convergence to fixed-point attractors
(survival or apoptosis), making the system less sensitive to update
order. Models with oscillation-prone topology (tight negative
feedback loops) may show larger sync/async divergence.


### P6: Clamp-to-0 vs clamp-to-1 asymmetry (EXPLORATORY)

For each model, compute the L1 distance between the clamp-to-0 and
clamp-to-1 energy spectra.

**Prediction:** G2 has the largest asymmetry among the original
models because it has both activator and repressor pathways to the
output. Clamping a repressor ON vs OFF has opposite effects. G1 and
G4 have moderate asymmetry because their output rules use only
negated inputs.

For new models (G5-G10), asymmetry is exploratory. Models with
mixed positive/negative regulation of the output node are expected
to show larger asymmetry than models with uniform regulation sign.


### P7: Wiring diagram faithfulness (secondary ground truth, EXPLORATORY)

AUROC for recovering true regulatory edges from Walsh pairwise
coefficients, using the wiring diagram (edge list) as ground truth:
- For each unordered pair (i, j), score = |w_{i,j}|
- Label = 1 if the edge list contains edge i->j OR j->i, else 0
- Compute AUROC over all C(n, 2) pairs

**Prediction:** AUROC > 0.55 for at least 4 of 9 models. Edge-based
AUROC should be lower than rules-based Spearman rho (P3) because the
edge list includes non-interacting regulatory connections.


### P8: Cross-regime rank stability (deferred)

Deferred until all regimes (transformer, protein, GRN) have coalition
tables and at least 4 methods are scored on all of them. Prediction:
method rankings correlate at Spearman rho > 0.5 between GRN and
transformer circuit regimes.


### P9: Cross-network composition hypothesis (SPLIT: EXPLORATORY + CONFIRMATORY)

**Hypothesis.** Networks with a higher fraction of positive feedback
loops (among all feedback loops up to length 4) show greater
order-3+ energy creation (positive delta_o3+ = global minus local
order-3+ energy fraction).

**Rationale (from exploratory pilot, G1/G2/G4).**
- G2 (tournier_apoptosis): 3 positive / 8 total = 38% positive
  feedback fraction. delta_o3+ = +8.5pp (creates).
- G1 (faure_cellcycle): 13 positive / 25 total = 52% positive
  feedback fraction. delta_o3+ = -3.0pp (destroys).
- G4 (davidich_yeast): 13 positive / 19 total = 68% positive
  feedback fraction. delta_o3+ = -0.4pp (neutral).

NOTE: The pilot data (G1/G2/G4) does NOT support the hypothesis —
the network with the lowest positive feedback fraction (G2) shows the
most creation. This is recorded honestly. The hypothesis is motivated
by the mechanistic argument (positive feedback creates bistability,
which generates multi-way threshold interactions), which the pilot
data may be too small to test. The positive feedback fraction may be
a poor predictor; an alternative is the number of positive feedback
loops involving the output node specifically.

**Exploratory set (already observed):** G1, G2, G4. These results
are pilot data that motivated the hypothesis. They cannot confirm
or falsify it.

**Confirmatory set (not yet observed):** Whichever of G5-G10 survive
the three-part quality gate. For these networks only, the feedback
loop counts below were computed from topology before running any
coalition sweep:

| Model | Pos FB | Neg FB | Mixed FB | Total | Pos fraction |
|-------|--------|--------|----------|-------|--------------|
| G5 drosophila_cellcycle | 5 | 2 | 29 | 36 | 14% (or 71% excl. mixed) |
| G6 myeloid_progenitors | 6 | 10 | 0 | 16 | 38% |
| G7 blood_stem_cell | 51 | 4 | 0 | 55 | 93% |
| G8 emt_switch | 15 | 25 | 0 | 40 | 38% |
| G9 fanconi_anemia | 48 | 37 | 13 | 98 | 49% (or 56% excl. mixed) |
| G10 arabidopsis_cellcycle | 16 | 43 | 22 | 81 | 20% (or 27% excl. mixed) |

**Test (on confirmatory set only):** Spearman correlation between
positive feedback fraction and delta_o3+. One-tailed (positive
direction). Alpha = 0.05. Report the point estimate and bootstrap
95% CI. With N <= 6, this is underpowered; a non-significant result
is uninformative (cannot distinguish "no relationship" from "too few
data points"). A significant result is suggestive, not definitive.

**Alternative predictor (exploratory):** Number of positive feedback
loops involving the output node, rather than total positive feedback
fraction. Reported alongside the primary predictor but not used for
the formal test.

**Framing commitment.** Regardless of the final N surviving the gate,
the paper presents per-network case studies with mechanistic
explanations as the primary structure. Any cross-network pattern
(including P9) is reported in a separate "exploratory cross-network"
section, explicitly labeled as suggestive. The N is stated and the
underpowered nature acknowledged. No regression is presented as a
main result.


## Experimental protocol

### Step 1: Generate coalition tables

Run `grn_coalition_sweep.py` for all 9 models, in both ablation modes
and both update schemes:

- **Ablation modes:** clamp-to-0 (knockout), clamp-to-1 (constitutive
  activation)
- **Update schemes:** synchronous, asynchronous (random-order
  sequential)
- **Input handling:** for models with input nodes (G1, G2, G4, G5),
  run with `--exclude-inputs` and sweep each input configuration
  (all-off, all-on). Input nodes are clamped and excluded from the
  player set.
- **Parameters:** `--n-init 512`, `--seed 42`, `--max-steps 200`
- **Limit cycle handling:** if the trajectory has not converged by
  max_steps, detect the cycle (compare state at step t to all states
  in a trailing window of length 50) and average the output over one
  full period. If no cycle is detected within the window, average the
  output over the last 50 steps.
- **Asynchronous update:** for each "step," update all nodes in a
  random permutation (one random permutation per step, seeded
  deterministically from the global seed + step number). Run 10
  independent async trajectories per initial state and average.
- **Total runs:** 9 models x 2 ablation types x 2 update schemes = 36
  (plus input-config variants for G1, G2, G4, G5)

Store results in `results/grn_v2/` within the epistasis-bench repo.
File naming: `{model}_{ablation}_{update}_coalition.npz`

### Step 2: Quality gate (three-part)

Before computing Walsh decompositions, apply all three gate criteria
to each coalition table. ALL THREE must pass; failure on any one
excludes the model from spectrum analysis.

1. Compute pooled-mean value function v(S) = mean over initial states
2. Apply gate criteria:

**Criterion 1 — Unique values >= 10.** Count unique values in
v(S) rounded to 6 decimal places. Models with fewer than 10
produce a value function too coarse for meaningful Walsh analysis.
(G3 was dropped at 9 unique; blood_stem_cell/GATA1 has 5.)

**Criterion 2 — Absolute order-3+ energy >= 1e-4.** Sum of squared
Walsh coefficients at order >= 3. Below this floor, order-3+
structure is indistinguishable from quantization noise.

**Criterion 3 — Not bimodal.** No single rounded value may occupy
> 90% of all 2^n coalitions. A bimodal value function (nearly all
coalitions produce 0 or 1) collapses Walsh analysis to a single
dominant coefficient. (emt_switch: 2048 zeros + 2032 ones out of
4096 coalitions.)

**Output-node rule.** Each model uses the output node from its
original publication. If the publication output fails the gate, the
model is dropped — no output-node substitution. The rationale: the
output node is part of the biological model, and swapping it to
improve gate passage is post-hoc selection on the dependent variable.

Models failing the gate are reported with the failing criterion.
Gate failure is itself a finding (the dynamics produce a degenerate
value function for that readout).

### Step 3: Walsh decomposition + interaction extraction

For each model passing the variance gate:
1. Compute normalized WHT coefficients w_T = (1/2^n) * WHT(v)[T]
2. Compute energy spectrum (fraction of variance at each order)
3. Extract k_99 sparsity (number of coefficients capturing 99% of
   energy)
4. Extract rule-derived interaction ground truth from the Boolean
   rules via exact WHT of each node's truth table (analytic, not
   data-derived). This produces continuous pairwise strengths, not
   binary labels.
5. Extract wiring diagram from the Boolean rules (already done during
   sweep -- stored in `*_wiring.json`)

### Step 4: Faithfulness scoring

For each model:
1. Compute global pairwise Walsh coefficients |w_{i,j}| for all
   C(n,2) pairs from the attractor-level value function
2. Compute local pairwise strengths from WHT of each node's truth
   table (primary ground truth)
3. Load wiring diagram edge labels (secondary ground truth)
4. Compute Spearman rho (local vs global pairwise magnitudes) with
   2000-bootstrap 95% CIs
5. Compute AUROC for wiring diagram edges
6. Apply spectrum gate (absolute order-3+ energy > 1e-4)
7. For models passing the gate, compute three-way motif enrichment
   (P4) against rule-derived three-way interactions
8. Score composition: compare local rule energy spectrum to global
   attractor energy spectrum. Report creation rate (global structure
   absent locally) and destruction rate (local structure absent
   globally).

### Step 5: Bootstrap confidence intervals

For each statistic (AUROC, precision@k, energy spectrum fractions):
1. Resample the initial states with replacement (512 draws from 512)
2. Recompute the pooled-mean value function from the resampled set
3. Recompute the statistic
4. Repeat 2000 times
5. Report the 2.5th and 97.5th percentiles as the 95% CI

This quantifies uncertainty due to initial-state sampling, which is
the only source of randomness in the protocol (the Boolean dynamics
are deterministic given initial states and update order).

### Step 6: Multiple comparison correction

Apply Benjamini-Yekutieli correction across all hypothesis tests
within each prediction family:
- P3 (faithfulness Spearman): 9 models x 2 ablation types = 18 tests
- P4 (motif enrichment): up to 9 models x 2 ablation types = 18 tests
- P5 (sync vs async): 9 models x 1 comparison = 9 tests
- P6 (ablation asymmetry): 36 pairwise model comparisons

BY correction is used rather than BH because the tests within each
family are positively dependent (they share the same underlying
Boolean models). Report both uncorrected and BY-corrected p-values.

### Step 7: Report

Table: G1/G2/G4/G5-G10 x ablation x update x order-1/order-2/order-3+/
Spearman-rho/AUROC-wiring/composition-creation/composition-destruction

Compare confirmatory results (P3, P4) to predictions.
Label exploratory results (P1, P2, P5, P6, P7) as such.


## Protocol fixes (not yet implemented)

The following fixes to `grn_coalition_sweep.py` are required before
running the registered experiment:

### Fix 1: Limit cycle detection and averaging

**Current behavior:** `simulate_to_attractor()` runs synchronous
dynamics for max_steps and returns whatever state the system is in.
If the trajectory is in a limit cycle, the readout is at an arbitrary
phase.

**Required fix:** After max_steps, check whether the final state
matches any state in the last 50 steps (cycle detection). If a cycle
is found, average the output node over one full period. If no cycle
is found, average the output node over the last 50 steps. This
eliminates phase-dependent noise in oscillating coalitions.

This fix is critical for models with tight negative feedback loops
that produce period-2 oscillations under synchronous update.

### Fix 2: Asynchronous update mode

**Current behavior:** Only synchronous update is implemented.

**Required fix:** Add an async update mode where nodes are updated
one at a time in a random permutation per step. G2 (tournier_apoptosis)
was published as an asynchronous model. Synchronous update may produce
spurious cyclic attractors. Run both update modes for all 9 models
and flag disagreements.

### Fix 3: Clamp-to-1 ablation mode

**Current behavior:** Only clamp-to-0 (knockout) is implemented.

**Required fix:** Add a `clamp_value` parameter (0 or 1). For
clamp-to-1, nodes not in the coalition are fixed to 1 at every step
instead of 0. This is the gain-of-function complement of knockout.

### Fix 4: Variance floor gate

**Current behavior:** No gate. Degenerate value functions (few unique
values, near-zero variance) produce numerically unstable energy
spectra.

**Required fix:** Before Walsh decomposition, check: (a) at least 20
unique values in the pooled-mean value function, (b) total variance
above 1e-6. Models failing either criterion are excluded from spectrum
analysis with an explicit flag.

### Fix 5: Rules-based ground truth extraction (IMPLEMENTED)

Replaced symbolic AND-clause parsing with exact WHT of each node's
truth table. `extract_rule_fourier()` in `grn_coalition_sweep.py`
computes, for each gene, the full Walsh spectrum of its truth table
and reports all pairwise and triple interactions with continuous
coefficients. This is representation-independent (does not depend on
DNF/CNF form). Output is stored in `*_wiring.json` under the
`rule_fourier` key.


## Confound controls

1. **Self-regulation.** CycD (G1), TNF (G2), Start (G4), Ago/CycD/
   Notch (G5) are external inputs (self-referential rules). Self-edges
   are excluded from both ground truths. UbcH10 (G1), p27 (G1), CycA
   (G1), and several nodes in other models have self-loops in their
   update rules. Self-loops contribute order-1 energy (main effect),
   not pairwise interaction.

2. **Input nodes.** External inputs persist indefinitely once set
   because their rule is f(x) = x. Knocking out an input node
   permanently removes its influence, whereas knocking out a regulated
   node removes its rule but other nodes can still produce the same
   downstream effect. The protocol uses `--exclude-inputs` to remove
   input nodes from the player set and clamp them to fixed
   configurations, running the sweep separately per input config.

3. **Attractor multiplicity.** Different initial states may converge
   to different attractors, including limit cycles. The value function
   v(S) averages over initial states, which blends multiple attractors
   into a continuous value. Report the fraction of coalition x
   initial-state combinations that converge to a fixed point vs a
   limit cycle, separately for sync and async update.

4. **Knockout asymmetry.** Clamping to 0 vs 1 produces different
   results because Boolean rules use both positive and negative
   regulation. The protocol runs both ablation types and reports
   discrepancies. This is the zero-vs-mean-ablation axis from the
   transformer setting.

5. **Update order sensitivity.** For asynchronous update, different
   random permutation orders produce different trajectories. The
   protocol averages over 10 independent async trajectories per
   initial state to reduce this variance. Report the standard
   deviation across async replicates.

6. **Functionally equivalent nodes.** G4's Ste9 and Rum1 have
   identical update rules. Functionally equivalent nodes inflate the
   apparent interaction order (a 2-way functional interaction appears
   as a higher-order interaction when measured over copies of the same
   signal). Report the number of functionally equivalent groups per
   model.


## What this pre-registration does NOT cover

- The budget-vs-accuracy experiment (requires method scoring framework)
- Specific detection methods scored against the GRN ground truth
- The Poelwijk 2019 protein landscape
- Transformer circuit regimes (separate pre-registration exists)
- Any cross-domain comparison (deferred to P8)
- Implementation of the protocol fixes (described above, not coded)

Those are separate pre-registrations that depend on this one producing
non-degenerate coalition tables with diverse interaction structure.


## Code freeze checklist

Before computing the SHA-256:

1. [ ] All 5 protocol fixes implemented and tested
2. [ ] grn_coalition_sweep.py runs for all 9 models in both ablation
       modes and both update schemes without error
3. [ ] Limit cycle detection produces correct output on a synthetic
       oscillating system (period-2 test case)
4. [ ] Async update produces different attractors than sync for at
       least one model
5. [ ] Variance floor gate correctly excludes degenerate value
       functions on a synthetic test case
6. [ ] Rules-based ground truth extraction matches hand-computed
       interaction sets for G1 (small enough to verify by hand)
7. [ ] Output npz files pass shape assertion (2^n coalitions)
8. [ ] Wiring JSON contains expected edge counts per model
9. [ ] data_utils.py Walsh functions verified against mobius_wht.py
       on a shared test vector
10. [ ] Bootstrap CI computation produces correct coverage on a
        known distribution

Seeds: initial-state seed 42; async permutation seeds derived
deterministically from 42 + step_number; bootstrap seed 2024.


## References

- Faure et al., Bioinformatics 22.14 (2006): e124-e131
- Faure et al., J Theor Biol 250.2 (2008): 332-343
- Tournier & Chaves, J Theor Biol 260.2 (2009): 196-209
- Davidich & Bornholdt, PLoS ONE 3.2 (2008): e1672
- Krumsiek et al., PLoS ONE 6.7 (2011): e22649
- Bonzanni et al., Bioinformatics 29.13 (2013): i80-i88
- Steinway et al., PLoS Comput Biol 10.8 (2014): e1003762
- Rodriguez et al., BMC Syst Biol 9 (2015): 79
- Ortiz-Gutierrez et al., PLoS Comput Biol 11.9 (2015): e1004486
- PyBoolNet: Klarner et al., Bioinformatics 33.5 (2017): 770-772
- biodivine-lib-param-bn: Benes et al., Bioinformatics 38.14 (2022)
- Benjamini & Yekutieli, Ann Statist 29.4 (2001): 1165-1188
