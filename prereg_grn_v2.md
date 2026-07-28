# Pre-registration: GRN Coalition Sweep v2

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
   This is the published update scheme for G2 (tournier_apoptosis)
   and G3 (saadatpour_guardcell). Under synchronous update, these
   models produce spurious limit cycles that do not appear in
   asynchronous dynamics.

Both update schemes are run for all four models. The PRIMARY result
uses each model's published scheme. Disagreements between sync and
async results are flagged and reported. For G2 and G3, synchronous
results are labeled SECONDARY and any predictions that differ between
schemes are noted.


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


### G3: saadatpour_guardcell (13 nodes)

**Output rule:** `KEV = Ca2`

A simple copy. The entire output function reduces to the dynamics
of Ca2.

**The Ca2 subsystem:**
- Ca2 = !Ca2ATP & CIS (2-input AND with negative self-feedback)
- Ca2ATP = Ca2 (copy: 1-step delay negative feedback)
- CIS = cGMP&cADPR | InsP3 (2-input AND OR'd with a bypass)

Because Ca2ATP = Ca2(t-1), the update becomes:
Ca2(t) = !Ca2(t-1) & CIS(t). When CIS is sustained at 1, Ca2
oscillates with period 2 under synchronous update.

**CIS depends on Ca2 through three redundant chains:**

Route A (to cGMP): Ca2 -> NOS -> NO -> GC -> cGMP (4 intermediate
copy nodes)

Route B (to cADPR): Ca2 -> NOS -> NO -> ADPRc -> cADPR (4
intermediate copy nodes, shares NOS and NO with Route A)

Route C (to InsP3): Ca2 -> PLC -> InsP3 (2 intermediate copy nodes)

CIS = (cGMP & cADPR) | InsP3. Routes A and B are AND'd together
(both required for the first term), then OR'd with Route C.

Since Routes A and B share NOS and NO, the effective requirement for
the first CIS term is: {NOS, NO, GC, ADPRc} all active. This is a
4-node AND chain. Route C requires only PLC.

**Simplified functional structure:**
CIS = (NOS & NO & GC & ADPRc all active) | (PLC active)

**Feedback loops (3):**
1. Ca2 -> Ca2ATP -> Ca2. Negative, period 1. Creates period-2
   oscillation.
2. Ca2 -> NOS -> NO -> GC -> cGMP -> CIS -> Ca2. Positive, delay 5.
   (Also through ADPRc -> cADPR, same delay.)
3. Ca2 -> PLC -> InsP3 -> CIS -> Ca2. Positive, delay 3.

**Critical dynamics issue: limit cycles under synchronous update.**
When CIS is sustained, Ca2 oscillates (Ca2 = !Ca2_prev & 1 toggles
every step). The original publication (Saadatpour et al. 2010) used
asynchronous update, which avoids this oscillation. Under synchronous
update, the readout at step 200 captures an arbitrary phase of the
oscillation, injecting noise into the value function.

**AND-gate analysis:** The dominant AND gate is the implicit 4-node
conjunction {NOS, NO, GC, ADPRc} required for Route A+B. A 4-input
AND function distributes Walsh energy uniformly across all orders
from 0 to 4:

For k variables in an AND gate, the fraction of total AND-gate
energy at order j is C(k,j) / 2^k. For k=4:
- order-0: 1/16 (6.25%)
- order-1: 4/16 (25%)
- order-2: 6/16 (37.5%)
- order-3: 4/16 (25%)
- order-4: 1/16 (6.25%)

This 4-node AND is the primary driver of high-order energy.

**Isolated nodes:** KAP = !Ca2 depends on Ca2 but is not upstream of
the output. Knocking out KAP has no effect on KEV. KAP contributes
only intercept energy. Similarly, Ca2ATP is a simple feedback node
whose knockout changes Ca2 dynamics but does not add pathway
complexity.

**Functionally equivalent chains:** Most intermediate nodes (NOS, NO,
GC, ADPRc, PLC, InsP3, cGMP, cADPR) are simple copies of their
single upstream regulator. This makes many nodes functionally
interchangeable within their chain. Knocking out any node in a chain
has the same downstream effect as knocking out any other node in that
chain.

**Prediction (clamp-to-0):**
- order-0: 5-15%
- order-1: 10-25%
- order-2: 25-35%
- order-3+: 30-55%

**Rationale:** The 4-node AND chain (NOS, NO, GC, ADPRc) distributes
energy across orders 0 through 4. The OR-redundancy between the
AND-chain pathway and the PLC bypass creates additional pairwise
interactions (knocking out PLC alone or any single AND-chain member
alone has limited effect, but knocking out PLC and any AND-chain
member eliminates both routes to CIS). The many copy nodes in chains
create additional multi-way dependencies: a signal must traverse
the entire chain, so removing any link kills the chain, generating
interactions among all chain members. The synchronous limit cycle
issue spreads energy into higher orders because oscillating
coalitions produce value ~0.5 while non-oscillating coalitions
produce 0 or 1, creating a complex step function.

G3 should have the highest order-3+ energy of all four models, driven
by the 4-node AND chain and the extensive chain topology.

**Key difference from G1/G4:** G3 has no "bottleneck AND" at the
output (KEV = Ca2 is a copy). Instead, the entire interaction
structure comes from the upstream feedback architecture. The deep
copy chains (length 4-5) create multi-node dependencies that G1 and
G4 (with shallower, more complex rules) do not have.


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

| Model | order-0 | order-1 | order-2 | order-3+ | Primary driver |
|-------|---------|---------|---------|----------|----------------|
| G1 faure_cellcycle | 25-40% | 40-55% | 15-30% | 2-10% | 2-input AND output, moderate feedback |
| G2 tournier_apoptosis | 10-25% | 25-40% | 25-40% | 10-25% | 3 positive feedback loops, competing pathways |
| G3 saadatpour_guardcell | 5-15% | 10-25% | 25-35% | 30-55% | 4-node AND chain, deep copy chains |
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

**G3:** KEV = Ca2 = !Ca2ATP & CIS. Clamping Ca2ATP to 1 permanently
kills Ca2 (and KEV). The asymmetry is sharp but narrow because most
nodes are in simple copy chains. Clamping any copy node to 1 does
not break the chain in the same way knocking it to 0 does (a node
clamped to 1 still passes signal downstream, since the next node in
a copy chain reads 1). Prediction: clamp-to-1 produces a SIMPLER
spectrum (more order-1, less order-3+) because clamping chain nodes
to 1 does not break the chain, eliminating the chain-AND structure
that drives G3's high-order energy under clamp-to-0.

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

**Interaction extraction protocol.** For each gene k and each AND
clause in k's Boolean rule, record the set of regulators that appear
in that clause. Every pair within the set is a functional pairwise
interaction for gene k. Every triple is a functional 3-way
interaction. The union over all clauses and all genes gives the full
ground-truth interaction set.

Example: IKKa = TNF & !C3a & !A20a. One clause with three regulators:
{TNF, C3a, A20a}. This generates 3 pairwise interactions and 1
three-way interaction.

The faithfulness test scores Walsh pairwise (and higher-order)
coefficients against these rule-derived interaction sets.

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

The four models span at least two distinct interaction regimes:
- At least one model with >50% order-1 energy (near-additive)
- At least one model with <30% order-1 energy (interactive)

**Falsification:** All four models have order-1 energy within 10
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
| G3 saadatpour_guardcell | 10-25% | 30-55% | 4-node AND chain, deep copy topology |
| G4 davidich_yeast | 40-55% | 3-12% | 3-input AND with redundant inhibitors |

**Falsification criteria:**
- G2 order-3+ < 5% would mean multiple positive feedback loops do
  not generate detectable higher-order epistasis in Boolean models,
  falsifying the "feedback drives higher-order" hypothesis.
- G3 order-1 > 40% would mean deep copy chains do not generate the
  multi-node AND dependencies predicted, falsifying the "chain depth
  drives order" hypothesis.
- G3 order-3+ < 15% would mean the 4-node AND chain does not
  dominate the spectrum as predicted.
- G4 and G1 differ by > 25 percentage points on order-1, despite
  modeling the same biological process (CDK-driven cell cycle entry)
  in different organisms.


### P3: Faithfulness to Boolean rules (CONFIRMATORY)

For each model, compute the AUROC for recovering rule-derived
pairwise interactions from Walsh pairwise coefficients:
- For each unordered pair (i, j), score = |w_{i,j}| (absolute Walsh
  pairwise coefficient from the pooled-mean value function)
- Label = 1 if genes i and j co-occur in any AND clause of any gene's
  Boolean rule (functional interaction), else 0
- Compute AUROC over all C(n, 2) pairs

**Prediction:** AUROC > 0.65 for at least 3 of 4 models. Walsh
pairwise coefficients should recover rule-based functional
interactions better than chance (0.50) because co-occurring AND
terms create multiplicative contributions to the value function,
which the Walsh pairwise coefficient directly measures.

Rules-based AUROC should exceed wiring-based AUROC for at least 3 of
4 models, because rules capture functional dependencies while the
edge list includes non-interacting regulatory connections (genes in
different OR clauses of the same rule).

**Falsification:** AUROC < 0.55 for all 4 models. This would mean
Walsh pairwise coefficients do not correspond to functional
dependencies in the Boolean rules. This is a genuine negative finding
(functional dependence in update rules does not predict
knockout-derived epistasis).

**Nuance on edge density and base rate.** The number of
rule-derived functional pairwise interactions varies by model. For
models with many AND clauses (G2, G4), the positive label rate is
higher, which compresses AUROC toward 0.50. Report the positive
label rate alongside AUROC, and also report precision@k for k = number
of true interactions (the biologically meaningful operating point).


### P4: Motif enrichment in higher-order coefficients (CONFIRMATORY)

For each model with detectable order-3+ energy (>0.5%):
- Take the top-10 three-way Walsh coefficients by magnitude
- For each triple (i, j, k), check whether the triple co-occurs in
  any single AND clause of any gene's Boolean rule (rule-derived
  three-way interaction)
- Compare the fraction of rule-matched triples to a permutation null
  (1000 random draws of 10 triples from the same order-3 set)

**Prediction:** Rule-matched triples are enriched in the top-10 vs
null (one-sided p < 0.05 after Benjamini-Yekutieli correction) for
models with order-3+ > 5%.

**Falsification:** No enrichment in any model. This would mean
higher-order Walsh structure does not correspond to higher-order
functional dependencies in the Boolean rules.


### P5: Sync vs async disagreement (EXPLORATORY)

For G2 and G3 (published asynchronous models):
- Compute the energy spectrum under both synchronous and asynchronous
  update
- Report the L1 distance between the two spectra

**Prediction:** For G3, sync and async spectra differ by > 10
percentage points in at least one order band. The Ca2 <-> Ca2ATP
oscillation under synchronous update produces period-2 limit cycles
that do not occur under asynchronous update. This changes which
coalitions produce v(S) near 0.5 vs exactly 0 or 1.

For G2, sync and async spectra differ by < 5 percentage points. The
positive feedback loops in G2 drive convergence to fixed-point
attractors (survival or apoptosis), making the system less sensitive
to update order.


### P6: Clamp-to-0 vs clamp-to-1 asymmetry (EXPLORATORY)

For each model, compute the L1 distance between the clamp-to-0 and
clamp-to-1 energy spectra.

**Prediction ranking (most to least asymmetry):** G2 > G1 > G4 > G3.

- G2 has the largest asymmetry because it has both activator and
  repressor pathways to the output. Clamping a repressor ON vs OFF
  has opposite effects.
- G1 and G4 have moderate asymmetry because their output rules use
  only negated inputs (both inhibitors).
- G3 has the smallest asymmetry because most nodes are in copy
  chains; clamping a copy node to 1 still passes signal downstream,
  unlike clamping to 0 which breaks the chain. However, the few
  critical nodes (Ca2ATP, CIS) may create sharp local asymmetries.

**Falsification:** G3 has the LARGEST asymmetry. This would falsify
the prediction that copy-chain topology reduces ablation-type
sensitivity.


### P7: Wiring diagram faithfulness (secondary ground truth, EXPLORATORY)

AUROC for recovering true regulatory edges from Walsh pairwise
coefficients, using the wiring diagram (edge list) as ground truth:
- For each unordered pair (i, j), score = |w_{i,j}|
- Label = 1 if the edge list contains edge i->j OR j->i, else 0
- Compute AUROC over all C(n, 2) pairs

**Prediction:** AUROC > 0.55 for at least 2 of 4 models. Edge-based
AUROC should be lower than rules-based AUROC (P3) because the edge
list includes non-interacting regulatory connections.


### P8: Cross-regime rank stability (deferred)

Deferred until all regimes (transformer, protein, GRN) have coalition
tables and at least 4 methods are scored on all of them. Prediction:
method rankings correlate at Spearman rho > 0.5 between GRN and
transformer circuit regimes.


## Experimental protocol

### Step 1: Generate coalition tables

Run `grn_coalition_sweep.py` for all 4 models, in both ablation modes
and both update schemes:

- **Ablation modes:** clamp-to-0 (knockout), clamp-to-1 (constitutive
  activation)
- **Update schemes:** synchronous, asynchronous (random-order
  sequential)
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
- **Total runs:** 4 models x 2 ablation types x 2 update schemes = 16

Store results in `results/grn_v2/` within the epistasis-bench repo.
File naming: `{model}_{ablation}_{update}_coalition.npz`

### Step 2: Variance floor gate

Before computing Walsh decompositions, apply a variance floor gate
to each coalition table:

1. Compute pooled-mean value function v(S) = mean over initial states
2. Count unique values: |{v(S) : S in 2^n}| rounded to 6 decimal
   places
3. Compute total variance: Var(v(S)) over all coalitions

**Exclusion criteria:**
- Fewer than 20 unique values across all coalitions
- Total variance below 1e-6

Models failing the gate are excluded from spectrum analysis. Report
which models pass/fail and under which update scheme and ablation
type. Failure is itself a finding (the dynamics produce a degenerate
value function).

### Step 3: Walsh decomposition + interaction extraction

For each model passing the variance gate:
1. Compute normalized WHT coefficients w_T = (1/2^n) * WHT(v)[T]
2. Compute energy spectrum (fraction of variance at each order)
3. Extract k_99 sparsity (number of coefficients capturing 99% of
   energy)
4. Extract rule-derived interaction ground truth from the Boolean
   rules (symbolic, not data-derived)
5. Extract wiring diagram from the Boolean rules (already done during
   sweep -- stored in `*_wiring.json`)

### Step 4: Faithfulness scoring

For each model:
1. Compute pairwise Walsh coefficients |w_{i,j}| for all C(n,2) pairs
2. Load rule-derived interaction labels (primary ground truth)
3. Load wiring diagram edge labels (secondary ground truth)
4. Compute AUROC, precision@k, recall@k for both ground truths
5. Compute three-way motif enrichment (P4) against rule-derived
   three-way interactions
6. Compute bootstrap 95% CIs over initial states (see below)

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
- P3 (faithfulness AUROC): 4 models x 2 ablation types = 8 tests
- P4 (motif enrichment): up to 4 models x 2 ablation types = 8 tests
- P5 (sync vs async): 4 models x 1 comparison = 4 tests
- P6 (ablation asymmetry): 6 pairwise model comparisons

BY correction is used rather than BH because the tests within each
family are positively dependent (they share the same underlying
Boolean models). Report both uncorrected and BY-corrected p-values.

### Step 7: Report

Table: G1/G2/G3/G4 x ablation x update x order-1/order-2/order-3+/
AUROC-rules/AUROC-wiring/precision@k

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

This fix is most critical for G3, where Ca2 oscillates with period 2
whenever CIS is sustained.

### Fix 2: Asynchronous update mode

**Current behavior:** Only synchronous update is implemented.

**Required fix:** Add an async update mode where nodes are updated
one at a time in a random permutation per step. G2 (tournier_apoptosis)
and G3 (saadatpour_guardcell) were published as asynchronous models.
Synchronous update produces spurious cyclic attractors in these models.
Run both update modes for all 4 models and flag disagreements.

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

### Fix 5: Rules-based ground truth extraction

**Current behavior:** Only the wiring diagram (edge list) is
extracted via `extract_interaction_graph()`.

**Required fix:** Add a function that extracts functional interactions
symbolically from the Boolean rules. For each gene k and each AND
clause in k's rule, record the set of co-occurring regulators. Every
pair within the set is a pairwise functional interaction. Every triple
is a 3-way interaction. This is the primary ground truth for
faithfulness scoring.


## Confound controls

1. **Self-regulation.** CycD (G1), TNF (G2), Start (G4) are external
   inputs (self-referential rules). Self-edges are excluded from both
   ground truths. UbcH10 (G1), p27 (G1), CycA (G1), and Ca2ATP (G3)
   have self-loops in their update rules. Self-loops contribute
   order-1 energy (main effect), not pairwise interaction.

2. **Input nodes.** CycD, TNF, and Start persist indefinitely once
   set because their rule is f(x) = x. This creates a different kind
   of "knockout" than nodes whose state depends on the network:
   knocking out an input node permanently removes its influence,
   whereas knocking out a regulated node removes its rule but other
   nodes can still produce the same downstream effect. Report results
   both with and without input nodes in the player set.

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

6. **Functionally equivalent nodes.** G3 has many copy nodes (NOS,
   NO, GC, ADPRc are all functionally equivalent to Ca2 at
   equilibrium). G4's Ste9 and Rum1 are identical. Functionally
   equivalent nodes inflate the apparent interaction order (a 2-way
   functional interaction appears as a 4-way interaction when
   measured over 4 copies of the same signal). Report the number of
   functionally equivalent groups per model.


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
2. [ ] grn_coalition_sweep.py runs for all 4 models in both ablation
       modes and both update schemes without error
3. [ ] Limit cycle detection produces correct output on a synthetic
       oscillating system (period-2 test case)
4. [ ] Async update produces different attractors than sync for G3
       (the Ca2 oscillation should disappear under async)
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
- Tournier & Chaves, J Theor Biol 260.2 (2009): 196-209
- Saadatpour et al., J Theor Biol 266.4 (2010): 641-656
- Davidich & Bornholdt, PLoS ONE 3.2 (2008): e1672
- PyBoolNet: Klarner et al., Bioinformatics 33.5 (2017): 770-772
- Benjamini & Yekutieli, Ann Statist 29.4 (2001): 1165-1188
