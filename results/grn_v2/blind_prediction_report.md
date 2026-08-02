# Structurally-Blinded Prediction Report: Boolean GRN Composition Experiment

**Date**: 2026-07-28
**Method**: Structural analysis of Boolean rules only (topology, feedback loops,
AND-gate depth, pathway structure). No prior experimental results consulted.
Predictions written and saved to `blind_predictions.json` before any
experiment was run.

---

## 1. Per-Model Structural Analysis

### 1.1 Faure Cell Cycle (G1)

**Network**: 10 nodes, 1024 coalitions, output = CycB
**Output rule**: `!cdh1 & !Cdc20` (2 regulators, pure AND of negations)
**Input nodes**: CycD (self-referential)
**Feedback**: 25 loops <= length 4 (13 positive, 12 negative), fb_ratio = 0.52
**AND-gate arity**: 4 (CycA rule)
**Pathway depth**: 8

**Key structural features**:
- CycB-Cdc20 negative feedback (length 2): CycB activates Cdc20, Cdc20 inhibits CycB
- CycB-cdh1 positive feedback (length 2): double negative creates bistability
- Multiple length-4 negative loops through E2F, CycA, UbcH10
- CycA rule has arity-4 AND gates creating strong local interactions
- Deep competing pathways from CycD input through p27 and Rb branches
- Local order-3+ energy: 3.8%

**Predictions**:
- Rho direction: positive [0.2, 0.5]
- Creation/destruction: CREATION (moderate confidence)
- Global order-3+: 3-10%
- Reasoning: Balanced feedback with bistability from CycB-cdh1 double negative. CycA arity-4
  AND gates create local higher-order interactions that compound through pathways.

### 1.2 Tournier Apoptosis (G2)

**Network**: 12 nodes, 4096 coalitions, output = C3a
**Output rule**: `!IAP & C8a` (2 regulators)
**Input nodes**: TNF
**Feedback**: 6 loops <= length 4 (3 positive, 3 negative), fb_ratio = 0.50
**AND-gate arity**: 3
**Pathway depth**: 9

**Key structural features**:
- ALL 3 positive loops directly involve the output C3a (C3a-C8a, C3a-IAP, C3a-CARP-C8a)
- ALL 3 negative loops are in the NF-kB subsystem, separate from output
- Bistable apoptosis switch: once C3a activates, positive feedback locks it on
- Competing pro-apoptotic (TNF->T2->C8a) vs anti-apoptotic (TNF->NFkB->IAP/CARP) pathways
- NFkB = NFkBnuc = !IkB (functional near-redundancy)
- Local order-3+ energy: 4.2%

**Predictions**:
- Rho direction: positive [0.2, 0.5]
- Creation/destruction: CREATION (high confidence)
- Global order-3+: 5-15%
- Reasoning: Although overall fb_ratio is 0.50, ALL positive loops involve the output while
  ALL negative loops are in a separate subsystem. The output-relevant feedback is 100%
  positive. This creates a bistable switch requiring multiple simultaneous conditions to flip.

### 1.3 Davidich Yeast (G4)

**Network**: 10 nodes, 1024 coalitions, output = Cdc2_Cdc13
**Output rule**: `!Ste9 & !Rum1 & (!Slp1 | Cdc2_Cdc13_A)` (4 regulators)
**Input nodes**: Start
**Feedback**: 16 loops <= length 4 (13 positive, 3 negative), fb_ratio = 0.81
**AND-gate arity**: 3
**Pathway depth**: 5

**Key structural features**:
- Overwhelmingly positive feedback (13+ vs 3-): strongest positive bias of all models
- Ste9 and Rum1 have IDENTICAL update rules (exact functional redundancy)
- Start and SK have identical rules (SK = Start)
- Delayed negative feedback via Slp1 pathway (length 3) creates oscillatory potential
- Shallowest network (depth 5)
- Local order-3+ energy: 2.5%

**Predictions**:
- Rho direction: positive [0.3, 0.6]
- Creation/destruction: CREATION (high confidence)
- Global order-3+: 3-10%
- Reasoning: Overwhelmingly positive feedback creates strong bistability. Ste9/Rum1
  redundancy creates strong pairwise interaction. Expected highest rho.

### 1.4 Drosophila Cell Cycle (G5)

**Network**: 14 nodes (11 dynamic + 3 inputs), 16384 coalitions, output = CycB
**Output rule**: `!Fzy & !Fzr & (!Wee1 | Stg)` (4 regulators)
**Input nodes**: Ago, CycD, Notch
**Feedback**: 36 loops <= length 4 (22 positive, 14 negative), fb_ratio = 0.61
**AND-gate arity**: 5 (CycA rule)
**Pathway depth**: 9

**Key structural features**:
- 3 input nodes create independent pathway origins
- Multiple positive feedback around CycB: CycB-Wee1, CycB-Stg (both length 2)
- Negative feedback: CycB-Fzy (length 2), CycB-E2F-CycA-Fzr (length 4)
- CycA has arity-5 AND gate creating strong local interactions
- Largest number of feedback loops (36 total)
- Competing pathways from 3 inputs through different intermediaries to CycB
- Local order-3+ energy: 4.2%

**Predictions**:
- Rho direction: positive [0.1, 0.4]
- Creation/destruction: CREATION (moderate confidence)
- Global order-3+: 5-15%
- Reasoning: Positive feedback dominant (22+/14-). Three independent inputs create
  competing routes. But 91 pairs dilute individual strengths, weakening rho.

### 1.5 Fanconi Anemia (G9)

**Network**: 15 nodes, 32768 coalitions, output = CHKREC
**Output rule**: `(!DSB & (TLS|HRR2|FAHRR|NHEJ)) | !(ICL|NHEJ|FAHRR|ADD|CHKREC|DSB|HRR2|TLS)` (8 regulators)
**Input nodes**: none
**Feedback**: 98 loops <= length 4 (36 positive, 45 negative, 17 mixed), fb_ratio = 0.44
**AND-gate arity**: 5
**Pathway depth**: 14

**Key structural features**:
- NEGATIVE FEEDBACK DOMINANT (45- vs 36+ plus 17 mixed)
- Output rule has 8 regulators with complex AND/OR structure
- Null-state clause: CHKREC=1 when ALL damage response genes are off
- 4 competing repair pathways: TLS, FAHRR, NHEJ, HRR2 (mutually exclusive)
- CHKREC acts as global reset: inhibits nearly all signaling components
- Deepest pathway (14 steps) through DNA damage cascade
- ICL conditional self-loop: ICL = ICL & !DSB (persists until DSB forms)
- Local order-3+ energy: 6.5%

**Predictions**:
- Rho direction: positive [0.0, 0.3]
- Creation/destruction: DESTRUCTION (moderate confidence)
- Global order-3+: 3-10%
- Reasoning: Negative feedback dominant. CHKREC global reset creates binary outcome
  that compresses higher-order structure. Null-state degeneracy flattens value function.
  Despite high local order-3+ from complex output rule, negative feedback linearization
  should destroy this structure.

### 1.6 Arabidopsis Cell Cycle (G10)

**Network**: 14 nodes, 16384 coalitions, output = CYCB1_1
**Output rule**: `!APC_C & (MYB77 | MYB3R1_4 | E2Fb&(!RBR|CYCD3_1&!KRP1)&!E2Fc)` (8 regulators)
**Input nodes**: none
**Feedback**: 81 loops <= length 4 (39 positive, 42 negative), fb_ratio = 0.48
**AND-gate arity**: 6 (E2Fe rule)
**Pathway depth**: 13

**Key structural features**:
- Near-balanced feedback with slight negative dominance (39+/42-)
- Output rule has 8 regulators with deeply nested AND/OR structure
- CYCB1_1-MYB3R1_4 positive feedback (length 2): only positive loop involving output
- Two negative loops from output: via APC_C (length 3) and via KRP1 (length 3)
- E2Fe has enormously complex rule (6 regulators, deeply nested)
- Very dense network: 63 edges for 14 nodes (4.5 edges/node)
- Deep pathways (13) through E2F family cascade
- Local order-3+ energy: 4.0%

**Predictions**:
- Rho direction: positive [0.0, 0.3]
- Creation/destruction: CREATION (low confidence)
- Global order-3+: 3-10%
- Reasoning: Borderline case. Near-balanced feedback. Slight creation expected from
  CYCB1_1-MYB3R1_4 positive feedback and deep compounding pathways.

---

## 2. Cross-Network Hypothesis

**Hypothesis**: The positive feedback loop ratio (fb_ratio = n_positive /
(n_positive + n_negative), ignoring mixed-sign loops) predicts the direction
of the composition gap. Networks with fb_ratio > 0.5 show creation of
higher-order interactions; networks with fb_ratio < 0.5 show destruction.

**Predictor**: fb_ratio (computed from structural topology)
**Outcome**: composition_gap = global_order3plus - local_order3plus
**Expected direction**: Positive correlation

**Predicted ordering** (strongest creation to strongest destruction):
1. davidich_yeast (fb=0.81): strong creation
2. drosophila_cellcycle (fb=0.61): moderate creation
3. faure_cellcycle (fb=0.52): mild creation
4. tournier_apoptosis (fb=0.50): borderline positive
5. arabidopsis_cellcycle (fb=0.48): borderline negative or neutral
6. fanconi_anemia (fb=0.44): destruction

**Caveats stated in advance**:
- Loops directly involving the output may matter more than subsystem loops
- Mixed-sign loops (fanconi has 17) are excluded but may contribute
- Network size confounds
- Tournier's output-relevant feedback is 100% positive despite overall ratio of 0.50

---

## 3. Results Table

| Model | n | rho | 95% CI | p | global o3+ | local o3+ | delta | Direction | cycling% |
|---|---|---|---|---|---|---|---|---|---|
| faure_cellcycle | 10 | 0.409 | [0.129, 0.637] | 0.005 | 0.0082 | 0.0385 | -0.030 | DESTRUCTION | 19.6% |
| tournier_apoptosis | 12 | 0.247 | [0.026, 0.450] | 0.046 | 0.1263 | 0.0417 | +0.085 | CREATION | 7.4% |
| davidich_yeast | 10 | 0.154 | [-0.148, 0.423] | 0.312 | 0.0219 | 0.0256 | -0.004 | DESTRUCTION | 17.7% |
| drosophila_cellcycle | 14 | 0.120 | [-0.089, 0.317] | 0.257 | 0.1042 | 0.0427 | +0.062 | CREATION | 14.0% |
| fanconi_anemia | 15 | 0.305 | [0.118, 0.476] | 0.002 | 0.0103 | 0.0652 | -0.055 | DESTRUCTION | 47.5% |
| arabidopsis_cellcycle | 14 | 0.205 | [-0.015, 0.409] | 0.052 | 0.3828 | 0.0398 | +0.343 | CREATION | 16.4% |

---

## 4. Prediction Accuracy

### 4.1 Spearman rho direction (6/6 correct)

All six models showed positive Spearman rho between local and global pairwise
interaction strengths, as predicted. However, only 3 of 6 reached statistical
significance at alpha=0.05: faure (p=0.005), tournier (p=0.046), fanconi (p=0.002).
Arabidopsis was borderline (p=0.052). Davidich (p=0.312) and drosophila (p=0.257)
were not significant.

| Model | Predicted direction | Actual rho | p-value | Correct? |
|---|---|---|---|---|
| faure_cellcycle | positive | +0.409 | 0.005 | YES |
| tournier_apoptosis | positive | +0.247 | 0.046 | YES |
| davidich_yeast | positive | +0.154 | 0.312 | YES (but n.s.) |
| drosophila_cellcycle | positive | +0.120 | 0.257 | YES (but n.s.) |
| fanconi_anemia | positive | +0.305 | 0.002 | YES |
| arabidopsis_cellcycle | positive | +0.205 | 0.052 | YES (borderline) |

### 4.2 Spearman rho range (4/6 correct)

| Model | Predicted range | Actual rho | In range? |
|---|---|---|---|
| faure_cellcycle | [0.2, 0.5] | 0.409 | YES |
| tournier_apoptosis | [0.2, 0.5] | 0.247 | YES |
| davidich_yeast | [0.3, 0.6] | 0.154 | NO (below) |
| drosophila_cellcycle | [0.1, 0.4] | 0.120 | YES |
| fanconi_anemia | [0.0, 0.3] | 0.305 | BORDERLINE (just above) |
| arabidopsis_cellcycle | [0.0, 0.3] | 0.205 | YES |

Davidich was predicted to have the highest rho due to the Ste9/Rum1 redundancy
creating strong pairwise structure. It actually had the second-lowest rho (0.154,
not significant). Faure had the highest rho (0.409). The redundancy argument was
wrong: identical update rules do not guarantee that local pairwise strengths
predict global ones.

### 4.3 Creation vs. destruction (4/6 correct)

| Model | Predicted | Actual | delta o3+ | Correct? |
|---|---|---|---|---|
| faure_cellcycle | CREATION | DESTRUCTION | -0.030 | NO |
| tournier_apoptosis | CREATION | CREATION | +0.085 | YES |
| davidich_yeast | CREATION | DESTRUCTION | -0.004 | NO |
| drosophila_cellcycle | CREATION | CREATION | +0.062 | YES |
| fanconi_anemia | DESTRUCTION | DESTRUCTION | -0.055 | YES |
| arabidopsis_cellcycle | CREATION | CREATION | +0.343 | YES |

The two failures (faure, davidich) were networks where destruction was mild
(delta = -0.030 and -0.004 respectively). The correct predictions included the
three most extreme outcomes: arabidopsis (+0.343, massive creation), tournier
(+0.085, clear creation), and fanconi (-0.055, clear destruction).

### 4.4 Global order-3+ energy spectrum (2/6 correct)

| Model | Predicted range | Actual global o3+ | Correct? |
|---|---|---|---|
| faure_cellcycle | 3-10% | 0.82% | NO (below) |
| tournier_apoptosis | 5-15% | 12.63% | YES |
| davidich_yeast | 3-10% | 2.19% | NO (below) |
| drosophila_cellcycle | 5-15% | 10.42% | YES |
| fanconi_anemia | 3-10% | 1.03% | NO (below) |
| arabidopsis_cellcycle | 3-10% | 38.28% | NO (enormously above!) |

The energy spectrum predictions were the weakest. The structural analysis
systematically overestimated order-3+ energy for destruction networks
(where dynamics compress to low orders) and massively underestimated it
for arabidopsis (where dynamics amplify higher-order interactions by ~10x).

### 4.5 Cross-network hypothesis (FAILED: 3/6 correct)

| Model | fb_ratio | Predicted gap | Actual gap | Correct? |
|---|---|---|---|---|
| davidich_yeast | 0.81 | creation | -0.004 (destruction) | NO |
| drosophila_cellcycle | 0.61 | creation | +0.062 (creation) | YES |
| faure_cellcycle | 0.52 | creation | -0.030 (destruction) | NO |
| tournier_apoptosis | 0.50 | borderline+ | +0.085 (creation) | YES |
| arabidopsis_cellcycle | 0.48 | borderline- | +0.343 (creation) | NO |
| fanconi_anemia | 0.44 | destruction | -0.055 (destruction) | YES |

The fb_ratio hypothesis does not predict the composition gap direction.
The network with the highest fb_ratio (davidich, 0.81) shows destruction,
while a network below the threshold (arabidopsis, 0.48) shows the strongest
creation. The Spearman rank correlation between fb_ratio and delta_o3+ is
approximately zero.

---

## 5. What the Structural Analysis Got Right and Wrong

### 5.1 What worked

**Rho direction was universally positive.** All six networks showed positive
rank correlation between local and global pairwise interaction strengths.
This was the easiest prediction: local Fourier structure provides a meaningful
signal about global structure regardless of topology, suggesting that
dynamical composition preserves the relative ordering of pairwise interactions
even when it changes their magnitudes.

**Fanconi destruction was correctly predicted.** The reasoning about negative
feedback dominance and the CHKREC global reset was confirmed. Fanconi's cycling
fraction (47.5%) was by far the highest, consistent with the prediction that
complex negative feedback creates limit cycles that linearize the value function.

**Tournier creation was correctly predicted with the right mechanism.** The
observation that all positive feedback loops directly involve the output was
validated. The bistable apoptosis switch creates strong higher-order interactions
in the global value function (order-3+ went from 4.2% local to 12.6% global).
Tournier also had the lowest cycling fraction (7.4%), consistent with the
positive-feedback-locks-to-fixed-point picture.

**Drosophila creation and diluted rho.** Correctly predicted creation and
that the large network (91 pairs) would dilute rho. Actual rho = 0.12 was
in the predicted range but not significant.

### 5.2 What failed

**The fb_ratio hypothesis is wrong.** Positive feedback loop ratio does
not predict creation vs. destruction. The correlation between fb_ratio
and composition gap is essentially zero. This was the central cross-network
hypothesis and it failed.

**Davidich yeast: highest fb_ratio (0.81) but destruction.** The prediction
that overwhelming positive feedback creates strong bistability and therefore
higher-order interactions was wrong. Despite having the most lopsided positive
feedback and identical redundant nodes (Ste9=Rum1), dynamics slightly
destroyed higher-order structure. The cycling fraction (17.7%) may explain
this: the Slp1-mediated delayed negative feedback creates oscillatory dynamics
that smooth out the higher-order structure despite the positive feedback
dominance.

**Faure: mild positive bias (fb=0.52) but destruction.** The prediction
that balanced feedback with slight positive dominance leads to mild creation
was wrong. The cycling fraction (19.6%) again correlates with destruction.

**Arabidopsis: the biggest surprise.** Predicted borderline creation with
low confidence and global order-3+ of 3-10%. Actual: 38.3% global order-3+,
the most dramatic creation in the entire set. Dynamics amplified higher-order
interactions by nearly 10x relative to local rules. The order-0 energy
collapsed from 55% (local) to 10% (global), meaning the mean-field
approximation is nearly useless for this network. This was a spectacular
prediction failure in magnitude, though the direction (creation) was correct.

### 5.3 The real predictor: cycling fraction

An unplanned finding that emerged from the results: **cycling fraction perfectly
separates creation from destruction** in this sample. Ordered by cycling
fraction:

| Model | Cycling % | Direction | delta o3+ |
|---|---|---|---|
| tournier_apoptosis | 7.4% | CREATION | +0.085 |
| drosophila_cellcycle | 14.0% | CREATION | +0.062 |
| arabidopsis_cellcycle | 16.4% | CREATION | +0.343 |
| davidich_yeast | 17.7% | DESTRUCTION | -0.004 |
| faure_cellcycle | 19.6% | DESTRUCTION | -0.030 |
| fanconi_anemia | 47.5% | DESTRUCTION | -0.055 |

The three networks with cycling fraction < 17% all show creation. The three
with cycling fraction >= 17.7% all show destruction. The threshold falls
between 16.4% and 17.7%.

**Mechanistic interpretation**: Limit cycles time-average the output node over
the cycle period. This averaging is a smoothing operation that compresses
higher-order Walsh coefficients toward zero. Fixed-point attractors preserve
the sharp threshold behavior of Boolean logic, which is inherently high-order
(AND gates require all inputs to be active simultaneously). When most initial
conditions converge to fixed points, the global value function inherits the
sharp transitions from the local rules and dynamics amplify them through
feedback. When many initial conditions cycle, the time-averaging destroys
the sharp transitions.

This predictor was not anticipated from structural analysis because cycling
fraction is a dynamical property, not a topological one. It cannot be computed
from the Boolean rules alone without running simulations. The structural
analysis identified feedback loops and AND-gate depth as the relevant features,
but missed that the propensity for limit cycles (a global dynamical property)
dominates.

**Caveat**: This is a post-hoc observation on 6 networks. The threshold at
~17% and the perfect separation could be coincidental. The observation should
be tested on additional Boolean GRN models to assess robustness.

### 5.4 Why arabidopsis was so extreme

Arabidopsis showed order-3+ energy of 38.3% -- 10x its local value and 3x
larger than any other network's. Several factors may explain this:

1. **Dense interconnectivity** (63 edges, 4.5 per node) creates many pathways
   for interactions to compound through dynamics.
2. **No input nodes**: all 14 nodes participate in the dynamics, so there are
   no "anchors" that constrain the system. Every node can be knocked out and
   every combination explored.
3. **Complex nested Boolean logic**: the output rule has 8 regulators with
   deeply nested AND/OR structure, and E2Fe has an arity-6 rule. These create
   a high-dimensional threshold surface that dynamics can amplify.
4. **Moderate cycling fraction** (16.4%): low enough that most initial states
   converge to fixed points, preserving the threshold behavior, but with
   enough variety in the attractor landscape to create structure.

The structural analysis identified all of these features but failed to predict
their combined magnitude. This suggests that interaction creation is a
nonlinear function of structural properties -- the whole is greater than the
sum of the parts.

### 5.5 Summary scorecard

| Prediction category | Score | Notes |
|---|---|---|
| Rho direction | 6/6 | All positive, as predicted |
| Rho range | 4/6 | Davidich below, fanconi borderline above |
| Creation vs destruction | 4/6 | Faure and davidich wrong |
| Global order-3+ range | 2/6 | Systematic bias; arabidopsis 10x over |
| Cross-network hypothesis | 3/6 | fb_ratio does not predict gap direction |
| **Overall** | **19/30** | **63%** |

The structural analysis was best at qualitative predictions (rho direction,
creation/destruction for extreme cases) and worst at quantitative ones
(energy spectrum magnitudes). The cross-network hypothesis based on
feedback loop ratio failed because cycling fraction -- a dynamical rather
than topological property -- is the dominant predictor of whether composition
creates or destroys higher-order interactions.

---

## Appendix: Experimental Details

- **n_init**: 512 random initial states per coalition
- **max_steps**: 200 synchronous update steps
- **seed**: 42
- **update_scheme**: synchronous
- **clamp_value**: 0 (absent genes clamped to 0)
- **Computation**: 6-worker multiprocessing for fanconi_anemia (32768 coalitions)
  and arabidopsis_cellcycle (16384 coalitions)
- **Predictions file**: `results/grn_v2/blind_predictions.json`
  (written before any experiment ran)
- **Structural analysis**: `results/grn_v2/structural_analysis.json`
- **Per-model results**: `results/grn_v2/{model}_composition_blind.json`
- **Full summary**: `results/grn_v2/blind_experiment_summary.json`
