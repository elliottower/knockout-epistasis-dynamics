# SPEC v2: revising the composition gap paper

**Absorbs `SPEC_measured_extension_v1.md`.** That file remains valid on the yeast
data specifics; this one is the whole plan. Nothing here has been run.

---

## 0. What the paper claims after revision

> **Dynamical composition creates higher-order epistasis that is absent from
> molecular wiring.** Established exactly across 28 Boolean gene regulatory
> networks, preserved under continuous ODE dynamics, absent in size- and
> degree-matched random networks, unpredictable from network topology,
> sign-reversing between loss- and gain-of-function, **surviving graded rather
> than complete perturbation**, and **reproducing in measured combinatorial
> perturbations**. It does *not* occur in protein fitness landscapes, where
> higher-order structure sits on the pairwise contact graph -- so the gap requires
> composition through **time**, not through structure.

One claim, six controls, one boundary condition.

---

## 1. Cut

**The sample-budget arm, and the title claim that came with it.** The v12 title
asserts "recoverable from 1% of knockout combinations"; no recovery experiment
exists in the body, and Faure et al. (*PLOS Comput Biol* 20(5):e1012132, 2024)
already subsampled a complete landscape to exactly 1% and recovered coefficients
at $r = 0.99$.

Revert to the v11 title. The sample-budget work moves to its own repository and is
reframed as **testing whether Faure's result generalises** across 28 systems
spanning $n = 7$--18, rather than claiming it. If recovery degrades with $n$, or
with gap size, or holds only for the near-additive networks, that is a finding
they cannot produce and we can.

---

## 2. Narrow

**C1/C2** -- "higher-order extension has not been quantified" is refuted by
Weinreich et al. 2013 and Sailer & Harms 2017. Replace with:

> Higher-order epistasis has been quantified spectrally in empirical fitness
> landscapes [Weinreich 2013; Sailer & Harms 2017] and Boolean update rules have
> been decomposed independently [Shmulevich & Kauffman 2004; Manicka 2023], but
> the two spectra have not been compared **within the same system**.

**P1** in `PLAN.md` -- "head-to-head benchmarking is rare" is refuted four ways
(Russ 2022; Herder 2015; Puy 2022; Muschalik 2024). Narrow to order-3+ recovery
against exhaustively enumerated truth under controlled budget.

---

## 3. Add (cheap: citation and framing only)

### 3.1 The boundary condition -- the single most valuable addition

**Brookes, Aghazadeh & Listgarten, *PNAS* 119(1):e2109649118 (2022)** found that in
proteins, higher-order epistasis is **localised on the pairwise contact graph**:
"structure-based GNK models correctly identify many of the important higher-order
epistatic interactions in the corresponding empirical fitness functions, despite
using only pairwise structural contact information."

So in proteins the wiring *does* predict the higher-order structure. In these 28
GRNs it does not. **That contrast is the paper's strongest framing move**, because
it converts "a gap exists" into "here is when a gap exists". The proposed
mechanism: protein fitness is close to a static energy function over a structure,
while a GRN phenotype is an attractor of a dynamical system. Composition through
time generates structure that composition through space does not.

Add as a Discussion subsection. It costs one paragraph and it pre-empts "so what".

### 3.2 Global epistasis theory as an ally

A referee will say global epistasis explains the higher-order terms away. The
primary theory says the opposite:

- **Reddy & Desai, *eLife* 10:e64740 (2021):** "diminishing-returns and
  increasing-costs epistasis emerge generically as a consequence of pervasive
  microscopic epistasis." Their derivation *requires* combinatorially many nonzero
  higher-order terms.
- **Bakerlee et al., *Science* 376:630 (2022):** all combinations of 10 mutations,
  six environments -- global trends "emerge from specific idiosyncratic
  interactions."
- **Lyons, Zou, Xu & Zhang, *Nat Ecol Evol* 4:1685 (2020):** same conclusion
  independently.

Cite all three where global epistasis is first mentioned.

### 3.3 An explicit additive baseline

Mandatory since the MULTI-evolve dispute (bioRxiv 10.64898/2026.04.23.719915
against *Science* 10.1126/science.aea1820): any claim of interaction structure
must beat explicit additivity. Likely already present -- make it visible and
labelled.

### 3.4 The cluster-expansion bridge

Effective cluster interactions **are** Walsh coefficients. Sanchez, *Phys Rev B*
81:224202 (2010): "the cluster basis... is a multidimensional discrete Fourier
transform." Barroso-Luque & Ceder, *npj Comput Mater* 10:158 (2024): "we identify
the cluster decomposition as an invariant ANOVA decomposition."

**Nobody in the alloy literature writes "Walsh--Hadamard"; nobody in biology
writes "cluster expansion."** A citation sweep of 100 papers citing Sanchez 2010
returned zero hits for Walsh, Hadamard, Boolean or Möbius.

One paragraph naming the correspondence pre-empts the prior-art attack and is a
genuine cross-domain observation. It also imports five sparse-recovery methods
(ARDR, split-Bregman $\ell_1$, $\ell_0\ell_2$, RFE-OLS, adaptive LASSO) for the
budget work.

### 3.5 Foreground the sign reversal

Creation dominance is specific to loss-of-function and **reverses under
gain-of-function clamping**. That is the strongest result in the paper and it is
currently buried under the magnitude numbers.

It is *directional*, so no measurement artifact reproduces it: a scale
nonlinearity inflates or deflates both arms, it does not flip the ordering
between them. Promote it toward the abstract.

---

## 4. Add (experiments)

### 4.1 Graded perturbation -- the continuous arm

**The question.** Every result so far uses a complete knockout, a binary clamp.
Real perturbations are graded: CRISPRi knocks down partially, drugs inhibit
partially. Is the composition gap an artifact of complete removal?

**Why it matters beyond realism.** The project's own measurement is that
zero-ablation inflates pairwise epistasis **10--50x** relative to mean ablation.
A binary clamp is the most extreme perturbation available, and the most extreme
regime is where an artifact would live.

**Design.** In the Hill-function ODE conversions already built, set each gene's
activity to a fraction $f$ of wild-type rather than to zero. Sweep
$f \in \{0, 0.25, 0.5, 0.75, 1\}$. The decomposition is no longer binary, so use
**Sobol indices** -- the continuous analogue, computable to 4th order with
`sensobol` (Puy et al., *JSS* 102(5), 2022). The Walsh spectrum is the special
case at $f = 0$.

**Predictions.**
- **H-C1:** the gap sign is preserved at every $f > 0$. If it flips at partial
  knockdown, the result is specific to complete removal and must be restated as
  such.
- **H-C2:** gap magnitude increases monotonically as $f \to 0$, consistent with
  the zero-versus-mean ablation measurement. A non-monotone dose-response would
  indicate something other than perturbation strength is driving it.
- **H-C3:** the loss-of-function / gain-of-function reversal (§3.5) interpolates
  smoothly through $f$, i.e. the two regimes are ends of one axis rather than
  distinct phenomena.

**Cost:** low. The ODE models exist; only the perturbation operator changes.

### 4.2 Measured combinatorial perturbations -- the real-world arm

**The objection it closes**, verbatim from the current abstract: the experimental
arm "tests an additive baseline rather than the local-rule comparison of the
Boolean analysis." Both sides of the Boolean comparison derive from one curated
network, so a bench biologist can ask whether the gap is biology or modelling.

**The requirement, so a weaker dataset is not accepted.** Two things from the same
organism: (a) measured combinatorial perturbations at order $\geq 3$, and (b) a
molecular wiring model built **without** using (a). A dataset satisfying only (a)
reproduces the existing arm and closes nothing.

**Candidates** (all *unverified*, counts and access routes must be checked):

| source | supplies | note |
|---|---|---|
| Kuzmin et al., *Science* 2018 | trigenic yeast interactions -- **order 3 measured** | the critical one; local density unknown |
| Costanzo et al., *Science* 2016 | digenic at very large scale | pairwise only |
| **BioGRID, physical evidence codes only** | the wiring side | **must exclude genetic-interaction evidence or the comparison is circular** -- the easiest way to invalidate this by accident |
| Norman et al., *Science* 2019 Perturb-seq | double perturbations, human, **CRISPRa** | the gain-of-function arm, and graded |

**Predictions.**
- **H-R1:** the gap reproduces in the same direction (creation), $|\Delta_{3+}|
  \geq 0.5$ pp.
- **H-R2:** measured gaps fall inside the Boolean range ($<1$ to $+83$ pp). A
  measured gap far outside it suggests a measurement-scale artifact.
- **H-R3, the sharpest:** the loss/gain-of-function **sign reversal reproduces
  across modalities** -- deletion and CRISPRi on one side, CRISPRa on the other.
  Directional, so artifact-resistant. **If only one prediction here is run, run
  this one.**

**Feasibility gate, before touching any biological data.** Real datasets are
combinatorially incomplete. Compute each candidate's effective coverage as a
fraction of its $2^n$, look it up on the sample-budget recovery curves, and
proceed only where recovery is adequate. If the curves say recovery at Kuzmin's
triple coverage is near zero, that is a two-hour finding instead of a two-month
one. **This is also what keeps the sample-budget work useful after it leaves the
paper: it becomes the feasibility argument for the biology.**

**Traps** (detail in `SPEC_measured_extension_v1.md` §5): global versus specific
epistasis and the scale problem; ablation-type mismatch, since deletion is closer
to zero-ablation and CRISPRi to attenuation; coverage bias, since screens
over-sample well-studied genes.

---

## 5. Order of work

1. **Title revert and the two narrowings** (§1, §2). Minutes. Do first -- the
   current title claims a result the paper does not contain.
2. **The four framing additions** (§3). Citations and paragraphs, no compute.
3. **Graded perturbation** (§4.1). ODE models exist; cheapest experiment.
4. **Feasibility gate** for the measured arm (§4.2). No new data.
5. **Measured perturbations** (§4.2), H-R3 first.
6. Pre-register 3 and 5 before computing any spectrum, in the style of
   `prereg_grn_v4.md`, SHA stamped.

---

## 6. Manual checks that could change the paper

1. **Do Sailer & Harms' (2017) seven combinatorially complete maps overlap the
   empirical landscapes in Table 2?** They report 2.2--31.0% higher-order energy,
   mean 12.7%. If they overlap, that section partly reproduces a 2017 *Genetics*
   result **on the same data** and must be reframed as replication. Highest
   priority.
2. Can BioGRID evidence codes cleanly separate physical from genetic interaction?
   If not, the real-world arm has no valid wiring side.
3. Does Kuzmin's trigenic set have any locally dense region -- a small gene set
   with all pairs *and* all triples measured -- where an exact low-$n$
   decomposition needs no recovery method at all? Even $n = 6$ complete is worth
   more than a large sparse sample.

---

## 7. Open questions

1. Does the graded arm belong in this paper or does it dilute it? It is a
   different perturbation regime, and a referee may read it as a second paper.
2. If the measured arm fails the feasibility gate, does the paper ship without it
   and concede the modelling objection explicitly, or wait?
3. Sobol at 4th order versus Walsh at all orders -- is the truncation acceptable,
   and does it make the continuous and binary results incommensurable?

---

## 8. What each arm tells us, regardless of publication

The §7 questions are publication questions. These are the scientific ones, and
they are the reason to run the arms even if some never appear in a manuscript.
Ordered by information per unit cost.

### 8.1 Is the gap a property of composition, or of complete removal? (cheap)

The graded arm answers this directly and it is the most basic mechanistic
question the project has not asked. A complete knockout is the most violent
perturbation available; the project's own measurement says zero-ablation inflates
pairwise epistasis 10--50x over mean ablation. If the gap dissolves at $f = 0.5$,
then "dynamics create higher-order epistasis" is really "removing a node entirely
creates higher-order epistasis", which is a much weaker and more specific claim --
and one worth knowing is true.

The **dose-response shape** is the informative part, not the endpoint. Monotone
in $f$ means perturbation strength drives it. A threshold means the system has a
regime boundary. Non-monotone means something else is going on and the current
account is incomplete.

### 8.2 Does the gap actually require dynamics? (cheap, and currently only asserted)

§3.1 cites Brookes et al. to argue proteins show no gap because protein fitness
is a static energy function over a structure. **That is a cross-literature
comparison, not a test**, and the two systems differ in a dozen ways besides
dynamics.

The controlled version is available and cheap: **build matched static and
dynamic systems over the same graph.** Take each of the 28 networks' wiring and
define (a) a static energy-function phenotype -- an Ising-style sum over the same
edges -- and (b) the existing attractor phenotype. Same topology, same nodes, same
perturbations. Compute both spectra.

If the static version shows no composition gap and the dynamic one does, the
mechanism is established **within one system** rather than inferred across two
literatures. That is a much stronger claim than §3.1 currently makes, and it is
the same cost as one more sweep.

**This is the experiment I would most want to run**, and it is the one that was
not in any earlier plan.

### 8.3 Is the gap biology or modelling? (expensive, fundamental)

The measured arm. Nothing else settles it, and no amount of additional simulation
substitutes. Worth doing even if the answer arrives after the paper ships, because
a negative here would mean the Boolean result describes a modelling formalism
rather than a biological phenomenon -- and we would want to know that.

### 8.4 Does Faure's 1% recovery generalise? (cheap, practical)

Their result is on one simulated landscape. Twenty-eight exact systems spanning
$n = 7$--18 and three orders of magnitude of complexity say whether it is a
general property of sparse landscapes or a property of their simulation. The
useful outcome is the **failure mode**: which systems break it, and whether the
breaking correlates with gap size, $n$, or order-3+ energy.

If recovery degrades exactly where the composition gap is large, that is a
connection between the two projects that neither predicted.

### 8.5 Where does the LoF/GoF reversal come from? (cheap)

The reversal is the paper's most artifact-resistant result and the project has no
mechanism for it. The graded arm gives one for free: sweeping $f$ from full
knockdown through wild-type to over-expression traces the whole axis. If the sign
flips at $f = 1$ exactly, that is a statement about the reference point rather
than about biology; if it flips elsewhere, the location is informative.

### 8.6 Reading order for the science, not the paper

1. graded sweep (§8.1, §8.5) -- one perturbation operator, models already exist
2. matched static-versus-dynamic (§8.2) -- the mechanism test
3. Faure generalisation (§8.4) -- reuses existing coalition tables
4. feasibility gate, then measured perturbations (§8.3)

The first three need no new data and no new systems. They are sweeps over assets
already on disk.
