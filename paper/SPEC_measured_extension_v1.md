# SPEC: extending the composition gap to measured combinatorial perturbations

**For Joe. Draft v1 — nothing here has been run, and the dataset facts marked
*unverified* must be checked before anyone commits time.**

---

## 1. The exact objection this closes

The composition gap paper shows, across 28 Boolean gene regulatory networks with
exhaustive $2^n$ knockout sweeps, that the interaction spectrum of the individual
gene rules does not predict the interaction spectrum of the attractor-level
knockout phenotype. Dynamics create higher-order epistasis. Structural predictors
score 48%, indistinguishable from chance. Hill-function ODE conversion preserves
the gap sign in 27 of 28, and the gap is absent in size- and degree-matched random
networks.

Two attacks are already answered by those controls: *"it's a Boolean
idealisation"* (ODE) and *"it's a size or connectivity artifact"* (matched random
networks).

**The attack that remains is that both sides of the comparison are derived from
the same curated network.** The "molecular wiring" and the "knockout phenotype"
are two readings of one model. A bench biologist can ask whether the gap is a
property of biology or of the modelling step, and nothing in the paper settles
it.

The current abstract concedes precisely this. The experimental arm — three
fitness landscapes showing higher-order epistasis — is described as testing "an
additive baseline rather than the local-rule comparison of the Boolean analysis".
That is the honest statement, and it is the hole.

**What closes it:** a system where the composition gap can be computed with
*measured* perturbation phenotypes on one side and an *independently determined*
molecular wiring model on the other.

---

## 2. The requirement, stated so we do not accept a weaker dataset

The local-rule comparison needs **two things from the same organism**:

| | what it supplies | why the paper needs it |
|---|---|---|
| **(a)** measured combinatorial perturbations, order $\geq 3$ | the phenotype-level interaction spectrum | without order-3 measurements there is no higher-order term to compare |
| **(b)** a molecular wiring model built **without** using (a) | the local-rule spectrum | if the wiring was inferred from genetic interactions the comparison is circular |

Datasets satisfying only (a) reproduce the existing experimental arm — higher-order
epistasis exists — and do not close the objection. **Do not spend time on a
dataset that fails (b).**

---

## 3. Candidate systems

### 3.1 Yeast — the only system likely to satisfy both

**(a) Measured combinatorial perturbations.**
- *Costanzo et al., Science 2016*, global genetic interaction network. Digenic, at
  very large scale. **Unverified: exact pair count (~23 million is the figure in
  circulation), coverage, and the current download route.**
- *Kuzmin et al., Science 2018*, systematic analysis of complex genetic
  interactions — **trigenic** interactions, i.e. order-3 measured rather than
  simulated. **Unverified: number of triples, how trigenic scores are defined, and
  whether the triple set is dense enough anywhere to support a local spectrum.**

**(b) Independent wiring.** BioGRID physical interactions, or a curated
regulatory network. **The critical check: BioGRID contains genetic interactions
as well as physical ones. The wiring side must be restricted to physical or
regulatory evidence codes, with genetic-interaction evidence excluded**, or the
comparison is circular in exactly the way §2(b) forbids. This is the single
easiest way to invalidate the whole extension by accident.

### 3.2 Human — Perturb-seq

*Norman et al., Science 2019* — double perturbations with rich single-cell
readout. **Unverified: number of pairs, whether any triples exist, and the
current access route.**

Note this is **CRISPRa (activation)**, which makes it the gain-of-function arm
rather than a replication of the loss-of-function result. See §4, H3.

### 3.3 What to reject

Any dataset that is pairwise-only satisfies (a) at order 2 and cannot address
higher-order creation. Any dataset whose network model was inferred from the
same perturbation data fails (b).

---

## 4. Analyses, and what each would establish

### H1 — the gap reproduces where interactions were measured

Compute the phenotype-level interaction spectrum from measured combinatorial
perturbations, and the local-rule spectrum from the independent wiring model.
Predict $|\Delta_{3+}| \geq 0.5$ pp in the same direction as the Boolean result,
i.e. **creation**: measured phenotypes carry more higher-order spectral energy
than the wiring predicts.

*If it holds:* the composition gap is a property of biology, not of Boolean
modelling, and the paper's central claim survives its strongest attack.
*If it fails:* the gap is a modelling artifact, and that is the finding. It must
be reported, not buried.

### H2 — the magnitude ordering survives

The Boolean result spans $<1$ to $+83$ pp. Predict that measured systems fall
inside that range rather than at either extreme. A measured gap far larger than
any simulated one would suggest a measurement-scale artifact rather than
biological composition (see §5.2).

### H3 — the loss-of-function / gain-of-function asymmetry replicates

The paper reports that **creation dominance is specific to loss-of-function
knockouts, and gain-of-function clamping reverses the ratio.** This maps onto
real perturbation modalities directly:

| paper | measured analogue |
|---|---|
| loss-of-function knockout | gene deletion (Costanzo, Kuzmin), CRISPRi |
| gain-of-function clamping | **CRISPRa** (Norman Perturb-seq) |

**Predict the sign reversal reproduces across modalities.** This is the sharpest
available test because it is a *directional* prediction that no measurement
artifact should reproduce: scale nonlinearity would inflate or deflate both arms,
not flip the ordering between them.

If only one prediction from this spec is run, run this one.

---

## 5. Methodological traps, all of which will be raised

### 5.1 The additive baseline is now mandatory

A live dispute in the protein literature — "Additive baselines furnish no
evidence for epistasis learning" (bioRxiv 10.64898/2026.04.23.719915) against
MULTI-evolve (*Science*, 10.1126/science.aea1820) — has established that **any
claim of interaction structure must beat an explicit additive baseline.**

The existing experimental arm already does this. The extension must keep it and
report it alongside the local-rule comparison, so that the two are not confused.

### 5.2 Global against specific epistasis — the most dangerous confound

Raw non-additivity conflates *specific* epistasis with *global* epistasis, the
latter arising from a nonlinear measurement scale rather than from any
interaction between loci. Fitness assays saturate; a log or logit transform
changes the apparent interaction spectrum.

The Boolean work is immune, because $v(S)$ is exact and the scale is chosen. The
measured extension is **not immune**, and a spectrum computed on the raw scale is
not comparable to one computed on Boolean attractors.

Requirements:
- report the spectrum on at least two scales (raw and a variance-stabilising
  transform) and show the gap sign is invariant;
- where a thermodynamic or explicit global-epistasis model is available for the
  system, fit it and compute the spectrum on the *residual* specific epistasis;
- state which scale the headline number uses.

**If the gap sign flips under a monotone rescaling, there is no result.**

### 5.3 Incompleteness — and why the sample-budget work is the enabler

Every measured dataset is combinatorially incomplete. Exact Walsh–Hadamard
decomposition needs $v(S)$ for all $2^n$ coalitions; Costanzo has a dense pairwise
slice, Kuzmin a sparse triple slice, Perturb-seq a sparse pair slice.

**This is exactly what `epistasis-sample-budget` was built to answer.** Its
recovery curves — what fraction of true order-3+ structure each method recovers
at budget $b$, calibrated against 28 systems with exact answers — tell us in
advance whether the coverage a real dataset actually has can support a spectral
claim at all.

**Do this before touching the biological data.** Compute the effective coverage
of each candidate dataset as a fraction of its $2^n$, look up the recovery curve
at that budget, and only proceed where recovery is adequate. If the curves say
recovery at Kuzmin's triple coverage is near zero, that is a two-hour finding
instead of a two-month one.

This also converts a limitation into a contribution: the budget paper stops being
a methods note and becomes the feasibility argument for the biology.

### 5.4 Ablation type

The project's own measurement is that **zero-ablation inflates pairwise epistasis
10–50× relative to mean ablation**. Gene deletion is closer to zero-ablation;
CRISPRi knockdown is partial and closer to an attenuation. The Boolean knockouts
are clamps.

These are not interchangeable, and any cross-system magnitude comparison must say
which regime it is in. Magnitudes should be compared within modality; only signs
and orderings should be compared across.

### 5.5 Coverage bias is not random

Genetic interaction screens are biased toward genes that are already
well-studied, and toward pairs someone had reason to test. A spectrum computed
over a biased sample is a spectrum of the sampling process as much as of the
biology. At minimum, report degree and study-bias distributions of the covered
set, and check the gap on a bias-matched subsample.

---

## 6. Work order

1. **Feasibility from the budget curves (§5.3).** No new data. Decides whether
   anything below is worth doing.
2. **Verify the dataset facts in §3** — counts, coverage, access routes, and
   above all whether BioGRID evidence codes can cleanly separate physical from
   genetic interactions. Any bulk download runs **on Modal against a volume**,
   never on a home connection.
3. **H3 first**, the loss-of-function / gain-of-function reversal, because it is
   directional and artifact-resistant.
4. **H1 and H2** on yeast, with both scales from §5.2.
5. Pre-register before computing any spectrum on real data, in the same style as
   `prereg_grn_v4.md`, with the decision rules fixed and the SHA stamped.

---

## 7. Open questions for Joe

1. Is there a yeast wiring model you would accept as independent of the genetic
   interaction data — and if BioGRID is the source, is filtering by evidence code
   sufficient, or is there contamination we would not see?
2. Does Kuzmin's trigenic set have any locally dense region — a small set of genes
   with all pairs *and* all triples measured — where an exact low-$n$ Walsh
   decomposition is possible without any recovery method? Even $n = 6$ complete
   would be worth more than a large sparse sample.
3. For the global-epistasis correction, is there an accepted transform in the
   yeast genetic-interaction literature we should adopt rather than inventing one?
4. Is the CRISPRa arm (H3) credible to that community as a gain-of-function
   analogue of the Boolean clamping, or is that mapping doing too much work?

---

## 8. Unverified — check before relying on any of it

Costanzo pair count, coverage and access route; Kuzmin triple count, trigenic
score definition, and local density; Norman Perturb-seq pair count and whether
triples exist; whether BioGRID evidence codes cleanly separate physical from
genetic interaction; and whether any of these datasets has a licence restricting
redistribution. None of this was fetched — it is written from background
knowledge and must be confirmed.

---

## 9. Other real-data epistasis systems, ranked

Added after the question "what about physics or chemistry". The criterion
throughout is **real measured or physically computed values, combinatorially
complete or nearly so** -- not a simulation of a curated model, and not
interactions between components of a neural network.

### 9.1 Protein fitness landscapes — real, measured, some genuinely complete

The strongest existing arm, and the one already partly in the paper.

| system | completeness | note |
|---|---|---|
| Poelwijk et al. 2019, fluorescent protein | **complete $2^{13}$ = 8,192, all measured** | already used; weighted Walsh--Hadamard, 260 significant terms to **7th order** |
| Weinreich et al. 2006, $\beta$-lactamase | **complete $2^5$**, 5 mutations, all 120 pathways | the classic; tiny but exact and famous |
| Khan et al. 2011, *E. coli* | **complete $2^5$**, 5 beneficial mutations | second independent small-complete landscape |
| Wu et al. 2016, GB1 | **complete over 4 sites $\times$ 20 amino acids** (~149k of 160k measured) | not $2^n$ but combinatorially complete on a 20-state alphabet |
| Olson et al. 2014, GB1 | all singles and all doubles, 55 positions | pairwise-complete, no order-3 |

**Why this matters more than it looks:** the small complete landscapes
($2^5$) let you compute an **exact** Walsh spectrum with no recovery method at
all. They are tiny, but they are exact and measured, which is the combination the
whole paper is short of. Several independent small-complete landscapes are worth
more here than one large sparse one.

**Unverified:** all counts above are from background knowledge. Weinreich and
Khan in particular need checking for exact completeness and current data access.

### 9.2 Chemistry — DFT as an exact, physics-grounded oracle

The most interesting new option, and the one that answers the circularity problem
most cleanly.

Pick a molecular scaffold with $n$ substitution positions and two choices at each.
That is $2^n$ molecules. Compute any property by density functional theory:
HOMO--LUMO gap, binding free energy, dipole, reaction barrier. **Non-additivity of
substituent effects is epistasis, and DFT gives it exactly.**

Three properties nothing else on this list has together:

1. **The oracle is physics, not a corpus.** No training distribution, no curated
   network, no measurement noise, no phylogeny.
2. **It is exhaustive at useful $n$.** $2^{10}$ = 1,024 small-molecule DFT
   calculations is routine; $2^{15}$ is a cluster job, not a research programme.
3. **The oracle evaluates anything.** Unlike a measured dataset, which can only
   score combinations someone already made, DFT scores whatever configuration you
   ask about — including ones in no dataset.

Substituent-effect non-additivity is also a real chemistry question with its own
literature (Hammett-type additivity and its failures), so the result has an
audience beyond the benchmark.

**This is the recommended second domain.** It is cheap, exact, and the cross-domain
claim — "the same method ranking holds on a Boolean GRN, a measured protein
landscape, and a DFT-computed molecular series" — is far stronger with a physical
oracle in it.

### 9.3 Materials — and a prior-art warning that must be checked first

Alloy configurational energetics: choose $n$ lattice sites, two species each,
compute energy by DFT for all $2^n$ configurations.

**Before anyone writes a line of this: the alloy community has been doing exactly
our decomposition since the 1980s.** The **cluster expansion** formalism (Sanchez,
Ducastelle and Gratias, 1984) expands configurational energy in *effective cluster
interactions* — pair, triplet, quadruplet terms — over a basis of site-occupation
variables. That is the same object as a Walsh--Hadamard or Möbius decomposition,
under a different name, with forty years of sparse-recovery methodology attached
(compressive sensing cluster expansion, genetic-algorithm cluster selection).

Two consequences, and both are load-bearing:

- **Opportunity.** Their sparse-recovery methods are directly importable to
  `epistasis-sample-budget`, and they have been stress-tested on physical systems
  for decades. If those methods beat LASSO-Walsh on our answer keys, that is a
  result in itself.
- **Risk.** Any claim of the form "nobody benchmarks higher-order interaction
  recovery against exact ground truth" **must be checked against the
  cluster-expansion literature before it is made.** They have exact DFT answer
  keys and they do compare recovery methods. The claim may need narrowing to
  biology, or reframing as a cross-domain transfer.

**Action: this check comes before the chemistry work, not after.** It is a
literature question, cheap to answer, and it can invalidate a framing claim.

### 9.4 What to reject, and why

- **Spin glasses and Ising models.** The Hamiltonian *is* the interaction
  spectrum, so recovering it is a self-consistency check rather than a discovery.
  Useful as a sanity test for a method, worthless as a scientific claim.
- **Ecology and microbial community assembly.** Higher-order interactions are a
  live question there, but combinatorially complete data is rare and the
  phenotype is noisy. Revisit only if a complete community dataset surfaces.
- **Anything involving interactions between neural-network components.** That is
  the cross-domain arm of `epistasis-bench` and it is a *different* paper. Keep
  the two separate: mixing them invites the reading that the biology is a
  metaphor for the interpretability, which weakens both.

### 9.5 Recommended shape

Three domains, one method ranking:

1. **Boolean GRNs** — exhaustive, exact, 28 systems. *Already done.*
2. **Measured protein landscapes** — real, several exactly complete at small $n$.
3. **DFT-computed molecular series** — physical oracle, exhaustive, extends
   off-distribution.

If a method ranks the same across a simulated network, a measured organism and a
physics calculation, that is as close to external validity as this kind of
benchmark can get. And each domain answers a different objection: the GRNs answer
scale, the proteins answer "is it real", the DFT answers "is your ground truth
just a summary of your data".
