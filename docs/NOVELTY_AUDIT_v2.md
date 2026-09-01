# Novelty audit v2: `composition_gap_v17.tex`

Supersedes `NOVELTY_AUDIT_v1.md`, which was written against v12. Scope here is the
composition-gap manuscript only. The estimator-benchmarking claims in v1 (P1--P4,
S1--S3) belong to the sample-budget arm and are not audited here.

**One claim refuted, two need narrowing, four citations are mandatory before
submission. The core result survives.**

## The check v1 flagged and did not perform, now performed

v1 ended on a question:

> Do the seven combinatorially complete maps in Sailer & Harms, *Genetics*
> 205(3):1079 (2017) include Weinreich 2006 TEM-1 and Hall 2010 yeast?

**Yes. Three of their seven maps are two of our three empirical datasets.**

Read from the full text (PMC5340324, EuropePMC), Table 1, *Data sets used in this
study*:

| ID | Genotype | Phenotype | L | Reference |
|---|---|---|---|---|
| I | Scattered genomic mutations | *E. coli* fitness | 5 | Khan et al. (2011) |
| II | Chromosomes in asexual fungi | *A. niger* fitness | 5 | de Visser et al. (2009) |
| **III** | **Protein point mutants** | **Bacterial fitness** | **5** | **Weinreich et al. (2006)** |
| IV | DNA/protein point mutants | binding affinity | 5 | Anderson et al. (2015) |
| V | Chromosomes in asexual fungi | *A. niger* fitness | 5 | de Visser et al. (2009) |
| **VI** | **Alleles in biosynthetic network** | ***S. cerevisiae* haploid growth** | **6** | **Hall et al. (2010)** |
| **VII** | **Alleles in biosynthetic network** | ***S. cerevisiae* diploid growth** | **6** | **Hall et al. (2010)** |

They report higher-order contributions of **2.2--31.0% of total variation, mean
12.7%**, after linearising for nonlinear scale.

Franke et al. (2011) is the one dataset of ours not among them; their *A. niger* maps
are de Visser et al. (2009), a different and smaller (L = 5) dataset, though de Visser
is an author on both.

## A0. The strongest attack is not priority, it is scale — and the answer is already in the paper

Recorded first because it outranks everything below.

**The attack.** *"Your order-3+ terms are an artifact of measurement scale. The value
function $v(S)$ is a bounded, saturating quantity in $[0,1]$. Push an additive
underlying effect through any saturating map and higher-order Walsh terms appear for
free. That is global (nonspecific) epistasis, a known phenomenon
\citep{otwinowski2018inferring, reddy2021global}, not composition."*

This is the objection a referee from the Lehner, Desai, McCandlish or Harms orbit
raises first, and it threatens the **core Boolean result**, not just the empirical arm.
It was researched in this project once already, against the protein-epistasis paper,
where it was written down as "the obvious attack."

**What the paper currently says.** Two sentences, in *Relationship to prior work*:

> The global epistasis literature \citep{otwinowski2018inferring, reddy2021global}
> establishes that nonlinear genotype-phenotype maps reshape interaction spectra. Our
> local-vs-global comparison operationalizes this reshaping as a quantitative spectral
> transformation, decomposed by interaction order.

That concedes the connection and does not rule out the deflationary reading. As
written, a referee can answer it with "yes, and that reshaping is global epistasis,
which is already understood."

**The paper already contains the control that defeats it, and does not use it that
way.** The matched nulls run the *same machinery* — same $2^n$ sweep, same 512 initial
conditions, same time-averaged bounded output, same saturating dynamics. If the gap
were a scale artifact of a bounded readout, the nulls would show it too. They do not:

| ensemble | mean $\Delta_{3+}$ |
|---|---|
| real networks | **+18.6 pp** |
| Kauffman NK | +3.5 pp |
| degree-preserving rewire | +3.8 pp |
| rule-preserving rewire | +10.3 pp |

Saturation is held constant across all four rows. What varies is the regulatory logic
and its arrangement. That is precisely the "residual specific epistasis after removing
the global component" control the earlier research went looking for and did not find in
the literature — and it is sitting in Section~\ref{sec:null-results} labelled as a
test of size and connectivity.

The Hill-function ODE arm strengthens it further: a Hill sigmoid is an *explicitly*
saturating readout with a different shape from the Boolean time-average, and the gap's
sign survives in every non-null network. A scale artifact would not be expected to
track across two different nonlinearities while random networks under both stay near
zero.

**Action.** Promote this. The nulls answer the global-epistasis objection and the
paper should say so where the objection is raised, not four subsections later under a
different heading. One paragraph, no new computation.

## Verdicts

| # | claim in v17 | verdict | threat |
|---|---|---|---|
| **A0** | **composition gap is a mechanism, not a scale artifact** | **survives, but under-defended** | **global epistasis \citep{otwinowski2018inferring, reddy2021global}** |
| A1 | abstract: "with what magnitude and sign, has not been quantified" | **narrow** | Sailer & Harms 2017; Weinreich et al. 2013 |
| A2 | intro l.115: "no known magnitude, no known sign structure" | **narrow** | same |
| A3 | §"Empirical fitness landscapes show higher-order creation" | **REFUTED as novel** | Sailer & Harms 2017, on our own data |
| A4 | Walsh--Hadamard framing, uncited prior art in the target journal | **citation mandatory** | Faure et al. 2024 |
| A5 | local rule spectra treated as unexamined | **narrow** | Shmulevich & Kauffman 2004 |
| A6 | the composition gap: rule spectrum vs attractor spectrum, same system | **survives** | -- |
| A7 | ODE replication, null-model decomposition, graded inversion, GNK control | **survives** | -- |

## A3, the one that must change

The subsection is titled *"Empirical fitness landscapes show higher-order creation"*
and the manuscript's own Methods contradict it:

> The empirical arm therefore tests whether higher-order structure exists beyond
> additivity, not whether composition created it.

Two defects, one presentational and one about priority. **Creation** cannot be
measured on these datasets because no local-rule model exists for them, which the
Methods already concede. And the measurement that remains -- higher-order energy in
Weinreich 2006 and Hall 2010 -- was published in 2017 on the same data.

**Action.** Retitle the subsection so it does not say *creation*. Cite Sailer & Harms
and state the relationship: these are three of their seven maps, our decomposition
recovers higher-order structure in the same data, and this section corroborates rather
than measures something new.

Citing them also does work the paper needs for A0. Sailer & Harms linearise each map
before extracting epistasis, which is the accepted correction for exactly the scale
objection, and their conclusion is that it survives: *"even after accounting for
nonlinearity, we found statistically significant high-order epistasis in all seven
maps,"* contributing 2.2--31.0% of variation, mean 12.7%. On two of our three
datasets, the scale-corrected answer is already published and it is favourable. That
is a stronger reason to cite them than priority.

## A4, the citation whose absence is conspicuous

**Faure, Lehner, Miró Pina, Serrano Colome and Weghorn, "An extension of the
Walsh--Hadamard transform to calculate and model epistasis in genetic landscapes of
arbitrary shape and complexity," *PLOS Computational Biology* 20(5):e1012132 (2024).**
DOI 10.1371/journal.pcbi.1012132, PMC11161127. Verified against EuropePMC.

Their extension is *multiallelic* -- it lifts the transform beyond two states per
position -- so it does not refute a binary knockout decomposition. That is not the
point. It is the current state of the art for Walsh--Hadamard epistasis, from the
Lehner group, published in the journal `NOVELTY_AUDIT_v1.md` names as this work's
target. A WHT-epistasis manuscript arriving at that journal without it reads as
unfamiliarity with the field.

## A1/A2, and why the narrowing is small

The claim as written is false: higher-order epistasis has a known magnitude in
empirical landscapes, and Sailer & Harms put a number on it. What is unmeasured is
narrower and is still ours -- nobody has compared a system's **measured local-rule
spectrum** against its **attractor phenotype spectrum**. Sailer & Harms decompose
phenotypes; they have no rule-level model to compare against. Shmulevich & Kauffman
(*Phys Rev Lett* 93:048701, 2004) decompose Boolean rules by activity and average
sensitivity, which is the first-order weight and the total influence, not the
higher-order spectrum, and they never carry it to an attractor phenotype.

The gap is real and unclaimed. Say that instead of the stronger thing.

**Proposed wording, adapted from v1:**

> Whether dynamical composition creates or destroys *higher-order* interactions, and
> with what sign, has not been measured. Higher-order epistasis has been quantified
> spectrally in empirical fitness landscapes \citep{weinreich2013,sailer2017} and
> Boolean update rules analysed by activity and average sensitivity
> \citep{shmulevich2004}, but the two spectra have not been compared within the same
> system.

## Mandatory before submission

1. Add `faure2024`, `sailer2017`, `weinreich2013`, `shmulevich2004` to `refs.bib`
   via `citations add --doi`. None is currently present.
2. Retitle §"Empirical fitness landscapes show higher-order creation" and reframe it
   as corroboration, naming the Sailer & Harms overlap explicitly.
3. Replace the abstract and intro sentences per A1/A2.
4. Re-run `citations audit --bib paper/refs.bib` after the additions.

## Unverified, do not quote

Manicka et al. 2023, named in v1 as decomposing Boolean update rules, did not resolve
in EuropePMC under that description and is omitted here rather than asserted. Whether
Weinreich et al. 2013 covers roughly fourteen landscapes is taken from v1 and was not
re-checked. Sailer & Harms's per-map percentages are reported as a range and mean from
their abstract; the per-map values for maps III, VI and VII were not extracted.

---

## Applied in v18 (2026-09-01)

- `faure2024walsh`, `sailer2017detecting`, `weinreich2013higherorder` and
  `shmulevich2004activities` added with `citations add --doi`; all four agree with
  their registry records on `citations audit`.
- A1/A2 narrowed in both the abstract and the introduction: higher-order epistasis
  has been measured spectrally in empirical landscapes and Boolean rules analyzed by
  activity, and what is unmeasured is the comparison of the two spectra in one system.
- A3 resolved by demotion. The empirical Results subsection is removed, its table
  moved to the supplement, and the point kept as two sentences in *Generalization and
  selection bias* naming the Sailer & Harms overlap and that higher-order epistasis
  survived their scale correction.
- A0 applied. The null ensembles are now given where the global-epistasis objection is
  raised, stating that they share the sweep, the initial conditions and the bounded
  readout, so a readout artifact would appear in them too.
- Three numeric defects fixed against `results/paper_number_reconciliation.json`:
  the degree-preserving null is +3.8 pp not +3.7; ODE sign preservation is 25 of 25
  non-null networks, not 28 of 28 (davidich_yeast flips inside the null band, and the
  ODE non-null count is 25 not 24, so 17/24 becomes 17/25); and the discussion's
  "18/6 ratio" is 19/6.
- Six British spellings corrected, two of them introduced by the v18 edits.

Everything else in the manuscript reconciles exactly. `results coverage` on v17 stood
at 24% of numbers bound to a run, up from 17%.
