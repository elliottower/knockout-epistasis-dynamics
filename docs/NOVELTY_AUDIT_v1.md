# Novelty audit: which priority claims survive

Adversarial literature check across cluster expansion, learning theory, global
sensitivity analysis, and genetics benchmarking. **Four claims are refuted and
five need narrowing.**

## The urgent one: the v12 title

**C3, the title of `composition_gap_v12.tex`:**
> "Higher-order knockout epistasis is unpredictable from molecular wiring **but
> recoverable from 1% of knockout combinations**"

**Refuted, and separately unsupported.**

*Refuted by:* **Faure, Lehner, Miró Pina, Serrano Colome and Weghorn, "An extension
of the Walsh--Hadamard transform to calculate and model epistasis in genetic
landscapes of arbitrary shape and complexity," *PLOS Comput Biol* 20(5):e1012132
(2024).** They subsample a combinatorially complete landscape **to exactly 1%**,
run LASSO, and recover epistatic coefficients at **r = 0.99** against known truth.
Verbatim: *"We sub-sampled the simulated phenotype (fitness) values to obtain
training dataset sizes ranging from 64% to 1% of all variants."*

Biology, Walsh--Hadamard framing, a journal this work would target, and **already
cited in `background.md` footnote 32** without the budget sweep being noticed.

*Unsupported:* the recovery result appears nowhere in v12's body, methods,
results or supplement. **The only difference between v11 and v12 is that title.**

**Action:** revert to the v11 title, or if the sample-budget arm is folded in,
avoid "1%" entirely. Suggested: *"Higher-order knockout epistasis is created by
dynamics and unpredictable from molecular wiring"*.

## Verdicts

| # | claim | verdict | threat |
|---|---|---|---|
| C1 | higher-order extension "has not been quantified" | **narrow** | Sailer & Harms 2017; Weinreich et al. 2013 |
| C2 | "no known magnitude, no known sign structure" | **narrow** | same |
| **C3** | **title: "recoverable from 1%"** | **REFUTED** | **Faure et al. 2024** |
| C4 | "extend this in three directions" | survives | -- |
| C5 | "local-vs-global comparison operationalizes this" | survives | -- |
| C6 | Table 2 empirical-landscape decomposition | **REFUTED as novel** | Sailer & Harms 2017 |
| P1 | "head-to-head benchmarking is rare" | **REFUTED** | Russ 2022; Herder 2015; Puy 2022; shapiq 2024 |
| P2 | "validated their own methods, not the installed base" | survives | -- |
| P3 | "nobody has scored [list] on the same exact answer key" | **narrow** | Russ 2022; Listopad 2025 (simulated) |
| P4 | "we extend this to a systematic head-to-head" | **narrow** | shapiq NeurIPS 2024 D&B |
| S1 | "nobody else has 28 systems with exact 2^n tables" | **narrow** | shapiq: 2,042 exact configurations, n ≤ 16 |
| S2 | "Sobol indices (the continuous analogue of WHT)" | **REFUTED, and backwards** | O'Donnell Ch. 8 |
| S3 | "grounded in exact truth, not simulation" | survives | -- |

## Cluster expansion: the correspondence is exact, and published as such

Not analogical. The alloy literature says it in its own abstracts:

- Sanchez, *Phys Rev B* **81**, 224202 (2010): *"It is shown that the cluster basis
  developed by Sanchez et al. [Physica A 128, 334 (1984)] is a multidimensional
  discrete Fourier transform."*
- Barroso-Luque and Ceder, *npj Comput Mater* **10**, 158 (2024): *"we identify the
  cluster decomposition as an invariant ANOVA decomposition."*

Effective cluster interaction $J_\alpha$ = Walsh coefficient $\hat f(S)$; cluster
function $\Gamma_\alpha$ = Walsh character $\chi_S$; empty cluster = $\hat f(\emptyset)$.

**And they benchmarked recovery against planted truth in 2013.** Nelson, Zhou,
Hart and Ozoliņš, *PRB* **87**, 035125: *"we choose a set of sparse coefficients and
then use them to compute the energies... knowing the exact solution a priori
allows us to easily determine the accuracy of the solution found."* Figure 2(c):
$\lVert J_\text{exact} - J_\text{fit}\rVert_1$ against number of fitting
structures and noise level, averaged over ~100 subsets. That is the sample-budget
sweep, thirteen years old.

Herder, Bray and Schneider, *Surface Science* **640**:104 (2015) is the head-to-head:
three fitting algorithms against synthetic data with hidden cluster expansions,
scored on support recovery.

**Symmetry averaging does not save us** -- it makes cluster expansion a *special
case* (the Walsh decomposition restricted to the trivial isotypic component of the
space-group action). Arguing we lack symmetry tying argues our object is more
generic, i.e. the textbook one.

**What does survive:** no lattice metric, so no cluster-diameter truncation --
knockout landscapes can only be truncated by order, a strictly weaker inductive
bias, and standard cluster-expansion model selection does not transfer. Also,
nobody in the alloy literature writes "Walsh--Hadamard": a citation sweep of 100
papers citing Sanchez 2010 returned zero hits for Walsh, Hadamard, Boolean or
Möbius. **The bridge is unwritten; the object is not novel in either field.**

## The seam that survives, and it is worth more than the refuted claims

Every global-sensitivity-analysis benchmark scores **first-order and total-order
indices only**. Puy, Becker, Lo Piano and Saltelli, *IJUQ* 12(2) (2022) compare
eight estimators with randomized sample budgets -- and **generate controlled
third-order structure** ($k_3 \sim U(0.1,0.3)$) **then discard it by scoring only
total-order.** The flagship benchmark had order-3 ground truth in hand and did
not score it.

Likewise every cluster-expansion software comparison scores *held-out prediction
error*, never *recovery of known true effective cluster interactions*.

> **No multi-estimator head-to-head on order $\geq 3$ coefficient recovery, with
> error-against-budget curves scored on exact truth, exists in any of the six
> literatures searched.**

That is the contribution. It is narrower than "nobody benchmarks epistasis
detection" and it is defensible.

## Proposed narrowed wordings

**C1/C2** → "Whether dynamical composition creates or destroys *higher-order*
interactions, and with what sign, has not been measured. Higher-order epistasis
has been quantified spectrally in empirical fitness landscapes [Weinreich et al.
2013; Sailer & Harms 2017] and Boolean update rules decomposed independently
[Shmulevich & Kauffman 2004; Manicka et al. 2023], but the two spectra have not
been compared within the same system."

**P1** → "Existing comparisons benchmark epistasis detectors on simulated or
semi-simulated genotypes with planted interactions [Shang 2011; Chatelain 2018;
Slim 2020; Russ 2022; Listopad 2025]. None scores detectors against an
exhaustively enumerated interaction spectrum."

**P3/P4** → "Multi-method benchmarks exist in genetics [Russ et al. 2022, twelve
detectors] and in explainable ML [Muschalik et al. 2024, seven interaction
approximators against exhaustively computed ground truth]. The genetics
benchmarks use simulated truth; the ML benchmark uses the Shapley-interaction
basis on non-biological games. Neither scores recovery of Walsh--Hadamard
coefficients at order $\geq 3$ on biological systems under a controlled sample
budget."

**S1** → "Twenty-eight *biological* systems with exact $2^n$ knockout tables,
spanning $n = 7$--18."

**S2** → "The Walsh--Hadamard decomposition is the Hoeffding/Efron--Stein ANOVA
decomposition specialised to the uniform Boolean cube; Sobol indices and Walsh
coefficients are the same object read in two bases."

## Methods now mandatory as competitors

From cluster expansion, all reachable via `sparse-lm` or scikit-learn:
**ARDR / Bayesian compressive sensing** (Nelson et al., *PRB* 88, 155105, 2013 --
beats OLS/RFE precisely in the underdetermined low-budget regime, which is our
regime); **split-Bregman $\ell_1$** (Nelson et al., *PRB* 87, 035125, 2013);
**$\ell_0\ell_2$ best subset + ridge** (Zhong et al., *PRB* 106, 024203, 2022);
**RFE-OLS** and **adaptive LASSO** (Ångqvist et al. 2019; Fransson et al. 2020).

From sparse-WHT and learning theory: **q-SFT / GFast**;
**ProxySPEX / SPEX** (NeurIPS 2025; ICML 2025) -- ProxySPEX explicitly recovers
interactions between attention heads, directly on-target for the circuit tables;
and **shapiq's own approximators** (KernelSHAP-IQ, SVARM-IQ, SHAP-IQ). The last is
now **mandatory** -- omitting the published benchmark's own methods is the fastest
route to rejection.

For a continuous / graded-knockdown arm: **sensobol** 4th-order estimator (Puy et
al., *JSS* 102(5), 2022) and **sparse PCE via LAR** (Blatman & Sudret 2011).

## Highest-priority manual check

**Do the seven combinatorially complete maps in Sailer & Harms, *Genetics*
205(3):1079 (2017) include Weinreich 2006 TEM-1 and Hall 2010 yeast?** They report
2.2--31.0% higher-order energy, mean 12.7%. If their maps overlap v12's Table 2,
that section partly reproduces a 2017 result **on the same data** and must be
reframed as replication.

## Bearing on `protein-epistasis`

Two of these threats apply directly to the new repo and its SPEC must be updated:
**Weinreich, Lan, Wylie and Heckendorn, *Curr Opin Genet Dev* 23(6):700 (2013)**
already Walsh-transformed ~14 combinatorially complete landscapes; **Sailer and
Harms 2017** already quantified higher-order energy across seven. The pairwise-ceiling
claim must be positioned against both, and the novel part is likely the *ceiling*
framing -- what a pairwise model structurally cannot reach -- rather than the
measurement of higher-order energy itself.

## Unverified -- do not quote

Sanchez 2010's body sentence naming "Walsh functions" (APS 403s; the
abstract-level "multidimensional discrete Fourier transform" is verified and
sufficient); the 1984 Sanchez--Ducastelle--Gratias abstract verbatim; shapiq's
budget-sweep sentences and figure numbers; Poelwijk 2019 Figure 4 legend; whether
Herder 2015 plots error explicitly against fit-set size; whether sensobol's JSS
paper plots 3rd/4th-order error against $N$.
