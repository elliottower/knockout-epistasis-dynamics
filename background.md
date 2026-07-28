# EpistasisBench Background Research: Cross-Domain Epistasis Detection Methods

## Genetics-Native Epistasis Detection Methods

### Epi-MEIF (mixed-effect conditional inference forests)

Epi-MEIF, published in *Nucleic Acids Research* in 2022 by Saha, Perrin, Röder, Brun, and Spinelli, is a method for detecting higher-order epistatic interactions using mixed-effect conditional inference forests fitted on a curated group of candidate SNPs. The tree structure of the forest is used to identify n-way interactions among the candidate SNPs, and the authors add supplementary testing strategies to improve robustness against false positives. The method was validated through extensive simulation on cross-sectional and longitudinal synthetic datasets and then applied to detect epistatic interactions underlying cardiac trait variation in *Drosophila*. Because it operates on tree splits rather than exhaustive combinatorial search, Epi-MEIF is explicitly designed to scale to "any GWAS data," trading exactness for tractability — precisely the compromise the EpistasisBench pitch highlights as universal in genetics.[^1][^2]

### DFIM (Deep Feature Interaction Maps)

DFIM, published in *Bioinformatics* (2018) by the Kundaje lab, estimates pairwise interactions between features (nucleotides or motifs) in a deep neural network's input by extending feature attribution methods to interaction attribution. It efficiently identifies all-pairs interactions in an input DNA sequence and was validated against embedded ground-truth motif interactions in simulated regulatory sequences, then applied to real transcription-factor binding models (GATA1–TAL1 synergy) and chromatin accessibility models. DFIM is the genetics method most structurally analogous to mech-interp circuit analysis, since it operates directly on a trained neural network's internal attributions rather than raw genotype-phenotype data — DFIM code is open-sourced at github.com/kundajelab/dfim.[^3]

### MDR and its descendants (GMDR, MB-MDR)

Multifactor Dimensionality Reduction (MDR), introduced by Ritchie et al. in 2001, is the foundational combinatorial approach in this space: it reduces high-dimensional multilocus genotype combinations into a single binary "high-risk/low-risk" variable via a case/control ratio threshold, then evaluates predictive power through 10-fold cross-validation and permutation testing. MDR is nonparametric and model-free, avoiding the curse-of-dimensionality problems that logistic regression faces with high-order interaction terms.[^4][^5]

Two major extensions emerged. Generalized MDR (GMDR), from Lou et al. (2007), generalizes MDR to allow covariate adjustment and continuous (not just dichotomous) phenotypes by using a residual score from a generalized linear model instead of a simple case/control ratio. Model-Based MDR (MB-MDR), from Cattaert et al., unifies parametric and nonparametric approaches, adjusts for lower-order effects and confounders, and generally shows higher empirical power than classic MDR, especially under genetic heterogeneity or low minor allele frequency. A pedigree-based extension (PGMDR) further adapts GMDR for family-based study designs.[^6][^7][^8]

### BOOST (Boolean Operation-based Screening and Testing)

BOOST, from Wan et al. (*American Journal of Human Genetics*, 2010), is an exhaustive pairwise SNP-SNP interaction screening method optimized for speed via bitwise Boolean operations on contingency tables. It uses a two-stage design: a fast approximate screen via the Kirkwood superposition approximation of the log-likelihood ratio, followed by rigorous re-evaluation of the surviving candidate pairs with a full log-likelihood ratio test. BOOST completed exhaustive pairwise analysis of ~360,000 SNPs on a standard desktop in under 60 hours per dataset and has since been incorporated into PLINK (`--fast-epistasis boost`). Its key limitation is sparse contingency tables at low minor allele frequency, which reduces detection power.[^9][^10][^11][^12][^13]

### AntEpiSeeker (ant colony optimization)

AntEpiSeeker, from Wang, Liu, Robbins, and Rekaya (2010), is a two-stage ant colony optimization (ACO) algorithm for detecting epistasis in case-control designs. Artificial "ants" represent candidate SNP sets; a shared pheromone-concentration probability distribution biases which SNPs are picked in each iteration based on prior chi-square association scores, and a second exhaustive-search stage refines the highest-pheromone candidate sets. This heuristic search strategy is representative of the broader combinatorial-optimization family in genetics (later extended in MACOED and epiACO), which trades completeness for tractable search over large SNP spaces — again illustrating that genetics methods systematically approximate rather than exhaustively verify.[^10][^14][^15][^9]

### Information-theoretic / entropy-based detection (GAIN and relatives)

A distinct family of methods detects epistasis via information theory rather than classification accuracy. The core quantity is information gain: \( IG(A;B;C) = I(A,B;C) - I(A;C) - I(B;C) \), where \( I \) denotes mutual information between SNPs and phenotype \( C \). This isolates "pure" interaction effects (synergy) from the additive main effects of each SNP individually. An extension to three-way interactions (a fast, nonparametric, model-free three-way epistasis measure) was applied to tuberculosis susceptibility data in a West African population and found a pure three-way interaction stronger than any lower-order association. Gene-based extensions (GBIGM) and Statistical Epistasis Networks (SEN) build on the same information-gain foundation to map genome-wide interaction structure as weighted networks.[^16][^17][^18][^19]

## Shapley Interaction Indices (the formal game-theoretic line)

### shapiq library and unified interaction computation

`shapiq` is an open-source Python package (Muschalik et al., NeurIPS 2024 Datasets and Benchmarks Track) that unifies state-of-the-art algorithms for computing Shapley values (SVs) and any-order Shapley interactions (SIs) in a single application-agnostic framework. Critically for EpistasisBench's design, `shapiq` ships with `shapiq.ExactComputer`, which computes 20 different interaction indices exactly for small "games" (feasible for roughly \(n \le 15\)–20 players before combinatorial blowup), plus a benchmarking suite of 11 ML applications with pre-computed ground truth. The library also includes multiple sampling-based approximators — KernelSHAP-IQ, SVARM-IQ, Permutation Sampling, Owen Sampling — explicitly benchmarked against exact ground truth as a function of "model evaluations" (forward-pass budget), which is directly the subsampling-vs-accuracy protocol proposed for EpistasisBench.[^20][^21]

### Faith-Shap (Faithful Shapley Interaction Index)

Faith-Shap, from Tsai, Yeh, and Ravikumar (*JMLR*, 2023), derives Shapley interactions from a different first-principles axiomatization: rather than starting from combinatorial game-theoretic axioms directly, it treats Shapley values as coefficients of the best linear approximation to a pseudo-Boolean value function, then generalizes this to higher-order polynomial approximations. Requiring the resulting "faithful interaction indices" to also satisfy interaction-extended versions of the standard Shapley axioms (dummy, symmetry, linearity, efficiency) yields a unique index, Faith-Shap. In the `shapiq` benchmark, Faith-Shap (FSII) is explicitly optimized for approximation faithfulness and shows the largest and most consistent gains over plain SHAP as interaction order increases.[^22][^23][^20]

### Shapley-Taylor Interaction Index (STII)

The Shapley-Taylor Interaction Index, introduced by Sundararajan, Dhamdhere, and Agarwal (ICML 2020), generalizes the Shapley value using a Taylor-series-style expansion over subsets of features rather than singletons, producing an interaction index that reduces to the ordinary Shapley value at order one. It is one of the four indices (alongside k-SII, FSII, and FBII) directly compared in the `shapiq` benchmarking analysis, where STII shows different faithfulness/order tradeoffs than Faith-Shap depending on application domain.[^24][^20]

### SHAP-IQ and k-SII

SHAP-IQ (2022 preprint, later folded into the `shapiq` package) provides a unified sampling-based estimator for the Shapley Interaction Index (SII) at arbitrary order, as well as the Shapley-Taylor Index (STI) and Shapley-Faith Index (FSI), using pairing tricks and stratified sampling kernels to improve budget efficiency. The k-order Shapley Interaction Index (k-SII) is the most commonly used any-order variant in practice and is one of the primary indices `shapiq`'s benchmarking suite reports results for across 100 unique benchmark games.[^21][^25][^20]

## Ground Truth and Combinatorially Complete Landscapes

### Weinreich and the combinatorially complete landscape concept

The concept of a "combinatorially complete fitness landscape" — measuring all \(2^n\) combinations of a fixed set of \(n\) mutations rather than sampling sparsely — traces to Weinreich's line of work, most directly discussed in a widely cited 2013 review on higher-order epistasis in evolutionary genetics. This approach "helps provide a thorough picture of the multidimensional patterns of epistasis," but has historically been limited to very small \(n\) (four or five loci) because experimental combinatorial coverage becomes intractable past that scale. This is the direct genetics-side analog of what the EpistasisBench pitch describes as "genetics can never run \(2^n\)."[^26][^27][^28]

A concrete example: a 2024 PNAS paper by Johnston et al. reports a combinatorially complete, 160,000-variant fitness landscape across four residues of an enzyme active site, explicitly framed as enabling "exact calculation of all epistasis in the landscape" because it is one of only a few studies extending combinatorial completeness beyond two or three sites. Even at four sites (\(2^4\) states per residue times mutational depth), this required a dedicated high-throughput selection assay — reinforcing how rare true combinatorial completeness is in wet-lab genetics, and how much smaller in scale than a 15-node exhaustive Shapley computation.[^29]

### Walsh-Hadamard transform as the formal epistasis decomposition

The Walsh-Hadamard transform (WHT) is the standard mathematical formalism unifying different quantitative definitions of epistasis in genetics, decomposing a fitness or phenotype function over a Boolean (or more generally multiallelic) genotype space into additive, pairwise, and higher-order interaction coefficients. Poelwijk's chapter on this method notes the WHT "unifies a number of different definitions of epistasis" and requires an explicit null hypothesis for independent mutational effects before higher-order coefficients can be interpreted as true epistasis. A 2023 extension generalizes the two-allele WHT formalism to arbitrary numbers of states per position (e.g., 20 amino acids or 4 nucleotides), which is directly relevant if EpistasisBench wants to extend beyond binary "head present/absent" circuit states. This Walsh-Hadamard/Fourier-over-the-hypercube view is also the natural bridge to Shapley interaction indices, since both decompose a set function over \(2^n\) subsets into an orthogonal or additive interaction hierarchy.[^30][^31][^32]

### Sparse/sublinear sample recovery of Walsh-Hadamard coefficients

The subsampling angle for EpistasisBench — recovering the interaction structure from a fraction of the full \(2^n\) coalition table — has a well-developed theoretical basis in sparse Fourier/Walsh-Hadamard transform recovery. The SPRIGHT algorithm (Li, Bradley, Pawar, and Ramchandran, building on Scheibler's SparseFHT) computes a \(K\)-sparse \(N\)-point Walsh-Hadamard transform using only \(O(K \log N)\) samples and \(O(K \log^2 N)\) time, and is provably robust to measurement noise at essentially no extra sample cost. A related line from theoretical computer science (Hassanieh-style sparse Fourier transform algorithms) achieves near-linear-in-\(K\) sample complexity for the Boolean-cube (Walsh-Hadamard) case specifically. On the compressed-sensing side, Boolean function learning via linear programming relaxations of Boolean group-testing has been shown to recover certain conjunctive/disjunctive rule structures exactly under conditions on the "sensing matrix," directly paralleling exact-recovery guarantees in classical compressed sensing. These results collectively suggest a theoretical ceiling for how efficiently *any* subsampling-based epistasis detector could recover sparse higher-order structure — a natural benchmark reference point for the compute-budget experiment.[^33][^34][^35][^36]

### Perturb-seq as combinatorial genetic-interaction ground truth in cells

Outside sequence-fitness landscapes, Perturb-seq (Dixit et al., *Cell*, 2016) is the closest genetics analog to exhaustive combinatorial perturbation at the level of gene regulatory circuits rather than sequence mutations: it combines pooled CRISPR perturbations with single-cell RNA-seq readout to measure transcriptional effects of many perturbations (and, in later extensions, perturbation pairs) simultaneously in the same pooled experiment. It is presented as accurately identifying "individual gene targets, gene signatures, and cell states affected by individual perturbations and their genetic interactions," but even Perturb-seq's combinatorial screens cover only a sparse subset of all possible perturbation combinations, not an exhaustive \(2^n\) grid — reinforcing that even modern single-cell functional genomics has not achieved the completeness Elliot's 15-head brute-force provides.[^37][^38]

## Mechanistic Interpretability Circuit Discovery Methods

### ACDC (Automated Circuit Discovery)

ACDC, introduced by Conmy, Mavor-Parker, Lynch, Heimersheim, and Garriga-Alonso at NeurIPS 2023, formalizes and automates the third step of the standard mech-interp workflow (isolate behavior, decompose network, discover circuit) via a greedy iterative pruning algorithm. It iterates from output nodes toward inputs, permanently removing an edge if doing so changes the KL divergence between the pruned and full model's output distributions by less than a threshold \( \tau \). Edges are removed via interchange (activation-patching) ablations rather than zero- or mean-ablation, since interchange ablations keep replaced activations within realistic ranges. ACDC rediscovered all 5 known component types of GPT-2 Small's Greater-Than circuit, selecting only 68 of 32,000 candidate edges, all previously identified manually. However, across a broader benchmark of five circuits, ACDC's overall AUC (0.596) trailed a simpler gradient-based baseline, Subnetwork Probing (0.692), and ACDC's results are known to be sensitive to parent iteration order, indicating the greedy local search is not globally optimal.[^39][^40][^41]

### Edge Attribution Patching (EAP)

EAP, from Syed, Rager, and Conmy (2023), replaces ACDC's per-edge forward pass with a first-order Taylor approximation of activation patching's effect: \( L(x_{clean} \mid do(E = e_{corr})) \approx L(x_{clean}) + (e_{corr} - e_{clean})^{T} \frac{\partial L(x_{clean} \mid do(E = e_{clean}))}{\partial e_{clean}} \). This requires only two forward passes and one backward pass to estimate importance for *all* edges simultaneously, versus one forward pass per edge for ACDC — a dramatic scalability improvement. Averaged across benchmark tasks, EAP achieves higher AUC for circuit recovery than ACDC despite its attribution scores being poorly correlated with true activation-patching scores; running ACDC as a refinement step after EAP performs best of all. A follow-on method, Edge Pruning (Bhaskar et al., 2024), formulates circuit discovery as gradient-based optimization over a continuous relaxation of a binary edge mask, outperforming both ACDC and EAP and scaling to 13B-parameter models, though it is slower on small datasets.[^42][^41][^39]

### Relevance to EpistasisBench's cross-domain framing

Both ACDC and EAP already use "ground truth" circuits from prior manual mech-interp work (e.g., the IOI circuit) as evaluation targets, framing circuit discovery evaluation as an edge-classification problem scored via ROC/AUC — structurally identical to how EpistasisBench proposes scoring epistasis-detection methods against exact ground truth. The key distinction the pitch identifies is that ACDC/EAP's "ground truth" circuits are themselves derived from earlier manual, potentially error-prone activation-patching work, not from an exhaustive combinatorial answer key — this is exactly the gap a 2^15 exact-Shapley dataset over attention heads would close, since it would be a mathematically exact ground truth rather than a consensus of prior noisy discovery methods.[^40][^41][^39]

## Cross-Domain Synthesis

The table below summarizes how each candidate method's original validation strategy compares to what EpistasisBench would newly provide.

| Method family | Origin field | Native validation approach | What EpistasisBench adds |
|---|---|---|---|
| Epi-MEIF (conditional inference forests) | Genetics | Simulation + Drosophila cardiac trait application, no exhaustive ground truth[^1] | Exact 2^15 Shapley ground truth for direct accuracy scoring |
| DFIM (deep feature interaction maps) | Genomics/NN interpretability | Embedded ground-truth motifs in simulated sequences[^3] | Real mech-interp circuit with exact higher-order ground truth |
| MDR / GMDR / MB-MDR | Genetics | Cross-validation, permutation testing, empirical power studies[^4][^6][^7] | Direct comparison against Walsh coefficients, not just prediction accuracy |
| BOOST | Genetics (GWAS) | Exhaustive pairwise only, no higher-order, contingency-table power limits[^10][^12] | Tests whether pairwise-only exhaustive search misses real higher-order structure |
| AntEpiSeeker (ACO) | Genetics | Chi-square scoring on simulated epistatic models[^9][^10] | Scored against exact answer rather than simulated approximate models |
| Information gain / GAIN | Genetics | Permutation p-values on case-only or case-control cohorts[^16][^18] | Exact mutual-information decomposition possible at n equal 15 |
| shapiq (k-SII, STII, FSII, FBII) | Explainable ML | Own benchmark suite with pre-computed games, but domain-general, not genetics or mech-interp specific[^20][^21] | Genetics-relevant framing plus mech-interp application domain |
| ACDC / EAP | Mechanistic interpretability | ROC/AUC against manually-discovered circuits (not exhaustive)[^40][^42] | Exact combinatorial ground truth replaces noisy manual "ground truth" |

The consistent pattern across every genetics-native method surveyed is validation against either simulated data with known planted effects, cross-validation/permutation testing, or comparison to other approximate methods — never exhaustive combinatorial ground truth beyond a handful of loci, because wet-lab combinatorial completeness tops out around four to five sites even in the most ambitious published landscapes. Mech-interp's own automated circuit discovery literature (ACDC, EAP) independently converged on the same evaluation paradigm — treating prior manually-discovered circuits as "ground truth" and scoring new methods via ROC/AUC — which inherits the same weakness of not being provably exhaustive or exact. The shapiq benchmarking framework is the only surveyed system that already computes exact ground truth via `ExactComputer`, but its 11 benchmark applications are general ML tasks (data valuation, local explanation, uncertainty explanation), not genetics or transformer-circuit specific. This confirms the pitch's central claim: no existing benchmark combines mechanistic-interpretability-scale exact ground truth with genetics' epistasis-detection method inventory, and the sparse Walsh-Hadamard recovery theory (SPRIGHT, sparse Fourier-over-hypercube algorithms) provides an established theoretical framework for the proposed compute-budget-vs-accuracy experiment.[^27][^28][^34][^35][^41][^20][^21][^33][^29][^40]

---

## References

1. [Epi-MEIF: detecting higher order epistatic interactions for ...](https://pubmed.ncbi.nlm.nih.gov/36107776/) - by S Saha · 2022 · Cited by 15 — In this article, we propose a novel method for higher-order epistas...

2. [Epi-MEIF, a flexible and efficient method for detection of ...](https://www.biorxiv.org/content/biorxiv/early/2021/12/22/2021.12.21.473474.full.pdf) - In this article, we propose a novel method for high-order epistasis detection using mixed effect con...

3. [Discovering epistatic feature interactions from neural network ...](https://pubmed.ncbi.nlm.nih.gov/30423062/) - Supplementary data are available at Bioinformatics online.

4. [Multifactor dimensionality reduction for detecting gene-gene and gene-environment interactions in pharmacogenomics studies - PubMed](https://pubmed.ncbi.nlm.nih.gov/16296945/) - In the quest for discovering disease susceptibility genes, the reality of gene-gene and gene-environ...

5. [Multifactor Dimensionality Reduction (MDR)](https://ritchielab.org/research/research-areas/genetic-architecture-of-complex-traits/mdr) - A multifactor dimensionality reduction (MDR) approach to detecting and characterizing high-order gen...

6. [Model-based multifactor dimensionality reduction for detecting epistasis in case-control data in the presence of noise - PubMed](https://pubmed.ncbi.nlm.nih.gov/21158747/) - Analyzing the combined effects of genes and/or environmental factors on the development of complex d...

7. [A combinatorial approach to detecting gene-gene and gene-environment interactions in family studies - PubMed](https://pubmed.ncbi.nlm.nih.gov/18834969/) - Widespread multifactor interactions present a significant challenge in determining risk factors of c...

8. [Multivariate generalized multifactor dimensionality reduction to detect gene-gene interactions](https://pmc.ncbi.nlm.nih.gov/articles/PMC4029529/) - Recently, one of the greatest challenges in genome-wide association studies is to detect gene-gene a...

9. [AntEpiSeeker: detecting epistatic interactions for case-control studies using a two-stage ant colony optimization algorithm - PubMed](https://pubmed.ncbi.nlm.nih.gov/20426808/) - AntEpiSeeker is a powerful and efficient tool for large-scale association studies and can be downloa...

10. [A survey about methods dedicated to epistasis detection](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2015.00285/full) - by C Niel · 2015 · Cited by 214 — BOOST runs an exhaustive analysis of all potential pairwise SNP-SN...

11. [macyang - BOOST](https://sites.google.com/site/eeyangc/software/boost) - X. Wan*, C. Yang*, Q. Yang, H. Xue, X. Fan, N. Tang and W. Yu. BOOST: a fast approach to detecting g...

12. [BOOST: A fast approach to detecting gene-gene interactions in genome- ...](https://pubmed.ncbi.nlm.nih.gov/20817139/) - by X Wan · 2010 · Cited by 641 — BOOST allows examination of all pairwise interactions in genome-wid...

13. [BOOST: A Fast Approach to Detecting Gene-Gene Interactions in Genome ...](https://pmc.ncbi.nlm.nih.gov/articles/PMC2933337/) - by X Wan · 2010 · Cited by 641 — BOOST allows examination of all pairwise interactions in genome-wid...

14. [MACOED: a multi-objective ant colony optimization algorithm for SNP epistasis detection in genome-wide association studies](https://academic.oup.com/bioinformatics/article/31/5/634/2748185)

15. [epiACO - a method for identifying epistasis based on ant Colony optimization algorithm](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5500974/) - Identifying epistasis or epistatic interactions, which refer to nonlinear interaction effects of sin...

16. [An information-gain approach to detecting three-way epistatic interactions in genetic association studies - PubMed](https://pubmed.ncbi.nlm.nih.gov/23396514/) - Our study provides a methodological basis for detecting and characterizing high-order gene-gene inte...

17. [A gene-based information gain method for detecting gene–gene interactions in case–control studies - European Journal of Human Genetics](https://www.nature.com/articles/ejhg201516) - Currently, most methods for detecting gene–gene interactions (GGIs) in genome-wide association studi...

18. [Characterizing gene-gene interactions in a statistical epistasis network of twelve candidate genes for obesity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4693412/) - Recent findings have reemphasized the importance of epistasis, or gene-gene interactions, as a contr...

19. [Entropy-Based Information Gain Approaches to Detect and ...](https://pmc.ncbi.nlm.nih.gov/articles/PMC3384547/) - For complex diseases, the relationship between genotypes, environment factors and phenotype is usual...

20. [shapiq: Shapley Interactions for Machine Learning](https://proceedings.neurips.cc/paper_files/paper/2024/file/eb3a9313405e2d4175a5a3cfcd49999b-Paper-Datasets_and_Benchmarks_Track.pdf)

21. [Benchmarking approximators of Shapley interactions#](https://shapiq.readthedocs.io/en/latest/notebooks/benchmark_notebooks/benchmark_approximators.html)

22. [Faith-Shap: The Faithful Shapley Interaction Index](https://jmlr.org/papers/v24/22-0202.html)

23. [the Faithful Shapley Interaction Index](https://icml.cc/media/icml-2023/Slides/25670.pdf)

24. [The Shapley Taylor Interaction Index](https://proceedings.mlr.press/v119/sundararajan20a/sundararajan20a.pdf) - by K Dhamdhere · Cited by 57 — The Shapley value is a commonly used method to attribute a model's pr...

25. [GitHub - FFmgll/shapiq](https://github.com/FFmgll/shapiq) - Contribute to FFmgll/shapiq development by creating an account on GitHub.

26. [Should evolutionary geneticists worry about higher-order ...](https://gwern.net/doc/genetics/heritable/2013-weinrich.pdf)

27. [Quantitative Description of a Protein Fitness Landscape Based on ...](https://academic.oup.com/mbe/article/32/7/1774/1017516) - Abstract. Understanding the driving forces behind protein evolution requires the ability to correlat...

28. [Next Level Challenges And...](https://pmc.ncbi.nlm.nih.gov/articles/PMC4254422/) - A combinatorially complete data set consists of studies of all possible combinations of a set of mut...

29. [A combinatorially complete epistatic fitness landscape in ...](https://pmc.ncbi.nlm.nih.gov/articles/PMC11317637/) - Predictive models for protein engineering seek to capture the relationship between protein sequence ...

30. [A Scalable Walsh-Hadamard Regularizer to Overcome the Low ...](https://proceedings.mlr.press/v216/gorji23a/gorji23a.pdf) - by A Gorji · 2023 · Cited by 10 — Higher-order epistasis shapes the fitness landscape of a xenobioti...

31. [An extension of the Walsh-Hadamard transform to calculate and model epistasis in genetic landscapes of arbitrary shape and complexity](https://www.biorxiv.org/content/10.1101/2023.03.06.531391v5) - Accurate models describing the relationship between genotype and phenotype are necessary in order to...

32. [Context-Dependent Mutation Effects in Proteins - PubMed](https://pubmed.ncbi.nlm.nih.gov/30298395/) - Defining the extent of epistasis-the nonindependence of the effects of mutations-is essential for un...

33. [A Fast and Robust Framework for Sparse Walsh-Hadamard Transform](https://arxiv.org/abs/1508.06336) - by X Li · 2015 · Cited by 20 — Abstract:We consider the problem of computing the Walsh-Hadamard Tran...

34. [Sparse fourier transform in any constant dimension with nearly-optimal sample complexity in sublinear time | Proceedings of the forty-eighth annual ACM symposium on Theory of Computing](https://dl.acm.org/doi/10.1145/2897518.2897650)

35. [LCAV/SparseFHT: A Fast Hadamard Transform for Signals with Sub- ...](https://github.com/LCAV/SparseFHT) - In this paper, we design a new iterative low-complexity algorithm for computing the Walsh-Hadamard t...

36. [Exact rule learning via Boolean compressed sensing | Proceedings of the 30th International Conference on International Conference on Machine Learning - Volume 28](https://dl.acm.org/doi/10.5555/3042817.3043022)

37. [Perturb-Seq: Dissecting Molecular Circuits with Scalable ...](https://www.sciencedirect.com/science/article/pii/S0092867416316105) - by A Dixit · 2016 · Cited by 2301 — We demonstrate Perturb-seq by analyzing 200,000 cells in immune ...

38. [Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens - PubMed](https://pubmed.ncbi.nlm.nih.gov/27984732/) - Genetic screens help infer gene function in mammalian cells, but it has remained difficult to assay ...

39. [Automated Circuit Discovery for Mechanistic Interpretability](https://erdogan.dev/sias.pdf) - Conmy, Arthur et al. (2023). “Towards Automated Circuit. Discovery for Mechanistic Interpretability”...

40. [Towards Automated Circuit Discovery for Mechanistic Interpretability](https://proceedings.neurips.cc/paper_files/paper/2023/hash/34e1dbe95d34d7ebaf99b9bcaeb5b2be-Abstract-Conference.html) - by A Conmy · 2023 · Cited by 809 — This paper systematizes the mechanistic interpretability process ...

41. [Automatic Circuit Discovery](https://sidn.baulab.info/autocircuits/) - Towards Automated Circuit Discovery for Mechanistic Interpretability (Conmy et al., 2023), proposes ...

42. [Attribution Patching Outperforms Automated Circuit Discovery](https://arxiv.org/abs/2310.10348) - by A Syed · 2023 · Cited by 205 — In this work, we show that a simple method based on attribution pa...

