# EpistasisBench: Benchmarking Epistasis Detection Against Exact Ground Truth

## Pitch

Genetics has dozens of methods for detecting higher-order epistasis, but
head-to-head benchmarking against exact ground truth is rare. Combinatorially
complete fitness landscapes exist (Poelwijk et al. 2019 measured all 2^13
mutants of a fluorescent protein; Johnston et al. 2024 covered 4 residues),
but these studies validated their own sparse-recovery methods, not the
installed base of epistasis detectors. Nobody has scored Epi-MEIF, MDR,
BOOST, GAIN, Shapley interaction indices, and information-theoretic methods
against each other on the same exact answer key.

We have two independent exact ground truths from different domains:
1. **Neural-network circuits**: exhaustive 2^15 coalition tables over
   attention heads in GPT-2-small (32,768 coalitions x 512 prompts).
   The players are computational components, not sequence mutations.
2. **Protein fitness**: Poelwijk et al.'s published 2^13 landscape
   (8,192 measurements). Different domain, same mathematical structure.

The benchmark scores every method against both answer keys. If a method
ranks the same on a 15-head transformer circuit and a 13-site protein
landscape, that's cross-domain external validity -- the strongest possible
evidence that the method detects real structure rather than domain artifacts.

Transfer goes both ways: import genetics algorithms, benchmark where truth
is known, export back which ones actually work.

## Ground truth assets (from weight paper)

Three exhaustive 2^15 coalition tables exist on GPT-2-small IOI task.
A fourth (C2) is planned but not yet computed.

| ID | Circuit | Regime | Status | File |
|----|---------|--------|--------|------|
| C1 | Weight-derived IOI (15 heads) | Interactive (38% order-1, 8.3% order-3+) | DONE (zero-abl); RUNNING (mean-abl) | `weight_ioi_zero_v2_coalition_values.npz` |
| C3 | Canonical IOI, Wang et al. 2023 (15 heads) | Near-additive (72% order-1, 2.8% order-3+) | DONE (zero-abl); RUNNING (mean-abl) | `ioi_zero_v2_coalition_values.npz` |
| C4 | Random-15 baseline | Null (99% intercept) | DONE (zero-abl); RUNNING (mean-abl) | `random15_zero_v2_coalition_values.npz` |
| C2 | Activation-based IOI (EAP or ACDC, 15 heads) | Predicted near-additive | NOT STARTED | TBD |

Each table: 32,768 coalitions x 512 prompts x (target_logit, foil_logit).
All live in `weight-circuit-discovery/experiments_batch2/genetics/`.

Exact Walsh-Hadamard coefficients, exact Shapley values, exact Mobius
inversion -- all derivable from these tables.

### Data integrity note

A 4-coalition smoke-test stub (`weight_ioi_zero_v2_coalition_values.npz`,
184 KB) exists at the REPO ROOT with the same filename as the real 256 MB
file in `experiments_batch2/genetics/`. Any script resolving a relative path
from the wrong cwd will load (4, 512) data and produce plausible garbage.
All benchmark loaders must assert `target_logits.shape[0] == 2**n_players`.

### Ablation-type dependency

Zero-ablation inflates pairwise epistasis 10-50x vs mean ablation (measured
via BOOST). The benchmark answer key MUST use mean-ablation tables once
available. Scoring methods against zero-ablation truth is scoring against
an inflated answer key. EpistasisBench is blocked on mean-ablation data
for its interactive regime (C1).

## Methods inventory

### Already implemented and run (in weight-circuit-discovery repo)

| Method | Type | What it measures | Script | Results |
|--------|------|-----------------|--------|---------|
| BOOST | Pairwise exhaustive | Deviation from additivity | `boost_pairwise_screen.py` | All 3 circuits, zero/mean/resample |
| Walsh-Hadamard | Spectral | Orthogonal order decomposition | `mobius_wht.py`, `walsh_test.py` | All 3 circuits, v2 |
| MDR | Classification | Pairwise prediction power | `mdr_circuit_discovery.py` | All 3 circuits |
| EDCF | Spectral correlation | Cross-head spectral co-variation | `edcf_spectral_discovery.py` | All 3 circuits |
| SPRIGHT | Sparse recovery | Sublinear WHT approximation | `sparse_wht.py` | All 3 circuits |
| Cross-method agreement | Meta-analysis | Rank correlation + top-k overlap | `cross_method_agreement.py` | 4 methods x 3 circuits |

### Already implemented, not yet run on v2 data

| Method | Type | Script | Notes |
|--------|------|--------|-------|
| shapiq k-SII | Shapley interaction index | `analyze_shapiq.py` | ExactComputer wrapper written; shapiq==1.6.0 installed; prereg exists (`prereg_shapiq_interaction_spectrum.md`). Only k-SII index; STII/FSII are one-line changes (`index="STII"` etc). |

### Implemented in epistasis-bench (new)

| Method | Role | Script | Status |
|--------|------|--------|--------|
| LASSO-Walsh | **Matched sparse-recovery baseline** (not a competitor — the dashed line on Pareto plots) | `lasso_walsh_oracle.py` | TESTED — intercept excluded from L1 penalty (fit_intercept=True), scored with order-0-excluded NMSE. Alpha never hits grid boundary. |
| iRF (iterative Random Forest) | Competitor | `irf_interaction_discovery.py` | TESTED — Basu et al. PNAS 2018. Hyperparameters tuned via nested 3-fold CV (min_leaf, max_depth) on training set only. CV saturates at min_leaf=1, max_depth=None (40/40 trials), suggesting the gap is capacity type (axis-aligned splits vs global Walsh rotation), not capacity amount. |
| Budget sweep harness | Infrastructure | `budget_sweep_harness.py` | Paired seeds: both methods get identical train/heldout splits. Reports median + 95% interval. |

### To implement

| Method | Origin | Key paper | Priority | Notes |
|--------|--------|-----------|----------|-------|
| shapiq STII, FSII, FBII | Explainable ML | Muschalik et al. NeurIPS 2024 | HIGH | Change `index=` param in existing script. NOTE: these indices disagree by construction -- FSII/FBII optimize faithfulness, k-SII/STII enforce structural axioms, STII projects higher-order onto top order. Which index counts as "ground truth" is itself a benchmark design decision to defend. |
| shapiq approximators | Explainable ML | same | HIGH | KernelSHAP-IQ, SVARM-IQ etc. Built into shapiq library. THE budget-vs-accuracy experiment. |
| Faith-Shap | Explainable ML | Tsai et al. JMLR 2023 | HIGH | Available as FSII in shapiq |
| Epi-MEIF | Genetics | Saha et al. NAR 2022 | MEDIUM | Trivially runnable on 15 binary features; headline prior-art, needed for citation credibility. |

### Full sweep results (7 budgets × 20 seeds × 3 circuits)

Results in `results/full_sweep_v1.json` (zero-ablation, 2.5h runtime).
All six method-circuit combinations are monotone.

**weight-IOI (interactive, k99=241)**

| Budget | LASSO R² [95% CI] | iRF R² [95% CI] | Error ratio |
|--------|-------------------|-----------------|-------------|
| 1% | 0.806 [0.726, 0.862] | 0.796 [0.721, 0.832] | 1.1x |
| 2% | 0.940 [0.923, 0.953] | 0.869 [0.786, 0.910] | 2.2x |
| 3% | 0.966 [0.957, 0.972] | 0.895 [0.836, 0.931] | 3.1x |
| 5% | 0.982 [0.980, 0.984] | 0.926 [0.887, 0.940] | 4.2x |
| 10% | 0.989 [0.987, 0.989] | 0.938 [0.904, 0.961] | 5.4x |
| 20% | 0.991 [0.991, 0.992] | 0.945 [0.912, 0.958] | 6.1x |
| 40% | 0.992 [0.992, 0.992] | 0.950 [0.918, 0.961] | 6.3x |

**Canonical IOI (near-additive, k99=142)**

| Budget | LASSO R² [95% CI] | iRF R² [95% CI] | Error ratio |
|--------|-------------------|-----------------|-------------|
| 1% | 0.880 [0.842, 0.923] | 0.845 [0.725, 0.869] | 1.3x |
| 2% | 0.971 [0.962, 0.980] | 0.846 [0.800, 0.891] | 5.3x |
| 3% | 0.989 [0.985, 0.991] | 0.878 [0.794, 0.910] | 10.6x |
| 5% | 0.996 [0.996, 0.997] | 0.887 [0.837, 0.926] | 30.7x |
| 10% | 0.998 [0.998, 0.999] | 0.896 [0.835, 0.937] | 63.6x |
| 20% | 0.999 [0.999, 0.999] | 0.904 [0.870, 0.934] | 83.0x |
| 40% | 0.999 [0.999, 0.999] | 0.911 [0.880, 0.946] | 88.6x |

**Random-15 (null, k99=34)**

| Budget | LASSO R² [95% CI] | iRF R² [95% CI] | Error ratio |
|--------|-------------------|-----------------|-------------|
| 1% | 0.987 [0.979, 0.991] | 0.935 [0.612, 0.954] | 5.0x |
| 5% | 1.000 [1.000, 1.000] | 0.960 [0.669, 0.969] | 113.6x |
| 40% | 1.000 [1.000, 1.000] | 0.965 [0.416, 0.968] | 169.5x |

#### Key findings

1. **iRF has an absolute error floor near 0.05-0.11 that is invariant to
   budget, circuit type, and hyperparameters.** CV selects min_leaf=1,
   depth=None on 120/120 trials across all circuits. The floor is a property
   of axis-aligned splits (can only threshold individual features) vs the
   global Walsh rotation that LASSO parameterizes directly. LASSO's absolute
   error keeps dropping toward zero while iRF's stays fixed.

2. **Error ratios in the table are dominated by LASSO approaching zero,
   not iRF getting worse.** On canonical IOI, iRF's absolute error (0.09-0.11)
   is similar to weight-IOI (0.05-0.07), but LASSO reaches R²=0.999 on the
   easier function, so the ratio explodes to 89x. This is a property of
   division by a near-zero denominator — the same trap as the earlier
   NMSE-over-total-variance bug. The substantive finding is the invariant
   iRF error floor; the 89x and 170x numbers are consequences of that floor
   combined with LASSO's vanishing error, not independent findings.

3. **The prediction was instructively wrong.** We predicted the gap would
   narrow on canonical IOI because LASSO's advantage is in sparse higher-order
   recovery. The gap ratio widened because LASSO's denominator shrank faster
   on the easier function. The failure illustrates why absolute error, not
   error ratios, should frame the comparison.

4. **Per-order decomposition at 5% (weight-IOI):**
   - LASSO: order-1 NMSE=0.0006, order-2=0.006, order-3=0.036, order-4=0.20
   - iRF: order-1 NMSE=0.029, order-2=0.053, order-3=0.174, order-4=0.45
   - iRF places ~0.5% of total variance at orders 5-7 where true energy is
     <0.4%. The per-order NMSE exceeds 1 at those orders, but this is the
     denominator trap again (true energy ~10^-4). In absolute terms, the
     excess energy is negligible (~0.005 of total variance on weight-IOI,
     ~0.017 on canonical IOI).

5. **Random-15 iRF: rare reweighting instability under null signal.**
   Extended to 50 seeds at 5% budget. One seed (43) crashes (R²=0.43) while
   LASSO on the same seed is R²=0.9998 — crash rate 1/50 = 2%. Seed 43's
   iRF reweighting locks onto spurious feature correlations and amplifies
   them across iterations; more data deepens the attractor (train R² degrades
   0.92→0.71 from 1%→40% budget). The remaining 49 seeds form a right-skewed
   distribution (median 0.961, CI [0.910, 0.969]) with a continuous left tail
   (18% of seeds below 0.93), not bimodal. Depth selection (None vs 10) does
   not correlate with R² — both clusters select both depths in similar
   proportions. The apparent CI widening with budget in the original 20-seed
   run was entirely driven by seed 43's increasing influence with fewer total
   seeds at higher budgets (10 vs 20).

6. **Budget for R²≥0.90:** LASSO reaches 0.90 by 2% on all circuits.
   iRF never reaches 0.90 median on canonical IOI (0.911 at 40%).
   On weight-IOI, iRF reaches 0.90 at ~3% but with wide CIs.

### Cut from v1

| Method | Why |
|--------|-----|
| PID (Partial Information Decomposition) | Williams-Beer redundancy lattice is astronomically large for n=15. Tractable for 2-3 sources; doing it triple-wise produces local quantities incomparable to global Walsh spectrum. |
| DFIM | Gradient-on-input-features method; doesn't naturally apply to a coalition table over components. Forcing the analogy reads as padding. |
| GMDR / AntEpiSeeker / GAIN | Weak variants or heuristics that won't add a distinct capability. |

## Scoring protocol

### Metrics

1. **Pairwise recovery**: precision@k and recall@k of top pairwise
   interactions vs exact Walsh pairwise coefficients
2. **Order spectrum correlation**: Spearman rho between estimated and exact
   energy-by-order profile
3. **Order-3+ detection sensitivity**: AUROC for detecting heads with >30%
   order-3+ spectral energy
4. **Ranking agreement**: Spearman rho of full head-importance ranking vs
   exact Shapley values

### The compute-budget experiment

For each method, vary the sample budget from 1% to 100% of the coalition
table. At each budget level, score recovery accuracy against exact truth.
Plot accuracy vs compute Pareto frontier for all methods.

Poelwijk et al. 2019 demonstrated that their fluorescent-protein landscape
was extraordinarily sparse, enabling prediction from limited measurements.
They validated one sparse-recovery approach on one landscape. We extend
this to a systematic head-to-head comparison of the full method inventory
across two ground-truth regimes from different domains.

The existing SPRIGHT sparse recovery work (k_99, blind NMSE at 10% samples)
is the prototype for this experiment. If both landscapes are sparse (our
k_99 < 250 out of 32,768; Poelwijk reported extraordinary sparsity),
sparse-recovery methods should dominate the Pareto frontier -- and we can
say why, not just observe it.

Frame as **detection power under a sampling budget**, not scalability. At
n=15 the search space is trivially small for methods built for 500k SNPs.
The question is: given B evaluations instead of 2^n, how much true
structure does each method recover?

### Regime analysis

Test each method across regimes with different interaction structure:

**Neural-network circuits (GPT-2-small IOI):**
- Near-additive (C3: 72% order-1) -- easy, most methods should work
- Interactive (C1: 38% order-1, 8.3% order-3+) -- harder, tests higher-order detection
- Null (C4: 99% intercept) -- false positive control

**Protein fitness landscape (Poelwijk et al. 2019):**
- Fluorescent protein, n=13, all 2^13 measured
- Published data; binary genotype space like our head on/off
- Different domain, same Walsh-Hadamard decomposition applies
- Cross-domain agreement = external validity for method rankings

## Relation to weight paper

Weight paper: "circuits found by different methods have different interaction
structure." Ships first, uses the coalition tables as evidence for the
selective-pressure thesis.

EpistasisBench: "here are N ways to detect that interaction structure, scored
against exhaustive ground truth." Cites weight paper for the ground truth
data. Two papers, two theses, zero overlap.

## What exists vs what to build

### Reuse from weight-circuit-discovery (already done)
- Coalition table infrastructure (sweep scripts, Modal wrappers)
- BOOST, Walsh, MDR, EDCF, SPRIGHT implementations + results
- Cross-method agreement framework
- shapiq k-SII exact computation wrapper + prereg
- Synthetic validation suite (`synthetic_validation_gate.py`, `validate_a5_a6.py`)

### Build new
- Run shapiq k-SII on v2 tables (script exists, just needs execution)
- Add STII/FSII/FBII indices (one-line change per index)
- shapiq approximator budget sweep (built-in to shapiq library)
- Epi-MEIF implementation
- DFIM adaptation
- Information gain computation
- Unified scoring framework (all metrics, all methods, all regimes)
- Paper

## Dependencies and blockers

1. **Mean-ablation coalition tables** (currently running on Modal). Required
   for credible ground truth. ETA: ~25h per circuit.
2. **Weight paper ships first.** EpistasisBench cites it for ground truth
   provenance.
3. **Probability-scale spectrum.** Same tables, apply softmax before
   subtracting. Needed to confirm ground truth is not a scale artifact.
4. **Poelwijk 2019 data.** Need to locate and download the published 2^13
   fitness landscape. Check supplementary data at Nature Communications
   10:4213 or contact authors. Binary genotype format should map directly
   to the coalition-table structure.
5. **Shapley index choice.** Must decide (and defend) which interaction
   index is the "correct" decomposition for scoring purposes. Options:
   use Walsh-Hadamard coefficients as the ground truth (basis-independent,
   orthogonal, unique) and score Shapley indices as approximations of it;
   or report all indices and let the benchmark reveal where they disagree.

## Key references

See `background.md` for full annotated bibliography with 42 references.

Core papers:
- Muschalik et al. NeurIPS 2024 (shapiq library + benchmark)
- Saha et al. NAR 2022 (Epi-MEIF)
- Kundaje lab Bioinformatics 2018 (DFIM)
- Ritchie et al. 2001 (MDR)
- Wan et al. AJHG 2010 (BOOST)
- Poelwijk et al. Nat Commun 2019 (2^13 combinatorially complete landscape + sparsity)
- Johnston et al. PNAS 2024 (combinatorially complete enzyme landscape)
- Li et al. 2015 (SPRIGHT sparse WHT)
- Tsai et al. JMLR 2023 (Faith-Shap)
