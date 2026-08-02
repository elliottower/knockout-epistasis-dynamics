# Epistasis Methods Paper — Project Scope

## One-sentence thesis

28 exhaustive Boolean GRN decompositions provide exact ground truth for answering how many combinatorial measurements each interaction-detection method needs to recover true higher-order structure.

## Why this paper exists

Paper 1 (composition gap) establishes 28 exhaustive Walsh decompositions (n=7 to n=18, up to 262144 coalitions each) and proves dynamics create higher-order epistasis. This paper asks: given that we *know* the true interaction structure, how efficiently can detection methods recover it from partial samples?

Nobody else has 28 systems with exact, complete 2^n coalition tables spanning three orders of magnitude in combinatorial complexity. That's the comparative advantage — no simulation assumptions, no held-out-split approximations, just exact answers.

## Phasing

**Phase 1 (this paper):** Sample-budget recovery + continuous perturbation, purely on the 28 GRNs. Self-contained, all data in hand.

**Phase 2 (follow-up, if warranted):** Cross-domain validation on real drug data (Brennan, Cokol). Folds in if results are clean; becomes a separate paper if not.

## Two contributions (Phase 1)

### Contribution 1: Sample-budget recovery curves

**Question:** If you can only measure k% of the 2^n coalitions, how much of the true order-3+ structure can you recover?

**Method:** For each of the 28 GRNs, subsample the coalition table at 5%, 10%, 25%, 50%, and reconstruct the Walsh spectrum from the subsample. Score by correlation with the true (exhaustive) spectrum and by correctly classifying the sign of the composition gap.

**Why it matters:** Real combinatorial screens can never measure everything. This gives practitioners a concrete answer to "how many combos do I need?" grounded in exact truth, not simulation.

**Detection methods to compare (capped at 4):**
1. Sparse Walsh recovery (SPRIGHT or compressed sensing) — designed for this
2. Additive/linear baseline — the adversarial floor (Ahlmann-Eltze result)
3. Shapley/k-SII interaction indices (shapiq) — the game-theory family
4. Random-subset averaging — naive strawman

Key: score by held-out coalition prediction (can the method predict v(S) for unseen coalitions), NOT by self-referential Walsh correlation. Perplexity's self-critique was right about the circularity.

### Contribution 2: Continuous perturbation extension

**Question:** Does the composition gap survive when knockouts are graded (30%, 50%, 70% expression) rather than binary (on/off)?

**Method:** Replace the binary coalition mask with continuous knockdown levels. Use the ODE Hill-function framework (already built) with intermediate clamping values. Compute Sobol indices (the continuous analogue of Walsh-Hadamard). Compare order-3+ energy with the binary result.

**Possible outcomes:**
- Boring (gap persists, magnitudes scale smoothly) → robustness result, one section
- Surprising (gap magnitude depends nonlinearly on knockdown depth, or new interaction orders emerge) → headline finding, changes the paper's weight distribution

Start with 2-3 networks to find out which case it is before scaling.

### Contribution 3: Cross-domain validation on real drug data

**Question:** Do the sample-budget curves from Boolean GRNs predict recovery performance on real experimental data?

**Datasets (in build order):**
1. COMBSecretomics (n=3, pipeline shakedown, zero friction)
2. Brennan 8-antibiotic (n=8, 256 combos, public, primary validation)
3. Cokol/DiaMOND TB (n=5, 32 combos, needs author request)

**Method:** Run the same 4 detection methods on the drug data. Compare the recovery curves against the GRN predictions. If GRN curves predict drug-data performance, that's a transferability result. If they don't, that's informative too — biological vs. pharmacological interaction structure differs.

**Explicitly NOT doing:** Norman Perturb-seq manifold-CI. The order-3+ extension is undefined and solving that is its own paper. The CRISPR data (Zhou) is post-hoc validation, not core.

## What this paper is NOT

- Not a method leaderboard (4 methods, pre-registered comparisons, Holm-Bonferroni)
- Not an extension of MechVal (one-sentence citation in discussion, not framework import)
- Not a deep learning benchmark (GEARS/DeepDrugs are a follow-up paper if anyone cares)
- Not a replication of Ahlmann-Eltze ("additive wins" is the expected baseline result, not the headline)

## Noise-floor check (from Perplexity's self-critique, correctly identified)

For Brennan and any dataset with replicates: compute replicate-based noise ceiling. Report all method scores relative to this ceiling, not in absolute terms. A method scoring 0.6 correlation means nothing if the noise ceiling is 0.65.

## Implementation-variant check (from MechVal E1, correctly identified)

For each of the 4 methods, run at least 2 defensible implementation variants (e.g., different reference points for "no interaction") and report agreement. Prevents cherry-picking the favorable variant.

## Build order

1. Sample-budget sweep on 3 small GRNs (n=7-10) — validate pipeline
2. Continuous perturbation on 2-3 GRNs — determine boring vs. surprising
3. Scale sample-budget to all 28 GRNs
4. COMBSecretomics shakedown (n=3)
5. Brennan primary analysis (n=8)
6. Cokol if data arrives (parallel author request)
7. Write-up

## Relationship to other papers

- **Cites:** composition gap paper (Paper 1) for the 28 decompositions
- **Cited by:** MechVal (as operationalization of criterion I6)
- **Independent of:** direction instability papers (different question — cross-context transport vs. combination structure)
- **Explicitly differentiates from:** Ajmal et al. 2025 (validates against curated gold standards, not exhaustive truth), GRAPE/Synulator 2026 (simulated truth, not real), Ahlmann-Eltze 2025 (same finding expected but with exact ground truth and budget curves)
