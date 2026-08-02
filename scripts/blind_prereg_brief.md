# Brief for structurally-blinded prediction + experiment agent

## What you are doing

You are writing structural predictions and running an experiment on
Boolean gene regulatory networks (GRNs). You have NOT seen any
experimental results. Your predictions come from analyzing the Boolean
rules (topology, feedback loops, AND-gate depth, pathway structure)
and NOTHING ELSE.

## The scientific question

Each gene's Boolean update rule is a function of its regulators.
The Walsh-Hadamard Transform (WHT) of each gene's truth table gives
exact, representation-independent Fourier coefficients at every
interaction order. These are the "local rule" interaction structure.

When you knock out genes (clamp their update rule to 0) and run
the full network dynamics to attractor convergence, you get a
different object: the "attractor-level" value function, whose WHT
gives the global interaction structure.

**The question**: How much of the local rule Fourier structure
survives dynamical composition into attractor-level epistasis? Does
composition create new interactions, destroy existing ones, or
reshuffle them? Does the answer depend on network topology?

## The method

For each network:
1. **Local Fourier extraction**: WHT of each gene's truth table
   (at most 2^6 = 64 entries). Gives exact local interaction
   coefficients at every order.
2. **Coalition sweep**: For each of 2^n coalitions S (subsets of
   genes), clamp genes not in S to 0 (knockout), run synchronous
   Boolean dynamics from 512 random initial states until attractor
   convergence, record output node value averaged over initial states.
3. **Global WHT**: Walsh-Hadamard transform of the 2^n-length
   value function gives global interaction coefficients.
4. **Composition scoring**: Spearman rank correlation between
   local |w_{i,j}| (summed across genes where pair co-occurs)
   and global |w_{i,j}| for all C(n,2) pairs. Bootstrap 95% CIs.
   Energy spectrum comparison (fraction of variance at each
   interaction order) for local rules vs global attractor.

## The models

The following models survived a three-part quality gate (unique
values >= 10, absolute order-3+ energy >= 1e-4, no single value
occupying > 90% of coalitions). Three models failed the gate and
are excluded (myeloid_progenitors, blood_stem_cell, emt_switch).

Surviving models are defined in `grn_coalition_sweep.py` in the
`BUILTIN_MODELS` dict:

1. faure_cellcycle (G1) — 10 nodes, Mammalian cell cycle
2. tournier_apoptosis (G2) — 12 nodes, Apoptosis signaling
3. davidich_yeast (G4) — 10 nodes, Fission yeast cell cycle
4. drosophila_cellcycle (G5) — 14 nodes (11 dynamic + 3 inputs), Drosophila cell cycle
5. fanconi_anemia (G9) — 15 nodes, Fanconi anemia and checkpoint recovery
6. arabidopsis_cellcycle (G10) — 14 nodes, Arabidopsis thaliana cell cycle

## What you must do

### Phase 1: Structural analysis and predictions

For EACH surviving model:

1. Read the Boolean rules from `grn_coalition_sweep.py`
2. Identify:
   - Output node and its rule
   - All feedback loops up to length 4 (classify as positive,
     negative, or mixed sign)
   - AND-gate arity (max number of inputs in any product term)
   - Pathway depth (longest non-feedback path from any input to output)
   - Competing pathways (if multiple independent paths reach the output)
   - Functionally redundant nodes (identical update rules)
   - Input nodes (self-referential rules f(x) = x)
3. Based on the structural analysis, predict:
   - **Energy spectrum**: What fraction of total energy at order-1,
     order-2, and order-3+? Give ranges (e.g., "30-45%").
   - **Spearman rho direction**: Will local pairwise structure
     positively correlate with global pairwise structure?
     (positive/zero/negative, with reasoning)
   - **Creation vs destruction**: Will dynamics create new
     higher-order interactions not present in any single rule
     (order-3+ energy increases), or destroy/compress them
     (order-3+ energy decreases)? Why?
4. For each prediction, state the structural reasoning (which
   specific feedback loops, AND gates, or pathway features drive
   the prediction).

### Phase 2: Cross-network hypothesis

After analyzing all models individually, state ONE cross-network
hypothesis about what structural feature predicts the composition
gap direction (creation vs destruction). Base this ONLY on the
structural analysis, not on any prior results. State:
- The predictor variable (computed from topology)
- The outcome variable (computed from experiment)
- The expected direction
- The structural reasoning

### Phase 3: Run the experiment

Run the composition experiment on ALL surviving models at n_init=512.
Use the existing code:

```bash
cd /Users/elliottower/Documents/GitHub/epistasis-bench
```

For each model:
1. Run `sweep_coalitions()` with n_init=512, max_steps=200, seed=42,
   update_scheme="sync", clamp_value=0
2. Run `score_composition()` on the result
3. Save coalition table (.npz) and composition scores (.json) to
   `results/grn_v2/`

### Phase 4: Compare predictions to results

For each model, compare your structural predictions to the actual
results. Report:
- Whether rho direction matched prediction
- Whether creation/destruction matched prediction  
- Whether energy spectrum fell within predicted ranges
- Any surprises and what they reveal about the structural analysis

For the cross-network hypothesis, test it on the full set.

### Phase 5: Write up

Save a complete report to:
`/Users/elliottower/Documents/GitHub/epistasis-bench/results/grn_v2/blind_prediction_report.md`

Structure:
1. Per-model structural analysis (rules, loops, predictions)
2. Cross-network hypothesis
3. Results table
4. Prediction accuracy
5. What the structural analysis got right and wrong

## Rules

- Use `uv run python` for all Python execution (not python3)
- Always save results to files, never rely on stdout
- Use tqdm for progress bars
- Use absolute paths
- Do NOT read any existing results files in results/grn_v2/ — those
  contain prior experimental results that would contaminate your
  predictions. Read only source code files.
- Do NOT read `drafts/composition_gap_draft_v2.md` — it contains
  interpreted results.
- Do NOT read `prereg_grn_v3.md` or `prereg_grn_v4.md` — they
  contain predictions influenced by prior results.
- You MAY read: `grn_coalition_sweep.py`, `composition_scorer.py`,
  `data_utils.py`, `scripts/run_composition_experiment.py`
