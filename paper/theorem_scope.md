# Theorem Scope: Feedback Necessity for Dynamical Composition Gap

## The Decomposition

The composition gap Δ₃₊ = (global order-3+ energy) - (local order-3+ energy)
can be decomposed into two components:

**Algebraic component (Δ_alg):**
Higher-order Walsh terms created by iterated substitution of local rules along
paths in the regulatory graph. Present even in DAGs. Determined entirely by
the local rules and the network topology.

**Dynamical component (Δ_dyn):**
Additional higher-order structure arising from attractor basin geometry —
the fact that knockout-dependent basin sizes weight different fixed points
and limit cycles differently. Requires feedback.

Δ₃₊ = Δ_alg + Δ_dyn

## Proposition 1 (DAG Case) — the easy theorem

**Statement:** Let G be a Boolean network whose regulatory graph Γ is acyclic
(after removing identity-function input nodes). Then:

(a) For every knockout set S and initial condition x₀, the synchronous dynamics
    converge to a fixed point within depth(Γ) steps.

(b) The fixed point is unique for each configuration of input nodes.

(c) No limit cycles exist.

(d) The value function v(S) = E_x₀[output at attractor] depends only on the
    composed Boolean function obtained by propagating local rules through the DAG
    and averaging over input configurations. It does not depend on attractor
    basin structure.

(e) Therefore Δ_dyn = 0, and Δ₃₊ = Δ_alg.

**Proof sketch for (a):**
Label nodes in topological order 1, ..., n. At step t=1 under synchronous
update, node 1 (source) reaches its fixed value (its rule depends only on
inputs and constants). At step t=2, node 2 sees node 1's correct value and
reaches its fixed value. After depth(Γ) steps, all nodes have stabilized.
Since the state is self-consistent, it is a fixed point. ∎

**Proof sketch for (b):**
In a DAG, each node's steady-state value is a deterministic function of the
source values (which are either clamped or input nodes). Given fixed input
configuration, the fixed point is unique. ∎

## Proposition 2 (DAGs still create higher-order epistasis) — the counterexample

**Statement:** There exist acyclic Boolean networks with only pairwise local
rules (each f_i has Walsh order ≤ 2) whose composed value function has
nonzero order-k Walsh terms for k > 2.

**Proof:** Explicit construction.
Let n = 4 with nodes x₁, x₂, h, y arranged as a depth-2 DAG:
- h = x₁ AND x₂  (order-2 local rule)
- y = h AND x₃   (order-2 local rule)
Output = y.

The composed function is y = (x₁ AND x₂) AND x₃ = x₁ AND x₂ AND x₃.

In ±1 encoding: y = ½(1 + z₁)(1 + z₂)(1 + z₃) where z_i = 1 - 2x_i.
Expanding: this has a z₁z₂z₃ term, which is Walsh order 3.

All local rules had order ≤ 2, but composition created order 3. ∎

## Proposition 3 (Feedback enables dynamical gap) — the distinguishing claim

**Statement:** For Boolean networks with feedback loops, the composition gap
can depend on attractor basin structure in ways that cannot be predicted from
local rules and topology alone.

**Proof approach:**
Construct two networks with identical local rules and identical regulatory
graph topology, but different knockout-dependent attractor structures (due to
different update orderings or parameter choices that don't change the graph).
Show they have different composition gaps. This proves that the dynamical
component Δ_dyn ≠ 0 in general for feedback networks.

**Alternative (empirical) approach:**
For each of the 27 networks:
1. Compute the "DAG skeleton" by removing back-edges (edges that create cycles)
2. Compute Δ_alg by running the value function on the acyclic skeleton
3. Δ_dyn = Δ₃₊ - Δ_alg
4. Test: does Δ_dyn correlate with feedback loop count?

## What needs to be computed

### Already have:
- Δ₃₊ for all 27 networks (merged_all_27_analysis.json)
- Feedback loop counts for 6 networks (structural_analysis.json)
- Cycling fraction for all 27 networks

### Need to compute:
1. **Feedback loop counts for all 27 networks** (extend structural_analysis.py)
   - Currently only 6 models have this
   - Need to load WEB_MODELS from run_batch2_blind_sweep.py + remaining models
2. **DAG skeleton computation** for each network
   - Remove minimum set of edges to make graph acyclic
   - Run value function on the acyclic version
   - Compute Δ_alg
3. **Correlation test**: feedback_loops vs |Δ_dyn| across all 27 networks

### Optional but strong:
4. **Depth × max-order bound**: In a depth-d DAG with max local Walsh order k,
   the maximum global Walsh order is bounded by min(n, d × k). Test empirically.

## What goes in the paper

### Theory section (~1 page):
- Propositions 1-2 with proofs (short, each ≤ half page)
- Proposition 3 stated with proof sketch
- Key insight sentence: "Composition of pairwise rules creates higher-order
  epistasis algebraically; feedback adds a dynamical component through
  attractor basin geometry."

### Empirical test (~1 figure):
- Scatter plot: feedback loop count vs Δ_dyn (or Δ₃₊ if Δ_alg is hard to compute)
- Alternative: cycling fraction (already have) vs Δ₃₊ — simpler but less clean

## Risks and mitigation

1. **Risk: Proposition 3 proof is hard.** Mitigation: The constructive proof
   (two networks, same topology, different basins) should work. If not, the
   empirical correlation is still publishable as a conjecture.

2. **Risk: DAG skeleton computation is ambiguous** (multiple valid minimum
   feedback arc sets). Mitigation: Use cycling fraction as a simpler proxy
   for "dynamical influence." Already have this data.

3. **Risk: Feedback count doesn't correlate.** Mitigation: Try other structural
   metrics — cycle length distribution, spectral radius of adjacency matrix,
   number of strongly connected components.

## Timeline estimate

- Propositions 1-2 proofs: 1 day (algebra only, no code)
- Extend feedback analysis to all 27 models: 1 day
- Correlation analysis: half day
- Write theory section in LaTeX: 1-2 days
- Total: ~4 days

## Decision: cycling fraction as proxy vs full DAG decomposition

The cycling fraction is already computed for all 27 models and captures
"how much dynamics matter" without needing the DAG skeleton computation.
The correlation between cycling fraction and Δ₃₊ is the simplest empirical
test. If it's strong, we don't need the full Δ_alg/Δ_dyn decomposition
for the paper — just cite it as future work.

Quick check needed: what's the correlation between cycling_fraction and
delta_o3plus across the 27 models?
