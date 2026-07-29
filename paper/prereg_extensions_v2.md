# Extensions to strengthen paper — v2

Revised per review. A is replication (not pre-reg). B is the only
genuine pre-registration. C is dropped. D is a sanity check.

---

## Extension A: Replication on additional empirical landscapes

NOT a pre-registration. These are published datasets analyzed for
epistasis in the literature. We apply the same Walsh pipeline and
report what we find. Framing: "we apply the same decomposition to
N additional published combinatorially complete landscapes."

### Candidate datasets
| Dataset | n | 2^n | Phenotype | Citation |
|---|---|---|---|---|
| Franke Aspergillus niger | 8 | 256 | Sporulation | PLOS Genetics 2011 |
| Khan E. coli | 5 | 32 | Growth rate | Science 2011 |
| de Visser Aspergillus | 8 | 256 | Sporulation | Evolution 1997 |
| Tan HIV-1 protease | 5 | 32 | Drug resistance | Genetics 2012 |

Priority: Franke (n=8, different organism, richer spectrum) then Khan.

### Effort
~2-4 hours per dataset. Download, parse, run normalized_wht().

---

## Extension B: Asynchronous update dynamics (PRE-REGISTERED)

This is the one genuine pre-registration. Nobody has computed the
composition gap under async update for these networks.

### Scope
Implement async random-order update in the coalition sweep. For each
time step, update genes in a uniformly random permutation. Average
over 100 random update-order sequences per initial condition. Run
on 6 batch-1 networks: faure_cellcycle (n=10), tournier_apoptosis
(n=12), davidich_yeast (n=10), drosophila_cellcycle (n=14),
fanconi_anemia (n=15), arabidopsis_cellcycle (n=14).

### Pre-registered predictions

**B1. Sign preservation**
Prediction: gap sign preserved in >= 4/6 networks.
If >= 3/6 flip sign, the composition gap is update-scheme-dependent
and the paper's claims need heavy qualification. That would itself
be a finding worth reporting.

**B2. Magnitude**
Prediction: median |delta_async| < median |delta_sync|. Async
dynamics introduce update-order noise that smooths basin boundaries.

### What would falsify the main paper
>= 3/6 sign flips. This is the test with teeth.

### Effort
~4-6 hours (implement, validate, run). drosophila (n=14), fanconi
(n=15), arabidopsis (n=14) are expensive (2^14 or 2^15 coalitions
x 100 update sequences). May need Modal for those three.

---

## Extension C: Alternative epistasis measures — DROPPED

Walsh and ANOVA are both orthogonal decompositions of the same
function. Showing they agree is a sanity check, not a robustness
test. One sentence citing Proposition 1 does the job. Skip unless
a reviewer specifically asks.

---

## Extension D: Sensitivity analysis

Sanity check, not pre-registration. Two parameters:

### D1. Convergence with N (initial conditions)
Run batch-1 networks at N = 64, 128, 256, 512 (current), 1024.
Report convergence of delta_3+.
Expected: converges by N=256.

### D2. Null threshold robustness
Tabulate classification at thresholds 0.1, 0.25, 0.5 (current),
1.0, 2.0 pp across all 27 networks.
Expected: only networks near boundary (davidich, remy, mendoza)
change classification.

### Effort
D1: ~2 hours (rerun 6 networks at 5 N values).
D2: ~30 minutes (just re-threshold existing data).

---

## Execution order

1. **D2** (threshold robustness) — 30 min, no computation needed
2. **D1** (N convergence) — 2h, rerun batch-1 sweeps
3. **A** (empirical landscapes) — 2-4h per dataset
4. **B** (async) — 4-6h, SHA-freeze predictions first
