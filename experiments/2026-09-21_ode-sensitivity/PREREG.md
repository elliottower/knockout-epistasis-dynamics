# Does class preservation depend on the Hill coefficient, the threshold, or the syntax of the conversion?

**Status:** FROZEN at `78071efa918c`
**Plan sha256:** `7a094857a455332c4fb6cc5884939f214c28cf0aa71fdd6d585e1cdf7becc907`
**Frozen:** 2026-09-22
**Log:** 3 entries, head `ec381430`

Sections use the [OSF Preregistration](https://osf.io/prereg/) question titles verbatim, so
this maps onto a registration without being rewritten. A question that does not apply is
answered **N/A** with the reason, never deleted.

## Research questions or hypotheses

The primary arm (`experiments/2026-09-21_ode-primary-rerun/PREREG.md`) tests class preservation
at one setting of the conversion: normalized HillCube, n_H = 10, K = 0.5. This arm asks whether
the answer depends on that setting, on the 17 networks with at most 12 nodes. There are five
settings:

| setting | construction | n_H | K |
|---|---|---|---|
| a | normalized HillCube | 2 | 0.5 |
| b | normalized HillCube | 4 | 0.5 |
| c | normalized HillCube | 10 | 0.3 |
| d | normalized HillCube | 10 | 0.7 |
| e | normalized operator-wise | 10 | 0.5 |

The reference setting, HillCube at n_H = 10 and K = 0.5, is the primary arm's run on these 17
networks. Classes, the null band, the estimators and their admissibility are as in
`experiments/ode_v2/CONTEXT.md`. All Δ3+ values are all-pairs, the primary estimator, unless
stated otherwise, and non-null means non-null under the Boolean all-pairs estimator.

**S1.** At each HillCube setting a–d, every scored non-null network among the 17 has the same
class as under Boolean dynamics.

**S2.** Under the normalized operator-wise conversion (e), every scored non-null network among
the 17 has the same class as under Boolean dynamics.

S1 carries the design. It decides whether class preservation is a property of the Boolean
rules or of one choice of Hill coefficient and threshold. S1 holding would strengthen the
primary arm's result for these networks. S2 asks the same question of the conversion's syntax,
since the operator-wise and HillCube conversions differ in 97 rules. If S2 fails where S1
holds, class preservation depends on how the rules are written.

n_H = 1 is not a setting. The arm tests switch-like interpolations of Boolean rules, and a Hill
function with n_H = 1 is z(1 + K)/(z + K), which is not switch-like. The decision rests on that
estimand, although calibration also showed that some trajectories at n_H = 1 had not settled by
t = 960 (`results/audit/classifier_calibration_rtol1e-8.json`).

## Foreknowledge of data or evidence

Common to all three registrations: `experiments/ode_v2/CONTEXT.md`. No Δ3+ exists at any of
settings a–e, and the paper reports no Hill-coefficient sweep. Calibration integrated
trajectories at n_H = 1, 2 and 10 on five networks and computed no Walsh quantity. It showed
that damped oscillations near a Hopf point settle slowly at n_H = 2.

## Explanation of foreknowledge and managing unintended influences

The settings were fixed before any sensitivity run. Settings c and d bracket the primary K by
0.2 on each side. Settings a and b are the coefficients at which calibration found the engine
accepts every trajectory. The engine, configuration and networks are frozen by digest, and
`scripts/analyze_ode_rerun.py` is frozen with this file.

## Study type

Computational experiment: deterministic simulation of fixed, published network models.

## Intention for causal interpretation

Within the models only: coalition values are outcomes of simulated knockouts. No claim about the
organisms is intended.

## Blinding of experimental treatments

N/A — no participants or assigned treatments; every run is a deterministic computation.

## Additional blinding during research or analysis

No sensitivity Δ3+ is inspected until `scripts/analyze_ode_rerun.py sensitivity` has run on
every record at every setting.

## Study design

Within-network comparison across the reference and five alternative settings. The clamp value
is 0 and there are 32 initial states, seed 42, with the frozen solver and classifier. The
networks are the 17 of the paper's 28 with at most 12 nodes: `albert_segment_polarity`,
`arellano_rootstem`, `asymmetric_cell_division`, `blood_stem_cell`,
`calzone_cellfate_reduced`, `cell_cycle_transcription`, `davidich_yeast`, `emt_switch`,
`faure_cellcycle`, `lambda_phage`, `li_budding_yeast`, `morphogenetic_checkpoint`,
`myeloid_progenitors`, `pair_rule_module`, `remy_p53_mdm2`, `tournier_apoptosis` and
`zanudo_tlgl`. The unit of analysis is the network at a setting.

## Randomization

N/A — nothing is assigned; the initial states and the replicate sample are fixed as in the
primary arm.

## Data collection procedures

The freeze and every launch follow `experiments/ode_v2/CONTEXT.md` (Freeze and launch). A
launch runs only from a checkout of the tagged commit, after
`uv run python -m scripts.verify_freeze` passes there. Five Modal launches, one per setting, with
run names `sensitivity-a` to `sensitivity-e`, each with `--keep-states`; for example setting c:

```
uv run --with modal==1.4.3 modal run --detach -m scripts.modal_ode_arm::main \
    --run-name sensitivity-c --networks <the 17 networks> \
    --construction hillcube_normalized --hill-n 10 --hill-k 0.3 \
    --solver experiments/ode_v2/solver_primary.json \
    --classifier experiments/ode_v2/classifier_primary.json \
    --n-init 32 --seed 42 --keep-states
```

Each run's records, and nothing else, are copied from the Modal volume into
`results/<run name>/` next to this file, for example:

```
uv run --with modal==1.4.3 modal run -m scripts.modal_ode_arm::fetch --run-name sensitivity-c \
    --out experiments/2026-09-21_ode-sensitivity/results/sensitivity-c
```

## Data collection procedures - File upload

N/A — the procedure is the code named above, frozen by digest.

## Sample size

17 networks at five settings. Under the published Boolean classes, 15 of the 17 are non-null;
`davidich_yeast` and `remy_p53_mdm2` are in the null band. S1 and S2 are censuses of these
networks, and say nothing about the 11 larger ones or about settings between the ones run.

## Sample size rationale

A sweep of the 17 networks costs about 12.5 core-hours, against about 270 for all 28, because
the three largest networks account for most of the cost. Restricting the arm to at most 12 nodes
buys five settings.

## Starting and stopping rules

Each network at each setting is run once to completion, with the primary arm's resume and rerun
rules.

## Manipulated variables

The conversion's Hill coefficient, threshold and syntax, per the table above. Within each
setting, the coalition of genes left unclamped: all 2^n.

## Measured variables

As in the primary arm, at each setting.

## Measured variables - File upload

N/A — the record format is defined by `scripts/run_ode_arm.py`.

## Indices

Δ3+ and its class at each setting.

## Indices - File upload

N/A — defined in `scripts/walsh_estimators.py`.

## Statistical models

N/A — S1 and S2 are counts.

## Statistical models - File upload

N/A — `scripts/analyze_ode_rerun.py` (`sensitivity`), frozen with this file.

## Transformations

N/A — Δ3+ in percentage points as computed.

## Inference criteria

| hypothesis | holds when |
|---|---|
| S1 | at every setting a–d, every scored non-null network has its Boolean class |
| S2 | at setting e, every scored non-null network has its Boolean class |

**S1 and S2 are void at a setting where more than 2 of the 17 networks cannot be scored.**

If S1 holds, the paper may say that class preservation holds at n_H ∈ {2, 4, 10} and
K ∈ {0.3, 0.5, 0.7} on the 17 networks with at most 12 nodes. It does not extend the statement
to the 11 larger networks or to settings that were not run. If S1 fails, it names the settings
and networks where a class changes, and makes no robustness claim beyond the settings where
every class held. The paper reports S2's outcome either way. If S2 fails and S1 holds, it says
that class preservation depends on the conversion's syntax for the networks that change.

## Data inclusion and exclusion

The 17 networks at all five settings. A network is excluded at a setting only if it cannot be
scored, is withheld for a replicate mismatch there, or has an inadmissible spectrum under the
estimator in use (`experiments/ode_v2/CONTEXT.md`, Admissibility), which counts toward the void
rule as a network that cannot be scored. Every exclusion is reported with its cause. Null-band
networks are excluded from S1 and S2.

## Missing data

As in the primary arm: a network with a coalition lacking a value is not scored at that setting.
If every failure is unclassified, a labeled sensitivity analysis scores it with provisional
outputs.

## Other planned analysis

1. **Secondary estimators.** S1 and S2 recomputed with the split-half and plain estimators.
2. **Admissibility gate.** S1 and S2 recomputed from the raw energies with the total-energy gate
   at 1e-7 and at 1e-5, with the networks each excludes.
3. **Correlation.** Boolean against ODE Δ3+ at each setting (Pearson and Spearman), descriptive.
4. **Magnitude.** For each network, the range of Δ3+ across the six settings, descriptive.

Anything else is exploratory and labeled so.

## Context and additional information

**Maximum claim under this registration.** Confirmatory: whether class preservation among the
17 networks with at most 12 nodes survives Hill coefficients 2 and 4 and thresholds 0.3 and 0.7
(S1), and the operator-wise conversion (S2). No claim about the 11 larger networks, about
settings between or beyond the ones run, or about n_H = 1.

---

## Log

Append only. Never edit above the line.

The last column is what distinguishes an amendment from a deviation, so you do not have to
decide which word to use: `nothing run`, `no results seen`, `results not opened`, `results seen`.

```
2026-09-21  created                              nothing run
2026-09-22  frozen at 851b7a96154f                nothing run  ·09adf740
2026-09-22  frozen at 78071efa918c                nothing run  ·ec381430
```
