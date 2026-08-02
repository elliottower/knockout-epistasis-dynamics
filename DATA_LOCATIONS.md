# Data Locations

Where everything lives across repos, local disk, and Modal volumes.

## Repos

| Repo | Path | What it does |
|------|------|-------------|
| epistasis-bench | `~/Documents/GitHub/epistasis-bench` | Coalition sweeps (Boolean + ODE), canalization, composition gap analysis |
| epistasis-sample-budget | `~/Documents/GitHub/epistasis-sample-budget` | Sample-budget recovery experiment (subsample coalitions, fit methods, score AUROC) |
| factorization-circuits | `~/Documents/GitHub/factorization-circuits` | Main repo (paper lives here at `paper/`) |

## Coalition data (inputs to everything else)

All 28 networks have exhaustive coalition sweeps (2^n knockout combinations).

| Network | n | Coalition NPZ location |
|---------|---|----------------------|
| albert_segment_polarity | 8 | LOCAL `epistasis-bench/results/grn_v2/` |
| arabidopsis_cellcycle | 10 | LOCAL `epistasis-bench/results/grn_v2/` |
| arellano_rootstem | 10 | LOCAL `epistasis-bench/results/grn_v2/` |
| asymmetric_cell_division | 6 | LOCAL `epistasis-bench/results/grn_v2/` |
| blood_stem_cell | 11 | LOCAL `epistasis-bench/results/grn_v2/` |
| calzone_cell_fate | 17 | LOCAL `epistasis-bench/results/grn_v2/` |
| calzone_cellfate_reduced | 9 | LOCAL `epistasis-bench/results/grn_v2/` |
| cell_cycle_transcription | 7 | LOCAL `epistasis-bench/results/grn_v2/` |
| davidich_yeast | 12 | LOCAL `epistasis-bench/results/grn_v2/` |
| drosophila_cellcycle | 10 | LOCAL `epistasis-bench/results/grn_v2/` |
| emt_switch | 8 | LOCAL `epistasis-bench/results/grn_v2/` |
| fanconi_anemia | 15 | LOCAL `epistasis-bench/results/grn_v2/` |
| faure_cellcycle | 10 | LOCAL `epistasis-bench/results/grn_v2/` |
| fumia_cellcycle | 14 | LOCAL `epistasis-bench/results/grn_v2/` |
| grieco_bladder | 18 | MODAL VOLUME `epistasis-bench-results` at `grn_v2/grieco_bladder_coalition_blind.npz` |
| hematopoiesis_aging | 15 | LOCAL `epistasis-bench/results/grn_v2/` |
| irons_cardiac | 15 | LOCAL `epistasis-bench/results/grn_v2/` |
| lac_operon | 10 | LOCAL `epistasis-bench/results/grn_v2/` |
| lambda_phage | 4 | LOCAL `epistasis-bench/results/grn_v2/` |
| li_budding_yeast | 11 | LOCAL `epistasis-bench/results/grn_v2/` |
| mendoza_thelper | 12 | LOCAL `epistasis-bench/results/grn_v2/` |
| morphogenetic_checkpoint | 7 | LOCAL `epistasis-bench/results/grn_v2/` |
| myeloid_progenitors | 11 | LOCAL `epistasis-bench/results/grn_v2/` |
| pair_rule_module | 6 | LOCAL `epistasis-bench/results/grn_v2/` |
| remy_p53_mdm2 | 4 | LOCAL `epistasis-bench/results/grn_v2/` |
| saadatpour_guardcell | 13 | LOCAL `epistasis-bench/results/grn_v2/` |
| tournier_apoptosis | 13 | LOCAL `epistasis-bench/results/grn_v2/` |
| zanudo_tlgl | 9 | LOCAL `epistasis-bench/results/grn_v2/` |

Why grieco is different: n=18 means 2^18 = 262,144 coalitions. Too large to store locally (shards were 1.1GB). Merged on Modal and kept on the volume.

## Composition gap results (epistasis-bench)

| Data | Location |
|------|----------|
| Boolean WHT + composition gap | LOCAL `epistasis-bench/results/grn_v2/` (per-network JSON) |
| ODE robustness results | LOCAL `epistasis-bench/results/ode_full/` (per-network JSON) |
| Canalization correlations (28 networks) | LOCAL `epistasis-bench/results/canalization/canalization_correlations_28.json` |
| Grieco analysis (delta, ODE, cycling) | LOCAL `epistasis-bench/results/grn_v2/grieco_bladder_analysis.json` AND MODAL VOLUME `epistasis-bench-results` at `grn_v2/grieco_bladder_analysis.json` |

## Sample-budget results (epistasis-sample-budget)

| Data | Location |
|------|----------|
| 22/28 completed budget results | LOCAL `epistasis-sample-budget/results/full_sweep/` (per-network JSON) |
| 22/28 budget results (backup) | MODAL VOLUME `sample-budget-results` at `budget_sweep/` |
| 22-network summary | LOCAL `epistasis-sample-budget/results/full_sweep/summary.json` |
| 6 remaining (in progress) | MODAL VOLUME `sample-budget-results` at `budget_sweep/` (will appear as they complete) |

Missing 6 (launched 2026-08-01):
- fumia_cellcycle (n=14) — script: `modal_remaining_5.py`
- fanconi_anemia (n=15) — script: `modal_remaining_5.py`
- hematopoiesis_aging (n=15) — script: `modal_remaining_5.py`
- irons_cardiac (n=15) — script: `modal_remaining_5.py`
- calzone_cell_fate (n=17) — script: `modal_remaining_5.py`
- grieco_bladder (n=18) — script: `modal_grieco_sweep.py`

## Modal volumes

| Volume name | What's in it |
|-------------|-------------|
| `epistasis-bench-results` | Coalition NPZs (grieco only), ODE shards, Boolean shards, canalization, shapiq sweeps |
| `sample-budget-results` | Budget sweep per-network JSONs (22 complete + 6 in progress) |

## Paper

| File | Location |
|------|----------|
| Current draft | `epistasis-bench/paper/composition_gap_v11.tex` |
| References | `epistasis-bench/paper/refs.bib` |
| Submission checklist | `epistasis-bench/paper/SUBMISSION_TODO.md` |
