# Extension Results Status

## D2: Threshold Sensitivity — COMPLETE

Classification is stable across null thresholds from 0.1 to 2.0 pp.
24/27 (89%) networks maintain their classification at all thresholds.

Only 3 near-boundary networks change:
- davidich_yeast (delta = -0.38 pp): destruction → null above 0.5 pp
- hematopoiesis_aging (delta = +0.93 pp): creation → null above 1.0 pp
- mendoza_thelper (delta = -0.34 pp): destruction → null above 0.5 pp

At default 0.5 pp threshold: 18 creation, 6 destruction, 3 null.
Results: `results/sensitivity/threshold_sensitivity.json`

## D1: N Convergence — PARTIAL (3/6 networks complete)

All 3 completed networks show delta_3+ converges by N=256.

| Network (n) | N=64 | N=128 | N=256 | N=512 | N=1024 | Max range |
|---|---|---|---|---|---|---|
| faure (10) | -3.09 | -3.04 | -3.00 | -3.02 | -3.01 | 0.09 pp |
| davidich (10) | +0.14 | -0.22 | -0.33 | -0.38 | -0.55 | 0.69 pp |
| tournier (12) | +8.98 | +9.48 | +8.51 | +8.46 | +8.56 | 1.02 pp |

Large networks (drosophila n=14, arabidopsis n=14, fanconi n=15) running
but computationally expensive. Script: `scripts/sensitivity_n_init.py`
Results will save to: `results/sensitivity/n_init_convergence.json`

## A: Empirical Landscapes — COMPLETE (Franke 2011)

Extracted Franke et al. 2011 A. niger fitness landscape (PLOS Comp Bio).
8 loci, 186/256 genotypes measured (70 missing = non-viable).

Walsh spectrum order-3+ energy across all empirical landscapes:
| Dataset | Organism | n | Order 3+ |
|---|---|---|---|
| Weinreich TEM-1 (log MIC) | E. coli | 5 | 2.9% |
| Weinreich TEM-1 (raw MIC) | E. coli | 5 | 9.4% |
| Hall haploid growth | S. cerevisiae | 6 | 0.1% |
| Hall mating efficiency | S. cerevisiae | 6 | 4.7% |
| Franke growth rate | A. niger | 8 | 11.4% |

Data: `data/empirical_landscapes/franke2011_aspergillus.{csv,json}`
Figure: `paper/figures/fig_empirical_spectra_v2.{pdf,png}`
Results: `results/empirical/franke2011_spectrum.json`

Khan 2011 E. coli data not accessible (Science paywall).

## B: Async Update — RUNNING (small networks)

Running on n<=12 networks (faure, davidich, tournier) with 100 async
replicates per initial condition. Large networks need Modal.

Script: `scripts/async_update_sweep.py`
Results will save to: `results/async_update/`
