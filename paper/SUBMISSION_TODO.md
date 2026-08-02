# Submission TODO — composition gap paper

Target: Genetics (GSA), fallback J Theor Bio or J Math Bio

## Before submission (required)

- [ ] Fill Zenodo DOI placeholders in paper (2 locations: conclusion + data availability)
- [ ] Upload data release to Zenodo (coalition values, Walsh spectra, prereg artifacts)

## Strengthening (recommended for Genetics)

### Sensitivity analysis — COMPLETE
- [x] D2: Threshold robustness — 24/27 stable across 0.1-2.0 pp thresholds
- [x] D1: N convergence — 3/6 networks done, all converge by N=256. Large networks running.
- Results: `results/sensitivity/threshold_sensitivity.json`, `results/sensitivity/n_init_convergence.json`

### More empirical fitness landscapes — PARTIAL
- [x] Franke 2011 A. niger (n=8, 186/256 measured): 11.4% order-3+ (1.1% excluding lethality)
- [ ] Khan 2011 E. coli (n=5, 32 genotypes) — behind Science paywall
- [ ] Expanded TEM-1 (bioRxiv July 2025, 55K genotypes, 18 mutations)
- Data: `data/empirical_landscapes/franke2011_aspergillus.{csv,json}`
- Figure: `paper/figures/fig_empirical_spectra_v2.{pdf,png}` (5 panels)

### Asynchronous update — IN PROGRESS
- [ ] Small networks (n<=12) running locally: `scripts/async_update_sweep.py --small-only`
- [ ] Large networks (n=14-15) need Modal
- Pre-reg: sign preserved >= 4/6, magnitude decreases (prereg_extensions_v2.md SHA 55cba8a)

### Full ODE sweep (nice to have)
- [ ] Complete remaining 15/27 models on Modal (`modal run --detach scripts/modal_ode_sweep.py`)
- [ ] Shard calzone_cell_fate n=17 across 16 containers (`scripts/modal_ode_shard_calzone.py`)
- [ ] Download results from Modal volume, update Section 3.5 table
- [ ] Rebuild fig_ode_comparison.pdf with all 27 networks
- 12/27 already saved on Modal volume epistasis-bench-results/ode_full/
- 4 local pilots at results/grn_v2/ode_pilot/
- Paper currently reports 4-network pilot; upgrade to "all 27" when done

### Alternative epistasis measures — DROPPED
- Walsh and ANOVA are equivalent (Proposition 1). One sentence suffices.

## Formatting for Genetics (GSA)

- [ ] Check Genetics "Investigations" format requirements
- [ ] Line numbers (already have)
- [ ] Double spacing (already have)
- [ ] Reference style (already unsrtnat — check if Genetics wants different)
- [ ] Data availability statement format

## Figures

- [x] fig1_delta_ranked.pdf — composition gap ranked bar chart
- [x] fig3_spectrum_comparison.pdf — local vs global spectra (3 examples)
- [x] fig_ode_comparison.pdf — Boolean vs ODE (4 pilots)
- [x] fig_empirical_spectra.pdf — Weinreich + Hall Walsh spectra
- [x] fig_empirical_spectra_v2.pdf — All 5 empirical landscapes including Franke
- [ ] Reference fig_ode_comparison and fig_empirical_spectra_v2 in the tex
- [ ] D1 convergence plot (once results complete)
- [ ] D2 threshold sensitivity table (for supplementary)
