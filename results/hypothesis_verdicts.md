# Hypothesis Verdicts — EpistasisBench Sweep v2

Results file: results/shapiq_sweep_v2.json
Harness version: 2.1-modal-auroc
Trials: 360, Completed: 360, Failed: 0
Runtime versions: {'shapiq': '1.6.0', 'numpy': '2.2.6', 'scipy': '1.15.3', 'scikit-learn': '1.6.1'}

---

# Shapiq Pre-registration (H1-H4)

## H1: Estimation efficiency (shapiq reaches 95% of ceiling at >=5% budget)


**Verdict: CONFIRMED**

- weight_ioi @ 5%: R²=0.7446, ceiling=0.757, efficiency=0.984 [PASS]
- weight_ioi @ 10%: R²=0.7516, ceiling=0.757, efficiency=0.993 [PASS]
- weight_ioi @ 20%: R²=0.7574, ceiling=0.757, efficiency=1.001 [PASS]
- weight_ioi @ 40%: R²=0.7574, ceiling=0.757, efficiency=1.001 [PASS]
- canonical_ioi @ 5%: R²=0.9514, ceiling=0.955, efficiency=0.996 [PASS]
- canonical_ioi @ 10%: R²=0.9534, ceiling=0.955, efficiency=0.998 [PASS]
- canonical_ioi @ 20%: R²=0.9544, ceiling=0.955, efficiency=0.999 [PASS]
- canonical_ioi @ 40%: R²=0.9552, ceiling=0.955, efficiency=1.000 [PASS]
- random15 @ 5%: R²=0.9941, ceiling=0.995, efficiency=0.999 [PASS]
- random15 @ 10%: R²=0.9945, ceiling=0.995, efficiency=0.999 [PASS]
- random15 @ 20%: R²=0.9946, ceiling=0.995, efficiency=1.000 [PASS]
- random15 @ 40%: R²=0.9947, ceiling=0.995, efficiency=1.000 [PASS]

## H2: Order-3 ceiling behavior per circuit


**Verdict: CONFIRMED**

- weight_ioi @ 5%: best_shapiq_o3=0.7645, irf=0.9268
- weight_ioi @ 10%: best_shapiq_o3=0.8297, irf=0.9369
- weight_ioi @ 20%: best_shapiq_o3=0.8600, irf=0.9335
- weight_ioi @ 40%: best_shapiq_o3=0.8684, irf=0.9261
- canonical_ioi @ 5%: best_shapiq_o3=0.9604, irf=0.8970
- canonical_ioi @ 10%: best_shapiq_o3=0.9793, irf=0.9015
- canonical_ioi @ 20%: best_shapiq_o3=0.9861, irf=0.9123
- canonical_ioi @ 40%: best_shapiq_o3=0.9877, irf=0.9111
- random15 @ 5%: best_shapiq_o3=0.9968, irf=0.9619
- random15 @ 10%: best_shapiq_o3=0.9987, irf=0.9641
- random15 @ 20%: best_shapiq_o3=0.9994, irf=0.9638
- random15 @ 40%: best_shapiq_o3=0.9996, irf=0.9348

weight_ioi order-3 vs order-2 gains:
  5%: o2=0.7446, o3=0.7645, gap=0.0198
  10%: o2=0.7516, o3=0.8297, gap=0.0782
  20%: o2=0.7578, o3=0.8600, gap=0.1022
  40%: o2=0.7574, o3=0.8684, gap=0.1110

Prediction: order-3 exceeds order-2 ceiling by >=5 R² points on weight_ioi: YES
canonical_ioi: all gaps <3 R² points: NO

## H3: Method ordering (KernelSHAP-IQ > SVARM-IQ > SHAPIQ-MC)


**Verdict: CONFIRMED**

- weight_ioi @ 3% o2: K=0.7322 SV=0.7249 SH=0.4242 OK
- weight_ioi @ 5% o2: K=0.7446 SV=0.7223 SH=0.5880 OK
- weight_ioi @ 10% o2: K=0.7516 SV=0.7427 SH=0.6871 OK
- canonical_ioi @ 3% o2: K=0.9479 SV=0.6185 SH=0.3243 OK
- canonical_ioi @ 5% o2: K=0.9514 SV=0.9194 SH=0.6362 OK
- canonical_ioi @ 10% o2: K=0.9534 SV=0.9421 SH=0.8163 OK
- random15 @ 3% o2: K=0.9937 SV=0.8358 SH=0.7313 OK
- random15 @ 5% o2: K=0.9941 SV=0.9531 SH=0.8635 OK
- random15 @ 10% o2: K=0.9945 SV=0.9781 SH=0.9408 OK

Low-to-mid budget: 0/9 reversals (qualifying cells with R²>0.3)
- weight_ioi @ 20% o2: K=0.7574 SV=0.7578 SH=0.7373 gap=0.0205
- weight_ioi @ 20% o3: K=0.8600 SV=0.8337 SH=0.7073 gap=0.1527
- weight_ioi @ 40% o2: K=0.7574 SV=0.7572 SH=0.7522 gap=0.0052
- weight_ioi @ 40% o3: K=0.8684 SV=0.8597 SH=0.8275 gap=0.0409
- canonical_ioi @ 20% o2: K=0.9544 SV=0.9503 SH=0.9165 gap=0.0379
- canonical_ioi @ 20% o3: K=0.9861 SV=0.9610 SH=0.6924 gap=0.2937
- canonical_ioi @ 40% o2: K=0.9552 SV=0.9542 SH=0.9459 gap=0.0093
- canonical_ioi @ 40% o3: K=0.9877 SV=0.9811 SH=0.9182 gap=0.0694
- random15 @ 20% o2: K=0.9946 SV=0.9894 SH=0.9803 gap=0.0142
- random15 @ 20% o3: K=0.9994 SV=0.9668 SH=0.8817 gap=0.1178
- random15 @ 40% o2: K=0.9947 SV=0.9931 SH=0.9901 gap=0.0045
- random15 @ 40% o3: K=0.9996 SV=0.9904 SH=0.9708 gap=0.0288
High budget convergence: 4/12 cells with gap > 0.05

## H4: Cross-circuit consistency of method ordering


**Verdict: CONFIRMED**

- 3% o2: consistent ('K', 'SV', 'SH')
- 5% o2: consistent ('K', 'SV', 'SH')
- 10% o2: consistent ('K', 'SV', 'SH')
- 20% o2: INCONSISTENT {'weight_ioi': ('SV', 'K', 'SH'), 'ioi': ('K', 'SV', 'SH'), 'random15': ('K', 'SV', 'SH')}
- 20% o3: consistent ('K', 'SV', 'SH')
- 40% o2: consistent ('K', 'SV', 'SH')
- 40% o3: consistent ('K', 'SV', 'SH')

1/7 inconsistent cells

---

# AUROC Addendum Pre-registration (HA1-HA7)

## HA1: LASSO achieves highest pairwise AUROC at >=10% budget


**Verdict: CONFIRMED**

- weight_ioi @ 10%: LASSO=0.9984, best_other=0.9283 (kernelshapiq_order2) [PASS]
- weight_ioi @ 20%: LASSO=0.9995, best_other=0.9508 (kernelshapiq_order2) [PASS]
- weight_ioi @ 40%: LASSO=1.0000, best_other=0.9733 (kernelshapiq_order2) [PASS]
- canonical_ioi @ 10%: LASSO=1.0000, best_other=0.9811 (kernelshapiq_order2) [PASS]
- canonical_ioi @ 20%: LASSO=1.0000, best_other=0.9849 (kernelshapiq_order2) [PASS]
- canonical_ioi @ 40%: LASSO=1.0000, best_other=0.9900 (kernelshapiq_order2) [PASS]

## HA2: Shapiq order-2 competitive with LASSO at low budgets


**Verdict: NOT FALSIFIED**

- weight_ioi @ 1%: LASSO=0.6941, best_shapiq_o2=0.7052, gap=-0.0111 [within 0.05]
- weight_ioi @ 2%: LASSO=0.9068, best_shapiq_o2=0.7863, gap=0.1206 [EXCEEDS 0.10]
- weight_ioi @ 3%: LASSO=0.9797, best_shapiq_o2=0.8278, gap=0.1519 [EXCEEDS 0.10]
- canonical_ioi @ 1%: LASSO=0.7941, best_shapiq_o2=0.8278, gap=-0.0337 [within 0.05]
- canonical_ioi @ 2%: LASSO=0.9730, best_shapiq_o2=0.9236, gap=0.0493 [within 0.05]
- canonical_ioi @ 3%: LASSO=0.9955, best_shapiq_o2=0.9392, gap=0.0562 [within 0.10]

## HA3: Per-circuit ordering of pairwise AUROC


**Verdict: CONFIRMED (random15 threshold FALSIFIED)**

- 5%: ioi=0.8365, weight_ioi=0.8270, random15=0.7661 OK
- 10%: ioi=0.8759, weight_ioi=0.8659, random15=0.8019 OK
- 20%: ioi=0.9138, weight_ioi=0.8982, random15=0.8499 OK
- 40%: ioi=0.9345, weight_ioi=0.9290, random15=0.8854 OK

random15 threshold check (all methods AUROC < 0.65 at 1-3%):
  VIOLATED: lasso_walsh @ 1% = 0.8197 > 0.75
  VIOLATED: kernelshapiq_order2 @ 1% = 0.8895 > 0.75
  VIOLATED: lasso_walsh @ 2% = 0.9927 > 0.75
  VIOLATED: kernelshapiq_order2 @ 2% = 0.9601 > 0.75
  VIOLATED: lasso_walsh @ 3% = 0.9987 > 0.75
  VIOLATED: kernelshapiq_order2 @ 3% = 0.9911 > 0.75

LASSO on random15 @ 1%: 0.8197

## HA5: Max-order effect on shapiq pairwise AUROC


**Verdict: HA5a=FALSIFIED, HA5b=FALSIFIED, HA5c=FALSIFIED**

### HA5a: weight_ioi, high budget (order-3 > order-2)
  kernelshapiq @ 10%: o2=0.9283, o3=0.9283 o2 wins
  shapiq_mc @ 10%: o2=0.7179, o3=0.7179 o2 wins
  svarmiq @ 10%: o2=0.8610, o3=0.8610 o2 wins
  kernelshapiq @ 20%: o2=0.9508, o3=0.9508 o2 wins
  shapiq_mc @ 20%: o2=0.8191, o3=0.8191 o2 wins
  svarmiq @ 20%: o2=0.9100, o3=0.9100 o2 wins
  kernelshapiq @ 40%: o2=0.9733, o3=0.9733 o2 wins
  shapiq_mc @ 40%: o2=0.9107, o3=0.9107 o2 wins
  svarmiq @ 40%: o2=0.9635, o3=0.9635 o2 wins

  HA5a violations: 9/9, verdict: FALSIFIED

### HA5b: weight_ioi, low budget (order-2 > order-3)
  kernelshapiq @ 1%: o2=0.7052, o3=0.6239 o2 wins
  shapiq_mc @ 1%: o2=0.5622, o3=0.5165 o2 wins
  svarmiq @ 1%: o2=0.5816, o3=0.5504 o2 wins
  kernelshapiq @ 2%: o2=0.7863, o3=0.7745 o2 wins
  shapiq_mc @ 2%: o2=0.5835, o3=0.5923 o3 wins
  svarmiq @ 2%: o2=0.6526, o3=0.6326 o2 wins
  kernelshapiq @ 3%: o2=0.8278, o3=0.8278 o3 wins
  shapiq_mc @ 3%: o2=0.6116, o3=0.6116 o3 wins
  svarmiq @ 3%: o2=0.7045, o3=0.7045 o3 wins

  HA5b violations: 4/9, verdict: FALSIFIED

### HA5c: canonical_ioi and random15 (order gap <= 0.03)
  VIOLATED: canonical_ioi kernelshapiq @ 1%: gap=0.0626
  VIOLATED: canonical_ioi svarmiq @ 2%: gap=0.0849
  VIOLATED: random15 kernelshapiq @ 1%: gap=0.0671

  HA5c violations: 3, verdict: FALSIFIED

## HA6: Shapiq inter-method ordering (K >= SV >= SH for AUROC)


**Verdict: CONFIRMED**

- REVERSAL: random15 @ 1% o2: K=0.8895 SV=0.4790 SH=0.5361
- REVERSAL: random15 @ 2% o2: K=0.9601 SV=0.5218 SH=0.5296

2/42 reversals in qualifying cells (threshold: >3 falsifies)

## HA7: iRF pairwise AUROC profile


**Verdict: CONFIRMED**

- weight_ioi @ 10%: LASSO=0.9984, iRF=0.8238 [LASSO wins]
- weight_ioi @ 20%: LASSO=0.9995, iRF=0.8115 [LASSO wins]
- weight_ioi @ 40%: LASSO=1.0000, iRF=0.7972 [LASSO wins]
- canonical_ioi @ 10%: LASSO=1.0000, iRF=0.7959 [LASSO wins]
- canonical_ioi @ 20%: LASSO=1.0000, iRF=0.7972 [LASSO wins]
- canonical_ioi @ 40%: LASSO=1.0000, iRF=0.7723 [LASSO wins]

iRF >= LASSO at high budget on both circuits: NO

iRF vs worst shapiq (order=2) at >=5%:
  weight_ioi @ 5%: iRF=0.8242, worst_shapiq_o2=0.6519 [iRF wins]
  weight_ioi @ 10%: iRF=0.8238, worst_shapiq_o2=0.7179 [iRF wins]
  weight_ioi @ 20%: iRF=0.8115, worst_shapiq_o2=0.8191 [shapiq wins]
  weight_ioi @ 40%: iRF=0.7972, worst_shapiq_o2=0.9107 [shapiq wins]
  canonical_ioi @ 5%: iRF=0.7935, worst_shapiq_o2=0.5893 [iRF wins]
  canonical_ioi @ 10%: iRF=0.7959, worst_shapiq_o2=0.6851 [iRF wins]
  canonical_ioi @ 20%: iRF=0.7972, worst_shapiq_o2=0.8237 [shapiq wins]
  canonical_ioi @ 40%: iRF=0.7723, worst_shapiq_o2=0.9300 [shapiq wins]

---

## Summary: median pairwise AUROC by method, circuit, budget


### weight_ioi

| Budget | lasso_walsh | irf | kernelshapiq_order2 | kernelshapiq_order3 | shapiq_mc_order2 | shapiq_mc_order3 | svarmiq_order2 | svarmiq_order3 |
|--------|------|------|------|------|------|------|------|------|
| 1% | 0.694 | 0.789 | 0.705 | 0.624 | 0.562 | 0.517 | 0.582 | 0.550 |
| 2% | 0.907 | 0.786 | 0.786 | 0.774 | 0.583 | 0.592 | 0.653 | 0.633 |
| 3% | 0.980 | 0.777 | 0.828 | 0.828 | 0.612 | 0.612 | 0.704 | 0.704 |
| 5% | 0.995 | 0.824 | 0.880 | 0.880 | 0.652 | 0.652 | 0.784 | 0.784 |
| 10% | 0.998 | 0.824 | 0.928 | 0.928 | 0.718 | 0.718 | 0.861 | 0.861 |
| 20% | 0.999 | 0.812 | 0.951 | 0.951 | 0.819 | 0.819 | 0.910 | 0.910 |
| 40% | 1.000 | 0.797 | 0.973 | 0.973 | 0.911 | 0.911 | 0.964 | 0.964 |

### canonical_ioi

| Budget | lasso_walsh | irf | kernelshapiq_order2 | kernelshapiq_order3 | shapiq_mc_order2 | shapiq_mc_order3 | svarmiq_order2 | svarmiq_order3 |
|--------|------|------|------|------|------|------|------|------|
| 1% | 0.794 | 0.732 | 0.828 | 0.765 | 0.490 | 0.493 | 0.507 | 0.525 |
| 2% | 0.973 | 0.788 | 0.924 | 0.910 | 0.540 | 0.519 | 0.550 | 0.634 |
| 3% | 0.995 | 0.786 | 0.939 | 0.939 | 0.555 | 0.555 | 0.743 | 0.743 |
| 5% | 0.999 | 0.794 | 0.969 | 0.969 | 0.589 | 0.589 | 0.831 | 0.831 |
| 10% | 1.000 | 0.796 | 0.981 | 0.981 | 0.685 | 0.685 | 0.917 | 0.917 |
| 20% | 1.000 | 0.797 | 0.985 | 0.985 | 0.824 | 0.824 | 0.963 | 0.963 |
| 40% | 1.000 | 0.772 | 0.990 | 0.990 | 0.930 | 0.930 | 0.980 | 0.980 |

### random15

| Budget | lasso_walsh | irf | kernelshapiq_order2 | kernelshapiq_order3 | shapiq_mc_order2 | shapiq_mc_order3 | svarmiq_order2 | svarmiq_order3 |
|--------|------|------|------|------|------|------|------|------|
| 1% | 0.820 | 0.600 | 0.890 | 0.822 | 0.536 | 0.520 | 0.479 | 0.524 |
| 2% | 0.993 | 0.597 | 0.960 | 0.968 | 0.530 | 0.505 | 0.522 | 0.531 |
| 3% | 0.999 | 0.599 | 0.991 | 0.991 | 0.541 | 0.541 | 0.598 | 0.598 |
| 5% | 1.000 | 0.598 | 0.996 | 0.996 | 0.577 | 0.577 | 0.659 | 0.659 |
| 10% | 1.000 | 0.599 | 0.998 | 0.998 | 0.636 | 0.636 | 0.776 | 0.776 |
| 20% | 1.000 | 0.600 | 0.999 | 0.999 | 0.775 | 0.775 | 0.876 | 0.876 |
| 40% | 1.000 | 0.591 | 0.999 | 0.999 | 0.879 | 0.879 | 0.958 | 0.958 |

## Summary: median held-out R² by method, circuit, budget


### weight_ioi

| Budget | lasso_walsh | irf | kernelshapiq_order2 | kernelshapiq_order3 | shapiq_mc_order2 | shapiq_mc_order3 | svarmiq_order2 | svarmiq_order3 |
|--------|------|------|------|------|------|------|------|------|
| 1% | 0.798 | 0.806 | 0.615 | -1.249 | -0.910 | -27.405 | -46.564 | -44.584 |
| 2% | 0.938 | 0.869 | 0.701 | -0.978 | -0.019 | -4.842 | -5.881 | -4.824 |
| 3% | 0.966 | 0.888 | 0.732 | 0.574 | 0.424 | -2.239 | 0.725 | -4.292 |
| 5% | 0.982 | 0.927 | 0.745 | 0.764 | 0.588 | -0.747 | 0.722 | -2.260 |
| 10% | 0.989 | 0.937 | 0.752 | 0.830 | 0.687 | 0.175 | 0.743 | 0.128 |
| 20% | 0.991 | 0.934 | 0.757 | 0.860 | 0.737 | 0.707 | 0.758 | 0.834 |
| 40% | 0.992 | 0.926 | 0.757 | 0.868 | 0.752 | 0.827 | 0.757 | 0.860 |

### canonical_ioi

| Budget | lasso_walsh | irf | kernelshapiq_order2 | kernelshapiq_order3 | shapiq_mc_order2 | shapiq_mc_order3 | svarmiq_order2 | svarmiq_order3 |
|--------|------|------|------|------|------|------|------|------|
| 1% | 0.884 | 0.821 | 0.912 | 0.572 | -2.869 | -56.900 | -125.855 | -115.560 |
| 2% | 0.971 | 0.869 | 0.937 | 0.473 | -0.758 | -10.670 | -19.400 | -11.583 |
| 3% | 0.989 | 0.880 | 0.948 | 0.913 | 0.324 | -5.449 | 0.619 | -11.345 |
| 5% | 0.996 | 0.897 | 0.951 | 0.960 | 0.636 | -2.276 | 0.919 | -6.364 |
| 10% | 0.998 | 0.901 | 0.953 | 0.979 | 0.816 | -0.388 | 0.942 | -0.710 |
| 20% | 0.999 | 0.912 | 0.954 | 0.986 | 0.917 | 0.692 | 0.950 | 0.961 |
| 40% | 0.999 | 0.911 | 0.955 | 0.988 | 0.946 | 0.918 | 0.954 | 0.981 |

### random15

| Budget | lasso_walsh | irf | kernelshapiq_order2 | kernelshapiq_order3 | shapiq_mc_order2 | shapiq_mc_order3 | svarmiq_order2 | svarmiq_order3 |
|--------|------|------|------|------|------|------|------|------|
| 1% | 0.985 | 0.943 | 0.988 | 0.926 | -0.470 | -24.995 | -38.421 | -34.830 |
| 2% | 0.998 | 0.958 | 0.992 | 0.941 | 0.330 | -3.573 | -5.663 | -4.924 |
| 3% | 0.999 | 0.959 | 0.994 | 0.992 | 0.731 | -1.505 | 0.836 | -4.679 |
| 5% | 1.000 | 0.962 | 0.994 | 0.997 | 0.863 | -0.297 | 0.953 | -2.301 |
| 10% | 1.000 | 0.964 | 0.994 | 0.999 | 0.941 | 0.456 | 0.978 | 0.232 |
| 20% | 1.000 | 0.964 | 0.995 | 0.999 | 0.980 | 0.882 | 0.989 | 0.967 |
| 40% | 1.000 | 0.935 | 0.995 | 1.000 | 0.990 | 0.971 | 0.993 | 0.990 |

---

## Verdict Summary

- **H1**: CONFIRMED
- **H2**: CONFIRMED
- **H3**: CONFIRMED
- **H4**: CONFIRMED
- **HA1**: CONFIRMED
- **HA2**: NOT FALSIFIED
- **HA3**: CONFIRMED (random15 threshold FALSIFIED)
- **HA5**: HA5a=FALSIFIED, HA5b=FALSIFIED, HA5c=FALSIFIED
- **HA6**: CONFIRMED
- **HA7**: CONFIRMED