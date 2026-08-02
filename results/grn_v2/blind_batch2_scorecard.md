# Blind Batch 2 Combined Scorecard

## Prediction Provenance

- `blind_predictions_batch2.json`: SHA-256 `2248fdeff1cda59e...`
- `blind_predictions_batch2b.json`: SHA-256 `2d4129fee179a67b...`

## Predictions vs Actuals

| Model | n | Pred Dir | Actual Dir | Pred rho Range | Actual rho | Pred o3+ | Actual o3+ | rho Hit | Dir Hit |
|-------|---|---------|-----------|---------------|------------|----------|-----------|---------|---------|
| lambda_phage | 7 | destruction | destruction | [0.00, 0.40] | 0.606 | 0.0279 | 0.0215 | NO | YES |
| arellano_rootstem | 9 | destruction | creation | [0.07, 0.47] | 0.603 | 0.0111 | 0.5245 | NO | NO |
| asymmetric_cell_division | 9 | destruction | creation | [-0.05, 0.35] | -0.104 | 0.0010 | 0.4409 | NO | NO |
| cell_cycle_transcription | 9 | creation | creation | [-0.09, 0.31] | -0.383 | 0.0380 | 0.8445 | NO | YES |
| remy_p53_mdm2 | 10 | creation | creation | [-0.06, 0.34] | 0.505 | 0.0372 | 0.0123 | NO | YES |
| albert_segment_polarity | 11 | destruction | creation | [0.07, 0.47] | 0.209 | 0.0026 | 0.2028 | YES | NO |
| blood_stem_cell | 11 | creation | creation | [0.00, 0.40] | 0.595 | 0.0691 | 0.0805 | NO | YES |
| calzone_cellfate_reduced | 11 | destruction | creation | [0.05, 0.45] | -0.075 | 0.0162 | 0.3974 | NO | NO |
| li_budding_yeast | 11 | creation | creation | [0.03, 0.43] | 0.528 | 0.0661 | 0.0943 | NO | YES |
| myeloid_progenitors | 11 | destruction | destruction | [-0.03, 0.37] | 0.514 | 0.0145 | 0.0014 | NO | YES |
| pair_rule_module | 11 | creation | destruction | [0.00, 0.40] | 0.660 | 0.1145 | 0.0356 | NO | NO |
| emt_switch | 12 | creation | destruction | [0.00, 0.40] | 0.377 | 0.0788 | 0.0010 | YES | NO |
| morphogenetic_checkpoint | 12 | destruction | creation | [0.02, 0.42] | 0.400 | 0.0174 | 0.1480 | YES | NO |
| zanudo_tlgl | 12 | destruction | creation | [-0.04, 0.36] | -0.193 | 0.0017 | 0.3353 | NO | NO |
| lac_operon | 13 | destruction | creation | [0.13, 0.53] | 0.079 | 0.0073 | 0.4866 | NO | NO |
| mendoza_thelper | 13 | creation | destruction | [-0.00, 0.40] | 0.215 | 0.0379 | 0.0085 | YES | NO |
| saadatpour_guardcell | 13 | destruction | creation | [-0.00, 0.40] | 0.065 | 0.0015 | 0.5164 | YES | NO |
| fumia_cellcycle | 14 | creation | creation | [0.03, 0.43] | 0.266 | 0.0872 | 0.5672 | YES | YES |
| hematopoiesis_aging | 15 | creation | creation | [-0.08, 0.32] | 0.631 | 0.0421 | 0.0241 | NO | YES |
| irons_cardiac | 15 | creation | creation | [-0.05, 0.35] | 0.634 | 0.0505 | 0.1296 | NO | YES |
| calzone_cell_fate | 17 | creation | creation | [0.02, 0.42] | -0.016 | 0.0432 | 0.5976 | NO | YES |

## Summary Statistics

- **Models scored**: 21
- **Direction hits**: 10/21 (48%)
- **Rho in range**: 6/21 (29%)
- **Mean rho prediction error**: 0.291

## All Results (sorted by n)

| Model | n | rho | 95% CI | p-value | Global o3+ | Local o3+ | Delta | Direction | Cycling |
|-------|---|-----|--------|---------|-----------|----------|-------|-----------|---------|
| lambda_phage | 7 | 0.606 | [0.140, 0.851] | 3.61e-03 | 0.0215 | 0.0465 | -0.0250 | destruction | 16.0% |
| arellano_rootstem | 9 | 0.603 | [0.320, 0.783] | 9.82e-05 | 0.5245 | 0.0186 | +0.5059 | creation | 0.0% |
| asymmetric_cell_division | 9 | -0.104 | [-0.339, 0.184] | 5.47e-01 | 0.4409 | 0.0000 | +0.4409 | creation | 6.3% |
| cell_cycle_transcription | 9 | -0.383 | [-0.708, -0.037] | 2.13e-02 | 0.8445 | 0.0120 | +0.8325 | creation | 0.4% |
| remy_p53_mdm2 | 10 | 0.505 | [0.190, 0.775] | 4.09e-04 | 0.0123 | 0.0115 | +0.0008 | creation | 13.4% |
| albert_segment_polarity | 11 | 0.209 | [-0.094, 0.498] | 1.26e-01 | 0.2028 | 0.0044 | +0.1984 | creation | 0.0% |
| blood_stem_cell | 11 | 0.595 | [0.381, 0.750] | 1.67e-06 | 0.0805 | 0.0327 | +0.0478 | creation | 24.5% |
| calzone_cellfate_reduced | 11 | -0.075 | [-0.328, 0.175] | 5.86e-01 | 0.3974 | 0.0271 | +0.3704 | creation | 11.1% |
| li_budding_yeast | 11 | 0.528 | [0.319, 0.692] | 3.44e-05 | 0.0943 | 0.0308 | +0.0635 | creation | 12.5% |
| myeloid_progenitors | 11 | 0.514 | [0.278, 0.704] | 5.87e-05 | 0.0014 | 0.0242 | -0.0228 | destruction | 6.6% |
| pair_rule_module | 11 | 0.660 | [0.455, 0.821] | 4.21e-08 | 0.0356 | 0.0630 | -0.0274 | destruction | 17.8% |
| emt_switch | 12 | 0.377 | [0.151, 0.575] | 1.80e-03 | 0.0010 | 0.0392 | -0.0382 | destruction | 0.1% |
| morphogenetic_checkpoint | 12 | 0.400 | [0.162, 0.602] | 8.63e-04 | 0.1480 | 0.0289 | +0.1191 | creation | 0.4% |
| zanudo_tlgl | 12 | -0.193 | [-0.382, 0.019] | 1.20e-01 | 0.3353 | 0.0029 | +0.3323 | creation | 10.3% |
| lac_operon | 13 | 0.079 | [-0.145, 0.289] | 4.92e-01 | 0.4866 | 0.0122 | +0.4743 | creation | 0.3% |
| mendoza_thelper | 13 | 0.215 | [-0.056, 0.449] | 5.90e-02 | 0.0085 | 0.0119 | -0.0034 | destruction | 1.7% |
| saadatpour_guardcell | 13 | 0.065 | [-0.122, 0.241] | 5.72e-01 | 0.5164 | 0.0025 | +0.5140 | creation | 5.7% |
| fumia_cellcycle | 14 | 0.266 | [0.067, 0.451] | 1.08e-02 | 0.5672 | 0.0448 | +0.5224 | creation | 12.5% |
| hematopoiesis_aging | 15 | 0.631 | [0.513, 0.728] | 5.55e-13 | 0.0241 | 0.0147 | +0.0093 | creation | 10.5% |
| irons_cardiac | 15 | 0.634 | [0.467, 0.770] | 3.68e-13 | 0.1296 | 0.0204 | +0.1092 | creation | 9.3% |
| calzone_cell_fate | 17 | -0.016 | [-0.193, 0.151] | 8.55e-01 | 0.5976 | 0.0155 | +0.5821 | creation | 0.5% |
