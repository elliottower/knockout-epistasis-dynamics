# Prediction 4 Result: Hall Yeast Composition Gap

## Prediction
- **Sign**: POSITIVE (creation), but small in magnitude
- **Confidence**: Medium-low
- **Reasoning**: Biosynthetic pathways are mostly independent locally; higher-order terms arise from whole-cell flux coupling and import competition.

## Result (by fitness component)

| Component             | Global o3+ | Sign at threshold 0.005 |
|-----------------------|-----------|------------------------|
| Haploid growth rate   | 0.0012    | null                   |
| Diploid growth rate   | 0.0019    | null                   |
| Mating efficiency     | 0.0469    | creation               |
| Sporulation efficiency| 0.0121    | creation               |

- **Primary component (haploid growth rate)**: 0.0012, below threshold — effectively null
- **Verdict**: **PARTIALLY CONFIRMED** — direction positive for all components, but primary fitness component is near-zero as predicted ("small"). Secondary components (mating, sporulation) show clear higher-order epistasis.

## Details

Haploid growth rate energy spectrum:
| Order | Energy fraction |
|-------|----------------|
| 0     | 0.9913         |
| 1     | 0.0072         |
| 2     | 0.0003         |
| 3     | 0.0005         |
| 4     | 0.0006         |
| 5     | 0.0002         |
| 6     | 0.0000         |

Growth rate is dominated by the intercept (99.1%) — most genotypes grow at similar rates regardless of auxotrophic marker combination. This is consistent with the prediction that biosynthetic knockouts have mainly additive effects on growth.

Mating efficiency spectrum shows richer interaction structure:
| Order | Energy fraction |
|-------|----------------|
| 0     | 0.8581         |
| 1     | 0.0630         |
| 2     | 0.0320         |
| 3     | 0.0187         |
| 4     | 0.0133         |
| 5     | 0.0093         |
| 6     | 0.0056         |

Mating efficiency has 4.7% order-3+ energy — auxotrophic markers interact through shared nutrient requirements during mating, a process that couples multiple biosynthetic deficiencies simultaneously.

## Interpretation
The prediction correctly identified the sign and expected small magnitude for growth rate. The surprise is how fitness-component-dependent higher-order epistasis is: the same genetic perturbations create negligible higher-order effects on growth but substantial effects on mating — consistent with mating imposing stronger multi-pathway coupling than vegetative growth.

## Data source
- Hall et al. (2010) J Heredity 101(suppl_1):S75-S84
- 6 biosynthetic gene knockouts, 64 genotypes, 4 fitness components
- Data via harmslab/notebooks-nonlinear-high-order-epistasis repository
- Results: results/empirical/hall_2010_walsh.json
