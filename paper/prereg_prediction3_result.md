# Prediction 3 Result: Weinreich TEM-1 Composition Gap

## Prediction
- **Sign**: POSITIVE (creation of higher-order epistasis)
- **Confidence**: High
- **Reasoning**: Structural contacts capture pairwise interactions; fitness integrates stability, catalysis, and dynamics into higher-order dependencies.

## Result
- **Global o3+ (log MIC)**: 0.0286
- **Sign**: POSITIVE (creation) — threshold 0.005
- **Verdict**: **CONFIRMED**

## Details

Energy spectrum (log MIC, fraction by Walsh order):
| Order | Energy fraction |
|-------|----------------|
| 0     | 0.2079         |
| 1     | 0.7307         |
| 2     | 0.0328         |
| 3     | 0.0192         |
| 4     | 0.0094         |
| 5     | 0.0000         |

Local model: purely additive (no structural contacts between any mutation pair at 8A Ca-Ca distance; closest pair E104K-G238S at 12.67A).

Local o3+ = 0.0 (additive model has no higher-order terms).

Composition gap = global o3+ - local o3+ = 0.0286.

The landscape is dominated by main effects (73.1% of energy) as expected for mutations spread across a protein, but 2.9% of total energy resides in order-3+ interactions. All of this higher-order structure is "created" by the protein's global functional integration (catalysis, stability, folding) rather than direct pairwise structural contacts.

Raw MIC spectrum shows a similar pattern (global o3+ also positive) but log transform is the appropriate phenotype scale for MIC data.

## Data source
- Weinreich et al. (2006) Science 312:111-114
- 5 mutations in TEM-1 beta-lactamase, 32 genotypes
- Structural contacts from PDB 1BTL
- Results: results/empirical/weinreich_2006_walsh.json
