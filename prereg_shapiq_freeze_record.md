# Freeze record: shapiq budget sweep pre-registration

Document: prereg_shapiq_budget_sweep_v2.md
SHA-256: 2caaffff9590e906305eb4fc2c6ee1e8a797fcfec31d4bc86f6b7782b985a976
Frozen: 2026-07-28

## Version pins

- shapiq==1.6.0
- numpy==2.2.6
- scipy==1.15.3
- scikit-learn==1.6.1
- Python >= 3.12

## Verify

```bash
sha256sum prereg_shapiq_budget_sweep_v2.md
# Must match: 2caaffff9590e906305eb4fc2c6ee1e8a797fcfec31d4bc86f6b7782b985a976
```

## Pre-computed ceilings (not predictions, used as inputs)

Exact k-SII reconstruction ceilings via shapiq.ExactComputer v1.6.0:

| Circuit | k-SII order≤2 | k-SII order≤3 |
|---------|---------------|---------------|
| weight_ioi | 0.757317 | 0.868135 |
| canonical_ioi | 0.955018 | 0.988119 |
| random15 | 0.994730 | 0.999635 |
