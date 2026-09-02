# 7. Multi-evidence GO term ranking

The following script ranks GO biological-process terms by combined selection + expression evidence,
separately per EE regime (Early/Late).

**Note:** This script was refactored from the original Jupyter notebooks used for the analysis. Non-environment specific analysis parameters are unchanged from the original implementation and are provided in `config.yaml`. The original notebooks are available from the author upon request.

## Method summary

For each GO term and evidence layer, a hypergeometric Z-score is computed
(negative scores capped to 0). Within each evidence category (Selection,
Expression/DEG), Z-scores are aggregated with a normalized p-norm (p = 3).

## Files

| File | Language | Role |
|---|---|---|
| `config.yaml` | YAML | All input/output paths, evidence-layer column groups, and parameters (p-norm power, permutation count, alpha, GO size filter, seed). Loaded by the script. |
| `01_integrated_go_ranking.py` | Python | Incidence matrix -> Z-scores -> p-norm integration -> permutation null -> empirical p-values -> per-category BH FDR -> Early/Late ranked tables. |
