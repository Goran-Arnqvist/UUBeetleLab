# CNV calling pipeline

**Author:** Athanasios Toros
**Email:** [athanasiostoros@gmail.com](mailto:athanasiostoros@gmail.com)

The following scripts implement a depth-of-coverage (DOC) CNV-calling pipeline that contrasts Early (E) vs Late
(L) selection-regime pool-seq lines. 

**Note:** The scripts were refactored from the original Jupyter notebooks used for the analysis. Non-environment specific analysis parameters are unchanged from the original implementation and are provided in `config.yaml`. The original notebooks are available from the author upon request.


## Files

| File | Language | Role |
|---|---|---|
| `config.yaml` | YAML | All paths, parameters, chromosome accessions/lengths. Loaded by every script. |
| `calculate_depth_of_coverage.batch` | bash script | mosdepth per-sample coverage ]Illustrative, cluster-specific. |
| `cnv_pipeline.py` | Python | CNV calling pipeline |
| `cnv_segment_cbs.R` | R | Per-chromosome PSCBS Circular Binary Segmentation of the denoised ratio.
| `analyze_cnv_calls.py` | Python | Post-call analysis: per-regime gain/loss calling, CDS-gene intersection, GO-prep gene lists, one summary figure. |
| `cnv_pipeline_figures.py` | Python | Code for relevant plots. |
| `correlation_cnv_selective_sweeps.Rmd` | R Markdown | GenometriCorr spatial association between CNV calls and selective-sweep regions. |

