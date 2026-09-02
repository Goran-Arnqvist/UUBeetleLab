# 5. Functional annotation and GO enrichment

**Author:** Athanasios Toros
**Email:** [athanasiostoros@gmail.com](mailto:athanasiostoros@gmail.com)

The following scripts compile the *A. obtectus* GO annotation universe, packages it as an organism
annotation database, and runs over-representation analysis (ORA) on the genes
implicated by the **selective sweeps** and **CNVs**
for the Early (E) and Late (L) experimental-evolution regimes.

**Note:** The scripts were refactored from the original scripts, Jupyter and/or RMarkdown notebooks used for the analysis. Non-environment specific analysis parameters are unchanged from the original implementation and are provided in `config.yaml. The original notebooks are available from the author upon request.


## Files

| File | Language | Role |
|---|---|---|
| `config.yaml` | YAML | All input/output paths, the two chromosome-accession sets, and the data-prep parameters.|
| `01_build_go_annotation.py` | Python | GO universe (InterPro ∪ EggNOG → obsolete remap), per-gene chromosome + blacklist flags, and the sweep gene lists. |
| `02_build_org_db.R` | R | AnnotationForge db |
| `03_ora_selective_sweeps.R` | R | Sweep ORA|
| `04_overlap_stats_cnv_sweeps.py` | Python | Set-overlap stats (counts, %-in, Jaccard): consistently-selected sweep genes vs CNV deletion/gain gene lists, per regime.|
| `05_ora_cnvs.R` | R | CNV ORA |

