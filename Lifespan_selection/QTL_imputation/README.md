# 3. Genotype imputation based on offline segmentation and recombination map calculation

**Author:** Athanasios Toros
**Email:** [athanasiostoros@gmail.com](mailto:athanasiostoros@gmail.com)

These scripts implement genotype imputation for noisy F2 QTL genotype calls and construct fine-scale recombination maps using offline changepoint (PELT) segmentation of genotype calls in 100-kb bins.

**Note:** The scripts were refactored from the original Jupyter notebooks used for the analysis. Non-environment specific analysis parameters are unchanged from the original implementation and are provided in `config.yaml`. The original notebooks are available from the author upon request.

Interactive versions of the resulting recombination maps are available for the autosomes and chromosome X at:

`https://arnqvist.org/myQTL/recombination_mapping_rerun.html` (autosomes) and

`https://arnqvist.org/myQTL/chromosome_X_recombination_rate.html` (X).


## Methods

Each F2 chromosome is modelled as a noisy 1-D 0/1/2 signal (EE=0, EL=1, LL=2).
PELT changepoint segmentation (`ruptures`, L1 cost, `jump=1`, `penalty=2`) finds
recombination events; biologically-informed heuristics drop implausible
transitions (interference distance, EL-mediated transitions, double crossovers).
Each segment is imputed to its modal genotype. Recombination counts per region →
Kosambi map function → cumulative cM → linear interpolation → window smoothing,
cross-checked against linkage disequilibrium.

The X script reuses the same helpers with `chromosome_X=True`, which (a) disallows
recombination on the hemizygous paternal X (only transitions through
heterozygosity are valid) and (b) halves the X gamete pool in the map.


## Files

| File                                                    | Language | Role                                                                                                                                                                                                                                              |
| ------------------------------------------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `config.yaml`                                           | YAML     | Paths, analysis parameters, chromosome/scaffold layout, and pedigree-specific thresholds used by both analysis scripts.                                                                                                                           |
| `qtl_common.py`                                         | Python   | Shared functions for data import and formatting, PELT segmentation, recombination-event validation, Kosambi map construction, segment-modal imputation, and the `Marker` and `Segment` data models. Also provides `load_config()`.                |
| `01_imputation_and_recombination_map_autosomes.py`      | Python   | Analysis of autosomes (chromosomes 1–9). Constructs separate recombination maps for the AB and CD F2 pedigrees, generates the combined map, performs the LD comparison using `scikit-allel`, and performs segmentation-based genotype imputation. |
| `02_imputation_and_recombination_map_chrX_scaffolds.py` | Python   | Analysis of chromosome X (chromosome 10) and scaffolds (chromosomes 11–54). Performs X-specific imputation and mapping, segmentation-based imputation of scaffold 12, and modal-genotype imputation of smaller scaffolds.                         |

## Configuration

Both analysis scripts call `qtl_common.load_config()`, which reads `config.yaml` from the script directory and exposes its contents as a nested configuration namespace, for example:

`cfg.segmentation.penalty`
`cfg.autosomes.input_calls`

Paths in the distributed configuration use generic placeholders such as `data/` and `results/`. Analysis parameter values correspond to those used for the paper and are unchanged from the original notebooks.
