# QTL mapping of body weight and lifespan: figures and tables

This folder holds the code and data that regenerate QTL-mapping figures and tables in the manuscript and its Supplementary Information, starting from the raw genotype and phenotype file:

- **Main text:** Fig. 4, Fig. 5, Fig. 6
- **Supplementary:** Fig. S13–S27, Table S1, Table S2

The data come from an F2 intercross of *Acanthoscelides obtectus* seed beetles from the E and L selection lines: 1330 F2 individuals phenotyped for body weight (mg), lifespan (days) and sex.

By default the script loads the saved results of its three slow steps from `saved_state/` and takes about 1–2 minutes on a laptop. To run the whole analysis from the raw data instead, see *Saved state* below. Figures are written to `figures/` (PDF), tables to `tables/` (CSV), and the R session details to `sessionInfo.txt`.

## Requirements

```r
install.packages(c("qtl", "scales"))
```

## Folder contents

```
Submission/
├── make_figures.R        master script (the only script)
├── README.md
├── data/                 raw input
├── saved_state/          saved results of the slow steps, as used for the manuscript
├── figures/              output PDFs
├── tables/               output CSVs
```

### Raw input (`data/`)

`pheno_input.csv`: R/qtl "csv" format: weight (`wgt`), lifespan (`life`), `sex` and `pgm` for 1330 F2 individuals, and their genotypes at 3487 markers (100-kb windows) on chromosomes 1–9 and X. Genotype codes: E, L (homozygous), H (heterozygous), C (not L), D (not E); `-` is missing. 
`retained_markers.txt`: The 2839 markers kept after reducing runs of identical adjacent markers to one marker each (see *Cross preparation*).

### Saved state (`saved_state/`)

| `seedle_cross.Rds`: The prepared cross: data from `pheno_input.csv` reduced to the 2839 markers, with the estimated genetic map, genotype probabilities and the 128 genotype imputations that all results are based on.
`wgt_add_perm.Rds`, `life_add_perm.Rds`: 1000 permutations for the sex-averaged scans, with sex as a covariate (weight, lifespan). 
`wgt_male_perm.Rds`, `wgt_female_perm.Rds`, `life_male_perm.Rds`, `life_female_perm.Rds`: 1000 permutations for the scans run in each sex separately.
`mqm_wgt_int.Rds`, `mqm_life_int.Rds`: The fitted six-QTL models for weight and lifespan, including each QTL's LOD profile.

All permutations used `perm.Xsp = TRUE`, so the X chromosome has its own threshold.

## Methods in brief

**Notation.** Alleles are E and L. At each position the additive effect of the L allele is *a* = (mean of LL − mean of EE) / 2, and the dominance deviation is *d* = mean of EL − (mean of EE + mean of LL) / 2. Units are mg for weight and days for lifespan. Sex is coded 1 = male, 0 = female.

**Cross preparation** The raw file is read with `read.cross(format = "csv", genotypes = c("E","H","L","C","D"), alleles = c("E","L"), convertXdata = TRUE)`. Runs of adjacent markers with identical genotypes are found with `findDupMarkers(adjacent.only = TRUE)` and reduced to one marker each, leaving 2839 markers. `findDupMarkers` chooses at random which marker of a run to keep, so the markers kept in the manuscript analysis are listed in `retained_markers.txt`. The genotypes are the same whichever marker of a run is kept; only the marker name (window position) differs. The genetic map is estimated with `est.map(error.prob = 0.01)`. Genotype probabilities are computed with `calc.genoprob(step = 1, error.prob = 0.01)`, and 128 genotype imputations are drawn with `sim.geno(n.draws = 128, step = 1, error.prob = 0.01)`. Positions are in cM, on a 1-cM grid plus the markers.

**Genome scans** The scans use single-QTL interval mapping by the EM algorithm (`scanone`):
- **Sex-averaged scans** use both sexes, with sex as an additive covariate.
- **Per-sex scans** are run on males and females separately.

Genome-wide 5% LOD thresholds come from 1000 permutations, with X-chromosome thresholds estimated separately.

**Effect scans** `effectscan` estimates *a* and *d* at every position from the 128 imputations. It is run on both sexes together ("sex-averaged") and on each sex separately. On the X chromosome, R/qtl reports separate additive effects for males, which are hemizygous, and for females from the forward and reverse cross directions. Comparisons between males and females therefore use the autosomes only.

**Multiple-QTL models** QTL are added one at a time by forward selection, on the autosomes. The first QTL is the highest peak of a scan with sex as both an additive and an interactive covariate. At each later step, `addqtl` is run with sex as a covariate, and the highest peak is added if its LOD is above 3. After each addition, all positions are refined jointly with `refineqtl` (multiple imputation). Each QTL's interaction with sex is then tested in the final model (`addcovarint`), and a QTL × sex term is included where *p* < 0.05. 

## Saved state

Three steps are slow: preparing the cross (map estimation), the six sets of 1000 permutations, and the multiple-QTL model search. A complete run from the raw data took ~80 minutes on an 8-core laptop. The switch `USE_SAVED_STATE` at the top of `make_figures.R` loads the objects used for the manuscript from `saved_state/`. Everything else (genome scans, effect scans, figures, tables) is recomputed from them. This reproduces the manuscript exactly. If set to `FALSE`, the whole analysis runs from `data/pheno_input.csv` (about 80 minutes). The results of the slow steps are saved to `results/`, and later runs reuse them; delete `results/` to start over.

A run from scratch reproduces the cross exactly: same markers, same genetic map and same genotype probabilities, so the genome scans are identical too. Two steps involve random numbers and will differ slightly from the manuscript: the permutations (thresholds), and the genotype imputations, which feed the effect scans and the multiple-QTL models. In our test run from scratch, both multiple-QTL models came out with the same QTL at the same positions, every estimate in Tables S1–S2 agreed with the manuscript to within rounding, no significance call changed, and the permutation thresholds differed by at most 0.17 LOD without changing which peaks are significant.
