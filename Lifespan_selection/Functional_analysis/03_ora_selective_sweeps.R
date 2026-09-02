#!/usr/bin/env Rscript

# Author: Athanasios Toros

# Email: athanasiostoros@gmail.com

# Over-representation analysis (ORA) of selective-sweep candidate genes.


suppressPackageStartupMessages({
library(tidyverse)
library(org.Aobtectus.eg.db)
library(GOSemSim)
library(clusterProfiler)
})

set.seed(42)

config <- yaml::read_yaml("config.yaml")

# ---------------------------------------------------------------------------

# Universe

# ---------------------------------------------------------------------------

gene_locs <- AnnotationDbi::select(
  org.Aobtectus.eg.db,
  keys    = keys(org.Aobtectus.eg.db, keytype = "GID"),
  columns = c("GID", "CHROMOSOME", "LOW_COVERAGE", "LOW_QUALITY", "HIGH_COVERAGE", "NEAR_GAP"),
  keytype = "GID"
)

# limit genes to only assembled scaffolds, excluding X (CAVLJG010000002.1)
included_chromosomes <- unlist(config$genome$chromosome_accessions_sweeps)

universe <- unique(gene_locs$GID[ gene_locs$CHROMOSOME %in% included_chromosomes ])

# additionally drop genes flagged for any quality reason (the gene lists are
# already filtered the same way)
cols_to_convert <- c("LOW_COVERAGE", "LOW_QUALITY", "HIGH_COVERAGE", "NEAR_GAP")

gene_locs <- gene_locs %>%
  mutate(across(all_of(cols_to_convert), ~ as.numeric(.)))

gene_locs_filtered <- gene_locs[rowSums(gene_locs[cols_to_convert]) == 0, ]

universe_filtered <- unique(gene_locs_filtered$GID[ gene_locs_filtered$CHROMOSOME %in% included_chromosomes ])

# ---------------------------------------------------------------------------

# ORA on consistently-selected genes 
# ---------------------------------------------------------------------------

consistently.selected.late <- read_tsv(
  config$outputs$consistently_selected_late,
  col_names = FALSE
)$X1

consistently.selected.early <- read_tsv(
  config$outputs$consistently_selected_early,
  col_names = FALSE
)$X1

consistently.selected.early.filtered = intersect(consistently.selected.early, universe_filtered)

consistently.selected.late.filtered = intersect(consistently.selected.late, universe_filtered)

## Stats

cat("========\nUnfiltered Universe\n========")

universe_size <- length(universe)

cat("Total genes in universe:", universe_size, "\n\n")

report_set_stats <- function(name, vec) {
  n <- length(vec)
  pct <- round(100 * n / universe_size, 2)
  cat(sprintf("%-30s: %5d genes (%.2f%% of universe)\n", name, n, pct))
}

cat("---- Stats for affected genes (filtered) ----\n")

report_set_stats("Early Consistently Selected Genes", consistently.selected.early)

report_set_stats("Late Consistently Selected Genes", consistently.selected.late)

report_set_stats("Genes common in both regimes:", (intersect(consistently.selected.early, consistently.selected.late)))


cat("\n\n========\nFiltered Universe\n========")

cat("---- Gene Set Filtering Statistics ----\n")

universe_size <- length(universe_filtered)

cat("Total genes in filtered universe:", universe_size, "\n\n")

report_set_stats <- function(name, vec) {
  n <- length(vec)
  pct <- round(100 * n / universe_size, 2)
  cat(sprintf("%-30s: %5d genes (%.2f%% of universe)\n", name, n, pct))
}

cat("---- Stats for affected genes (filtered) ----\n")

report_set_stats("Early Consistently Selected Genes", consistently.selected.early.filtered)

report_set_stats("Late Consistently Selected Genes", consistently.selected.late.filtered)

report_set_stats("Genes common in both regimes:", (intersect(consistently.selected.early.filtered, consistently.selected.late.filtered)))

### Biological Process 

ego.late.bp = enrichGO(gene         = consistently.selected.late.filtered,
OrgDb         =org.Aobtectus.eg.db,
keyType       = 'GID',
ont           = "BP",
pAdjustMethod = "BH",
universe = universe_filtered,
qvalueCutoff  = 0.05)

ego.early.bp = enrichGO(gene         = consistently.selected.early.filtered,
OrgDb         =org.Aobtectus.eg.db,
keyType       = 'GID',
ont           = "BP",
pAdjustMethod = "BH",
universe = universe_filtered,
qvalueCutoff  = 0.05)

semData <- godata(org.Aobtectus.eg.db, ont="BP", keytype = "GID")

ego.late.bp.slim <- simplify(ego.late.bp,
cutoff      = 0.7,
by          = "p.adjust",
select_fun  = min,
measure     = "Wang",
semData     = semData)

ego.early.bp.slim <- simplify(ego.early.bp,
cutoff      = 0.7,
by          = "p.adjust",
select_fun  = min,
measure     = "Wang",
semData     = semData)

# GO-term overlap stats, Early vs Late (BP)

late_ids  <- ego.late.bp@result$ID

early_ids <- ego.early.bp@result$ID

late_slim_ids  <- ego.late.bp.slim@result$ID

early_slim_ids <- ego.early.bp.slim@result$ID

common_before_slim <- intersect(late_ids, early_ids)

common_after_slim  <- intersect(late_slim_ids, early_slim_ids)

common_late_before  <- ego.late.bp[ego.late.bp@result$ID %in% common_before_slim, ]

common_early_before <- ego.early.bp[ego.early.bp@result$ID %in% common_before_slim, ]

common_late_slim  <- ego.late.bp.slim[ego.late.bp.slim@result$ID %in% common_after_slim, ]

common_early_slim <- ego.early.bp.slim[ego.early.bp.slim@result$ID %in% common_after_slim, ]

cat("\nCommon GO terms before simplification (n =", length(common_before_slim), "):\n")

cat("\nCommon GO terms after simplification (n =", length(common_after_slim), "):\n")

cat("\nTotal GO terms before simplification: Late =", length(late_ids), ", Early =", length(early_ids), "\n")

cat("Percentage overlap =", round(100 * length(common_before_slim) / (length(late_ids) + length(early_ids))), "%\n")

cat("\nTotal GO terms after simplification: Late =", length(late_slim_ids), ", Early =", length(early_slim_ids), "\n")

cat("Percentage overlap =", round(100 * length(common_after_slim) / (length(late_slim_ids) +length(early_slim_ids))), "%\n")

### Molecular Function 

ego.late.mf = enrichGO(gene         = consistently.selected.late,
OrgDb         =org.Aobtectus.eg.db,
keyType       = 'GID',
ont           = "MF",
pAdjustMethod = "BH",
pvalueCutoff  = 0.01,
qvalueCutoff  = 0.05)

ego.early.mf = enrichGO(gene         = consistently.selected.early,
OrgDb         =org.Aobtectus.eg.db,
keyType       = 'GID',
ont           = "MF",
pAdjustMethod = "BH",
pvalueCutoff  = 0.01,
qvalueCutoff  = 0.05)

semData <- godata(org.Aobtectus.eg.db, ont="MF", keytype = "GID")

ego.late.mf.slim <- simplify(ego.late.mf,
cutoff      = 0.7,
by          = "p.adjust",
select_fun  = min,
measure     = "Wang",
semData     = semData)

ego.early.mf.slim <- simplify(ego.early.mf,
cutoff      = 0.7,
by          = "p.adjust",
select_fun  = min,
measure     = "Wang",
semData     = semData)

# GO-term overlap stats, Early vs Late (MF)

late_ids  <- ego.late.mf@result$ID

early_ids <- ego.early.mf@result$ID

late_slim_ids  <- ego.late.mf.slim@result$ID

early_slim_ids <- ego.early.mf.slim@result$ID

common_before_slim <- intersect(late_ids, early_ids)

common_after_slim  <- intersect(late_slim_ids, early_slim_ids)

common_late_before  <- ego.late.mf[ego.late.mf@result$ID %in% common_before_slim, ]

common_early_before <- ego.early.mf[ego.early.mf@result$ID %in% common_before_slim, ]

common_late_slim  <- ego.late.mf.slim[ego.late.mf.slim@result$ID %in% common_after_slim, ]

common_early_slim <- ego.early.mf.slim[ego.early.mf.slim@result$ID %in% common_after_slim, ]

cat("\nCommon GO terms before simplification (n =", length(common_before_slim), "):\n")

cat("\nCommon GO terms after simplification (n =", length(common_after_slim), "):\n")

### Cellular Component 
ego.late.cc = enrichGO(gene         = consistently.selected.late,
OrgDb         =org.Aobtectus.eg.db,
keyType       = 'GID',
ont           = "CC",
pAdjustMethod = "BH",
pvalueCutoff  = 0.01,
qvalueCutoff  = 0.05)

ego.early.cc = enrichGO(gene         = consistently.selected.early,
OrgDb         =org.Aobtectus.eg.db,
keyType       = 'GID',
ont           = "CC",
pAdjustMethod = "BH",
pvalueCutoff  = 0.01,
qvalueCutoff  = 0.05)

semData <- godata(org.Aobtectus.eg.db, ont="CC", keytype = "GID")

ego.late.cc.slim <- simplify(ego.late.cc,
cutoff      = 0.5,
by          = "p.adjust",
select_fun  = min,
measure     = "Wang",
semData     = semData)

ego.early.cc.slim <- simplify(ego.early.cc,
cutoff      = 0.5,
by          = "p.adjust",
select_fun  = min,
measure     = "Wang",
semData     = semData)