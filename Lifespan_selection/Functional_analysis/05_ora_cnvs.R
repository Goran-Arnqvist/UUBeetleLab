#!/usr/bin/env Rscript

# Author: Athanasios Toros
# Email: athanasiostoros@gmail.com

# Over-representation analysis (ORA) of CNV candidate genes.


suppressPackageStartupMessages({
  library(tidyverse)
  library(org.Aobtectus.eg.db)
  library(GOSemSim)
  library(clusterProfiler)
})

set.seed(42)
config <- yaml::read_yaml("config.yaml")

# ---------------------------------------------------------------------------
# Universe (all 10 assembled scaffolds)
# ---------------------------------------------------------------------------
gene_locs <- AnnotationDbi::select(
  org.Aobtectus.eg.db,
  keys    = keys(org.Aobtectus.eg.db, keytype = "GID"),
  columns = c("GID", "CHROMOSOME", "LOW_COVERAGE", "LOW_QUALITY", "HIGH_COVERAGE", "NEAR_GAP"),
  keytype = "GID"
)

# limit genes to only assembled scaffolds
included_chromosomes <- unlist(config$genome$chromosome_accessions_all)

universe <- unique(gene_locs$GID[ gene_locs$CHROMOSOME %in% included_chromosomes ])

# Filter out genes in bad regions
cols_to_convert <- c("LOW_COVERAGE", "LOW_QUALITY", "HIGH_COVERAGE", "NEAR_GAP")
gene_locs <- gene_locs %>%
  mutate(across(all_of(cols_to_convert), ~ as.numeric(.)))

gene_locs_filtered <- gene_locs[rowSums(gene_locs[cols_to_convert]) == 0, ]
universe_filtered <- unique(gene_locs_filtered$GID[ gene_locs_filtered$CHROMOSOME %in% included_chromosomes ])

# ---------------------------------------------------------------------------
# CNV candidate-gene lists, filtered to the universe
# ---------------------------------------------------------------------------
genes.late.gain = scan(config$inputs$cnv_gene_lists$late_gain, what = character())
genes.early.gain = scan(config$inputs$cnv_gene_lists$early_gain, what = character())

genes.late.deletion = scan(config$inputs$cnv_gene_lists$late_deletions, what = character())
genes.early.deletion = scan(config$inputs$cnv_gene_lists$early_deletions, what = character())

genes.late.deletion.filtered = intersect(genes.late.deletion, universe_filtered)
genes.early.deletion.filtered = intersect(genes.early.deletion, universe_filtered)
genes.late.gain.filtered     = intersect(genes.late.gain, universe_filtered)
genes.early.gain.filtered    = intersect(genes.early.gain, universe_filtered)

cat("---- Gene Set Filtering Statistics ----\n")

universe_size <- length(universe_filtered)
cat("Total genes in filtered universe:", universe_size, "\n\n")

report_set_stats <- function(name, vec) {
  n <- length(vec)
  pct <- round(100 * n / universe_size, 2)
  cat(sprintf("%-30s: %5d genes (%.2f%% of universe)\n", name, n, pct))
}

cat("---- Stats for affected genes (filtered) ----\n")

report_set_stats("Early Deletions", genes.early.deletion.filtered)
report_set_stats("Late Deletions", genes.late.deletion.filtered)
report_set_stats("Late Gains",      genes.late.gain.filtered)
report_set_stats("Early Gains",     genes.early.gain.filtered)

# ===========================================================================
# Deletions
# ===========================================================================

### Biological Process
late.deletion.bp = enrichGO(gene         = genes.late.deletion.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "BP",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

early.deletion.bp = enrichGO(gene         = genes.early.deletion.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "BP",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

# The default measure of Wang is left, because it might be
# more robust to sparse, or bad quality annotation
semData <- godata(org.Aobtectus.eg.db, ont="BP", keytype = "GID")
late.deletion.bp.slim = simplify(late.deletion.bp,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

early.deletion.bp.slim = simplify(early.deletion.bp,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

### Molecular Function
late.deletion.mf = enrichGO(gene         = genes.late.deletion.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "MF",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

early.deletion.mf = enrichGO(gene         = genes.early.deletion.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "MF",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

semData <- godata(org.Aobtectus.eg.db, ont="MF", keytype = "GID")
late.deletion.mf.slim = simplify(late.deletion.mf,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

early.deletion.mf.slim = simplify(early.deletion.mf,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

### Cellular Component
late.deletion.cc = enrichGO(gene         = genes.late.deletion.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "CC",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

early.deletion.cc = enrichGO(gene         = genes.early.deletion.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "CC",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

semData <- godata(org.Aobtectus.eg.db, ont="CC", keytype = "GID")
late.deletion.cc.slim = simplify(late.deletion.cc,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

early.deletion.cc.slim = simplify(early.deletion.cc,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

# ===========================================================================
# Gains
# ===========================================================================

### Biological Process
late.gain.bp = enrichGO(gene         = genes.late.gain.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "BP",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

early.gain.bp = enrichGO(gene         = genes.early.gain.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "BP",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

semData <- godata(org.Aobtectus.eg.db, ont="BP", keytype = "GID")
late.gain.bp.slim = simplify(late.gain.bp,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

early.gain.bp.slim = simplify(early.gain.bp,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

### Molecular Function
late.gain.mf = enrichGO(gene         = genes.late.gain.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "MF",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

early.gain.mf = enrichGO(gene         = genes.early.gain.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "MF",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

semData <- godata(org.Aobtectus.eg.db, ont="MF", keytype = "GID")
late.gain.mf.slim = simplify(late.gain.mf,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

early.gain.mf.slim = simplify(early.gain.mf,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

### Cellular Component
late.gain.cc = enrichGO(gene         = genes.late.gain.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "CC",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

early.gain.cc = enrichGO(gene         = genes.early.gain.filtered,
                OrgDb         =org.Aobtectus.eg.db,
                keyType       = 'GID',
                ont           = "CC",
                pAdjustMethod = "BH",
                universe = universe_filtered,
                qvalueCutoff  = 0.05)

semData <- godata(org.Aobtectus.eg.db, ont="CC", keytype = "GID")
late.gain.cc.slim = simplify(late.gain.cc,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

early.gain.cc.slim = simplify(early.gain.cc,
                     cutoff      = 0.7,
                     by          = "p.adjust",
                     select_fun  = min,
                     measure     = "Wang",
                     semData     = semData)

# ===========================================================================
# Common genes: gained in one regime and deleted in the other
# ===========================================================================
# Genes that are gains in Early, deletions in Late:
genes.late.deletion.early.gain.filtered = intersect(genes.late.deletion.filtered, genes.early.gain.filtered)

# The reverse
genes.early.deletion.late.gain.filtered = intersect(genes.early.deletion.filtered, genes.late.gain.filtered)

cat("---- Gene Set Filtering Statistics ----\n")

universe_size <- length(universe_filtered)
cat("Total genes in filtered universe:", universe_size, "\n\n")

report_set_stats <- function(name, vec) {
  n <- length(vec)
  pct <- round(100 * n / universe_size, 2)
  cat(sprintf("%-30s: %5d genes (%.2f%% of universe)\n", name, n, pct))
}

cat("---- Stats for affected genes (filtered) ----\n")

report_set_stats("Late Deletions Early Gains", genes.late.deletion.early.gain.filtered)
report_set_stats("Early Deletions Late Gains", genes.early.deletion.late.gain.filtered)
