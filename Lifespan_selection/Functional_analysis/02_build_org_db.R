#!/usr/bin/env Rscript
#
# Author: Athanasios Toros
# Email: athanasiostoros@gmail.com
# Build the organism annotation package org.Aobtectus.eg.db (AnnotationForge).


suppressPackageStartupMessages({
  library(GO.db)
  library(AnnotationDbi)
  library(tidyverse)
  library(AnnotationForge)
})

config <- yaml::read_yaml("config.yaml")

# ---------------------------------------------------------------------------
# GID -> GO table (+ IEA evidence) and GO term metadata
# ---------------------------------------------------------------------------
go_df <- read.delim(config$outputs$go_universe, header = TRUE)
colnames(go_df) <- c("GID", "GO")
go_df$EVIDENCE  <- "IEA"
go_df$ONTOLOGY  <- Ontology(go_df$GO)
go_df           <- go_df[, c("GID", "GO", "EVIDENCE")]

go_ids  <- unique(go_df$GO)
go_meta <- select(
  GO.db,
  keys    = go_ids,
  keytype = "GOID",
  columns = c("TERM", "DEFINITION", "ONTOLOGY")
)

final_go <- merge(go_df, go_meta, by.x = "GO", by.y = "GOID", all.x = TRUE)
go_for_pkg <- final_go[, c("GID", "GO", "EVIDENCE")]

# ---------------------------------------------------------------------------
# GID -> chromosome and GID -> quality-flag tables
# ---------------------------------------------------------------------------
gene2chrom_df <- read.delim(config$outputs$chromosome_annotation, sep = "\t")
gene2chrom_df <- gene2chrom_df %>% drop_na() %>% distinct()
colnames(gene2chrom_df) <- c("GID", "CHROMOSOME")

gene2exclusion_df <- read.delim(config$outputs$quality_annotations, sep = "\t")
gene2exclusion_df <- gene2exclusion_df %>% drop_na() %>% distinct()
colnames(gene2exclusion_df) <- c("GID", "LOW_COVERAGE", "LOW_QUALITY", "HIGH_COVERAGE", "NEAR_GAP")

go_for_pkg <- go_for_pkg %>% drop_na() %>% distinct()

# ---------------------------------------------------------------------------
# Compile the org package
# ---------------------------------------------------------------------------
makeOrgPackage(
  go        = go_for_pkg,
  gene2chrom = gene2chrom_df,
  gene2exclusion = gene2exclusion_df,
  version   = "0.2.0",
  maintainer= "Athanasios Toros athanasios.toros.0729@student.uu.se",
  author    = "Athanasios Toros",
  outputDir = ".",
  tax_id    = "200917",               # NCBI taxon for *A. obtectus*
  genus     = "Acanthoscelides",
  species   = "obtectus",
  goTable   = "go",
  verbose = TRUE
)
