#!/usr/bin/env Rscript

# Author: Athanasios Toros
# Email: athanasiostoros@gmail.com


suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(PSCBS)
  library(parallel)
})

## ---------- CLI ----------
ap <- ArgumentParser(description = "Parallel PSCBS segmentation per chromosome")
ap$add_argument("-i", "--input",   required = TRUE, help = "Input bedGraph")
ap$add_argument("-o", "--output",  required = TRUE, help = "Output bedGraph")
ap$add_argument("-t", "--threads", default = detectCores() - 1, type = "integer",
                help = "Workers [default: %(default)s]")
ap$add_argument("--undo",     default = 0,     type = "integer")
ap$add_argument("--alpha",    default = 1e-4,  type = "double")
ap$add_argument("--minwidth", default = 2,     type = "integer")
ap$add_argument("--verbose",  action = "store_true",
                help = "Print per-chrom debug info")
args <- ap$parse_args()

## ---------- I/O ----------
dt <- fread(args$input, col.names = c("chrom","start","end","logR"))
setorder(dt, chrom, start)
chroms  <- unique(dt$chrom)
workers <- max(1L, args$threads)

## ---------- Function to do one chromosome ----------
do_chr <- function(chr) {
  d0 <- dt[chrom == chr]

  # 1) call outliers WITHOUT dropping rows
  d0[, logR_clean := dropSegmentationOutliers(
         y          = logR,
         chromosome = 1L,
         x          = start)]

  # 2) build PSCBS input – all rows preserved, outliers are NA
  cna <- data.frame(
    chromosome = 1L,
    x          = d0$start,
    y          = d0$logR_clean
  )
  gaps      <- PSCBS::findLargeGaps(cna, minLength = 5e4)
  knownSegs <- PSCBS::gapsToSegments(gaps)

  # 3) segment
  fit <- PSCBS::segmentByCBS(
    y             = cna$y,
    chromosome    = cna$chromosome,
    x             = cna$x,
    knownSegments = knownSegs,
    undo          = args$undo,
    joinSegments  = TRUE,
    seed          = 0xA5EED,
    verbose       = FALSE,
    alpha         = args$alpha,
    min.width     = args$minwidth
  )

  # 4) remap & compute mean **in the same table**
  out <- fit$output
  for (i in seq_len(nrow(out))) {
    sr <- fit$segRows$startRow[i]
    er <- fit$segRows$endRow[i]
    if (is.na(sr) || is.na(er)) next
    out$start[i] <- d0$start[sr]
    out$end[i]   <- d0$end[er]
    out$mean[i]  <- mean(d0$logR_clean[sr:er], na.rm = TRUE)
  }
  data.table(chrom = chr, start = out$start, end = out$end, mean = out$mean)
}

## ---------- Parallel run with error capture ----------
seg_list <- mclapply(
  chroms, mc.cores = workers,
  FUN = function(chr) {
    tryCatch(
      do_chr(chr),
      error = function(e) {
        message(sprintf("ERROR on chr %s: %s", chr, e$message))
        NULL
      }
    )
  }
)

## drop any NULLs
seg_list <- Filter(Negate(is.null), seg_list)
if (length(seg_list)==0L) {
  stop("All chromosomes failed.  See errors above.")
}

## ---------- Collate & write ----------
# prevent scientific notation, otherwise fwrite will print
# 1e+05 instead of 100000
res <- rbindlist(seg_list)[order(factor(chrom, levels=chroms), start)]

options(scipen = 999)

# coerce to integer so fwrite prints them as plain numbers
res[, start := as.integer(start)]
res[, end   := as.integer(end)]

# Sort for bedGraph output
res <- res[order(chrom, start)]

fwrite(res, args$output, sep="\t", col.names=FALSE, quote=FALSE)

## ---------- Summary ----------
cat(sprintf(
  "Done: %d chromosomes → %d total segments → %s\n",
  length(seg_list), nrow(res), args$output
))

# ---------------------------------------------------------------------------
# Commands downstream of this script; produces output TSV that
# cnv_pipeline.py reads (pipeline.paths.segmentation_indices_file):
#
#   sort -k 1,1 -k 2,2n segmentation_calls_100bp_undo.bedgraph \
#     > segmentation_calls_100bp_undo_sorted.bedgraph
#
#   # Some of the segments are empty, so we need to remove them
#   awk 'BEGIN{OFS="\t"} NF>=4 {print $1,$2,$3,NR,$4}' \
#     segmentation_calls_100bp_undo_sorted.bedgraph \
#     > segmentation_calls_100bp_undo_index_filtered.bedgraph
#
#   bedtools intersect \
#     -a A_obtectus_window_100bp.bed \
#     -b segmentation_calls_100bp_undo_index_filtered.bedgraph \
#     -wa -wb \
#   | cut -f1-3,7 | sort -k1,1 -k2,2n > genome_windows_with_segmentation_indices.tsv
# ---------------------------------------------------------------------------
