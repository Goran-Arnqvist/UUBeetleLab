###############################################################################
## QTL mapping of body weight and lifespan in an E x L F2 intercross of the
## seed beetle Acanthoscelides obtectus.
##
## This master script runs the whole analysis from the raw genotype/phenotype
## file (data/pheno_input.csv) and produces:
##   Main text : Fig. 4, Fig. 5, Fig. 6
##   SI        : Fig. S13 - Fig. S27, Table S1, Table S2
##
## Run from the Submission/ folder:
##     Rscript make_figures.R
## Figures are written to figures/, tables to tables/. See README.md.
##
## Conventions used throughout
##   * Alleles: E and L, from the E and L selection lines respectively.
##     Additive effect  a = (mean_LL - mean_EE) / 2  ("effect of the L allele").
##   * sex is coded 1 = male, 0 = female.
##   * Trait units: weight in mg, lifespan in days.
###############################################################################

## ===========================================================================
## 0. Setup
## ===========================================================================

library(qtl)       # R/qtl 1.74 
library(scales)    
library(parallel) 

DATA_DIR <- "data"          # raw input
FIG_DIR  <- "figures"
TAB_DIR  <- "tables"
dir.create(FIG_DIR, showWarnings = FALSE)
dir.create(TAB_DIR, showWarnings = FALSE)

## --- Saved state -------------------------------------------------------------
## Three steps are slow: preparing the cross (~4 min), the permutation tests
## (~35 min) and the multiple-QTL model search (~40 min). Their results are
## kept as .Rds files.
##   USE_SAVED_STATE = TRUE  : load the objects used for the manuscript from
##                             saved_state/ (exact reproduction, ~1-2 min).
##   USE_SAVED_STATE = FALSE : compute everything from data/pheno_input.csv
##                             (~80 min) and save the results to results/. Later runs reuse
##                             what is in results/; delete it to start over.
USE_SAVED_STATE <- TRUE
STATE_DIR <- if (USE_SAVED_STATE) "saved_state" else "results"
dir.create(STATE_DIR, showWarnings = FALSE)

## Load <file> from the state folder if it is there; otherwise evaluate <expr>,
## save the result there, and return it.
cached <- function(file, expr) {
  path <- file.path(STATE_DIR, file)
  if (file.exists(path)) {
    return(readRDS(path))
  }
  if (USE_SAVED_STATE) stop(path, " is missing; set USE_SAVED_STATE <- FALSE to recompute it")
  obj <- force(expr)
  saveRDS(obj, path)
  obj
}

N_PERM    <- 1000                           # permutations per genome scan
N_CORES   <- max(1, detectCores() - 1)
SIM_SEED  <- 1                              # genotype imputations (sim.geno)
PERM_SEED <- 20250601                       # permutations

## --- Traits and plotting style ---------------------------------------------
TRAITS     <- c(wgt = 1L, life = 2L)                 # phenotype columns in the cross
TRAIT_NAME <- c(wgt = "Weight",   life = "Lifespan")
TRAIT_UNIT <- c(wgt = "mg",       life = "days")

COL_M     <- "royalblue"        # males
COL_F     <- "firebrick"        # females
COL_ALL   <- "grey40"           # sexes pooled (sex-averaged)
COL_XF    <- "#4daf4a"          # X chromosome, females of the forward cross
COL_XR    <- "#ff7f00"          # X chromosome, females of the reverse cross
COL_SHADE <- "grey93"           # alternate-chromosome shading
COL_THR   <- "black"            # permutation thresholds
ADD_RAMP  <- colorRampPalette(c("#762a83", "white", "#1b7837"))(101)  # additive colour scale

AUTOSOMES <- as.character(1:9)

## Open a PDF, run the plotting code, always close the device.
save_pdf <- function(file, width, height, code) {
  path <- file.path(FIG_DIR, file)
  pdf(path, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off())
  force(code)
}

## ===========================================================================
## 1. The F2 cross
## ===========================================================================
## pheno_input.csv: 1330 F2 individuals with weight, lifespan and sex, genotyped
## at 3487 markers (100-kb windows) on chromosomes 1-9 and X.
##  * Runs of adjacent markers with identical genotypes are reduced to one
##    marker each (2839 remain). findDupMarkers() picks the marker kept from
##    each run at random, so the markers kept in the manuscript analysis are
##    listed in retained_markers.txt; the genotype data are the same whichever
##    marker of a run is kept, only its name (window position) differs.
##  * Genetic map: est.map(error.prob = 0.01).
##  * Genotype probabilities on a 1-cM grid (calc.genoprob) and 128 genotype
##    imputations (sim.geno), both with error.prob = 0.01. Every result below
##    is based on these.
seed <- cached("seedle_cross.Rds", {
  x <- suppressWarnings(read.cross(format = "csv", dir = DATA_DIR, file = "pheno_input.csv",
                                   na.strings = "-", genotypes = c("E", "H", "L", "C", "D"),
                                   alleles = c("E", "L"), estimate.map = FALSE,
                                   convertXdata = TRUE))
  dups <- findDupMarkers(x, adjacent.only = TRUE)
  runs <- mapply(c, names(dups), dups, SIMPLIFY = FALSE)       # each run of duplicates
  keep <- readLines(file.path(DATA_DIR, "retained_markers.txt"))
  stopifnot(all(sapply(runs, function(r) sum(r %in% keep)) == 1))
  x <- drop.markers(x, setdiff(unlist(runs), keep))
  stopifnot(identical(markernames(x), keep))
  x <- replace.map(x, est.map(x, error.prob = 0.01, n.cluster = N_CORES))
  x <- calc.genoprob(x, step = 1, error.prob = 0.01)
  set.seed(SIM_SEED)
  sim.geno(x, n.draws = 128, step = 1, error.prob = 0.01)
})
sex <- as.numeric(pull.pheno(seed, "sex") == "M")
stopifnot(all(sex %in% c(0, 1)))

males   <- subset(seed, ind = sex == 1)
females <- subset(seed, ind = sex == 0)

## ===========================================================================
## 2. Single-QTL genome scans (interval mapping, EM algorithm)
## ===========================================================================
##   sex-averaged : both sexes, sex as an additive covariate     (Fig. S26, S27)
##   per sex      : males and females scanned separately          (Fig. 5)
scans <- list()
for (tr in names(TRAITS)) {
  pc <- TRAITS[[tr]]
  scans[[tr]] <- list(
    all = scanone(seed, addcovar = sex, pheno.col = pc),
    m   = scanone(males,   pheno.col = pc),
    f   = scanone(females, pheno.col = pc))
}

## --- Permutation thresholds (genome-wide 5 %; X chromosome separately) ----
perm_scanone <- function(cross, pheno.col, addcovar = NULL) {
  ## split the permutations across cores, each with its own seed
  n_each <- diff(round(seq(0, N_PERM, length.out = N_CORES + 1)))
  n_each <- n_each[n_each > 0]
  set.seed(PERM_SEED)
  seeds <- sample.int(1e8, length(n_each))
  res <- mclapply(seq_along(n_each), function(i) {
    set.seed(seeds[i])
    scanone(cross, pheno.col = pheno.col, addcovar = addcovar,
            n.perm = n_each[i], perm.Xsp = TRUE, verbose = FALSE)
  }, mc.cores = length(n_each))
  do.call(c, res)
}

PERM_FILE <- c(all = "%s_add_perm.Rds", m = "%s_male_perm.Rds", f = "%s_female_perm.Rds")
perms <- list()
for (tr in names(TRAITS)) {
  for (grp in c("all", "m", "f")) {
    cross <- switch(grp, all = seed, m = males, f = females)
    cov   <- if (grp == "all") sex else NULL
    perms[[paste(tr, grp, sep = "_")]] <- cached(sprintf(PERM_FILE[[grp]], tr),
      suppressWarnings(perm_scanone(cross, TRAITS[[tr]], addcovar = cov)))
  }
}

## 5 % genome-wide LOD thresholds, autosomes (A) and X separately.
thr <- lapply(perms, function(p) {
  s <- summary(p, alpha = 0.05)
  c(A = as.numeric(s$A), X = as.numeric(s$X))
})

## ===========================================================================
## 3. Single-QTL effect scans (additive and dominance effects at every position)
## ===========================================================================
## effectscan() estimates the genotype means at each position from the 128
## imputations. Sex-averaged = both sexes together; on the X chromosome R/qtl
## reports separate additive effects for males and for females of the forward
## and reverse cross directions. Per sex = males and females separately.

as_df <- function(x) { x <- as.data.frame(x); x$chr <- as.character(x$chr); x }
effects <- list()
for (tr in names(TRAITS)) {
  pc <- TRAITS[[tr]]
  effects[[tr]] <- list(
    all = as_df(effectscan(seed,    pheno.col = pc, draw = FALSE, get.se = TRUE)),
    m   = as_df(effectscan(males,   pheno.col = pc, draw = FALSE, get.se = TRUE)),
    f   = as_df(effectscan(females, pheno.col = pc, draw = FALSE, get.se = TRUE)))
  ## all three scans are on the same positions, in the same order
  with(effects[[tr]], stopifnot(identical(all$pos, m$pos), identical(all$pos, f$pos)))
}
CHRS <- unique(effects$wgt$all$chr)            # "1" ... "9", "X"

## ===========================================================================
## 4. Multiple-QTL models 
## ===========================================================================
## Forward selection, autosomes only, sex as a covariate throughout:
##  1. The first QTL is the highest peak of a single-QTL scan (Haley-Knott)
##     with sex as an additive and an interactive covariate.
##  2. addqtl() scans for one more QTL given those already in the model; the
##     highest peak is added if its LOD exceeds 3.
##  3. After each addition all QTL positions are refined jointly (refineqtl,
##     multiple imputation). Repeat from 2 until no peak exceeds LOD 3.
##  4. Each QTL's interaction with sex is tested in the final model
##     (addcovarint); a QTL x sex term is included where p < 0.05.
## The final model is refined once more, keeping each QTL's LOD profile
## (Fig. 6).
MQM <- list(
  wgt  = list(file = "mqm_wgt_int.Rds",  digits = 3),
  life = list(file = "mqm_life_int.Rds", digits = 2))

sexdf <- data.frame(sex = sex)
qtl_formula <- function(n, sex_int = integer(0), add_sex = FALSE)
  as.formula(paste("y ~", paste(c(if (add_sex) "sex", paste0("Q", seq_len(n)),
                                  if (length(sex_int)) paste0("Q", sex_int, ":sex")),
                                collapse = " + ")))

forward_select <- function(pheno.col, lod_threshold = 3, max_qtl = 10) {
  auto <- function(sc) sc[sc$chr != "X", ]
  top  <- max(auto(scanone(seed, pheno.col = pheno.col, method = "hk",
                           addcovar = sexdf, intcovar = sexdf)))
  qtl  <- makeqtl(seed, chr = as.character(top$chr), pos = top$pos)
  while (qtl$n.qtl < max_qtl) {
    top <- max(auto(addqtl(seed, pheno.col = pheno.col, qtl = qtl, covar = sexdf,
                           forceXcovar = TRUE)))
    if (top$lod < lod_threshold) break
    qtl <- addtoqtl(seed, qtl = qtl, chr = as.character(top$chr), pos = top$pos)
    qtl <- refineqtl(seed, pheno.col = pheno.col, qtl = qtl, covar = sexdf,
                     formula = qtl_formula(qtl$n.qtl), forceXcovar = TRUE,
                     verbose = FALSE)
    cat(sprintf("    %d QTL: %s\n", qtl$n.qtl, paste(qtl$name, collapse = ", ")))
  }
  qtl
}

mqm <- list()
for (tr in names(MQM)) {
  spec <- MQM[[tr]]
  mqm[[tr]] <- cached(spec$file, {
    qtl <- forward_select(TRAITS[[tr]])
    ci  <- addcovarint(seed, pheno.col = TRAITS[[tr]], qtl = qtl, covar = sexdf,
                       formula = qtl_formula(qtl$n.qtl), icovar = "sex")
    print(ci)
    sex_int <- which(ci[, "Pvalue(F)"] < 0.05)
    q <- refineqtl(seed, pheno.col = TRAITS[[tr]], qtl = qtl, covar = sexdf,
                   formula = qtl_formula(qtl$n.qtl, sex_int), forceXcovar = TRUE,
                   keeplodprofile = TRUE, verbose = FALSE)
    attr(q, "sex_int") <- sex_int
    q
  })
  ## The saved manuscript models predate the "sex_int" attribute; their
  ## QTL x sex terms were QTL 2-5 for lifespan and none for weight.
  sex_int <- attr(mqm[[tr]], "sex_int")
  if (is.null(sex_int)) sex_int <- list(wgt = integer(0), life = 2:5)[[tr]]
  MQM[[tr]]$fit_formula <- qtl_formula(mqm[[tr]]$n.qtl, sex_int,
                                       add_sex = length(sex_int) == 0)
  cat("  QTL:", paste(mqm[[tr]]$name, collapse = ", "), "\n")
}

## ===========================================================================
## 5. Plotting helpers
## ===========================================================================

## Concatenate chromosomes (cM) into one genome axis, gap cM between them.
genome_axis <- function(chr, pos, gap = 25) {
  chr  <- as.character(chr)
  chrs <- unique(chr)
  x <- numeric(length(pos)); st <- en <- numeric(length(chrs)); off <- 0
  for (i in seq_along(chrs)) {
    k <- chr == chrs[i]
    x[k]  <- pos[k] - min(pos[k]) + off
    st[i] <- off
    en[i] <- max(x[k])
    off   <- en[i] + gap
  }
  list(x = x, chr = chr, chrs = chrs, start = st, end = en, mid = (st + en) / 2)
}

## Empty genome-wide panel with alternate chromosomes shaded.
genome_frame <- function(g, ylim, ylab, main = "", yaxt = "s") {
  plot(NA, xlim = range(c(g$start, g$end)) + c(-12.5, 12.5), ylim = ylim, xaxs = "i",
       xaxt = "n", yaxt = yaxt, xlab = "", ylab = ylab, main = main)
  usr <- par("usr")
  for (i in seq_along(g$chrs))
    if (i %% 2 == 0)
      rect(g$start[i] - 12.5, usr[3], g$end[i] + 12.5, usr[4], col = COL_SHADE, border = NA)
  box()
  axis(1, at = g$mid, labels = g$chrs, tick = FALSE, line = -0.4)
}

## LOD curve, drawn chromosome by chromosome (sgn = -1 mirrors it downwards).
lod_curve <- function(g, lod, col, sgn = 1, lwd = 1.3) {
  for (ch in g$chrs) {
    k <- g$chr == ch
    lines(g$x[k], sgn * lod[k], col = col, lwd = lwd)
  }
}

## Threshold line over each chromosome: autosomal threshold over 1-9, the
## separately estimated X threshold over the X.
threshold_lines <- function(g, t, col = COL_THR, sgn = 1, lty = 1, lwd = 1.2) {
  for (i in seq_along(g$chrs)) {
    h <- if (g$chrs[i] == "X") t[["X"]] else t[["A"]]
    segments(g$start[i], sgn * h, g$end[i], sgn * h, col = col, lty = lty, lwd = lwd)
  }
}

## Genotyped marker positions (pseudomarkers are named "c<chr>.loc<pos>").
marker_rug <- function(g, scan, y, tick) {
  mk <- !grepl("^c[0-9X]+\\.loc", rownames(scan))
  segments(g$x[mk], y, g$x[mk], y + tick, lwd = 0.3, col = "grey20")
}

density_panel <- function(dens, cols, xlim, main) {
  dens <- lapply(dens, function(d) { d$y <- d$y / max(d$y); d })
  plot(NA, xlim = xlim, ylim = c(0, 1.35), yaxt = "n", xlab = "", ylab = "",
       main = main, yaxs = "i")
  for (i in seq_along(dens)) {
    polygon(c(dens[[i]]$x, rev(dens[[i]]$x)), c(dens[[i]]$y, rep(0, length(dens[[i]]$y))),
            col = alpha(cols[i], 0.3), border = NA)
    lines(dens[[i]], col = cols[i], lwd = 1.5)
  }
  abline(v = 0, lty = 2, lwd = 1.2)
}

scatter_panel <- function(x, y, xlim, ylim, main = "", bg = alpha(COL_ALL, 0.35),
                          border = alpha("black", 0.35), xlab = "", ylab = "",
                          cex = 0.6, add = FALSE) {
  if (!add) {
    plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)
    abline(h = 0, v = 0, col = "grey60", lty = 2)
  }
  points(x, y, pch = 21, col = border, bg = bg, cex = cex, lwd = 0.3)
}

panel_legend <- function(legend, cols, cex = 1.05) {
  legend("top", horiz = TRUE, bty = "n", legend = legend, cex = cex,
         fill = alpha(cols, 0.4), border = cols, x.intersp = 0.5, inset = c(0, 0.01))
}

r_label <- function(r, col = "black", line = 0) {
  legend("topleft", bty = "n", text.col = col, cex = 0.9,
         legend = sprintf("r = %.2f", r), inset = c(0, 0.08 * line))
}

lim_sym <- function(...) c(-1, 1) * max(abs(c(...)), na.rm = TRUE) * 1.05

hist_by_sex <- function(xm, xf, main, xlab) {
  xm <- xm[is.finite(xm)]; xf <- xf[is.finite(xf)]
  br <- pretty(range(c(xm, xf)), n = 60)
  hm <- hist(xm, breaks = br, plot = FALSE); hf <- hist(xf, breaks = br, plot = FALSE)
  dm <- density(xm); df <- density(xf)
  ymax <- max(hm$density, hf$density, dm$y, df$y)
  plot(NA, xlim = range(c(dm$x, df$x)), ylim = c(0, ymax * 1.05), yaxs = "i",
       xlab = xlab, ylab = "Density", main = main)
  plot(hm, freq = FALSE, add = TRUE, col = alpha(COL_M, 0.30), border = alpha(COL_M, 0.6))
  plot(hf, freq = FALSE, add = TRUE, col = alpha(COL_F, 0.30), border = alpha(COL_F, 0.6))
  lines(dm, col = COL_M, lwd = 2); lines(df, col = COL_F, lwd = 2)
  abline(v = 0, lty = 2)
  legend("topright", bty = "n", fill = alpha(c(COL_M, COL_F), 0.4),
         border = c(COL_M, COL_F), legend = c("Males", "Females"))
}

hist_pooled <- function(x, main, xlab) {
  x  <- x[is.finite(x)]
  h  <- hist(x, breaks = pretty(range(x), n = 50), plot = FALSE)
  d  <- density(x)
  plot(NA, xlim = range(d$x), ylim = c(0, max(h$density, d$y) * 1.05), yaxs = "i",
       xlab = xlab, ylab = "Density", main = main)
  plot(h, freq = FALSE, add = TRUE, col = "grey80", border = "grey45")
  lines(d, lwd = 2)
  abline(v = 0, lty = 2)
}

eff_lab <- function(tr) sprintf("Additive effect of the L allele (%s)", TRAIT_UNIT[[tr]])

## ===========================================================================
## 6. Figures
## ===========================================================================

## --- Fig. 4 and Fig. S18: per-sex additive effects, all markers -------------
## Additive effects estimated separately in males and in females.
save_pdf("Fig4_lifespan_additive_effects_by_sex.pdf", 9, 5.5, {
  par(mar = c(4.5, 4.5, 2.5, 1))
  hist_by_sex(effects$life$m$a, effects$life$f$a, "Lifespan", eff_lab("life"))
})
save_pdf("FigS18_weight_additive_effects_by_sex.pdf", 9, 5.5, {
  par(mar = c(4.5, 4.5, 2.5, 1))
  hist_by_sex(effects$wgt$m$a, effects$wgt$f$a, "Weight", eff_lab("wgt"))
})

## --- Fig. 5: per-sex LOD profiles, females up, males mirrored down ---------
save_pdf("Fig5_LOD_by_sex.pdf", 9, 10, {
  par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1))
  for (tr in names(TRAITS)) {
    s <- scans[[tr]]
    g <- genome_axis(s$f$chr, s$f$pos)
    tf <- thr[[paste0(tr, "_f")]]; tm <- thr[[paste0(tr, "_m")]]
    top <- max(s$f$lod, tf) * 1.08; bot <- max(s$m$lod, tm) * 1.08
    tick <- (top + bot) * 0.025
    genome_frame(g, ylim = c(-bot - tick, top), yaxt = "n",
                 ylab = sprintf("LOD score (%s)", TRAIT_NAME[[tr]]))
    at <- pretty(c(-bot, top))
    axis(2, at = at, labels = abs(at), las = 1)
    abline(h = 0, lwd = 1.5)
    threshold_lines(g, tf, col = COL_F, sgn = 1, lty = 2)
    threshold_lines(g, tm, col = COL_M, sgn = -1, lty = 2)
    lod_curve(g, s$f$lod, COL_F, sgn = 1)
    lod_curve(g, s$m$lod, COL_M, sgn = -1)
    marker_rug(g, s$f, par("usr")[3], tick)
    leg <- c("Females (upwards)", "Males (downwards)",
             "5% genome-wide threshold (permutations)")
    legend("bottom", inset = c(0, 1), xpd = NA, horiz = TRUE, bty = "n", cex = 0.85,
           text.width = strwidth(leg, cex = 0.85) * 1.1,
           lwd = c(1.5, 1.5, 1.2), lty = c(1, 1, 2), col = c(COL_F, COL_M, "grey30"),
           legend = leg)
    mtext("Chromosome", side = 1, line = 2.3)
  }
})

## --- Fig. 6: LOD profiles of the multiple-QTL models -----------------------
## Each QTL's profile compares the full model with the model omitting that QTL.
save_pdf("Fig6_MQM_LOD_profiles.pdf", 8, 10, {
  par(mfrow = c(2, 1), mar = c(3, 4.5, 2.5, 1), oma = c(2, 0, 0, 0))
  xmax <- c(wgt = 585, life = 875)
  for (tr in names(mqm)) {
    plotLodProfile(mqm[[tr]], main = TRAIT_NAME[[tr]], xlab = "",
                   ylab = "Profile LOD score", lwd = 1.5,
                   xlim = c(-20, xmax[[tr]]), las = 1)
  }
  mtext("Chromosome", side = 1, outer = TRUE, line = 0.5, cex = 1.2)
})

## --- Fig. S13, S14: sex-averaged additive effects, all markers -------------
save_pdf("FigS13_weight_additive_effects.pdf", 9, 5.5, {
  par(mar = c(4.5, 4.5, 2.5, 1))
  hist_pooled(effects$wgt$all$a, "Weight", eff_lab("wgt"))
})
save_pdf("FigS14_lifespan_additive_effects.pdf", 9, 5.5, {
  par(mar = c(4.5, 4.5, 2.5, 1))
  hist_pooled(effects$life$all$a, "Lifespan", eff_lab("life"))
})

## --- Fig. S15, S16 (sex-averaged) and S19, S20 (per sex): distributions of
##     additive effects per chromosome. X chromosome: males (hemizygous) and
##     females of the forward and reverse cross directions.
x_densities <- function(tr) {
  ex <- effects[[tr]]$all[effects[[tr]]$all$chr == "X", ]
  list(density(ex$a.male, na.rm = TRUE),
       density(ex$a.femaleforw, na.rm = TRUE),
       density(ex$a.femalerev, na.rm = TRUE))
}
per_chrom_density_fig <- function(tr, by_sex, file) {
  e <- effects[[tr]]
  dl <- list()
  for (ch in AUTOSOMES) {
    k <- e$all$chr == ch
    dl[[ch]] <- if (by_sex) list(density(e$m$a[k]), density(e$f$a[k]))
                else        list(density(e$all$a[k]))
  }
  dl[["X"]] <- x_densities(tr)
  xl <- lim_sym(unlist(lapply(dl, function(l) sapply(l, function(d) range(d$x)))))
  save_pdf(file, 8.5, 10.5, {
    par(mfrow = c(5, 2), mar = c(2.2, 1.5, 2, 1), oma = c(3.5, 2.5, 3, 0))
    for (ch in AUTOSOMES) {
      density_panel(dl[[ch]], if (by_sex) c(COL_M, COL_F) else COL_ALL, xl,
                    paste("Chromosome", ch))
      if (by_sex && ch == "1") panel_legend(c("Males", "Females"), c(COL_M, COL_F))
    }
    density_panel(dl[["X"]], c(COL_M, COL_XF, COL_XR), xl, "Chromosome X")
    panel_legend(c("Males", "Females (forward)", "Females (reverse)"),
                 c(COL_M, COL_XF, COL_XR))
    mtext(eff_lab(tr), side = 1, outer = TRUE, line = 1.8, cex = 1.1)
    mtext("Density (each curve scaled to the same height)", side = 2, outer = TRUE,
          line = 0.8, cex = 1)
    mtext(sprintf("%s: additive effect", TRAIT_NAME[[tr]]),
          side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2)
  })
}
per_chrom_density_fig("wgt",  FALSE, "FigS15_weight_additive_effects_by_chromosome.pdf")
per_chrom_density_fig("life", FALSE, "FigS16_lifespan_additive_effects_by_chromosome.pdf")
per_chrom_density_fig("wgt",  TRUE,  "FigS19_weight_additive_effects_by_chromosome_and_sex.pdf")
per_chrom_density_fig("life", TRUE,  "FigS20_lifespan_additive_effects_by_chromosome_and_sex.pdf")

## --- Fig. S17: weight vs lifespan additive effects per chromosome (sex-averaged)
## --- Fig. S25: the same, males and females separately ------------------------
au_rows <- effects$wgt$all$chr %in% AUTOSOMES
wl_lim <- function(groups) list(
  x = lim_sym(sapply(groups, function(g) effects$wgt[[g]]$a[au_rows])),
  y = lim_sym(sapply(groups, function(g) effects$life[[g]]$a[au_rows])))

save_pdf("FigS17_weight_vs_lifespan_effects_by_chromosome.pdf", 8, 8, {
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(3, 3, 1, 0), mgp = c(2, 0.6, 0))
  for (ch in AUTOSOMES) {
    k <- effects$wgt$all$chr == ch
    x <- effects$wgt$all$a[k]; y <- effects$life$all$a[k]
    lims <- wl_lim("all")
    scatter_panel(x, y, lims$x, lims$y, main = paste("Chromosome", ch))
    r_label(cor(x, y, use = "complete.obs"))
  }
  mtext(eff_lab("wgt"),  side = 1, outer = TRUE, line = 1.5, cex = 1.1)
  mtext(eff_lab("life"), side = 2, outer = TRUE, line = 1.2, cex = 1.1)
})

save_pdf("FigS25_weight_vs_lifespan_effects_by_chromosome_and_sex.pdf", 8, 8, {
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(3, 3, 1, 0), mgp = c(2, 0.6, 0))
  for (ch in AUTOSOMES) {
    k <- effects$wgt$all$chr == ch
    xm <- effects$wgt$m$a[k]; ym <- effects$life$m$a[k]
    xf <- effects$wgt$f$a[k]; yf <- effects$life$f$a[k]
    lims <- wl_lim(c("m", "f"))
    scatter_panel(xm, ym, lims$x, lims$y, main = paste("Chromosome", ch),
                  bg = alpha(COL_M, 0.45), border = alpha(COL_M, 0.7))
    scatter_panel(xf, yf, bg = alpha(COL_F, 0.45), border = alpha(COL_F, 0.7), add = TRUE)
    r_label(cor(xm, ym, use = "complete.obs"), COL_M, line = 0)
    r_label(cor(xf, yf, use = "complete.obs"), COL_F, line = 1)
    if (ch == "1")
      legend("bottomright", bty = "n", legend = c("Males", "Females"), pch = 21,
             pt.cex = 1.3, pt.lwd = 0.3, col = c(COL_M, COL_F),
             pt.bg = alpha(c(COL_M, COL_F), 0.6), text.col = c(COL_M, COL_F))
  }
  mtext(eff_lab("wgt"),  side = 1, outer = TRUE, line = 1.5, cex = 1.1)
  mtext(eff_lab("life"), side = 2, outer = TRUE, line = 1.2, cex = 1.1)
})

## --- Fig. S21, S22: male vs female additive effects, genome-wide -----------
## --- Fig. S23, S24: the same per chromosome ---------------------------------
## Autosomes only: the X chromosome has no single female additive effect
## (forward and reverse females differ) and males are hemizygous.
mf_lab <- function(tr, who) sprintf("Additive effect, %s (%s)", who, TRAIT_UNIT[[tr]])
for (tr in names(TRAITS)) {
  e  <- effects[[tr]]
  au <- e$all$chr %in% AUTOSOMES
  lim <- lim_sym(e$m$a[au], e$f$a[au])
  num <- c(wgt = "S21", life = "S22")[[tr]]
  save_pdf(sprintf("Fig%s_%s_male_vs_female_effects.pdf", num, tolower(TRAIT_NAME[[tr]])), 6.5, 6.5, {
    par(mar = c(4.5, 4.5, 2.5, 1))
    scatter_panel(e$f$a[au], e$m$a[au], lim, lim, main = TRAIT_NAME[[tr]],
                  xlab = mf_lab(tr, "females"), ylab = mf_lab(tr, "males"), cex = 0.7)
    r_label(cor(e$f$a[au], e$m$a[au]))
  })
  num <- c(wgt = "S23", life = "S24")[[tr]]
  save_pdf(sprintf("Fig%s_%s_male_vs_female_effects_by_chromosome.pdf", num, tolower(TRAIT_NAME[[tr]])), 8, 8.4, {
    par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(3, 3, 2.5, 0), mgp = c(2, 0.6, 0))
    for (ch in AUTOSOMES) {
      k <- e$all$chr == ch
      scatter_panel(e$f$a[k], e$m$a[k], lim, lim, main = paste("Chromosome", ch))
      r_label(cor(e$f$a[k], e$m$a[k], use = "complete.obs"))
    }
    mtext(TRAIT_NAME[[tr]], side = 3, outer = TRUE, line = 0.8, cex = 1.2, font = 2)
    mtext(mf_lab(tr, "females"), side = 1, outer = TRUE, line = 1.5, cex = 1.1)
    mtext(mf_lab(tr, "males"),   side = 2, outer = TRUE, line = 1.2, cex = 1.1)
  })
}

## --- Fig. S26: sex-averaged LOD profiles with permutation thresholds -------
save_pdf("FigS26_LOD_sex_averaged.pdf", 10, 6.5, {
  par(mfrow = c(2, 1), mar = c(3.5, 4.5, 1, 1), oma = c(1.5, 0, 0, 0))
  for (tr in names(TRAITS)) {
    s <- scans[[tr]]$all; t <- thr[[paste0(tr, "_all")]]
    g <- genome_axis(s$chr, s$pos)
    top <- max(s$lod, t) * 1.08; tick <- top * 0.03
    genome_frame(g, ylim = c(-tick, top), ylab = sprintf("LOD score (%s)", TRAIT_NAME[[tr]]))
    threshold_lines(g, t, col = COL_F, lwd = 1.5)
    lod_curve(g, s$lod, COL_ALL, lwd = 1.5)
    marker_rug(g, s, par("usr")[3], tick)
    legend("top", bty = "n", lwd = 1.5, col = COL_F, cex = 0.85,
           legend = "5% genome-wide threshold")
  }
  mtext("Chromosome", side = 1, outer = TRUE, line = 0.3)
})

## --- Fig. S27: sex-averaged LOD, each position coloured by its additive effect
## Autosomes only: the X chromosome has no sex-averaged additive effect.
save_pdf("FigS27_LOD_coloured_by_additive_effect.pdf", 10, 7.5, {
  layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(6, 0.7))
  for (tr in names(TRAITS)) {
    au <- scans[[tr]]$all$chr %in% AUTOSOMES
    s <- scans[[tr]]$all[au, ]; a <- effects[[tr]]$all$a[au]
    g <- genome_axis(s$chr, s$pos)
    amax <- max(abs(a))
    brks <- seq(-amax, amax, length.out = length(ADD_RAMP) + 1)
    fill <- ADD_RAMP[cut(a, brks, include.lowest = TRUE)]
    par(mar = c(4, 4.5, 2.5, 0.5))
    genome_frame(g, ylim = c(0, max(s$lod) * 1.05), ylab = "LOD score",
                 main = TRAIT_NAME[[tr]])
    points(g$x, s$lod, pch = 21, bg = fill, col = "grey20", lwd = 0.15, cex = 0.7)
    mtext("Chromosome", side = 1, line = 2.3, cex = 0.8)
    ## colour bar
    par(mar = c(4, 0.5, 2.5, 4))
    mids <- (brks[-1] + brks[-length(brks)]) / 2
    image(x = 1, y = mids, z = matrix(mids, nrow = 1), col = ADD_RAMP,
          axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1, cex.axis = 0.8); box()
    mtext(eff_lab(tr), side = 4, line = 2.8, cex = 0.7)
  }
})

## ===========================================================================
## 7. Tables S1 and S2: parameter estimates of the multiple-QTL models
## ===========================================================================
## Additive ('a') and dominance ('d') effects of each QTL, the sex effect (sex
## coded 1 = male, so negative = males lighter / shorter-lived) and, for
## lifespan, QTL x sex interactions. Significance: two-sided t-test on the
## residual degrees of freedom, p < 0.05.
mqm_table <- function(tr) {
  spec <- MQM[[tr]]
  fit  <- suppressWarnings(fitqtl(seed, pheno.col = TRAITS[[tr]], qtl = mqm[[tr]],
                                  covar = sexdf, formula = spec$fit_formula,
                                  method = "imp", get.ests = TRUE))
  est <- fit$ests$ests
  se  <- sqrt(diag(fit$ests$covar))
  t   <- est / se
  p   <- 2 * pt(-abs(t), df = fit$result.full["Error", "df"])
  out <- data.frame(Parameter = names(est),
                    Estimate  = round(est, spec$digits),
                    SE        = round(se, spec$digits),
                    t         = round(t, 3),
                    p         = signif(p, 3),
                    Significant = ifelse(p < 0.05, "*", ""),
                    row.names = NULL, stringsAsFactors = FALSE)
  names(out)[2] <- sprintf("Estimate (%s, %s)", TRAIT_NAME[[tr]], TRAIT_UNIT[[tr]])
  out
}
for (tr in names(MQM)) {
  tab  <- mqm_table(tr)
  num  <- c(wgt = "S1", life = "S2")[[tr]]
  file <- file.path(TAB_DIR, sprintf("Table%s_%s_MQM_parameters.csv", num, tolower(TRAIT_NAME[[tr]])))
  write.csv(tab, file, row.names = FALSE)
  print(tab, row.names = FALSE)
}

