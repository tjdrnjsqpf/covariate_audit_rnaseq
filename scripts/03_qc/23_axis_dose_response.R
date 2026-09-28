#!/usr/bin/env Rscript
## Are the alternative axes (quality, composition, globin, batch) confounders "in their own right": for each axis, the change when adjusting for that axis (relative to M1), by the axis's own imbalance class
suppressMessages(library(data.table)); med <- function(x) round(median(x, na.rm=TRUE), 3)
setwd("/var2/lsg/Claude_Code/covariate_audit_rnaseq")
for (org in c("human","mouse")) {
  d <- fread(sprintf("results/main_%s_summary_v2.tsv", org)); A <- d[study %in% readLines(sprintf("config/c2_analysis_set_%s.txt", org))]
  cat(sprintf("\n################ %s (analysis set %d) ################\n", org, nrow(A)))
  axes <- list(
    quality     = list(x = pmax(abs(A$qc_mito_smd), abs(A$qc_unique_smd), na.rm=TRUE), j = A$M4v1_jaccard, jp = if ("M4permv1_jaccard" %in% names(A)) A$M4permv1_jaccard else NA_real_, cuts = c(0, 0.5, 1, 99), lab = "|quality SMD| (max of mito, unique-map)"),
    composition = list(x = pmax(abs(A$cell_a_smd), abs(A$cell_b_smd), na.rm=TRUE), j = A$M5v1_jaccard, jp = if ("M5permv1_jaccard" %in% names(A)) A$M5permv1_jaccard else NA_real_, cuts = c(0, 0.5, 1, 99), lab = "|composition SMD|"),
    globin      = list(x = abs(A$globin_smd), j = A$M7v1_jaccard, jp = if ("M7permv1_jaccard" %in% names(A)) A$M7permv1_jaccard else NA_real_, cuts = c(0, 0.5, 1, 99), lab = "|globin SMD| (globin>=1% studies)"),
    batch       = list(x = A$batch_v_max, j = A$M6v1_jaccard, jp = if ("M6permv1_jaccard" %in% names(A)) A$M6permv1_jaccard else NA_real_, cuts = c(0, 0.3, 0.6, 1.01), lab = "batch Cramer V"))
  for (nm in names(axes)) { a <- axes[[nm]]; ok <- !is.na(a$x) & !is.na(a$j)
    b <- cut(a$x[ok], a$cuts, include.lowest=TRUE)
    cat(sprintf("\n== %s: by %s -> Jaccard(M1 vs M1+axis); lower = adjusting for the axis changes more ==\n", nm, a$lab))
    jp <- if (length(a$jp) == 1 && is.na(a$jp)) rep(NA_real_, sum(ok)) else a$jp[ok]
    print(data.table(bin = b, J = a$j[ok], Jp = jp, sexJ = A$M1_jaccard[ok])[, .(n = .N, J_axis_vs_M1 = med(J), J_perm_vs_M1 = med(Jp), dJ = med(J - Jp), J_sex_vs_M0 = med(sexJ)), by=bin][order(bin)])
    cat(sprintf("  Spearman(axis imbalance, J) = %.3f  (negative = more imbalance, more change)\n", cor(a$x[ok], a$j[ok], method="spearman")))
  }
}
