#!/usr/bin/env Rscript
## C3–C7 aggregation: how SVA (M3), quality (M4), cell composition (M5, tissue-specific panels), batch (M6) and globin (M7) overlap with sex adjustment (M1)
suppressMessages(library(data.table)); ORG <- commandArgs(trailingOnly=TRUE)[1]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
d <- fread(sprintf("results/main_%s_summary_v2.tsv", ORG)); A <- d[study %in% readLines(sprintf("config/c2_analysis_set_%s.txt", ORG))]
lv <- c("balanced","mildly_imbalanced","partially_confounded"); med <- function(x) round(median(x, na.rm=TRUE), 3)
cat(sprintf("########## Extra-model aggregation (%s, analysis set %d) ##########\n", ORG, nrow(A)))
if ("M3_jaccard" %in% names(A)) {
  cat("\n== C5 SVA (M3) ==  M3=M0 vs SVA; M3v1=M1 vs SVA (high means SVA≈sex adjustment); sv_sex_r2=max R² of sex among the SVs\n")
  print(A[!is.na(M3_jaccard), .(n = .N, n_sv = med(n_sv), sv_sex_r2 = med(sv_sex_r2), frac_r2_gt_0.5 = round(mean(sv_sex_r2 > 0.5, na.rm=TRUE), 2), J_M1 = med(M1_jaccard), J_perm = med(M0perm_jaccard), J_SVA = med(M3_jaccard), J_M1_vs_SVA = med(M3v1_jaccard)), by=design_sex][order(match(design_sex, lv))]) }
if ("M4_jaccard" %in% names(A)) {
  cat("\n== C4 RNA quality ==\n")
  print(A[, .(n = .N, mito_smd = med(abs(qc_mito_smd)), unique_smd = med(abs(qc_unique_smd)), intron_smd = med(abs(qc_intron_smd)), frac_any_qc_smd_gt_0.8 = round(mean(pmax(abs(qc_mito_smd), abs(qc_unique_smd), abs(qc_intron_smd), na.rm=TRUE) > 0.8, na.rm=TRUE), 2), J_M1_vs_M4 = med(M4v1_jaccard)), by=design_sex][order(match(design_sex, lv))])
  cat(sprintf("  sex_diff vs |mito SMD| Spearman: %.3f\n", cor(A$sex_diff, abs(A$qc_mito_smd), method="spearman", use="complete.obs"))) }
if ("cell_panel" %in% names(A)) {
  cat("\n== C3 cell composition (tissue-specific panels; a=parenchymal/specific cells, b=immune/stromal) ==\n")
  print(A[!is.na(M5_jaccard), .(n = .N, a_smd = med(abs(cell_a_smd)), b_smd = med(abs(cell_b_smd)), J_M1 = med(M1_jaccard), J_M0_vs_M5 = med(M5_jaccard), J_M1_vs_M5 = med(M5v1_jaccard)), by=cell_panel][order(-n)])
  cat("  by sex grade (M1 vs M5: how much composition changes on top of sex):\n")
  print(A[!is.na(M5_jaccard), .(n = .N, J_M1 = med(M1_jaccard), J_perm = med(M0perm_jaccard), J_M0_vs_M5 = med(M5_jaccard), J_M1_vs_M5 = med(M5v1_jaccard)), by=design_sex][order(match(design_sex, lv))])
  cat(sprintf("  sex_diff vs |cell_b SMD| Spearman: %.3f\n", cor(A$sex_diff, abs(A$cell_b_smd), method="spearman", use="complete.obs"))) }
if ("batch_v_max" %in% names(A)) {
  cat("\n== C6 technical batch (analysis set) ==\n")
  cat(sprintf("  batch structure present %d (%.0f%%), of which V≥0.5 %d, V=1 %d;  most entangled variable: %s\n", sum(!is.na(A$batch_v_max)), 100*mean(!is.na(A$batch_v_max)), sum(A$batch_v_max >= 0.5, na.rm=TRUE), sum(A$batch_v_max >= 0.999, na.rm=TRUE), paste(names(table(A$batch_var)), table(A$batch_var), collapse=", ")))
  print(A[!is.na(M6_jaccard), .(n = .N, batch_v = med(batch_v_max), J_M1 = med(M1_jaccard), J_perm = med(M0perm_jaccard), J_M0_vs_M6 = med(M6_jaccard), J_M1_vs_M6 = med(M6v1_jaccard)), by=design_sex][order(match(design_sex, lv))])
  cat(sprintf("  sex_diff vs batch_v_max Spearman: %.3f;  intron range (%%p, library heterogeneity) median %.1f, studies >10%%p %.0f%%\n", cor(A$sex_diff, A$batch_v_max, method="spearman", use="complete.obs"), med(A$intron_range), 100*mean(A$intron_range > 10, na.rm=TRUE))) }
if ("globin_max" %in% names(A)) {
  G <- A[globin_max >= 0.01]
  cat(sprintf("\n== C7 globin (studies with max fraction ≥1%%: %d, mostly blood) ==\n", nrow(G)))
  print(G[, .(n = .N, globin_max = med(globin_max), globin_range = med(globin_range), globin_smd = med(abs(globin_smd)), J_M1 = med(M1_jaccard), J_M0_vs_M7 = med(M7_jaccard), J_M1_vs_M7 = med(M7v1_jaccard)), by=design_sex][order(match(design_sex, lv))])
  cat(sprintf("  studies with globin range >10%%p (suspected mixed depletion status): %d / %d\n", sum(G$globin_range > 0.1), nrow(G))) }
