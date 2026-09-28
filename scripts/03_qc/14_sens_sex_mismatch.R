#!/usr/bin/env Rscript
## Sex-mismatch sensitivity: keep metadata sex (v2) vs replace with expression-inferred sex (expr) vs drop mismatched samples (drop)
suppressMessages(library(data.table)); ORG <- commandArgs(trailingOnly=TRUE)[1]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
grade <- function(x) fifelse(is.na(x), NA_character_, fifelse(x >= 0.99, "completely_confounded", fifelse(x >= 0.3, "partially_confounded", fifelse(x >= 0.15, "mildly_imbalanced", "balanced"))))
base <- fread(sprintf("results/main_%s_summary_v2.tsv", ORG))[, .(study, n0 = n, mism = n_sex_mismatch, sd0 = sex_diff, ds0 = design_sex, J0 = M1_jaccard, P0 = M0perm_jaccard, tS0 = M1_t_spearman)]
out <- base[mism > 0 & !is.na(J0)]
for (v in c("expr","drop")) {
  fs <- list.files(sprintf("results/sens_sexfix_%s_%s", v, ORG), "summary.tsv", recursive=TRUE, full.names=TRUE)
  if (!length(fs)) { cat(v, ": no results\n"); next }
  s <- rbindlist(lapply(fs, fread), fill=TRUE)[, .(study, n = n, sd = sex_diff, J = M1_jaccard, P = M0perm_jaccard, tS = M1_t_spearman)]
  setnames(s, c("n","sd","J","P","tS"), paste0(c("n","sd","J","P","tS"), "_", v)); out <- merge(out, s, by="study", all.x=TRUE)
}
fwrite(out, sprintf("results/qc/sens_sex_mismatch_%s.tsv", ORG), sep="\t")
cat(sprintf("########## Sex-mismatch sensitivity: %s — %d analysis-set studies with mismatch>0 ##########\n", ORG, nrow(out)))
for (v in c("expr","drop")) { J <- out[[paste0("J_", v)]]; P <- out[[paste0("P_", v)]]; sd <- out[[paste0("sd_", v)]]; if (is.null(J)) next
  ok <- !is.na(J) & !is.na(out$J0)
  cat(sprintf("\n[%s] %d comparable studies\n", v, sum(ok)))
  cat(sprintf("  M1 Jaccard median: v2 %.3f → %s %.3f   (Δ=J−perm: v2 %.3f → %.3f)\n", median(out$J0[ok]), v, median(J[ok]), median((out$J0 - out$P0)[ok], na.rm=TRUE), median((J - P)[ok], na.rm=TRUE)))
  cat(sprintf("  studies with |Jaccard change|>0.1: %d / %d;  studies whose design class changed: %d\n", sum(abs(J - out$J0)[ok] > 0.1), sum(ok), sum(grade(sd[ok]) != out$ds0[ok], na.rm=TRUE)))
  cat("  median Δ by class (v2 classes), v2 vs variant:\n")
  print(out[ok, .(n = .N, d_v2 = round(median(J0 - P0, na.rm=TRUE), 3), d_var = round(median(get(paste0("J_", v)) - get(paste0("P_", v)), na.rm=TRUE), 3)), by=ds0][order(ds0)])
}
cat("\nInterpretation: if Δ in both variants has the same direction and magnitude as in v2, metadata sex errors do not produce the result.\n")
