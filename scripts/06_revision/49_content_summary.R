#!/usr/bin/env Rscript
## Summary of 46 output: by class, (i) p-value of the real labels against the permutation distribution, (ii) content of the dropped set (real vs within-group permutation overlap, and its null), sex-chromosome fraction, (iii) Jaccard excluding XY, (iv) pathway-level p-values.
## Usage: LC_ALL=en_US.UTF-8 Rscript scripts/06_revision/49_content_summary.R [PROJ]
suppressMessages(library(data.table)); PROJ <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(PROJ)) PROJ <- getwd(); setwd(PROJ)
lv <- c("balanced","mildly_imbalanced","partially_confounded")
fs <- list.files("results/revision/content", pattern="^[A-Z]+[0-9]+\\.tsv$", full.names=TRUE); d <- rbindlist(lapply(fs, fread), fill=TRUE)
cls <- rbindlist(lapply(c("human","mouse"), function(o) fread(sprintf("results/main_%s_summary_v2.tsv", o))[, .(study, design_sex, sex_diff, tissue_category)]))
d <- merge(d, cls, by="study")[design_sex %in% lv]
## 2026-09-23 (review J A1 / L 1): instead of the p = k/B saved by 46, recompute p = (k+1)/(B+1) from the per-permutation files (k = number of permutations at or below the real labels); under the null, p ≤ 0.05 has probability exactly 1/(B+1)
pf <- rbindlist(lapply(list.files("results/revision/content", pattern="_perms[.]tsv$", full.names=TRUE), fread), fill=TRUE)
for (k in c("p_J_within","p_J_across","p_tsp_within","p_PJ_within","p_PJ_across","B_within","B_across")) d[, (k) := NULL]
pp <- merge(pf, d[, .(study, J_M1, tsp_M1, PJ_M1)], by="study")
pr <- pp[, .(B=.N, p_J=(1 + sum(J <= J_M1)) / (.N + 1), p_tsp=(1 + sum(t_sp <= tsp_M1)) / (.N + 1), p_PJ=if (all(is.na(PJ)) || is.na(PJ_M1[1])) NA_real_ else (1 + sum(PJ <= PJ_M1, na.rm=TRUE)) / (sum(!is.na(PJ)) + 1)), by=.(study, tag)]
prw <- dcast(pr, study ~ tag, value.var=c("B","p_J","p_tsp","p_PJ"))
d <- merge(d, prw, by="study", all.x=TRUE)
setnames(d, c("p_J_within","p_J_across","p_tsp_within","p_tsp_across","p_PJ_within","p_PJ_across"), c("p_J_within","p_J_across","p_tsp_within","p_tsp_across","p_PJ_within","p_PJ_across"), skip_absent=TRUE)
## 2026-09-23: if CONTENT_VALID_ONLY=1, only studies without interleaved groups (group runs <= 2) — provisional summary before rerunning after the 46 row-order bug fix; output suffix _valid
SUF <- ""; if (nzchar(Sys.getenv("CONTENT_VALID_ONLY")) && file.exists("results/revision/content_group_runs_all.tsv")) { gr <- fread("results/revision/content_group_runs_all.tsv", header=FALSE, col.names=c("study","org","kind","runs")); d <- d[study %in% gr[runs <= 2, study]]; SUF <- "_valid" }
bci <- function(x) { x <- x[!is.na(x)]; if (length(x) < 3) return(NA_character_); set.seed(1); b <- replicate(2000, median(sample(x, replace=TRUE))); sprintf("%.3f [%.3f, %.3f]", median(x), quantile(b,.025), quantile(b,.975)) }
S <- d[, .(n=.N, B=as.numeric(median(B_within)), B_across_lt10=sum(B_across < 10, na.rm=TRUE), k_p_within=sum(p_J_within <= 1/(B_within+1) + 1e-9, na.rm=TRUE), k_p_across=sum(p_J_across <= 1/(B_across+1) + 1e-9, na.rm=TRUE), null_rate=round(1/(median(B_within)+1), 3),
  J_M1=bci(J_M1), J_within=bci(J_within_mean), J_across=bci(J_across_mean),
  d_within=bci(J_M1 - J_within_mean), d_within_noxy=bci(J_M1_noxy - J_within_noxy_mean),
  frac_p_within_le05=round(mean(p_J_within <= 1/(B_within+1) + 1e-9, na.rm=TRUE), 3), frac_p_across_le05=round(mean(p_J_across <= 1/(B_across+1) + 1e-9, na.rm=TRUE), 3),
  tsp_M1=bci(tsp_M1), tsp_within=bci(tsp_within_mean), frac_p_tsp_le05=round(mean(p_tsp_within <= 1/(B_within+1) + 1e-9, na.rm=TRUE), 3), frac_p_tsp_across_le05=round(mean(p_tsp_across <= 1/(B_across+1) + 1e-9, na.rm=TRUE), 3), k_p_tsp=sum(p_tsp_within <= 1/(B_within+1) + 1e-9, na.rm=TRUE),
  lost_overlap_real_vs_within=bci(J_lost_M1_vs_within), lost_overlap_null=bci(J_lost_within_vs_within), lost_overlap_ratio=bci(J_lost_M1_vs_within / J_lost_within_vs_within),
  gained_overlap_real_vs_within=bci(J_gained_M1_vs_within), gained_overlap_null=bci(J_gained_within_vs_within),
  xy_lost_M1=bci(xy_lost_M1), xy_lost_within=bci(xy_lost_within_mean), xy_gained_M1=bci(xy_gained_M1), xy_gained_within=bci(xy_gained_within_mean),
  lfc_shift_M1=bci(lfc_shift_M1), lfc_shift_within=bci(lfc_shift_within_mean), xy_share_shift_M1=bci(xy_share_shift_M1), xy_share_shift_within=bci(xy_share_shift_within_mean),
  PJ_M1=bci(PJ_M1), PJ_within=bci(PJ_within_mean), frac_p_PJ_within_le05=round(mean(p_PJ_within <= 1/(B_within+1) + 1e-9, na.rm=TRUE), 3), frac_p_PJ_across_le05=round(mean(p_PJ_across <= 1/(B_across+1) + 1e-9, na.rm=TRUE), 3), k_p_PJ=sum(p_PJ_within <= 1/(B_within+1) + 1e-9, na.rm=TRUE), n_PJ=sum(!is.na(p_PJ_within)),
  top1_changed_M1=round(mean(top1_changed_M1, na.rm=TRUE), 3), top1_changed_within=round(mean(top1_changed_within_frac, na.rm=TRUE), 3), top1_changed_across=round(mean(top1_changed_across_frac, na.rm=TRUE), 3)),
  by=.(organism, design_sex)][order(organism, match(design_sex, lv))]
fwrite(S, sprintf("results/revision/C_content_by_class%s.tsv", SUF), sep="\t"); fwrite(d, sprintf("results/revision/C_content_per_study%s.tsv", SUF), sep="\t")
options(width=300); print(t(S)); cat("studies:", nrow(d), "\n")
