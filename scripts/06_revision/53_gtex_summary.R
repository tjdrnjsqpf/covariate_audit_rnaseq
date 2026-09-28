#!/usr/bin/env Rscript
## Summary of 52 output: by class, (i) enrichment (OR) of GTEx autosomal sex-biased genes in the unadjusted/adjusted DE sets, (ii) enrichment of the lost/gained sets — real vs within-group vs across-group,
## (iii) share of sex-biased genes in the logFC shift, paired comparisons. Usage: Rscript scripts/06_revision/53_gtex_summary.R → results/revision/G_gtex_by_class.tsv
suppressMessages(library(data.table)); set.seed(1)
R <- fread("results/revision/G_gtex_per_study.tsv"); lv <- c("balanced","mildly_imbalanced","partially_confounded")
bci <- function(x) { x <- x[is.finite(x)]; if (length(x) < 3) return(NA_character_); b <- replicate(2000, median(sample(x, replace=TRUE))); sprintf("%.2f [%.2f, %.2f]", median(x), quantile(b,.025), quantile(b,.975)) }
pw <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 5) return(NA_real_); signif(wilcox.test(a[ok], b[ok], paired=TRUE)$p.value, 2) }
S <- R[, .(n=.N, n_matched=sum(set_kind != "shared"), sb_frac_universe=bci(frac_U_sb_auto), n_sb_auto=median(n_sb_auto),
  or_de0=bci(or_de0_sb), or_de1=bci(or_de1_sb), p_de0_vs_de1=pw(or_de0_sb, or_de1_sb),
  lost_frac_real=bci(frac_lost_real), lost_frac_within=bci(frac_lost_w), lost_frac_across=bci(frac_lost_a),
  or_lost_real=bci(or_lost_real), or_lost_within=bci(or_lost_w), or_lost_across=bci(or_lost_a), p_lost_real_vs_within=pw(or_lost_real, or_lost_w),
  or_gain_real=bci(or_gain_real), or_gain_within=bci(or_gain_w), p_gain_real_vs_within=pw(or_gain_real, or_gain_w),
  shift_share_real=bci(shift_share_real), shift_share_within=bci(shift_share_w), shift_share_across=bci(shift_share_a), p_share_real_vs_within=pw(shift_share_real, shift_share_w),
  shift_ratio_real=bci(shift_ratio_real), shift_ratio_within=bci(shift_ratio_w), p_ratio_real_vs_within=pw(shift_ratio_real, shift_ratio_w)),
  by=design_sex][order(match(design_sex, lv))]
fwrite(S, "results/revision/G_gtex_by_class.tsv", sep="\t"); options(width=250); print(t(S))
## tissue-matched sets only (excluding shared), by class
cat("\n## tissue-matched only (set_kind != shared):\n"); M <- R[set_kind != "shared", .(n=.N, or_lost_real=bci(or_lost_real), or_lost_within=bci(or_lost_w), p=pw(or_lost_real, or_lost_w), shift_share_real=bci(shift_share_real), shift_share_within=bci(shift_share_w), p2=pw(shift_share_real, shift_share_w), or_de0=bci(or_de0_sb)), by=design_sex][order(match(design_sex, lv))]; print(t(M))
