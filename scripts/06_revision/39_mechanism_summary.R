#!/usr/bin/env Rscript
## E10 + E9b (local part). Summarises the per-study table produced on the server by 39a_server_extract_mechanism_leakage.R (frozen DE tables; no refitting).
##  PER_STUDY_HUMAN / PER_STUDY_MOUSE default to results/revision/E10_per_study_{org}.tsv
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
PSF <- c(human=Sys.getenv("PER_STUDY_HUMAN", file.path(OUTDIR, "E10_per_study_human.tsv")), mouse=Sys.getenv("PER_STUDY_MOUSE", file.path(OUTDIR, "E10_per_study_mouse.tsv")))
P <- rbindlist(lapply(ORGS, function(o) fread(PSF[[o]])), fill=TRUE); if ("error" %in% names(P)) { cat("studies with errors:", sum(!is.na(P$error) & P$error != ""), "\n"); P <- P[is.na(error) | error == ""] }
D <- load_summary(); X <- merge(D[, .(organism=as.character(organism), study, analysis, design_sex, sex_diff, M1_jaccard, M0perm_jaccard)], P, by=c("organism", "study"))
X <- X[design_sex %in% CLS]; X[, cls := factor(design_sex, levels=CLS)]; X[, organism := factor(organism, levels=ORGS)]; X[, sqrt_vif := sqrt(vif)]; X[, delta := M1_jaccard - M0perm_jaccard]
A <- X[analysis == TRUE]; cat("studies: all", nrow(X), " analysis set", nrow(A), "\n"); print(A[, .N, keyby=.(organism, cls)])
md <- function(x) as.double(median(x, na.rm=TRUE))
## (a) SE inflation vs sqrt(VIF)
SE <- A[, .(n_studies=.N, sqrt_vif_median=md(sqrt_vif), se_ratio_M1_median=md(se_ratio_M1_median), se_ratio_M1_q1=quantile(se_ratio_M1_median, .25), se_ratio_M1_q3=quantile(se_ratio_M1_median, .75), se_ratio_perm_median=md(se_ratio_perm_median),
            sigma_ratio_median=md(se_ratio_M1_median/sqrt_vif), pct_se_ratio_below_1=100*mean(se_ratio_M1_median < 1), pct_se_ratio_gt_1.1=100*mean(se_ratio_M1_median > 1.1), pct_se_ratio_gt_1.25=100*mean(se_ratio_M1_median > 1.25)), keyby=.(organism, cls)]; out(SE, "E10_se_inflation_by_class.tsv")
FT <- rbindlist(lapply(ORGS, function(o) { x <- A[organism == o & is.finite(vif)]; f <- summary(lm(log(se_ratio_M1_median) ~ log(sqrt_vif), x)); cf <- f$coefficients
  data.table(organism=o, n=nrow(x), spearman_se_ratio_vs_sqrtvif=cor(x$se_ratio_M1_median, x$sqrt_vif, method="spearman"), loglog_slope=cf[2, 1], slope_se=cf[2, 2], intercept=cf[1, 1], exp_intercept=exp(cf[1, 1]), r2=f$r.squared,
             spearman_delta_vs_se_ratio=cor(x$delta, x$se_ratio_M1_median, method="spearman", use="complete.obs"), spearman_delta_vs_dM1_med_abs_rel=cor(x$delta, x$dM1_med_abs/x$med_abs_l0, method="spearman", use="complete.obs")) })); out(FT, "E10_se_ratio_vs_sqrtvif_fit.tsv")
## (b) where does logFC(M1) − logFC(M0) sit? concentration metrics, real vs permuted sex
SH <- A[, .(n_studies=.N, rel_shift_M1=md(dM1_med_abs/med_abs_l0), rel_shift_perm=md(dPerm_med_abs/med_abs_l0), p99_rel_shift_M1=md(dM1_p99_abs/med_abs_l0), p99_rel_shift_perm=md(dPerm_p99_abs/med_abs_l0),
            top1pct_share_M1=md(dM1_top1pct_share_ss), top1pct_share_perm=md(dPerm_top1pct_share_ss), gini_M1=md(dM1_gini_abs), gini_perm=md(dPerm_gini_abs), kurt_M1=md(dM1_excess_kurtosis), kurt_perm=md(dPerm_excess_kurtosis),
            panel_over_bg_M1=md(dM1_panel_over_bg), panel_over_bg_perm=md(dPerm_panel_over_bg), panel_share_ss_M1=md(dM1_panel_share_ss), panel_share_ss_perm=md(dPerm_panel_share_ss), panel_pctile_M1=md(dM1_panel_mean_pctile), panel_pctile_perm=md(dPerm_panel_mean_pctile),
            cor_l0_l1=md(cor_l0_l1)), keyby=.(organism, cls)]; out(SH, "E10_logfc_shift_by_class.tsv")
## (b2) shift per unit of group–sex correlation: for OLS, logFC(M0) − logFC(M1) = δ × (estimated sex coefficient) and δ_perm/δ_real = r_perm/r_real (same margins), so
##      (median|shift| / |r|) is proportional to the typical estimated covariate effect; ratio real/perm > 1 means real sex carries more signal per gene than a random label (studies with |r|, |r_perm| > 0.05)
B2 <- A[abs(r_group_sex) > 0.05 & abs(r_group_sexperm) > 0.05, .(n_studies=.N, abs_r_real=md(abs(r_group_sex)), abs_r_perm=md(abs(r_group_sexperm)), shift_per_r_real=md(dM1_med_abs/abs(r_group_sex)), shift_per_r_perm=md(dPerm_med_abs/abs(r_group_sexperm)),
            ratio_real_over_perm=md((dM1_med_abs/abs(r_group_sex))/(dPerm_med_abs/abs(r_group_sexperm)))), keyby=.(organism, cls)]; out(B2, "E10_shift_per_unit_r.tsv")
## (c) why M0-DE genes are lost: logFC shrinkage vs SE inflation (studies with >= 5 lost genes)
LS <- A[!is.na(lost_med_log_t_ratio), .(n_studies=.N, lost_lfc_ratio=exp(md(lost_med_log_lfc_ratio)), lost_se_ratio=exp(md(lost_med_log_se_ratio)), lost_t_ratio=exp(md(lost_med_log_t_ratio)),
            share_of_log_t_drop_from_lfc=md(lost_med_log_lfc_ratio/lost_med_log_t_ratio), frac_genes_lfc_dominant=md(lost_frac_lfc_shrunk_gt_se_inflated)), keyby=.(organism, cls)]; out(LS, "E10_lost_gene_decomposition.tsv")
## (d) E9b: list genes among lost vs retained M0-DE genes; studies with >= 10 lost and >= 10 retained genes; random same-size gene-set control and permuted-model control
lk <- function(pre, minn=10) { v <- function(z) A[[paste0(pre, "_", z)]]; ok <- !is.na(v("n_lost")) & v("n_lost") >= minn & v("n_kept") >= minn & v("n_list_in_universe") > 0; if (!any(ok)) return(NULL)
  Y <- data.table(organism=A$organism, cls=A$cls, lost=v("share_lost"), kept=v("share_kept"), lostperm=v("share_lostperm"), rnd=v("random_share_lost_mean"), bg=v("bg_share"), k_lost=v("k_lost"), n_lost=v("n_lost"), k_kept=v("k_kept"), n_kept=v("n_kept"))[ok]
  Y[, { w <- suppressWarnings(wilcox.test(lost, kept, paired=TRUE)); .(list=pre, n_studies=.N, share_lost_median=100*md(lost), share_kept_median=100*md(kept), share_lostperm_median=100*md(lostperm), random_control_median=100*md(rnd), background_median=100*md(bg),
        pooled_share_lost=100*sum(k_lost)/sum(n_lost), pooled_share_kept=100*sum(k_kept)/sum(n_kept), pooled_ratio=(sum(k_lost)/sum(n_lost))/(sum(k_kept)/sum(n_kept)),
        pooled_ratio_lo=quantile({ set.seed(1); replicate(1000, { i <- sample.int(.N, replace=TRUE); (sum(k_lost[i])/sum(n_lost[i]))/(sum(k_kept[i])/sum(n_kept[i])) }) }, .025, na.rm=TRUE),
        pooled_ratio_hi=quantile({ set.seed(1); replicate(1000, { i <- sample.int(.N, replace=TRUE); (sum(k_lost[i])/sum(n_lost[i]))/(sum(k_kept[i])/sum(n_kept[i])) }) }, .975, na.rm=TRUE),
        studies_with_list_gene_lost=sum(k_lost > 0), studies_with_list_gene_kept=sum(k_kept > 0), enrichment_lost_vs_bg_median=md(lost/bg), enrichment_kept_vs_bg_median=md(kept/bg), pct_studies_lost_gt_kept=100*mean(lost > kept), paired_wilcoxon_p=w$p.value) }, keyby=.(organism, cls)] }
LK <- rbindlist(list(lk("xy"), if ("list_n_lost" %in% names(A)) lk("list")), fill=TRUE); out(LK, "E9_lost_vs_retained_genelist.tsv")
## (e) panel (X/Y) genes GAINED by the adjusted model (not M0-DE, M1-DE): the adjusted model itself calls sex-linked genes for the group contrast
GN <- A[, { g <- wilson(sum(xy_k_gained >= 1), .N); l <- wilson(sum(xy_k_lost >= 1), .N); .(n_studies=.N, k_gained=g$k, pct_studies_panel_gained_M1=g$pct, gained_lo=g$lo, gained_hi=g$hi, k_lost=l$k, pct_studies_panel_lost_M1=l$pct, lost_lo=l$lo, lost_hi=l$hi,
            se_ratio_panel_median=md(se_ratio_M1_panel_median), se_ratio_all_median=md(se_ratio_M1_median)) }, keyby=.(organism, cls)]; out(GN, "E9_xy_panel_gained_lost_by_class.tsv")
options(width=250); print(GN, digits=3); print(SE, digits=3); print(FT, digits=3); print(SH, digits=3); print(B2, digits=3); print(LS, digits=3); print(LK, digits=3)
