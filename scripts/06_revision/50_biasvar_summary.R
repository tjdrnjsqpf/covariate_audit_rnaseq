#!/usr/bin/env Rscript
## Summary of 40 v2 output: error relative to reference A (full-sample M1), reference C (out-of-sample) and reference B (balanced draw), plus logFC bias²/variance decomposition, per step.
## Usage: LC_ALL=en_US.UTF-8 Rscript scripts/06_revision/50_biasvar_summary.R [PROJ]
suppressMessages(library(data.table)); PROJ <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(PROJ)) PROJ <- getwd(); setwd(PROJ)
fs <- list.files("results/revision/synthetic", pattern="^[A-Z]+[0-9]+\\.tsv$", full.names=TRUE); d <- rbindlist(lapply(fs, fread), fill=TRUE)
bv <- rbindlist(lapply(list.files("results/revision/synthetic", pattern="_biasvar\\.tsv$", full.names=TRUE), fread), fill=TRUE)
## 2026-09-23 (review J B2): the bias² estimate from the mean of 10 draws includes variance/10, so corrected value bias2c = bias2 − variance/10; rms_bias = sqrt(max(bias2c,0))
bv[, bias2c := bias2 - variance/10]; bv[, rms_bias := sqrt(pmax(bias2c, 0))]; bv[, msec := bias2c + variance]
bci <- function(x) { x <- x[!is.na(x)]; if (length(x) < 3) return(NA_character_); set.seed(1); b <- replicate(2000, median(sample(x, replace=TRUE))); sprintf("%.3f [%.3f, %.3f]", median(x), quantile(b,.025), quantile(b,.975)) }
cols <- intersect(c("J_M1","J_M0perm","J_M0permw","M0_precision","M1_precision","M0_recall","M1_recall","M0_t_sp_ref","M1_t_sp_ref","M0_precision_C","M1_precision_C","M0_recall_C","M1_recall_C","M0_t_sp_C","M1_t_sp_C","J_M0_vs_balanced","J_M1_vs_balanced","ref_share","n_comp"), names(d))
s <- d[, lapply(.SD, mean, na.rm=TRUE), by=.(organism, study, target_df), .SDcols=cols]
E <- s[, c(list(n=.N), lapply(.SD, bci)), by=.(organism, target_df), .SDcols=cols][order(organism, target_df)]
Bv <- bv[, .(n=.N, bias2=bci(bias2), bias2_corrected=bci(bias2c), rms_bias=bci(rms_bias), variance=bci(variance), mse=bci(mse), signed_bias=bci(signed_bias), rmse_per_draw=bci(rmse_per_draw)), by=.(organism, target_df, model)][order(organism, target_df, model)]
## Paired comparison: M1 − M0 within the same study and step (mse, bias2, variance)
w <- dcast(bv, organism + study + target_df ~ model, value.var=c("bias2","bias2c","rms_bias","variance","mse","msec"))
Pd <- w[, .(n=.N, d_bias2=bci(bias2_M1 - bias2_M0), d_bias2c=bci(bias2c_M1 - bias2c_M0), rms_bias_M0=bci(rms_bias_M0), rms_bias_M1=bci(rms_bias_M1), d_msec=bci(msec_M1 - msec_M0), p_msec=signif(suppressWarnings(wilcox.test(msec_M1, msec_M0, paired=TRUE)$p.value), 3), d_variance=bci(variance_M1 - variance_M0), d_mse=bci(mse_M1 - mse_M0), p_mse=signif(suppressWarnings(wilcox.test(mse_M1, mse_M0, paired=TRUE)$p.value), 3), p_bias2=signif(suppressWarnings(wilcox.test(bias2_M1, bias2_M0, paired=TRUE)$p.value), 3)), by=.(organism, target_df)][order(organism, target_df)]
fwrite(E, "results/revision/S2_error_by_rung.tsv", sep="\t"); fwrite(Bv, "results/revision/S2_biasvar_by_rung.tsv", sep="\t"); fwrite(Pd, "results/revision/S2_biasvar_paired.tsv", sep="\t")
options(width=300); cat("studies:", uniqueN(d$study), "\n"); print(t(E)); print(Bv); print(Pd)
