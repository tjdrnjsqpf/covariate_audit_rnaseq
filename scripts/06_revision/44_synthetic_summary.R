#!/usr/bin/env Rscript
## S summary: within-study synthetic imbalance. Per step (target |Δf|): (1) change metrics J(M1), J(full permutation), J(within-group permutation); (2) error metrics vs the balanced reference (full-sample M1).
## Replicates are averaged within each study, then median across studies with bootstrap CI. Usage: Rscript 44_synthetic_summary.R [PROJ]
suppressMessages(library(data.table)); PROJ <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(PROJ)) PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
fs <- list.files("results/revision/synthetic", pattern="\\.tsv$", full.names=TRUE); d <- rbindlist(lapply(fs, fread), fill=TRUE)
cat("studies:", uniqueN(d$study), " (", paste(d[, uniqueN(study), by=organism][, paste(organism, V1)], collapse=", "), ") rows:", nrow(d), " median m per arm:", median(d[, m_per_arm[1], by=study]$V1), "\n")
d[, `:=`(d_across = J_M1 - J_M0perm, d_within = J_M1 - J_M0permw, fd_M0 = 1 - M0_precision, fd_M1 = 1 - M1_precision)]
s <- d[, lapply(.SD, mean, na.rm=TRUE), by=.(organism, study, target_df), .SDcols=c("J_M1","J_M0perm","J_M0permw","d_across","d_within","M0_n_de","M1_n_de","M0_precision","M1_precision","fd_M0","fd_M1","M0_recall","M1_recall","M0_t_sp_ref","M1_t_sp_ref","M0_top100_ref","M1_top100_ref","M0_sexgene_de","M1_sexgene_de","n_ref_de")]
bci <- function(x, R=2000) { x <- x[!is.na(x)]; if (length(x) < 3) return(NA_character_); set.seed(1); b <- replicate(R, median(sample(x, replace=TRUE))); sprintf("%.3f [%.3f, %.3f]", median(x), quantile(b,.025), quantile(b,.975)) }
chg <- s[, .(n=.N, J_M1=bci(J_M1), J_across=bci(J_M0perm), J_within=bci(J_M0permw), delta_across=bci(d_across), delta_within=bci(d_within)), by=.(organism, target_df)][order(organism, target_df)]
err <- s[, .(n=.N, nDE_M0=round(median(M0_n_de)), nDE_M1=round(median(M1_n_de)), falseDE_frac_M0=bci(fd_M0), falseDE_frac_M1=bci(fd_M1), recall_M0=bci(M0_recall), recall_M1=bci(M1_recall),
             tsp_ref_M0=bci(M0_t_sp_ref), tsp_ref_M1=bci(M1_t_sp_ref), top100_ref_M0=bci(M0_top100_ref), top100_ref_M1=bci(M1_top100_ref), sexgene_M0=round(mean(M0_sexgene_de > 0),3), sexgene_M1=round(mean(M1_sexgene_de > 0),3)), by=.(organism, target_df)][order(organism, target_df)]
## Within-study paired comparison: increase in the false-DE fraction of M0 at each step relative to step 0, and the M1 − M0 difference
w <- dcast(s, organism + study ~ target_df, value.var=c("fd_M0","fd_M1","M0_t_sp_ref","M1_t_sp_ref","M0_n_de","M1_n_de"))
pair <- rbindlist(lapply(c("0.2","0.4","0.6","0.8"), function(k) w[, .(target_df=k, n=.N, fd_M0_rise=bci(get(paste0("fd_M0_",k)) - fd_M0_0), fd_M1_rise=bci(get(paste0("fd_M1_",k)) - fd_M1_0),
   fd_M1_minus_M0=bci(get(paste0("fd_M1_",k)) - get(paste0("fd_M0_",k))), tsp_M1_minus_M0=bci(get(paste0("M1_t_sp_ref_",k)) - get(paste0("M0_t_sp_ref_",k))),
   p_fd=signif(suppressWarnings(wilcox.test(get(paste0("fd_M1_",k)), get(paste0("fd_M0_",k)), paired=TRUE)$p.value),3), p_tsp=signif(suppressWarnings(wilcox.test(get(paste0("M1_t_sp_ref_",k)), get(paste0("M0_t_sp_ref_",k)), paired=TRUE)$p.value),3)), by=organism]))
fwrite(s, "results/revision/S_per_study_rung.tsv", sep="\t"); fwrite(chg, "results/revision/S_change_by_rung.tsv", sep="\t"); fwrite(err, "results/revision/S_error_by_rung.tsv", sep="\t"); fwrite(pair, "results/revision/S_paired.tsv", sep="\t")
options(width=260); cat("\n== Change metrics ==\n"); print(chg); cat("\n== Error vs balanced reference ==\n"); print(err); cat("\n== Within-study paired comparison ==\n"); print(pair[order(organism, target_df)])
