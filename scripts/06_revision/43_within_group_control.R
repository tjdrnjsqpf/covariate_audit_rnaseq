#!/usr/bin/env Rscript
## W summary: decomposition of the sex-adjustment effect against the full-permutation (M0perm) vs within-group-permutation (M0permw) baselines.
##   Δ_across = J(M1) − J(M0perm)   : only the degrees-of-freedom cost removed (existing metric)
##   Δ_within = J(M1) − J(M0permw)  : df + collinearity cost removed → the share attributable to the sex–expression link itself
##   collinearity cost = J(M0permw) − J(M0perm)
## Usage: Rscript 43_within_group_control.R [PROJ]  → results/revision/W_*.tsv
suppressMessages(library(data.table)); PROJ <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(PROJ)) PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
dir.create("results/revision", showWarnings=FALSE, recursive=TRUE); lv <- c("balanced","mildly_imbalanced","partially_confounded")
bci <- function(x, R=2000) { x <- x[!is.na(x)]; if (length(x) < 3) return(sprintf("%.3f", median(x))); set.seed(1); b <- replicate(R, median(sample(x, replace=TRUE))); sprintf("%.3f [%.3f, %.3f]", median(x), quantile(b,.025), quantile(b,.975)) }
m <- function(x) as.numeric(median(x, na.rm=TRUE)); cls <- list(); bins <- list(); deg <- list()
for (org in c("human","mouse")) { d <- fread(sprintf("results/main_%s_summary_v2.tsv", org)); A <- readLines(sprintf("config/c2_analysis_set_%s.txt", org))
  a <- d[study %in% A & design_sex %in% lv & !is.na(M1_jaccard) & !is.na(M0perm_jaccard) & !is.na(M0permw_jaccard)]
  a[, `:=`(d_across = M1_jaccard - M0perm_jaccard, d_within = M1_jaccard - M0permw_jaccard, collin = M0permw_jaccard - M0perm_jaccard)]
  cls[[org]] <- a[, .(organism=org, n=.N, median_n=m(n), J_real=round(m(M1_jaccard),3), J_perm_across=round(m(M0perm_jaccard),3), J_perm_within=round(m(M0permw_jaccard),3),
                      delta_across=bci(d_across), delta_within=bci(d_within), collinearity_cost=bci(collin),
                      p_wilcox_within=signif(suppressWarnings(wilcox.test(d_within)$p.value), 3)), by=design_sex][order(match(design_sex, lv))]
  a[, bin := cut(sex_diff, c(-0.01,0.1,0.2,0.3,0.4,0.5,0.7,1), labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.7","0.7-1"))]
  bins[[org]] <- a[, .(organism=org, n=.N, delta_across=bci(d_across), delta_within=bci(d_within)), by=bin][order(bin)]
  deg[[org]] <- a[, .(organism=org, n=.N, deg_M0=m(M1_n_deg_ref), deg_M1=m(M1_n_deg_alt), deg_perm_across=m(M0perm_n_deg_alt), deg_perm_within=m(M0permw_n_deg_alt),
                      drop_real=round(1-m(M1_n_deg_alt)/m(M1_n_deg_ref),3), drop_across=round(1-m(M0perm_n_deg_alt)/m(M1_n_deg_ref),3), drop_within=round(1-m(M0permw_n_deg_alt)/m(M1_n_deg_ref),3),
                      tsp_real=round(m(M1_t_spearman),3), tsp_within=round(m(M0permw_t_spearman),3), top100_real=round(m(M1_top100_overlap),2), top100_within=round(m(M0permw_top100_overlap),2),
                      sexgene_M0_any=round(mean(M1_sexgene_deg_ref > 0),3), sexgene_within_any=round(mean(M0permw_sexgene_deg_alt > 0),3), sexgene_M1_any=round(mean(M1_sexgene_deg_alt > 0),3)), by=design_sex][order(match(design_sex, lv))] }
fwrite(rbindlist(cls), "results/revision/W_class.tsv", sep="\t"); fwrite(rbindlist(bins), "results/revision/W_bins.tsv", sep="\t"); fwrite(rbindlist(deg), "results/revision/W_deg_counts.tsv", sep="\t")
options(width=250); cat("== By class ==\n"); print(rbindlist(cls)); cat("\n== By |Δf| bin ==\n"); print(rbindlist(bins)); cat("\n== DE counts and other metrics ==\n"); print(rbindlist(deg))
