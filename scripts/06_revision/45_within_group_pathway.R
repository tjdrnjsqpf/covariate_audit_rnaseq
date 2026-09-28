#!/usr/bin/env Rscript
## W (pathway level): compare the Jaccard of significant Hallmark pathway sets and top-pathway changes across real sex / full permutation / within-group permutation.
## Usage: Rscript 45_within_group_pathway.R [PROJ] → results/revision/W_pathway_class.tsv
suppressMessages(library(data.table)); PROJ <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(PROJ)) PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
lv <- c("balanced","mildly_imbalanced","partially_confounded"); out <- list()
bci <- function(x, R=2000) { x <- x[!is.na(x)]; if (length(x) < 3) return(NA_character_); set.seed(1); b <- replicate(R, median(sample(x, replace=TRUE))); sprintf("%.3f [%.3f, %.3f]", median(x), quantile(b,.025), quantile(b,.975)) }
for (org in c("human","mouse")) { g <- fread(sprintf("results/c2/gsea_summary_%s.tsv", org)); d <- fread(sprintf("results/main_%s_summary_v2.tsv", org))[, .(study, design_sex, sex_diff)]
  a <- merge(g, d, by="study")[design_sex %in% lv & !is.na(P1_H_jaccard) & !is.na(P0perm_H_jaccard) & !is.na(P0permw_H_jaccard)]
  out[[org]] <- a[, .(organism=org, n=.N, PJ_real=round(median(P1_H_jaccard),3), PJ_across=round(median(P0perm_H_jaccard),3), PJ_within=round(median(P0permw_H_jaccard),3),
                      d_across=bci(P1_H_jaccard - P0perm_H_jaccard), d_within=bci(P1_H_jaccard - P0permw_H_jaccard),
                      changed_vs_across=round(mean(P1_H_jaccard < P0perm_H_jaccard - 0.2),3), changed_vs_within=round(mean(P1_H_jaccard < P0permw_H_jaccard - 0.2),3), within_changed_vs_across=round(mean(P0permw_H_jaccard < P0perm_H_jaccard - 0.2),3),
                      top1_changed_real=round(mean(!P1_H_top1_same, na.rm=TRUE),3), top1_changed_across=round(mean(!P0perm_H_top1_same, na.rm=TRUE),3), top1_changed_within=round(mean(!P0permw_H_top1_same, na.rm=TRUE),3), n_top1=sum(!is.na(P1_H_top1_same))), by=design_sex][order(match(design_sex, lv))] }
res <- rbindlist(out); fwrite(res, "results/revision/W_pathway_class.tsv", sep="\t"); options(width=250); print(res)
