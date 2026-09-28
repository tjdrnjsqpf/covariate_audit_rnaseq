#!/usr/bin/env Rscript
## Regression Δ(=M1 Jaccard − M0perm Jaccard) ~ |Δf| (review D 9.4: pin the regression reported in the manuscript Results to code)
## Covariates: log n, log(number of M0 DE genes). Stratum: studies whose sex labels all come from metadata (n_sex_meta_missing == 0).
## Usage: Rscript 25_delta_regression.R [PROJ]  → results/qc/delta_regression.tsv
suppressMessages(library(data.table))
PROJ <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(PROJ)) PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
lv <- c("balanced","mildly_imbalanced","partially_confounded"); out <- list()
for (org in c("human","mouse")) {
  d <- fread(sprintf("results/main_%s_summary_v2.tsv", org)); A <- readLines(sprintf("config/c2_analysis_set_%s.txt", org))
  a <- d[study %in% A & design_sex %in% lv & !is.na(M1_jaccard) & !is.na(M0perm_jaccard)]
  a[, delta := M1_jaccard - M0perm_jaccard]
  fits <- list(all = lm(delta ~ sex_diff + log(n) + log(M1_n_deg_ref), a),
               metadata_only = lm(delta ~ sex_diff + log(n) + log(M1_n_deg_ref), a[n_sex_meta_missing == 0]))
  if ("M0permw_jaccard" %in% names(a)) { a[, delta_w := M1_jaccard - M0permw_jaccard]
    fits$within_group_control <- lm(delta_w ~ sex_diff + log(n) + log(M1_n_deg_ref), a[!is.na(delta_w)]) }
  for (k in names(fits)) { cf <- summary(fits[[k]])$coefficients
    out[[length(out)+1]] <- data.table(organism=org, stratum=k, n=nobs(fits[[k]]), beta_sexdiff=round(cf["sex_diff",1],4), se=round(cf["sex_diff",2],4), p=signif(cf["sex_diff",4],3),
                                       beta_logn=round(cf["log(n)",1],4), beta_logdeg=round(cf["log(M1_n_deg_ref)",1],4)) }
}
res <- rbindlist(out); dir.create("results/qc", showWarnings=FALSE); fwrite(res, "results/qc/delta_regression.tsv", sep="\t"); print(res)
