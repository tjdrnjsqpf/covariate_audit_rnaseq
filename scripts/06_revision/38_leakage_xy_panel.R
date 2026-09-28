#!/usr/bin/env Rscript
## E9a. Sex-chromosome leakage from the summary table. Panel = config/sex_linked_genes.tsv (human: XIST + 9 Y-linked; mouse: 5-gene panel of the mouse DE script).
##  M1_sexgene_deg_ref = panel genes in the M0 (unadjusted) DE set; M1_sexgene_deg_alt = in the M1 (sex-adjusted) DE set; M0perm_sexgene_deg_alt = in the permuted-sex model's DE set.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
D <- load_summary(); V <- D[valid == TRUE & single_sex == FALSE & design_sex %in% CLS & !is.na(M1_sexgene_deg_ref)]; V[, cls := factor(design_sex, levels=CLS)]
lk <- function(d, set) d[, { a <- wilson(sum(M1_sexgene_deg_ref >= 1), .N); b <- wilson(sum(M1_sexgene_deg_alt >= 1), .N); ok <- !is.na(M0perm_sexgene_deg_alt); p <- wilson(sum(M0perm_sexgene_deg_alt[ok] >= 1), sum(ok)); a3 <- wilson(sum(M1_sexgene_deg_ref >= 3), .N)
  .(set=set, n_studies=.N, k_M0=a$k, pct_M0_ge1=a$pct, M0_lo=a$lo, M0_hi=a$hi, k_M1=b$k, pct_M1_ge1=b$pct, M1_lo=b$lo, M1_hi=b$hi, k_perm=p$k, pct_perm_ge1=p$pct, perm_lo=p$lo, perm_hi=p$hi,
    pct_M0_ge3=a3$pct, median_panel_genes_when_present=as.double(median(M1_sexgene_deg_ref[M1_sexgene_deg_ref >= 1])), removed_by_adjustment=sum(M1_sexgene_deg_ref >= 1 & M1_sexgene_deg_alt == 0)) }, keyby=.(organism, cls)]
L <- rbind(lk(V, "all valid mixed-sex, M1 estimable"), lk(V[analysis == TRUE], "analysis set"), lk(V[M1_n_deg_ref >= 1], "M0 DE >= 1")); out(L, "E9_xy_leakage_by_class.tsv")
## dose-response over |Δf| bins (analysis set) and trend test
B <- V[analysis == TRUE, { a <- wilson(sum(M1_sexgene_deg_ref >= 1), .N); .(n_studies=.N, k=a$k, pct_M0_ge1=a$pct, lo=a$lo, hi=a$hi) }, keyby=.(organism, bin=cut(sex_diff, c(-0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1)))]; out(B, "E9_xy_leakage_by_bin.tsv")
TR <- rbindlist(lapply(ORGS, function(o) { x <- V[organism == o & analysis == TRUE]; f <- summary(glm(I(M1_sexgene_deg_ref >= 1) ~ sex_diff + log(n) + log1p(M1_n_deg_ref), binomial, x))$coefficients
  data.table(organism=o, n=nrow(x), OR_per_0.1_sexdiff=exp(0.1*f["sex_diff", 1]), lo=exp(0.1*(f["sex_diff", 1] - 1.96*f["sex_diff", 2])), hi=exp(0.1*(f["sex_diff", 1] + 1.96*f["sex_diff", 2])), p=f["sex_diff", 4]) })); out(TR, "E9_xy_leakage_trend.tsv")
options(width=250); print(L[, .(organism, cls, set, n_studies, k_M0, pct_M0_ge1, M0_lo, M0_hi, k_M1, pct_M1_ge1, k_perm, pct_perm_ge1, removed_by_adjustment)], digits=3); print(B, digits=3); print(TR, digits=3)
