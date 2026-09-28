#!/usr/bin/env Rscript
## E5. SVA: n_sv distribution by class; share where the cap floor((n−4)/2) binds; capture (max sv_sex_r2 > thr) within n strata and at R² 0.3/0.5/0.7.
## DE script: SVA attempted only if n >= 10; nsv <- min(num.sv(be), 5, floor((n−4)/2)); sv_sex_r2 = max over SVs of R²(sv ~ sex); NA when n_sv = 0.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
D <- load_summary(); A <- D[analysis == TRUE & design_sex %in% CLS]; A[, cls := factor(design_sex, levels=CLS)]; A[, nstr := nstrata(n)]
A[, cap_n := floor((n - 4)/2)]; A[, sva_attempted := n >= 10]
A[, cap_state := fifelse(!sva_attempted, "not_attempted(n<10)", fifelse(n_sv == 0, "n_sv=0", fifelse(n_sv == cap_n & cap_n < 5, "n_cap_binds", fifelse(n_sv == 5, "cap5_binds", "num.sv_decides"))))]
NS <- A[, .(n_studies=.N, n_median=as.double(median(n)), sva_not_attempted=sum(!sva_attempted), nsv0=sum(sva_attempted & n_sv == 0), nsv1=sum(n_sv == 1), nsv2=sum(n_sv == 2), nsv3=sum(n_sv == 3), nsv4=sum(n_sv == 4), nsv5=sum(n_sv == 5),
            nsv_median_given_ge1=as.double(median(n_sv[n_sv >= 1])), nsv_mean_given_ge1=mean(n_sv[n_sv >= 1]), n_with_sv=sum(n_sv >= 1),
            pct_ncap_binds_given_sv=100*mean(cap_state[n_sv >= 1] == "n_cap_binds"), pct_cap5_binds_given_sv=100*mean(cap_state[n_sv >= 1] == "cap5_binds")), keyby=.(organism, cls)]; out(NS, "E5_nsv_by_class.tsv")
S <- A[!is.na(sv_sex_r2)]   # same set as Fig 2B
cap <- function(d, by) d[, { w3 <- wilson(sum(sv_sex_r2 > 0.3), .N); w5 <- wilson(sum(sv_sex_r2 > 0.5), .N); w7 <- wilson(sum(sv_sex_r2 > 0.7), .N)
  .(n_studies=.N, n_median=as.double(median(n)), nsv_median=as.double(median(n_sv)), r2_median=median(sv_sex_r2), r2_q1=quantile(sv_sex_r2, .25), r2_q3=quantile(sv_sex_r2, .75),
    capture_r2_0.3=w3$pct, lo_0.3=w3$lo, hi_0.3=w3$hi, capture_r2_0.5=w5$pct, lo_0.5=w5$lo, hi_0.5=w5$hi, capture_r2_0.7=w7$pct, lo_0.7=w7$lo, hi_0.7=w7$hi) }, keyby=by]
C1 <- cap(S, c("organism", "cls")); out(C1, "E5_capture_by_class_thresholds.tsv")
C2 <- cap(S, c("organism", "nstr", "cls")); out(C2, "E5_capture_by_class_nstrata.tsv")
C3 <- cap(S[, nsv_grp := fifelse(n_sv == 1, "1", fifelse(n_sv <= 3, "2-3", "4-5"))], c("organism", "nsv_grp", "cls")); out(C3, "E5_capture_by_class_nsv.tsv")
## does class still predict capture after n and n_sv? logistic, confounded vs balanced
LG <- rbindlist(lapply(ORGS, function(o) { x <- S[organism == o]; f0 <- glm(I(sv_sex_r2 > 0.5) ~ cls, binomial, x); f1 <- glm(I(sv_sex_r2 > 0.5) ~ cls + log(n) + n_sv, binomial, x)
  g <- function(f, lab) { cf <- summary(f)$coefficients; data.table(organism=o, model=lab, term=rownames(cf), OR=exp(cf[, 1]), lo=exp(cf[, 1] - 1.96*cf[, 2]), hi=exp(cf[, 1] + 1.96*cf[, 2]), p=cf[, 4]) }
  rbind(g(f0, "class_only"), g(f1, "class+log_n+n_sv")) }))[term != "(Intercept)"]; out(LG, "E5_capture_logistic.tsv")
options(width=250); print(NS, digits=3); print(C1, digits=3); print(C2[, .(organism, nstr, cls, n_studies, nsv_median, capture_r2_0.5, lo_0.5, hi_0.5)], digits=3); print(C3[, .(organism, nsv_grp, cls, n_studies, capture_r2_0.5, lo_0.5, hi_0.5)], digits=3); print(LG, digits=3)
