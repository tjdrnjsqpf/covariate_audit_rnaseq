#!/usr/bin/env Rscript
## E2. Sample size: median n by |Δf| class; Δ = M1_jaccard − M0perm_jaccard by class within n strata; n >= 20 subset; Δ vs log VIF.
## VIF = 1/(1 − r²), r = phi correlation of group and sex indicators from the 2×2 counts.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
D <- load_summary()
V <- D[valid == TRUE & single_sex == FALSE & design_sex %in% CLS]
A <- D[analysis == TRUE & design_sex %in% CLS & !is.na(M1_jaccard) & !is.na(M0perm_jaccard)]
A[, delta := M1_jaccard - M0perm_jaccard]; A[, nstr := nstrata(n)]
V[, cls := factor(design_sex, levels=CLS)]; A[, cls := factor(design_sex, levels=CLS)]
## (a) median n by class
qn <- function(d, set) d[, .(set=set, n_studies=.N, n_median=as.double(median(n)), n_q1=quantile(n, .25), n_q3=quantile(n, .75), n_min=min(n), n_max=max(n), pct_n_lt12=100*mean(n < 12), pct_n_ge20=100*mean(n >= 20), pct_n_ge40=100*mean(n >= 40)), keyby=.(organism, cls)]
NB <- rbind(qn(V, "all_valid_mixed_sex"), qn(A, "analysis_set")); out(NB, "E2_n_by_class.tsv")
## (b) Δ by class within n strata (+ all, + n>=20)
ds <- function(d, lab) d[, { z <- bci(delta); .(n_studies=.N, n_median=as.double(median(n)), delta_median=z$med, lo=z$lo, hi=z$hi, M1_jaccard_median=median(M1_jaccard), M0perm_jaccard_median=median(M0perm_jaccard),
                              pct_delta_below_m0.2=100*mean(delta < -0.2)) }, keyby=.(organism, cls)][, nstr := lab]
DS <- rbind(ds(A, "all"), rbindlist(lapply(levels(A$nstr), function(s) ds(A[nstr == s], s))), ds(A[n >= 20], "n>=20"))
DS[, nstr := factor(nstr, levels=c("all", "<12", "12-19", "20-39", ">=40", "n>=20"))]; setkey(DS, organism, nstr, cls); setcolorder(DS, c("organism", "nstr", "cls")); out(DS, "E2_delta_by_class_nstrata.tsv")
## stratified tests: confounded vs balanced within n strata (van Elteren via stratified Wilcoxon = sum of stratum rank statistics), and Δ ~ class + n-stratum
vt <- rbindlist(lapply(ORGS, function(o) { a <- A[organism == o & cls != "mildly_imbalanced"]
  s <- a[, { nb <- sum(cls == "balanced"); nc <- sum(cls == "partially_confounded"); if (nb < 1 || nc < 1) NULL else { r <- rank(delta); W <- sum(r[cls == "partially_confounded"]); N <- .N
             .(w=1/(N + 1), W=W, EW=nc*(N + 1)/2, VW=nb*nc*(N + 1)/12) } }, by=nstr]
  z <- sum(s$w*(s$W - s$EW))/sqrt(sum(s$w^2*s$VW)); f <- lm(delta ~ sex_diff + nstr, A[organism == o]); f20 <- lm(delta ~ sex_diff + log(n), A[organism == o & n >= 20]); cf <- summary(f)$coefficients; c20 <- summary(f20)$coefficients
  data.table(organism=o, van_elteren_z=z, van_elteren_p=2*pnorm(-abs(z)), slope_sexdiff_adj_nstrata=cf["sex_diff", 1], se=cf["sex_diff", 2], p=cf["sex_diff", 4],
             n_ge20=nobs(f20), slope_sexdiff_n_ge20_adj_logn=c20["sex_diff", 1], se_n_ge20=c20["sex_diff", 2], p_n_ge20=c20["sex_diff", 4]) })); out(vt, "E2_stratified_tests.tsv")
## (c) Δ vs log VIF
A[, `:=`(a=n_female_case, b=n_case - n_female_case, c=n_female_control, d=n_control - n_female_control)]
A[, r := (a*d - b*c)/sqrt((a + b)*(c + d)*(a + c)*(b + d))]; A[, vif := 1/(1 - r^2)]; A[, logvif := log(vif)]
VF <- rbindlist(lapply(ORGS, function(o) { x <- A[organism == o & is.finite(logvif)]; sp <- suppressWarnings(cor.test(x$delta, x$logvif, method="spearman")); f1 <- summary(lm(delta ~ logvif, x))$coefficients; f2 <- summary(lm(delta ~ logvif + log(n), x))$coefficients
  f3 <- summary(lm(delta ~ logvif + sex_diff + log(n), x))$coefficients
  data.table(organism=o, n_studies=nrow(x), vif_median=median(x$vif), vif_q3=quantile(x$vif, .75), vif_p95=quantile(x$vif, .95), vif_max=max(x$vif), spearman_delta_logvif=sp$estimate, spearman_p=sp$p.value,
             spearman_delta_sexdiff=cor(x$delta, x$sex_diff, method="spearman"), spearman_sexdiff_logvif=cor(x$sex_diff, x$logvif, method="spearman"),
             slope_logvif=f1["logvif", 1], se=f1["logvif", 2], slope_logvif_adj_logn=f2["logvif", 1], se_adj=f2["logvif", 2], p_adj=f2["logvif", 4],
             slope_logvif_adj_sexdiff_logn=f3["logvif", 1], p_logvif_given_sexdiff=f3["logvif", 4], slope_sexdiff_given_logvif=f3["sex_diff", 1], p_sexdiff_given_logvif=f3["sex_diff", 4]) })); out(VF, "E2_delta_vs_logvif.tsv")
VB <- A[is.finite(vif), { z <- bci(delta); .(n_studies=.N, vif_median=median(vif), sqrt_vif_median=median(sqrt(vif)), delta_median=z$med, lo=z$lo, hi=z$hi) }, keyby=.(organism, vif_bin=cut(vif, c(1, 1.02, 1.1, 1.25, 1.5, 2, Inf), right=FALSE))]; out(VB, "E2_delta_by_vif_bin.tsv")
VC <- A[is.finite(vif), .(vif_median=median(vif), vif_max=max(vif), sqrt_vif_median=median(sqrt(vif))), keyby=.(organism, cls)]; out(VC, "E2_vif_by_class.tsv")
out(A[, .(organism, study, n, design_sex, sex_diff, r, vif, delta)], "E2_per_study_vif.tsv")
options(width=220); print(NB, digits=3); print(DS, digits=3); print(vt, digits=3); print(VF, digits=3); print(VB, digits=3); print(VC, digits=3)
