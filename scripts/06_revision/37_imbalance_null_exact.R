#!/usr/bin/env Rscript
## E8. Exact test of the imbalance prevalence against the hypergeometric (random-allocation) null.
##  Under H0 the count of studies with |Δf| >= thr is Poisson-binomial with per-study probabilities p_i. We give the exact upper-tail P(X >= observed) by DP convolution,
##  a z-score, and a CI for observed/expected (Wilson / Clopper–Pearson interval of the observed proportion divided by the expected proportion; E is fixed given the margins).
##  (a) from the null_p_ge_* columns written by scripts/03_qc/15_prevalence_null.R (server results/qc/, md5-identical local copy);
##  (b) recomputed from the 2×2 counts in the summary table with one tolerance (1e-9) used for BOTH observed and null — 15_prevalence_null.R counts observed from a 3-d.p.-rounded
##      |Δf| but evaluates the null with an untolerated `dd >= 0.3`, so boundary atoms (e.g. exactly 0.3) can be treated asymmetrically.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
NULLF <- c(human=Sys.getenv("NULL_HUMAN", "results/qc/prevalence_null_sex_human.tsv"), mouse=Sys.getenv("NULL_MOUSE", "results/qc/prevalence_null_sex_mouse.tsv"))
pb_tail <- function(p, k) { d <- 1; for (pi in p) d <- c(d*(1 - pi), 0) + c(0, d*pi); list(upper=sum(d[(k + 1):length(d)]), lower=sum(d[1:(k + 1)])) }   # P(X>=k), P(X<=k)
res <- function(o, src, thr, obs, p) { N <- length(p); E <- sum(p); V <- sum(p*(1 - p)); t <- pb_tail(p, obs); w <- wilson(obs, N); cp <- binom.test(obs, N)$conf.int
  data.table(organism=o, source=src, threshold=thr, n_studies=N, observed=obs, observed_pct=100*obs/N, expected=E, expected_pct=100*E/N, ratio=obs/E, ratio_wilson_lo=(w$lo/100)/(E/N), ratio_wilson_hi=(w$hi/100)/(E/N),
             ratio_exact_lo=cp[1]/(E/N), ratio_exact_hi=cp[2]/(E/N), excess_pp=100*(obs - E)/N, z=(obs - E)/sqrt(V), p_upper_poisson_binomial=t$upper, p_lower_poisson_binomial=t$lower, p_two_sided=min(1, 2*min(t$upper, t$lower))) }
D <- load_summary(); R <- list(); PS <- list()
for (o in ORGS) { f <- fread(NULLF[[o]])
  for (thr in c(0.3, 0.15)) R[[length(R) + 1]] <- res(o, "15_prevalence_null.R columns", thr, sum(f$obs_diff >= thr), f[[sprintf("null_p_ge_%s", thr)]])
  d <- D[organism == o & valid == TRUE & single_sex == FALSE & !is.na(n_female_case)]
  x <- d[, { k <- 0:min(n_female, n_case); p <- dhyper(k, n_female, n - n_female, n_case); dd <- abs(k/n_case - (n_female - k)/(n - n_case)); ob <- abs(n_female_case/n_case - (n_female - n_female_case)/n_control)
             .(n=n, obs=ob, p03=sum(p[dd >= 0.3 - 1e-9]), p015=sum(p[dd >= 0.15 - 1e-9]), p_geq_obs=sum(p[dd >= ob - 1e-9])) }, by=study]
  x[, organism := o]; PS[[o]] <- x
  R[[length(R) + 1]] <- res(o, "recomputed, tolerance 1e-9 both sides", 0.3, sum(x$obs >= 0.3 - 1e-9), x$p03); R[[length(R) + 1]] <- res(o, "recomputed, tolerance 1e-9 both sides", 0.15, sum(x$obs >= 0.15 - 1e-9), x$p015)
  for (s in levels(nstrata(1))) { y <- x[nstrata(n) == s]; if (nrow(y) >= 10) R[[length(R) + 1]] <- res(o, paste0("recomputed, n ", s), 0.3, sum(y$obs >= 0.3 - 1e-9), y$p03) }
  ## excluding completely confounded designs (|Δf| >= 0.99), which are arguably deliberate sex comparisons mislabelled as group contrasts
  y <- x[obs < 0.99]; R[[length(R) + 1]] <- res(o, "recomputed, excluding |Δf| >= 0.99 (null prob. unchanged)", 0.3, sum(y$obs >= 0.3 - 1e-9), y$p03) }
R <- rbindlist(R); out(R, "E8_imbalance_null_exact.tsv"); out(rbindlist(PS), "E8_per_study_null_prob.tsv")
options(width=250); print(R[, .(organism, source, threshold, n_studies, observed, expected=round(expected, 1), ratio=round(ratio, 2), lo=round(ratio_wilson_lo, 2), hi=round(ratio_wilson_hi, 2), z=round(z, 2), p_upper=signif(p_upper_poisson_binomial, 3), p_lower=signif(p_lower_poisson_binomial, 3))])
