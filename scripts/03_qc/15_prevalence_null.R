#!/usr/bin/env Rscript
## C1-b: does the observed sex/age imbalance exceed what random assignment (hypergeometric / permutation null) would produce?
suppressMessages(library(data.table)); ORG <- commandArgs(trailingOnly=TRUE)[1]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ); dir.create("results/qc", showWarnings=FALSE)
d <- fread(sprintf("results/main_%s_summary_v2.tsv", ORG))[valid == TRUE & single_sex == FALSE & !is.na(n_female_case)]
## Sex: k of K females fall in the case group (nc) → k ~ Hypergeom(N, K, nc). |Δfemale fraction| = |k/nc − (K−k)/(N−nc)|
null_sex <- function(N, K, nc) { k <- 0:min(K, nc); p <- dhyper(k, K, N - K, nc); dd <- abs(k / nc - (K - k) / (N - nc)); list(k = k, p = p, dd = dd) }
r <- d[, { z <- null_sex(n, n_female, n_case); obs <- abs(n_female_case / n_case - (n_female - n_female_case) / n_control)
  .(obs_diff = round(obs, 3), p_geq_obs = round(sum(z$p[z$dd >= obs - 1e-9]), 4), null_p_ge_0.3 = round(sum(z$p[z$dd >= 0.3]), 4), null_p_ge_0.15 = round(sum(z$p[z$dd >= 0.15]), 4)) }, by=study]
r <- merge(d[, .(study, n, n_case, n_control, n_female, design_sex)], r, by="study")
fwrite(r, sprintf("results/qc/prevalence_null_sex_%s.tsv", ORG), sep="\t")
wil <- function(x, n) { p <- x / n; z <- 1.96; den <- 1 + z^2/n; c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/den, (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/den) }
cat(sprintf("########## Prevalence null comparison (%s): %d mixed-sex valid studies ##########\n", ORG, nrow(r)))
for (thr in c(0.15, 0.3)) { obs <- sum(r$obs_diff >= thr); exp <- sum(r[[if (thr == 0.3) "null_p_ge_0.3" else "null_p_ge_0.15"]]); ci <- wil(obs, nrow(r))
  cat(sprintf("|Δfemale fraction| ≥ %.2f: observed %d (%.1f%% [%.1f-%.1f]) vs expected under random assignment %.1f (%.1f%%) → observed/expected %.2fx, excess %.1f%%p\n", thr, obs, 100*obs/nrow(r), 100*ci[1], 100*ci[2], exp, 100*exp/nrow(r), obs/exp, 100*(obs-exp)/nrow(r))) }
cat(sprintf("Per-study p<0.05 (imbalance beyond chance): %d (%.1f%%);  Bonferroni p<0.05/n: %d\n", sum(r$p_geq_obs < 0.05), 100*mean(r$p_geq_obs < 0.05), sum(r$p_geq_obs < 0.05/nrow(r))))
cat("\n== By class: observed count vs null expectation (threshold ≥0.3) ==\n"); print(r[, .(n = .N, exp_ge_0.3 = round(sum(null_p_ge_0.3), 1), p_lt_0.05 = sum(p_geq_obs < 0.05)), by=design_sex][order(design_sex)])
## Age: SMD null by permuting group labels
da <- fread(sprintf("results/main_%s_summary_v2.tsv", ORG))[valid == TRUE & !is.na(age_smd)]
ra <- rbindlist(lapply(da$study, function(s) { f <- sprintf("results/main_%s_v2/%s/pheno.tsv", ORG, s); if (!file.exists(f)) return(NULL)
  ph <- fread(f)[!is.na(age)]; if (nrow(ph) < 6 || uniqueN(ph$group) < 2) return(NULL)
  smd <- function(g) { a <- ph$age[g == "case"]; b <- ph$age[g == "control"]; s <- sqrt((var(a) + var(b))/2); if (is.na(s) || s == 0) NA else abs(mean(a) - mean(b))/s }
  obs <- smd(ph$group); set.seed(1); nl <- replicate(2000, smd(sample(ph$group))); nl <- nl[!is.na(nl)]
  data.table(study = s, n_age = nrow(ph), obs_smd = round(obs, 3), p_geq_obs = round(mean(nl >= obs), 4), null_p_ge_0.8 = round(mean(nl >= 0.8), 4), null_p_ge_1.6 = round(mean(nl >= 1.6), 4)) }))
if (nrow(ra)) { fwrite(ra, sprintf("results/qc/prevalence_null_age_%s.tsv", ORG), sep="\t")
  cat(sprintf("\n== Age (n=%d) ==\n", nrow(ra)))
  for (thr in c(0.8, 1.6)) { obs <- sum(ra$obs_smd >= thr, na.rm=TRUE); exp <- sum(ra[[if (thr == 1.6) "null_p_ge_1.6" else "null_p_ge_0.8"]])
    cat(sprintf("SMD ≥ %.1f: observed %d (%.1f%%) vs permutation expectation %.1f (%.1f%%) → %.2fx\n", thr, obs, 100*obs/nrow(ra), exp, 100*exp/nrow(ra), obs/max(exp, 1e-9))) }
  cat(sprintf("Per-study p<0.05: %d (%.1f%%)\n", sum(ra$p_geq_obs < 0.05, na.rm=TRUE), 100*mean(ra$p_geq_obs < 0.05, na.rm=TRUE))) }
