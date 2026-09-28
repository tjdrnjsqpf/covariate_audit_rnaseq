#!/usr/bin/env Rscript
## Prevalence tables (Aim 1) rewritten: denominators stated explicitly, unresolved cases via Manski bounds.
suppressMessages(library(data.table)); PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ); dir.create("results/tables", showWarnings=FALSE)
out <- list()
for (org in c("human","mouse")) {
  ps <- fread(if (org == "human") "results/project_summary_all.tsv" else "results/mouse_project_summary_sexinferred.tsv")
  v2 <- fread(sprintf("results/main_%s_summary_v2.tsv", org))
  row <- function(k, v) out[[length(out) + 1]] <<- data.table(organism = org, item = k, value = as.character(v))
  row("recount3 SRA projects (parsed)", nrow(ps))
  if ("n_sex" %in% names(ps)) { row("  projects with sex in metadata (>=1 sample)", sum(ps$n_sex > 0, na.rm=TRUE)); row("  projects with sex in metadata (%)", round(100*mean(ps$n_sex > 0, na.rm=TRUE), 1)) }
  if ("n_age" %in% names(ps)) row("  projects with age in metadata (%)", round(100*mean(ps$n_age > 0, na.rm=TRUE), 1))
  V <- v2[valid == TRUE]; row("LLM-judged usable design → reanalysis valid", nrow(V))
  row("  single-sex (not adjustable) n", sum(V$single_sex)); row("  single-sex (%)", round(100*mean(V$single_sex), 1))
  M <- V[single_sex == FALSE]; row("  mixed-sex studies n", nrow(M))
  for (g in c("balanced","mildly_imbalanced","partially_confounded")) row(sprintf("    among mixed-sex, %s (%%)", g), round(100*mean(M$design_sex == g, na.rm=TRUE), 1))
  row("    among mixed-sex, |Δfemale fraction| ≥ 0.3 (%)", round(100*mean(M$sex_diff >= 0.3, na.rm=TRUE), 1))
  nall <- nrow(v2); k <- sum(M$sex_diff >= 0.3, na.rm=TRUE); unres <- nall - nrow(M)
  row("    |Δfemale fraction| ≥ 0.3 bounds (all studies with results as denominator, lower–upper %)", sprintf("%.1f – %.1f", 100*k/nall, 100*(k + unres)/nall))
  row("  samples with missing metadata sex → recovered from expression (number of samples)", sum(V$n_sex_meta_missing)); row("  samples with metadata–expression sex mismatch", sum(V$n_sex_mismatch))
  Ag <- V[!is.na(age_smd)]; row("  age available, SMD computable n (%)", sprintf("%d (%.1f)", nrow(Ag), 100*nrow(Ag)/nrow(V)))
  for (g in c("balanced","mildly_imbalanced","partially_confounded","completely_confounded")) row(sprintf("    age %s (%%)", g), round(100*mean(Ag$design_age == g, na.rm=TRUE), 1))
}
tab <- rbindlist(out); fwrite(tab, "results/tables/prevalence_v2.tsv", sep="\t"); print(tab, nrows=200)
