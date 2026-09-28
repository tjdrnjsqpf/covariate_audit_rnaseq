#!/usr/bin/env Rscript
## S eligible studies: analysis set ∩ balanced ∩ all four group×sex cells ≥6. Usage: Rscript 42_make_synthetic_list.R  → config/synthetic_list_{org}.txt
suppressMessages(library(data.table)); PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
for (org in c("human","mouse")) { d <- fread(sprintf("results/main_%s_summary_v2.tsv", org)); A <- readLines(sprintf("config/c2_analysis_set_%s.txt", org))
  a <- d[study %in% A & design_sex == "balanced"]; a[, mincell := pmin(n_female_case, n_female_control, n_case - n_female_case, n_control - n_female_control)]
  writeLines(a[mincell >= 6, study], sprintf("config/synthetic_list_%s.txt", org)); cat(org, sum(a$mincell >= 6), "eligible\n") }
