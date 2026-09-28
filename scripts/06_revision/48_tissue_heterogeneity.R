#!/usr/bin/env Rscript
## Tissue / comparison-type heterogeneity (review I v3 request 3): Δ_across, Δ_within and sex-chromosome leakage by tissue category and comparison type; tissue composition of the 85 synthetic-experiment studies.
## Usage: LC_ALL=en_US.UTF-8 Rscript scripts/06_revision/48_tissue_heterogeneity.R [PROJ]
suppressMessages(library(data.table)); PROJ <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(PROJ)) PROJ <- getwd(); setwd(PROJ)
lv <- c("balanced","mildly_imbalanced","partially_confounded"); bci <- function(x) { x <- x[!is.na(x)]; if (length(x) < 3) return(sprintf("%.3f (n<3)", median(x))); set.seed(1); b <- replicate(2000, median(sample(x, replace=TRUE))); sprintf("%.3f [%.3f, %.3f]", median(x), quantile(b,.025), quantile(b,.975)) }
tcl <- function(t) { t <- tolower(ifelse(is.na(t) | t == "", "unknown", t)); fifelse(grepl("blood|immune|pbmc", t), "blood/immune", fifelse(grepl("brain|cns|neur", t), "brain/CNS", fifelse(grepl("tumo", t), "tumour", fifelse(grepl("liver", t), "liver", "other"))) ) }
out <- list(); ctype <- list()
for (o in c("human","mouse")) { d <- fread(sprintf("results/main_%s_summary_v2.tsv", o)); A <- readLines(sprintf("config/c2_analysis_set_%s.txt", o))
  a <- d[study %in% A & design_sex %in% lv & !is.na(M1_jaccard) & !is.na(M0perm_jaccard) & !is.na(M0permw_jaccard)]
  a[, tissue := tcl(tissue_category)]; a[, ctype := fifelse(grepl("case_control", design_type), "case-control", fifelse(grepl("genotype", design_type), "genotype", fifelse(grepl("treatment|dose|time", design_type), "treatment/time", "other")))]
  a[, `:=`(d_across=M1_jaccard - M0perm_jaccard, d_within=M1_jaccard - M0permw_jaccard, leak0=M1_sexgene_deg_ref > 0)]
  for (grp in c("tissue","ctype")) { g <- a[, .(organism=o, stratum=grp, n=.N, n_confounded=sum(design_sex=="partially_confounded"),
      d_across_all=bci(d_across), d_within_all=bci(d_within),
      d_across_conf=bci(d_across[design_sex=="partially_confounded"]), d_within_conf=bci(d_within[design_sex=="partially_confounded"]),
      leak_M0_conf=round(mean(leak0[design_sex=="partially_confounded"]),3)), by=.(level=get(grp))][order(-n)]; out[[length(out)+1]] <- g }
  ## Tissue composition of the synthetic-experiment studies
  sl <- readLines(sprintf("config/synthetic_list_%s.txt", o)); ctype[[o]] <- a[study %in% sl, .(organism=o, n=.N), by=.(tissue)][order(-n)] }
res <- rbindlist(out); dir.create("results/revision", showWarnings=FALSE); fwrite(res, "results/revision/T_tissue_heterogeneity.tsv", sep="\t"); fwrite(rbindlist(ctype), "results/revision/T_synthetic_tissue_composition.tsv", sep="\t")
options(width=240); print(res[n >= 8]); print(rbindlist(ctype))
