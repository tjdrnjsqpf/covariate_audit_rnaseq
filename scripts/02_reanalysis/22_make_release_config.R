#!/usr/bin/env Rscript
## Unified public config (review D 9.2): append tissue:<category> to every row of main_{org}_runs.tsv so the public config matches the config actually run.
## Usage: Rscript 22_make_release_config.R → config/release_{org}_runs.tsv
suppressMessages(library(data.table)); PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
for (org in c("human","mouse")) {
  cfg <- fread(sprintf("config/main_%s_runs.tsv", org), colClasses="character")
  tc <- if (org == "human") fread("results/project_summary_all.tsv")[, .(study, tissue_category)] else fread("results/mouse_project_summary_sexinferred.tsv")[, .(study, tissue_category)]
  tc <- unique(tc, by="study"); cfg <- merge(cfg, tc, by="study", all.x=TRUE, sort=FALSE)
  cfg[is.na(tissue_category) | tissue_category == "", tissue_category := "unknown"]
  cfg[, exclude := gsub("(^|;)tissue:[^;]*", "", exclude)][, exclude := sub("^;", "", exclude)]   # strip any existing tissue: option, then reassign
  cfg[, exclude := ifelse(nzchar(exclude), paste0(exclude, ";tissue:", tissue_category), paste0("tissue:", tissue_category))]; cfg[, tissue_category := NULL]
  for (cc in names(cfg)) cfg[[cc]] <- gsub("[\t\r\n]+", " ", cfg[[cc]])   # 2026-09-20: tabs/newlines inside values (SRP159556) broke the TSV columns, so the collector stopped reading at row 2,192
  fwrite(cfg, sprintf("config/release_%s_runs.tsv", org), sep="\t", quote=FALSE); cat(org, nrow(cfg), "rows\n") }
