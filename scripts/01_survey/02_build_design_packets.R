#!/usr/bin/env Rscript
# Build input packets for LLM design classification: per-project title / abstract / tissue / n / attribute keys+levels / GEOMeta flags -> JSONL
# Usage: Rscript 02_build_design_packets.R <study_list.txt> <out.jsonl>
suppressMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE); LIST <- args[1]; OUTF <- args[2]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; PROC <- file.path(PROJ, "data/processed"); RES <- file.path(PROJ, "results")
studies <- readLines(LIST)
hum <- fread(file.path(RES, "project_summary_all.tsv"))[, .(study, organism="human", n, study_title, abstract, tissue_example, tissue_category, strategy, source, any_sc, any_small, any_cell_line, age_min, age_max, n_sex, n_age)]
mou <- fread(file.path(RES, "mouse_project_summary_sexinferred.tsv"))[, .(study, organism="mouse", n, study_title, abstract=NA_character_, tissue_example, tissue_category, strategy, source, any_sc, any_small, any_cell_line=FALSE, age_min=NA_real_, age_max=NA_real_, n_sex, n_age)]
## mouse abstracts come from the sample table
ms <- fread(file.path(PROC, "mouse_sample_table.tsv.gz"), select=c("study","abstract"))[, .(abstract = abstract[1]), by=study]
mou <- merge(mou[, -"abstract"], ms, by="study", all.x=TRUE)
ps <- rbind(hum, mou, fill=TRUE)[study %in% studies]
al_h <- fread(file.path(PROC, "sample_attributes_long.tsv.gz"), colClasses="character", na.strings=c("","NA"))[study %in% studies]
al_m <- fread(file.path(PROC, "mouse_sample_attributes_long.tsv.gz"), colClasses="character", na.strings=c("","NA"))[study %in% studies]
al <- rbind(al_h, al_m)[!is.na(key) & !is.na(value)][, value := tolower(trimws(value))]
gm <- tryCatch(fread(file.path(RES, "geometa_project_coverage.tsv"))[, .(study, geometa_n=n_gm, geometa_invitro_frac=round(invitro, 2), geometa_pert_frac=round(pert, 2))], error=function(e) NULL)
con <- file(OUTF, "w")
for (st in studies) {
  p <- ps[study == st]; if (!nrow(p)) next
  a <- al[study == st]
  keys <- a[, .(n_samples = uniqueN(external_id), n_levels = uniqueN(value)), by=key][order(-n_samples, n_levels)]
  keys <- keys[!key %in% c("sample name","sample_name","title","alias","sra accession","biosample","insdc status","insdc first public","insdc last update","insdc center name","insdc center alias","ena checklist","broker name","ena-checklist","description")]
  kv <- lapply(seq_len(min(25, nrow(keys))), function(i) { k <- keys$key[i]; tb <- sort(table(a[key == k, value]), decreasing=TRUE)
    list(key=k, n_samples=keys$n_samples[i], n_levels=keys$n_levels[i], levels=if (length(tb) <= 12) as.list(setNames(as.integer(tb), names(tb))) else c(as.list(setNames(as.integer(tb[1:10]), names(tb[1:10]))), list(`...`=paste0(length(tb)-10, " more"))))})
  pk <- list(study=st, organism=p$organism, n_samples=p$n, study_title=p$study_title, abstract=substr(gsub("\\s+", " ", ifelse(is.na(p$abstract), "", p$abstract)), 1, 2500),
             tissue_example=p$tissue_example, tissue_category=p$tissue_category, library_strategy=p$strategy, library_source=p$source,
             flags=list(single_cell_keywords=isTRUE(p$any_sc), small_rna_keywords=isTRUE(p$any_small), cell_line_attribute=isTRUE(p$any_cell_line)),
             metadata=list(n_with_sex=p$n_sex, n_with_age=p$n_age, age_min=p$age_min, age_max=p$age_max),
             sample_attributes=kv)
  if (!is.null(gm) && st %in% gm$study) { g <- gm[study == st]; pk$geometa <- list(n_annotated=g$geometa_n, in_vitro_fraction=g$geometa_invitro_frac, perturbation_fraction=g$geometa_pert_frac) }
  writeLines(toJSON(pk, auto_unbox=TRUE, null="null", na="null"), con)
}
close(con); cat("packets written:", length(readLines(OUTF)), "of", length(studies), "\n")
