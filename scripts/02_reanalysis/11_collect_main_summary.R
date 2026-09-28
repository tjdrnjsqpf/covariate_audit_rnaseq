#!/usr/bin/env Rscript
# Collect main reanalysis results + aggregate by design grade / tissue. Usage: Rscript 11_collect_main_summary.R <outdir> <config> <organism>
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE); OUTDIR <- args[1]; CFG <- args[2]; ORG <- args[3]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
fs <- list.files(file.path(PROJ, OUTDIR), pattern="summary.tsv", recursive=TRUE, full.names=TRUE)
d <- rbindlist(lapply(fs, fread), fill=TRUE)
cfg <- fread(file.path(PROJ, CFG))[, .(study, note)]
d <- merge(cfg, d, by="study", all.y=TRUE)
ud <- fread(file.path(PROJ, "results/design_llm/usable_designs_v1.tsv"))[organism == ORG, .(study, design_type, experimental_setting, sex_composition, scope_flags, confidence)]
d <- merge(d, ud, by="study", all.x=TRUE)
if (ORG == "human") { ps <- fread(file.path(PROJ, "results/project_summary_all.tsv"))[, .(study, tissue_category, first_published)]
} else ps <- fread(file.path(PROJ, "results/mouse_project_summary_sexinferred.tsv"))[, .(study, tissue_category, first_published = NA_character_)]
d <- merge(d, ps, by="study", all.x=TRUE)
## Validity: number of detected genes
d[, valid := is.na(n_genes_detected) | n_genes_detected >= 5000]
## Design grade
d[, design := fifelse(is.na(sex_diff), NA_character_, fifelse(sex_diff >= 0.99, "completely_confounded",
              fifelse(sex_diff >= 0.3 | (!is.na(age_smd) & age_smd >= 0.8), "partially_confounded",
              fifelse(sex_diff >= 0.15 | (!is.na(age_smd) & age_smd >= 0.4), "mildly_imbalanced", "balanced"))))]
fwrite(d, file.path(PROJ, sprintf("results/main_%s_summary.tsv", ORG)), sep="\t")
cat("projects with results:", nrow(d), " valid:", sum(d$valid), "\n")
cat("models fitted:\n"); print(table(d[valid == TRUE, models]))
v <- d[valid == TRUE]
cat("\n== Design grade distribution ==\n"); print(table(v$design, useNA="ifany"))
cat("\n== M0 vs M1(+sex) change by grade (median) ==\n")
print(v[!is.na(design) & !is.na(M1_jaccard), .(n = .N, deg_M0 = as.numeric(median(M1_n_deg_ref)), deg_M1 = as.numeric(median(M1_n_deg_alt)),
       jaccard = round(median(M1_jaccard[M1_n_deg_ref > 0]), 3), top100 = round(median(M1_top100_overlap), 3),
       tSpearman = round(median(M1_t_spearman), 3), sexgene_M0 = round(as.numeric(median(M1_sexgene_deg_ref)), 1)), by=design][order(design)])
cat("\n== M0 vs M2(+sex+age) change by grade ==\n")
print(v[!is.na(design) & !is.na(M2_jaccard), .(n = .N, jaccard = round(median(M2_jaccard[M2_n_deg_ref > 0]), 3),
       top100 = round(median(M2_top100_overlap), 3), tSpearman = round(median(M2_t_spearman), 3)), by=design][order(design)])
cat("\n== By tissue (partially_confounded or worse) ==\n")
print(v[design %in% c("partially_confounded","completely_confounded") & !is.na(M1_jaccard),
        .(n = .N, jaccard = round(median(M1_jaccard[M1_n_deg_ref > 0]), 2), tSp = round(median(M1_t_spearman), 3)), by=tissue_category][order(-n)][1:12])
cat("\n== Source of sex information ==\n"); print(v[, .(n = .N, sex_meta_missing = sum(n_sex_meta_missing), mismatch = sum(n_sex_mismatch))])
