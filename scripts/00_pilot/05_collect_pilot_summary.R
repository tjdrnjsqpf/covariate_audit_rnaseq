#!/usr/bin/env Rscript
# Collect results/pilot/*/summary.tsv into one table with batch (core/conditional/nonblood) and tissue category.
suppressMessages(library(data.table))
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
fs <- list.files(file.path(PROJ, "results/pilot"), pattern="summary.tsv", recursive=TRUE, full.names=TRUE)
d <- rbindlist(lapply(fs, fread), fill=TRUE)
cfg <- rbindlist(lapply(c(core="pilot_core_runs.tsv", conditional="pilot_conditional_runs.tsv", nonblood="pilot_nonblood_runs.tsv"),
                        function(f) fread(file.path(PROJ, "config", f))[, .(study, note)]), idcol="batch")
ps <- fread(file.path(PROJ, "results/project_summary_all.tsv"))[, .(study, tissue_category, tissue_example, study_title)]
d <- merge(cfg, d, by="study", all.x=TRUE); d <- merge(d, ps, by="study", all.x=TRUE)
## Invalid: targeted sequencing rather than whole transcriptome (detected genes < 5,000) or manually confirmed
d[, valid := !is.na(n)]
d[!is.na(n_genes_detected) & n_genes_detected < 5000, valid := FALSE]
d[study %in% c("SRP100652","SRP164913","SRP065758"), valid := FALSE]
d[study == "SRP100652", note := paste(note, "| invalid: ERAP targeted amplicon")]
d[study == "SRP164913", note := paste(note, "| invalid: TCR repertoire")]
d[study == "SRP065758", note := paste(note, "| invalid: 1 detected gene")]
d[, design := fifelse(is.na(sex_diff), NA_character_, fifelse(sex_diff >= 0.99, "completely_confounded",
              fifelse(sex_diff >= 0.3 | (!is.na(age_smd) & age_smd >= 0.8), "partially_confounded",
              fifelse(sex_diff >= 0.15 | (!is.na(age_smd) & age_smd >= 0.4), "mildly_imbalanced", "balanced"))))]
setcolorder(d, c("study","batch","tissue_category","valid","design","note"))
fwrite(d, file.path(PROJ, "results/pilot_all_summary.tsv"), sep="\t")
cat("runs:", nrow(d), " with results:", sum(!is.na(d$n)), " valid:", sum(d$valid, na.rm=TRUE), "\n")
print(table(d$batch, is.na(d$n)))
print(d[!is.na(n), .(study, batch, tissue_category, design, n, n_case, n_control, sex_diff, age_smd, n_sex_mismatch, n_genes_detected,
                     M2_ref = fifelse(is.na(M2_n_deg_ref), M1_n_deg_ref, M2_n_deg_ref), M2_alt = M2_n_deg_alt, M2_jaccard, M2_top100_overlap, M2_t_spearman)])
