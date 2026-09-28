#!/usr/bin/env Rscript
suppressMessages(library(data.table))
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
fs <- list.files(file.path(PROJ, "results/pilot_mouse"), pattern="summary.tsv", recursive=TRUE, full.names=TRUE)
d <- rbindlist(lapply(fs, fread), fill=TRUE)
cfg <- rbindlist(lapply(c(batch1="pilot_mouse_runs.tsv", batch2="pilot_mouse_runs2.tsv"), function(f) fread(file.path(PROJ, "config", f))[, .(study, note)]), idcol="batch")
ps <- fread(file.path(PROJ, "results/mouse_project_summary_sexinferred.tsv"))[, .(study, tissue_category, strain_top, design_final, study_title)]
d <- merge(cfg, d, by="study", all.x=TRUE); d <- merge(d, ps, by="study", all.x=TRUE)
d[, valid := !is.na(n) & (is.na(n_genes_detected) | n_genes_detected >= 5000)]
## Invalid studies revealed by design adjudication (2026-09-07): SRP059331 (4 biological samples; the 40 runs are technical replicates), SRP080852 (in vitro stimulation)
d[study %in% c("SRP059331","SRP080852"), valid := FALSE]
d[study == "SRP059331", note := paste(note, "| invalid: technical replicates (biological n=4)")]
d[study == "SRP080852", note := paste(note, "| invalid: in vitro stimulation")]
d[, design := fifelse(is.na(sex_diff), NA_character_, fifelse(sex_diff >= 0.99, "completely_confounded", fifelse(sex_diff >= 0.3, "partially_confounded", fifelse(sex_diff >= 0.15, "mildly_imbalanced", "balanced"))))]
setcolorder(d, c("study","batch","tissue_category","valid","design","note"))
cat("\n== by tissue x design: median M1 Jaccard / top100 / t-Spearman / M0 sex-gene leakage ==\n")
print(d[valid == TRUE & !is.na(M1_jaccard), .(n=.N, jaccard=round(median(M1_jaccard[M1_n_deg_ref > 0]), 2), top100=round(median(M1_top100_overlap), 2), tSp=round(median(M1_t_spearman), 3), sexg=round(median(M1_sexgene_deg_ref), 1)), by=.(tissue_category, design)][order(tissue_category, design)])
fwrite(d, file.path(PROJ, "results/pilot_mouse_summary.tsv"), sep="\t")
cat("runs:", nrow(d), " with results:", sum(!is.na(d$n)), " valid:", sum(d$valid), "\n")
print(d[!is.na(n), .(study, tissue_category, design, n, n_case, n_control, sex_diff, n_sex_mismatch, n_genes_detected, models,
                     M0_deg = fifelse(is.na(M2_n_deg_ref), M1_n_deg_ref, M2_n_deg_ref), M1_deg = M1_n_deg_alt, M1_jaccard, M1_top100 = M1_top100_overlap, M1_tSp = M1_t_spearman, sexg_M0 = M1_sexgene_deg_ref)])
