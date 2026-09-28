#!/usr/bin/env Rscript
# Collect main reanalysis results, v2 (revised). Usage: Rscript 11b_collect_main_summary_v2.R <outdir> <config> <organism>
# Changes from v1: B1 common denominator, B2 single_sex category, B5 separate sex/age grades,
#               B6 bootstrap CI + minimum n in the tissue table, C1 Δ reported relative to the negative-control (M0perm) baseline
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE); OUTDIR <- args[1]; CFG <- args[2]; ORG <- args[3]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
fs <- list.files(file.path(PROJ, OUTDIR), pattern="summary.tsv", recursive=TRUE, full.names=TRUE)
d <- rbindlist(lapply(fs, fread), fill=TRUE)
cfg <- fread(file.path(PROJ, CFG))[, .(study, note)]
n_dirs <- nrow(d)
d <- merge(cfg, d, by="study")   # config is authoritative: studies removed from the config are excluded even if a results folder exists
cat("Result folders:", n_dirs, "; excluded as not in config:", n_dirs - nrow(d), "\n")
ud <- fread(file.path(PROJ, "results/design_llm/usable_designs_v1.tsv"))[organism == ORG, .(study, design_type, experimental_setting, sex_composition, scope_flags, confidence)]
d <- merge(d, ud, by="study", all.x=TRUE)
if (ORG == "human") { ps <- fread(file.path(PROJ, "results/project_summary_all.tsv"))[, .(study, tissue_category, first_published)]
} else ps <- fread(file.path(PROJ, "results/mouse_project_summary_sexinferred.tsv"))[, .(study, tissue_category, first_published = NA_character_)]
d <- merge(d, ps, by="study", all.x=TRUE)

## ---- Validity ----
d[, valid := is.na(n_genes_detected) | n_genes_detected >= 5000]

## ---- B2: Separate single-sex studies ----
d[, single_sex := !is.na(n_female) & !is.na(n_male) & (n_female == 0 | n_male == 0)]

## ---- B5: Grade sex and age imbalance separately ----
grade <- function(x, cuts) fifelse(is.na(x), NA_character_,
           fifelse(x >= cuts[3], "completely_confounded",
           fifelse(x >= cuts[2], "partially_confounded",
           fifelse(x >= cuts[1], "mildly_imbalanced", "balanced"))))
d[, design_sex := grade(sex_diff, c(0.15, 0.30, 0.99))]
d[, design_age := grade(age_smd, c(0.40, 0.80, 1.60))]
d[single_sex == TRUE, design_sex := "single_sex"]
## Combined grade for backward compatibility (the worse of the two) — not used in reporting, kept only for comparison with v1
lv <- c("balanced","mildly_imbalanced","partially_confounded","completely_confounded")
d[, design := fifelse(single_sex == TRUE, "single_sex",
       fifelse(is.na(design_sex) & is.na(design_age), NA_character_,
         lv[pmax(match(design_sex, lv), match(design_age, lv), na.rm=TRUE)]))]

fwrite(d, file.path(PROJ, sprintf("results/main_%s_summary_v2.tsv", ORG)), sep="\t")

## ---- B1: Compute all metrics on the same analysis set ----
## Analysis set = studies that are valid & sex-adjustable (M1 fitted) & have at least one DEG in M0
## (with M0 DEG=0 the Jaccard is undefined, so no metric can use such a study)
A <- d[valid == TRUE & !is.na(M1_jaccard) & M1_n_deg_ref > 0]
cat("== Analysis set (B1: common denominator for all metrics) ==\n")
cat("Studies with results:", nrow(d), " / valid:", sum(d$valid),
    " / M1 fitted:", sum(d$valid & !is.na(d$M1_jaccard)),
    " / M0 DEG>0 (analysis set):", nrow(A), "\n")
cat("Single-sex (sex adjustment impossible):", sum(d$valid & d$single_sex, na.rm=TRUE), "\n\n")

boot_ci <- function(x, R = 2000) {
  x <- x[!is.na(x)]; if (length(x) < 3) return(c(NA_real_, NA_real_))
  set.seed(1); b <- replicate(R, median(sample(x, length(x), replace=TRUE)))
  round(quantile(b, c(0.025, 0.975)), 3)
}
med_ci <- function(x) { ci <- boot_ci(x); sprintf("%.3f [%.3f-%.3f]", median(x, na.rm=TRUE), ci[1], ci[2]) }

## ---- C1: Negative-control baseline ----
cat("== C1 negative control: change caused merely by adding a meaningless covariate (randomly permuted sex) ==\n")
print(A[!is.na(M0perm_jaccard), .(n = .N, perm_jaccard = med_ci(M0perm_jaccard),
        perm_top100 = med_ci(M0perm_top100_overlap), perm_tSp = med_ci(M0perm_t_spearman))])
cat("\n")

## ---- By sex grade (incl. Δ relative to baseline) ----
cat("== B5/C1: M0 vs M1(+sex) by sex-imbalance grade ==\n")
tab <- A[!is.na(design_sex), .(n = .N,
      deg_M0 = as.numeric(median(M1_n_deg_ref)), deg_M1 = as.numeric(median(M1_n_deg_alt)),
      jaccard = med_ci(M1_jaccard), jaccard_perm = med_ci(M0perm_jaccard),
      d_jaccard = round(median(M1_jaccard - M0perm_jaccard, na.rm=TRUE), 3),
      jaccard_permw = if ("M0permw_jaccard" %in% names(A)) med_ci(M0permw_jaccard) else NA_character_,
      d_jaccard_w = if ("M0permw_jaccard" %in% names(A)) round(median(M1_jaccard - M0permw_jaccard, na.rm=TRUE), 3) else NA_real_,
      top100 = med_ci(M1_top100_overlap), tSpearman = med_ci(M1_t_spearman)),
      by=design_sex][order(match(design_sex, c("balanced", lv[-1], "single_sex")))]
print(tab)

cat("\n== B5: M0 vs M2(+sex+age) by age-imbalance grade — separate from sex grade ==\n")
A2 <- d[valid == TRUE & !is.na(M2_jaccard) & M2_n_deg_ref > 0]
has_m2perm <- "M2perm_jaccard" %in% names(A2)
if (!has_m2perm) A2[, M2perm_jaccard := NA_real_]
print(A2[!is.na(design_age), .(n = .N, jaccard = med_ci(M2_jaccard),
      jaccard_perm_age = if (has_m2perm) med_ci(M2perm_jaccard) else NA_character_,
      d_jaccard = round(median(M2_jaccard - M2perm_jaccard, na.rm=TRUE), 3),
      top100 = med_ci(M2_top100_overlap), tSpearman = med_ci(M2_t_spearman)),
      by=design_age][order(match(design_age, lv))])
if (has_m2perm) { cat("\n== C1 negative control (age): change caused merely by adding permuted age (M2 analysis set) ==\n")
  print(A2[!is.na(M2perm_jaccard), .(n = .N, perm_jaccard = med_ci(M2perm_jaccard), perm_top100 = med_ci(M2perm_top100_overlap), perm_tSp = med_ci(M2perm_t_spearman))]) }
cat("\n== B5 cross-table: sex grade x age grade (M2 analysis set) ==\n")
print(table(A2$design_sex, A2$design_age, useNA="ifany"))

## ---- B6: By tissue — minimum n=5, with CI, exploratory ----
cat("\n== B6 (exploratory, n>=5 only): M1 Jaccard by tissue, sex imbalance partially or worse ==\n")
tt <- A[design_sex %in% c("partially_confounded","completely_confounded"),
        .(n = .N, jaccard = med_ci(M1_jaccard), tSp = round(median(M1_t_spearman), 3)), by=tissue_category]
tt <- tt[n >= 5][order(-n)]
if (nrow(tt) == 0) cat("  no tissue satisfies n>=5 — tissue comparison not reported\n") else print(tt)

## ---- B7 (mouse): By comparison type ----
cat("\n== B7: By experiment type (design_type) ==\n")
print(A[, .(n = .N, jaccard = round(median(M1_jaccard), 3),
        perm = round(median(M0perm_jaccard, na.rm=TRUE), 3)), by=design_type][order(-n)][1:8])

cat("\n== Source of sex information ==\n"); print(d[valid == TRUE, .(n = .N, sex_meta_missing = sum(n_sex_meta_missing), mismatch = sum(n_sex_mismatch))])
