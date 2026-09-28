#!/usr/bin/env Rscript
# Call sex for all mouse samples from the extracted Xist/Y counts -> compare with reported sex -> update prevalence figures (reported ∪ inferred)
suppressMessages(library(data.table))
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; PROC <- file.path(PROJ, "data/processed"); RES <- file.path(PROJ, "results")
fs <- list.files(file.path(PROC, "mouse_sexgenes"), pattern="tsv$", full.names=TRUE)
x <- rbindlist(lapply(fs, fread, header=FALSE, col.names=c("study","external_id","libsize","xist","ddx3y","eif2s3y","kdm5d","uty")))
cat("projects:", uniqueN(x$study), " samples:", nrow(x), "\n")
for (v in c("libsize","xist","ddx3y","eif2s3y","kdm5d","uty")) x[[v]] <- as.numeric(x[[v]])   # avoid integer64 overflow
lc <- function(v, L) log2(v / L * 1e6 + 1)
x[, `:=`(xist_l = lc(xist, libsize), y_l = (lc(ddx3y, libsize) + lc(eif2s3y, libsize) + lc(kdm5d, libsize) + lc(uty, libsize)) / 4)]
x[, low_cov := libsize < 1e7]   # flag samples with very small totals (low depth / targeted)
x[, sex_expr := fifelse(xist_l - y_l >= 2, "F", fifelse(y_l - xist_l >= 2, "M", "ambiguous"))]
x[, ambig_type := fifelse(sex_expr != "ambiguous", NA_character_, fifelse(xist_l > 2 & y_l > 2, "both_high", fifelse(xist_l <= 2 & y_l <= 2, "both_low", "intermediate")))]
cat("-- expression sex calls --\n"); print(round(100*prop.table(table(x$sex_expr)), 1)); print(table(x$ambig_type))

## Compare with reported sex
s <- fread(file.path(PROC, "mouse_sample_table.tsv.gz"), colClasses=list(character=c("study","external_id","sex")), na.strings=c("","NA"))
m <- merge(s, x[, .(study, external_id, libsize, xist_l, y_l, sex_expr, ambig_type, low_cov)], by=c("study","external_id"), all.x=TRUE)
both <- m[sex %in% c("F","M") & sex_expr %in% c("F","M")]
cat("\n== reported vs expression (confident) ==\n"); print(table(reported = both$sex, expr = both$sex_expr))
cat("agreement:", sprintf("%.2f%%", 100*mean(both$sex == both$sex_expr)), " n =", nrow(both), "\n")
## Diagnostics: are mismatches concentrated in specific projects, signal strength of mismatched samples, tissue distribution
both[, mism := sex != sex_expr]
pm <- both[, .(n = .N, n_mism = sum(mism), frac = mean(mism)), by=study]
cat("-- per-project mismatch fraction among projects with >=1 mismatch --\n"); print(summary(pm[n_mism > 0, frac]))
cat("projects with mismatch frac >= 0.5 (likely inverted/scrambled labels):", pm[n_mism > 0 & frac >= 0.5, .N], " samples in them:", pm[n_mism > 0 & frac >= 0.5, sum(n)], " their mismatches:", pm[frac >= 0.5, sum(n_mism)], "\n")
cat("mismatches from projects with frac < 0.2 (sporadic):", pm[frac < 0.2, sum(n_mism)], "\n")
cat("-- mismatched samples: xist_l / y_l quantiles --\n"); print(both[mism == TRUE, .(xist_q25 = quantile(xist_l, .25), xist_med = median(xist_l), xist_q75 = quantile(xist_l, .75), y_q25 = quantile(y_l, .25), y_med = median(y_l), y_q75 = quantile(y_l, .75)), by=.(reported = sex, expr = sex_expr)])
cat("-- concordant samples for reference --\n"); print(both[mism == FALSE, .(xist_med = median(xist_l), y_med = median(y_l)), by=.(reported = sex)])
cat("-- mismatch rate by tissue --\n"); print(both[, .(n = .N, mism = round(100*mean(mism), 1)), by=tissue_category][order(-n)])
cat("-- top 12 projects by mismatch count --\n"); print(merge(pm[order(-n_mism)][1:12], unique(s[, .(study, study_title = substr(study_title, 1, 70), tissue_category)])[!duplicated(study)], by="study")[order(-n_mism)])
## Suspected project-level label flip (frac>=0.5): do not trust reported sex, replace with the expression call
inv <- pm[frac >= 0.5, study]
m[study %in% inv & sex_expr %in% c("F","M"), sex := sex_expr]
cat("reported NA & expr confident:", m[is.na(sex) & sex_expr %in% c("F","M"), .N], " of reported NA:", m[is.na(sex) & !is.na(sex_expr), .N], "\n")
m[, sex_final := fifelse(sex %in% c("F","M"), sex, fifelse(sex_expr %in% c("F","M"), sex_expr, NA_character_))]
m[, sex_source := fifelse(sex %in% c("F","M"), "reported", fifelse(sex_expr %in% c("F","M"), "inferred", "none"))]
m[, sex_mismatch := sex %in% c("F","M") & sex_expr %in% c("F","M") & sex != sex_expr]
fwrite(m, file.path(PROC, "mouse_sample_table_sexinferred.tsv.gz"), sep="\t")

## Re-aggregate at project level
p0 <- fread(file.path(RES, "mouse_project_summary_all.tsv"))
pj <- m[, .(n_expr = sum(!is.na(sex_expr)), n_final = sum(!is.na(sex_final)), n_F_final = sum(sex_final == "F", na.rm=TRUE), n_M_final = sum(sex_final == "M", na.rm=TRUE),
            n_ambig = sum(sex_expr == "ambiguous", na.rm=TRUE), n_both_high = sum(ambig_type == "both_high", na.rm=TRUE), n_mismatch = sum(sex_mismatch), n_low_cov = sum(low_cov, na.rm=TRUE)), by=study]
p <- merge(p0, pj, by="study", all.x=TRUE)
p[, sex_design_final := fifelse(is.na(n_final) | n_final/n < 0.9, "unresolved", fifelse(n_F_final > 0 & n_M_final == 0, "female_only", fifelse(n_M_final > 0 & n_F_final == 0, "male_only", "mixed")))]
base <- p[study %in% readLines(file.path(PROC, "mouse_base_projects.txt"))]
cat("\n== MOUSE base projects:", nrow(base), "==\n")
cat("sex resolvable (reported ∪ inferred) >=90%:", sum(base$sex_design_final != "unresolved"), sprintf("(%.1f%%)  [reported only: %.1f%%]\n", 100*mean(base$sex_design_final != "unresolved"), 100*mean(base$sex_design == "unreported") * 0 + 100*mean(base$n_sex/base$n >= .9)))
cat("-- sex design, reported only (among reported) --\n"); print(round(100*prop.table(table(base[sex_design != "unreported", sex_design])), 1))
cat("-- sex design, reported ∪ inferred (among resolvable) --\n"); print(round(100*prop.table(table(base[sex_design_final != "unresolved", sex_design_final])), 1))
cat("-- sex design among previously UNREPORTED projects now resolved --\n"); print(round(100*prop.table(table(base[sex_design == "unreported" & sex_design_final != "unresolved", sex_design_final])), 1))
cat("-- projects with >=1 reported/expression mismatch:", sum(base$n_mismatch > 0, na.rm=TRUE), " of", sum(base$n_sex > 0), "with reported sex\n")
cat("-- ambiguous samples:", sum(base$n_ambig, na.rm=TRUE), " both_high:", sum(base$n_both_high, na.rm=TRUE), " low_cov:", sum(base$n_low_cov, na.rm=TRUE), "\n")

## Between-group imbalance in (final) mixed-sex studies
al <- fread(file.path(PROC, "mouse_sample_attributes_long.tsv.gz"), colClasses="character", na.strings=c("","NA"))
imb <- function(st, gkey) { a <- al[study == st & key == gkey, .(external_id, grp = value)]; z <- merge(a, m[study == st, .(external_id, sex_final)], by="external_id")[sex_final %in% c("F","M")]
  top2 <- names(sort(table(z$grp), dec=TRUE))[1:2]; z <- z[grp %in% top2]; if (uniqueN(z$grp) < 2) return(NA_real_); pf <- z[, mean(sex_final == "F"), by=grp]$V1; abs(diff(pf)) }
mixed <- base[sex_design_final == "mixed" & !is.na(group_key) & group_score > 0]
mixed[, sex_diff_final := mapply(imb, study, group_key)]
mixed[, design_final := fifelse(is.na(sex_diff_final), NA_character_, fifelse(sex_diff_final >= 0.99, "completely_confounded", fifelse(sex_diff_final >= 0.3, "partially_confounded", fifelse(sex_diff_final >= 0.15, "mildly_imbalanced", "balanced"))))]
cat("\n== mixed-sex (final) projects with group key:", nrow(mixed), "==\n"); print(table(mixed$design_final, useNA="ifany"))
p <- merge(p, mixed[, .(study, sex_diff_final, design_final)], by="study", all.x=TRUE)
fwrite(p, file.path(RES, "mouse_project_summary_sexinferred.tsv"), sep="\t")
bt <- base[, .(projects=.N, resolvable=round(100*mean(sex_design_final != "unresolved")), female_only=round(100*mean(sex_design_final=="female_only")), male_only=round(100*mean(sex_design_final=="male_only")), mixed=round(100*mean(sex_design_final=="mixed"))), by=tissue_category][order(-projects)]
cat("\n-- by tissue (reported ∪ inferred): projects / resolvable% / female_only% / male_only% / mixed% --\n"); print(bt)
fwrite(bt, file.path(RES, "mouse_sexdesign_by_tissue_inferred.tsv"), sep="\t")
