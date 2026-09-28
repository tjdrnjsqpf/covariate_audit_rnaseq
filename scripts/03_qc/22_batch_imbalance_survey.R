#!/usr/bin/env Rscript
## C6 survey: Cramér V between group and technical batch variables (instrument, layout, center, publication month, read length) across all valid studies — metadata only
suppressMessages(library(data.table)); ORG <- commandArgs(trailingOnly=TRUE)[1]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
MD <- if (ORG == "human") "data/raw/recount3_metadata" else "data/raw/recount3_metadata_mouse"
d <- fread(sprintf("results/main_%s_summary_v2.tsv", ORG))[valid == TRUE]
cramer <- function(x, g) { ok <- !is.na(x) & x != ""; if (sum(ok) < 4) return(NA_real_); tab <- table(x[ok], g[ok]); if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct=FALSE)$statistic); round(as.numeric(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1)))), 3) }
one <- function(s) {
  f <- file.path(MD, sprintf("sra.sra.%s.MD.gz", s)); pf <- sprintf("results/main_%s_v2/%s/pheno.tsv", ORG, s)
  if (!file.exists(f) || !file.exists(pf)) return(NULL)
  md <- fread(f, select=c("external_id","platform_model","library_layout","run_center_name","run_published","submission_acc"), na.strings=c("","NA"))
  ph <- fread(pf)[, .(external_id, group)]; m <- merge(ph, md, by="external_id")
  if (nrow(m) < 4) return(NULL)
  m[, pubmonth := substr(run_published, 1, 7)]
  v <- sapply(c(platform="platform_model", layout="library_layout", center="run_center_name", pubmonth="pubmonth", submission="submission_acc"), function(k) cramer(m[[k]], m$group))
  nl <- sapply(c("platform_model","library_layout","run_center_name","pubmonth","submission_acc"), function(k) uniqueN(m[[k]][!is.na(m[[k]])]))
  data.table(study = s, n = nrow(m), n_platform = nl[1], n_layout = nl[2], n_center = nl[3], n_pubmonth = nl[4], n_submission = nl[5],
             v_platform = v[["platform"]], v_layout = v[["layout"]], v_center = v[["center"]], v_pubmonth = v[["pubmonth"]], v_submission = v[["submission"]],
             batch_v_max = suppressWarnings(max(v, na.rm=TRUE)), batch_var = if (all(is.na(v))) NA_character_ else names(v)[which.max(v)])
}
r <- rbindlist(lapply(d$study, function(s) tryCatch(one(s), error=function(e) NULL)))
r[is.infinite(batch_v_max), batch_v_max := NA_real_]
r <- merge(r, d[, .(study, design_sex, single_sex, sex_diff)], by="study")
fwrite(r, sprintf("results/qc/batch_imbalance_%s.tsv", ORG), sep="\t")
cat(sprintf("########## Technical batch survey (%s, valid %d) ##########\n", ORG, nrow(r)))
cat(sprintf("Studies with batch structure (any variable with >=2 levels): %d (%.1f%%)\n", sum(!is.na(r$batch_v_max)), 100*mean(!is.na(r$batch_v_max))))
cat(sprintf("  >=2 publication months: %.1f%%  >=2 instruments: %.1f%%  >=2 submissions: %.1f%%\n", 100*mean(r$n_pubmonth > 1), 100*mean(r$n_platform > 1), 100*mean(r$n_submission > 1)))
b <- r[!is.na(batch_v_max)]
cat(sprintf("Of %d studies with batch structure, group–batch Cramér V ≥ 0.5: %d (%.1f%%);  V = 1 (complete confounding): %d (%.1f%%)\n", nrow(b), sum(b$batch_v_max >= 0.5), 100*mean(b$batch_v_max >= 0.5), sum(b$batch_v_max >= 0.999), 100*mean(b$batch_v_max >= 0.999)))
cat("With all valid studies as denominator: V≥0.5 ", round(100*sum(b$batch_v_max >= 0.5)/nrow(r), 1), "%,  V=1 ", round(100*sum(b$batch_v_max >= 0.999)/nrow(r), 1), "%\n")
cat("\n== Most entangled batch variable ==\n"); print(table(b$batch_var))
cat("\n== Batch confounding by sex class (fraction with V≥0.5) ==\n"); print(r[, .(n = .N, has_batch = round(mean(!is.na(batch_v_max)), 2), v_ge_0.5 = round(mean(batch_v_max >= 0.5, na.rm=TRUE), 2), v_eq_1 = round(mean(batch_v_max >= 0.999, na.rm=TRUE), 2)), by=design_sex][order(design_sex)])
cat(sprintf("\nsex_diff vs batch_v_max Spearman: %.3f\n", cor(r$sex_diff, r$batch_v_max, method="spearman", use="complete.obs")))
