#!/usr/bin/env Rscript
## Gate 2 summary: adjustment practice in the original papers, crossed with our design class and Δ
suppressMessages({library(data.table); library(jsonlite)}); ORG <- commandArgs(trailingOnly=TRUE)[1]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
ft <- rbindlist(lapply(readLines(sprintf("results/gate2/fulltext_%s.jsonl", ORG)), function(l) { x <- fromJSON(l); data.table(study = x$study, route = x$route %||% NA_character_, n_pmid = length(x$pmids), has_fulltext = !is.null(x$methods)) }), fill=TRUE)
ad <- rbindlist(lapply(readLines(sprintf("results/gate2/adjust_%s.jsonl", ORG)), function(l) { x <- fromJSON(l); data.table(study = x$study, text_relevant = x$text_relevant, de_method = x$de_method %||% NA_character_, sex_adjusted = x$sex_adjusted, age_adjusted = x$age_adjusted, sex_matched = x$sex_matched_design, age_matched = x$age_matched_design, batch_adj = x$batch_or_latent_adjusted, sex_reported = x$sex_reported_in_demographics, confidence = x$confidence, usd = x$usd) }), fill=TRUE)
m <- fread(sprintf("results/main_%s_summary_v2.tsv", ORG))[, .(study, valid, design_sex, design_age, single_sex, M1_jaccard, M0perm_jaccard, M2_jaccard, M2perm_jaccard, tissue_category)]
d <- merge(merge(m, ft, by="study", all.x=TRUE), ad, by="study", all.x=TRUE); fwrite(d, sprintf("results/gate2/gate2_merged_%s.tsv", ORG), sep="\t")
cat(sprintf("########## Gate 2 (%s) ##########\n", ORG))
cat(sprintf("valid studies %d → paper PMID found %d (%.0f%%) → OA full text %d (%.0f%%) → LLM-judged %d, relevant text %d\n", sum(d$valid), sum(d$n_pmid > 0, na.rm=TRUE), 100*mean(d$n_pmid > 0, na.rm=TRUE), sum(d$has_fulltext, na.rm=TRUE), 100*mean(d$has_fulltext, na.rm=TRUE), sum(!is.na(d$sex_adjusted)), sum(d$text_relevant, na.rm=TRUE)))
r <- d[text_relevant == TRUE]
cat("\n== Adjustment practice in the original papers (studies with relevant text) ==\n")
for (k in c("sex_adjusted","age_adjusted","sex_matched","age_matched","batch_adj","sex_reported")) { t <- table(r[[k]]); cat(sprintf("%-14s yes %3d (%4.1f%%)  no %3d  unclear %3d\n", k, t["yes"] %||% 0, 100*(t["yes"] %||% 0)/nrow(r), t["no"] %||% 0, t["unclear"] %||% 0)) }
cat("\n== Sex adjusted or not x our design class (mixed-sex studies) ==\n"); print(table(r[single_sex == FALSE, design_sex], r[single_sex == FALSE, sex_adjusted], useNA="ifany"))
cat("\n== Change we observed in studies that did not adjust (Δ = M1 − perm) ==\n")
print(r[!is.na(M1_jaccard) & !is.na(M0perm_jaccard), .(n = .N, dJ = round(median(M1_jaccard - M0perm_jaccard), 3)), by=.(sex_adjusted, design_sex)][order(design_sex, sex_adjusted)])
cat("\n== DE tools ==\n"); print(head(sort(table(tolower(r$de_method)), decreasing=TRUE), 8))
cat(sprintf("\nTotal LLM cost: $%.2f\n", sum(ad$usd, na.rm=TRUE)))
