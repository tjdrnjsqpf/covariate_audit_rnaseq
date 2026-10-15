#!/usr/bin/env Rscript
## Advisor comment Q4: do d_across and d_within depend on the DE engine (limma-voom / edgeR QL / DESeq2) when the same designs and stored permuted labels are refitted?
## Input: results/response/engine/<SRP>.tsv (54), class labels from results/main_{human,mouse}_summary_v2.tsv. Output: results/response/E_engine_by_class.tsv, E_engine_pairs.tsv
suppressMessages(library(data.table)); PROJ <- Sys.getenv("COV_AUDIT_PROJ", "/path/to/covariate_audit_rnaseq"); setwd(PROJ)
L <- rbind(fread("results/main_human_summary_v2.tsv", select=c("study","design_sex")), fread("results/main_mouse_summary_v2.tsv", select=c("study","design_sex"))); fl <- list.files("results/response/engine", "\\.tsv$", full.names=TRUE)
E <- rbindlist(lapply(fl, fread)); E <- merge(E, L[, .(study, design_sex)], by="study")
E[, engine := factor(engine, c("limma","edgeR","DESeq2"))]
med <- function(x) { x <- x[is.finite(x)]; if (length(x) < 3) return(sprintf("%.3f", median(x))); b <- replicate(2000, median(sample(x, replace=TRUE))); sprintf("%.3f [%.3f, %.3f]", median(x), quantile(b, .025), quantile(b, .975)) }
S <- E[, .(n=uniqueN(study), n_de_M0=as.numeric(median(n_de_M0)), J_M1=med(J_M1), J_perm=med(J_perm), J_permw=med(J_permw), d_across=med(d_across), d_within=med(d_within),
           collinearity_cost=med(J_permw - J_perm), frac_within_lt_across=mean(abs(d_within) < abs(d_across), na.rm=TRUE)), by=.(organism, design_sex, engine)][order(organism, design_sex, engine)]
fwrite(S, "results/response/E_engine_by_class.tsv", sep="\t"); print(S, nrows=100)
## Per-study agreement between engines: Spearman correlation of Δ with limma and median difference
W <- dcast(E, study + organism + design_sex ~ engine, value.var=c("d_across","d_within","J_M1"))
P <- rbindlist(lapply(c("edgeR","DESeq2"), function(e) W[, .(engine=e,
  rho_across=cor(d_across_limma, get(paste0("d_across_", e)), method="spearman", use="complete"), rho_within=cor(d_within_limma, get(paste0("d_within_", e)), method="spearman", use="complete"),
  med_diff_across=median(get(paste0("d_across_", e)) - d_across_limma, na.rm=TRUE), med_diff_within=median(get(paste0("d_within_", e)) - d_within_limma, na.rm=TRUE), n=.N), by=.(organism, design_sex)]))
fwrite(P, "results/response/E_engine_pairs.tsv", sep="\t"); print(P)
