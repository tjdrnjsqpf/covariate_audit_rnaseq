#!/usr/bin/env Rscript
## Source data for Keynote charts: histogram → counts per bin, boxplot → summary statistics, point/line+CI → value + error-bar length (±).
## Run: LC_ALL=en_US.UTF-8 COV_AUDIT_PROJ="$(pwd)" Rscript scripts/05_report/29_chart_ready_source.R  → results/figures/draft_v3/source_data_keynote/*.tsv
suppressMessages(library(data.table)); PROJ <- Sys.getenv("COV_AUDIT_PROJ", "/var2/lsg/Claude_Code/covariate_audit_rnaseq"); setwd(PROJ)
SD <- "results/figures/draft_v3/source_data"; KD <- "results/figures/draft_v3/source_data_keynote"; dir.create(KD, showWarnings=FALSE)
out <- function(d, f) fwrite(d, file.path(KD, paste0(f, ".tsv")), sep="\t")
## 1A: number of studies per 0.05 bin (columns = Human, Mouse)
## 1B: the per-bin count table was already produced by script 26 with the same computation as the figure, so only copy it
invisible(file.copy(file.path(SD, "Fig1B_histogram_bins.tsv"), file.path(KD, "Fig1B_histogram_bins.tsv"), overwrite=TRUE))
## 1B: box summary (median, quartiles, whiskers) — drawn in Keynote as bar (median) + error bars (Q1–Q3)
b <- melt(fread(file.path(SD, "Fig2A_jaccard_three_controls.tsv")), id.vars=c("study","organism","cls"), variable.name="series", value.name="J")
bs <- b[, { q <- quantile(J, c(.25,.5,.75), na.rm=TRUE); iqr <- q[3]-q[1]; .(n=.N, median=q[2], q1=q[1], q3=q[3], whisker_low=max(min(J, na.rm=TRUE), q[1]-1.5*iqr), whisker_high=min(max(J, na.rm=TRUE), q[3]+1.5*iqr), err_minus=q[2]-q[1], err_plus=q[3]-q[2]) }, by=.(organism, cls, series)]
out(bs[order(organism, cls, series)], "Fig2A_box_summary"); out(dcast(bs, organism + cls ~ series, value.var="median"), "Fig2A_medians_wide")
## Point/line + CI panels: values in wide format + error-bar lengths
wide_err <- function(f, xcol, scol, mcol="m", lo="lo", hi="hi") { d <- fread(file.path(SD, paste0(f, ".tsv"))); d$err_minus <- d[[mcol]] - d[[lo]]; d$err_plus <- d[[hi]] - d[[mcol]]
  v <- dcast(d, as.formula(sprintf("organism + %s ~ %s", xcol, scol)), value.var=c(mcol, "err_minus", "err_plus")); out(v, paste0(f, "_wide")) }
wide_err("Fig2B_dose_response_two_controls", "bin", "control")
for (f in c("Fig3A_synthetic_delta_two_controls","Fig3B_synthetic_false_de_fraction","Fig3C_synthetic_recall","Fig3D_synthetic_t_correlation")) wide_err(f, "target_df", "variable")
e <- fread(file.path(SD, "Fig3E_synthetic_leakage.tsv")); out(e, "Fig3E_synthetic_leakage_wide")
## Bar panels are already wide — add percentage columns
for (f in c("Fig4A_pathway_conclusion_changed","Fig4B_top_pathway_changed","FigS1_observed_leakage")) { d <- fread(file.path(SD, paste0(f, ".tsv"))); num <- setdiff(names(d), c("organism","cls","n")); d[, (paste0(num, "_pct")) := lapply(.SD, function(x) round(100*x, 1)), .SDcols=num]; out(d, paste0(f, "_wide")) }
s <- fread(file.path(SD, "Fig4C_sva_capture.tsv")); out(dcast(s[, pct := round(100*frac,1)], design_sex ~ organism, value.var="pct"), "Fig4C_sva_capture_wide")
g <- fread(file.path(SD, "Fig4D_publication_audit.tsv")); out(dcast(g[, pct := round(100*frac,1)], item ~ organism, value.var="pct"), "Fig4D_publication_audit_wide"); out(g, "Fig4D_publication_audit_counts")
invisible(file.copy(file.path(SD, "Fig1_schematic_counts.tsv"), file.path(KD, "Fig1_schematic_counts.tsv"), overwrite=TRUE))
cat("chart-ready tables:", length(list.files(KD)), "\n")
