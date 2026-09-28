#!/usr/bin/env Rscript
## A5: check of the expression-based sex-call threshold (|XIST − Y| ≥ 2 log2CPM)
suppressMessages({library(data.table); library(ggplot2)})
setwd("/var2/lsg/Claude_Code/covariate_audit_rnaseq"); dir.create("results/figures", showWarnings=FALSE)
th <- theme_bw(base_size=10) + theme(panel.grid.minor=element_blank(), legend.position="top")

## Mouse: all samples (including ambiguous)
m <- fread("data/processed/mouse_sample_table_sexinferred.tsv.gz", na.strings=c("","NA"))
m <- m[!is.na(xist_l) & !is.na(y_l)][, d := xist_l - y_l]
m[, call := fifelse(d >= 2, "F", fifelse(d <= -2, "M", "ambiguous"))]
## Human: analysis-set pheno (ambiguous already excluded)
hs <- readLines("config/c2_analysis_set_human.txt")
h <- rbindlist(lapply(hs, function(s) { f <- sprintf("results/main_human_v2/%s/pheno.tsv", s); if (file.exists(f)) fread(f)[, .(study = s, xist, y_score, sex)] else NULL }))
h[, d := xist - y_score]

cat("== Mouse, all samples", nrow(m), "==\n"); print(m[, .N, by=call])
cat("Threshold sensitivity (mouse): ambiguous fraction\n")
for (t in c(1, 1.5, 2, 2.5, 3)) cat(sprintf("  |d| < %.1f : %5.2f%%   (F %.1f%% / M %.1f%%)\n", t, 100*mean(abs(m$d) < t), 100*mean(m$d >= t), 100*mean(m$d <= -t)))
cat("\n== Human analysis-set samples", nrow(h), "(after excluding ambiguous) ==\n")
cat("Fraction that becomes newly ambiguous when the threshold is raised:\n")
for (t in c(2, 2.5, 3, 4)) cat(sprintf("  |d| < %.1f : %5.2f%%\n", t, 100*mean(abs(h$d) < t)))

## (a) Mouse scatter plot
pa <- ggplot(m[sample(.N, min(.N, 20000))], aes(y_l, xist_l, colour=call)) + geom_point(size=0.35, alpha=0.4) +
  geom_abline(intercept=c(-2, 2), slope=1, linetype=2, colour="grey30") +
  scale_colour_manual(values=c(F="#c0392b", M="#2E6F6C", ambiguous="#e67e22")) +
  labs(x="mean log2 CPM, Y-linked genes", y="log2 CPM, Xist", colour=NULL, title="Mouse: all samples (n=26,839)") + th
## (b) Human scatter plot
pb <- ggplot(h, aes(y_score, xist, colour=sex)) + geom_point(size=0.35, alpha=0.4) +
  geom_abline(intercept=c(-2, 2), slope=1, linetype=2, colour="grey30") +
  scale_colour_manual(values=c(F="#c0392b", M="#2E6F6C")) +
  labs(x="mean log2 CPM, Y-linked genes", y="log2 CPM, XIST", colour=NULL, title="Human: analysis-set samples (ambiguous already excluded)") + th
## (c) Distribution of the difference
dd <- rbind(m[, .(d, organism="mouse")], h[, .(d, organism="human")])
pc <- ggplot(dd, aes(d)) + geom_histogram(binwidth=0.25, fill="grey55", colour=NA) +
  geom_vline(xintercept=c(-2, 2), linetype=2, colour="#c0392b") + facet_wrap(~organism, scales="free_y") +
  labs(x="XIST − mean(Y)  [log2 CPM]", y="Samples", title="Separation between the two calls; dashed = ±2 threshold") + th
## (d) Ambiguous fraction per study
amb <- rbindlist(lapply(c("human","mouse"), function(o) { d <- fread(sprintf("results/main_%s_summary_v2.tsv", o))[valid == TRUE]
  d[, .(organism = o, frac = n_sex_ambiguous / pmax(n + n_sex_ambiguous, 1))] }))
pd <- ggplot(amb, aes(frac)) + geom_histogram(binwidth=0.02, fill="grey55", colour=NA) + facet_wrap(~organism, scales="free_y") +
  labs(x="Fraction of a study's samples called ambiguous", y="Studies", title="Per-study ambiguous rate") + th
cat("\nStudies with a per-study ambiguous fraction >10%:\n"); print(amb[, .(n = .N, frac_gt_10pct = round(100*mean(frac > 0.1), 1), frac_gt_25pct = round(100*mean(frac > 0.25), 1)), by=organism])

for (nm in c("a","b","c","d")) { p <- get(paste0("p", nm))
  ggsave(sprintf("results/figures/a5_sex_threshold_%s.png", nm), p, width=6.5, height=4, dpi=170) }
fwrite(dd, "results/qc/a5_sex_delta.tsv.gz", sep="\t")
cat("\nSaved: results/figures/a5_sex_threshold_{a,b,c,d}.png\n")
