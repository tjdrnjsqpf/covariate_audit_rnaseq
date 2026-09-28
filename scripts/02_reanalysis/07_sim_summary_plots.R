#!/usr/bin/env Rscript
# Collect simulation results + figures (base R). Usage: Rscript 07_sim_summary_plots.R
suppressMessages(library(data.table))
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; SIM <- file.path(PROJ, Sys.getenv("SIMDIR", "results/simulation"))
fs <- list.files(SIM, pattern="sim_results.tsv", recursive=TRUE, full.names=TRUE)
d <- rbindlist(lapply(fs, fread)); fwrite(d, file.path(SIM, "sim_results_all.tsv"), sep="\t")
cat("fits:", nrow(d), "\n"); print(d[, .N, by=.(study, n_per, age_shift)])
## Mean summary
sm <- d[, .(reps = .N, M0_precision = mean(M0_precision, na.rm=TRUE), M1_precision = mean(M1_precision, na.rm=TRUE), M2_precision = mean(M2_precision, na.rm=TRUE),
            M0_precision_strict = mean(M0_precision_strict, na.rm=TRUE), M2_precision_strict = mean(M2_precision_strict, na.rm=TRUE),
            M0_recall = mean(M0_recall, na.rm=TRUE), M1_recall = mean(M1_recall, na.rm=TRUE), M2_recall = mean(M2_recall, na.rm=TRUE),
            M0_top100 = mean(M0_top100_vs_truth, na.rm=TRUE), M1_top100 = mean(M1_top100_vs_truth, na.rm=TRUE), M2_top100 = mean(M2_top100_vs_truth, na.rm=TRUE),
            M0_ndeg = mean(M0_n_deg, na.rm=TRUE), M2_ndeg = mean(M2_n_deg, na.rm=TRUE), M0_sexgene = mean(M0_sexgene_deg, na.rm=TRUE), M2_sexgene = mean(M2_sexgene_deg, na.rm=TRUE),
            M0_fdp = 1 - mean(M0_precision, na.rm=TRUE), M2_fdp = 1 - mean(M2_precision, na.rm=TRUE),
            obs_jaccard = mean(obs_jaccard_M0_M2, na.rm=TRUE), obs_top100 = mean(obs_top100_M0_M2, na.rm=TRUE), obs_sex_diff = mean(obs_sex_diff), obs_age_smd = mean(obs_age_smd)),
        by=.(study, n_per, sex_delta, age_shift)][order(study, n_per, age_shift, sex_delta)]
fwrite(sm, file.path(SIM, "sim_summary.tsv"), sep="\t")
## Figure 1: precision / recall vs sex_delta, one line per model, panels = study × n_per × age_shift
cols <- c(M0="#d62728", M1="#ff7f0e", M2="#1f77b4")
studies <- unique(sm$study)
for (metric in c("precision", "recall", "top100")) {
  NPS <- sort(unique(sm$n_per)); ASS <- sort(unique(sm$age_shift)); ncol_p <- length(NPS) * length(ASS)
  png(file.path(SIM, paste0("fig_sim_", metric, ".png")), width=450 * ncol_p, height=450 * length(studies), res=150)
  par(mfrow=c(length(studies), ncol_p), mar=c(4, 4, 3, 1))
  for (st in studies) for (np in NPS) for (as in ASS) {
    x <- sm[study == st & n_per == np & age_shift == as]
    plot(NA, xlim=c(0, 1), ylim=c(0, 1), xlab="sex imbalance (Δ female fraction)", ylab=metric,
         main=sprintf("%s  n/group=%d  age shift=%s", st, np, c("0"="none", "0.5"="moderate", "1"="strong")[as.character(as)]), cex.main=0.9)
    for (m in c("M0", "M1", "M2")) { col <- paste0(m, "_", metric); if (!col %in% names(x)) next
      lines(x$sex_delta, x[[col]], col=cols[m], lwd=2, type="b", pch=16) }
    legend("bottomleft", legend=c("M0 ~group", "M1 +sex", "M2 +sex+age"), col=cols, lwd=2, bty="n", cex=0.8)
  }
  dev.off()
}
## Figure 2: observable change metric (M0 vs M2 Jaccard) vs true error (M0 precision) — how well does the observable metric predict the error?
png(file.path(SIM, "fig_sim_obs_vs_truth.png"), width=500 * length(studies), height=500, res=150); par(mfrow=c(1, length(studies)), mar=c(4, 4, 3, 1))
for (st in studies) { x <- d[study == st]
  plot(x$obs_jaccard_M0_M2, x$M0_precision, col=c("0"="#1f77b4", "0.5"="#ff7f0e", "1"="#d62728")[as.character(x$age_shift)], pch=ifelse(x$n_per == 30, 16, 1),
       xlab="observed Jaccard (M0 vs M2 DEG)", ylab="M0 precision vs truth", main=st, xlim=c(0, 1), ylim=c(0, 1))
  legend("bottomright", legend=c("age none", "age moderate", "age strong", "n=30", "n=15"), col=c("#1f77b4", "#ff7f0e", "#d62728", "black", "black"), pch=c(15, 15, 15, 16, 1), bty="n", cex=0.8)
  cat(st, " cor(obs_jaccard, M0_precision) =", round(cor(x$obs_jaccard_M0_M2, x$M0_precision, use="complete.obs", method="spearman"), 3), "\n") }
dev.off()
print(sm[, .(study, n_per, sex_delta, age_shift, M0_precision = round(M0_precision, 2), M2_precision = round(M2_precision, 2), M0_recall = round(M0_recall, 3), M2_recall = round(M2_recall, 3), M0_sexgene = round(M0_sexgene, 1), obs_jaccard = round(obs_jaccard, 2))])
