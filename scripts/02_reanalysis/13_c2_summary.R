#!/usr/bin/env Rscript
## C2 aggregation: report pathway-level change (Hallmark/Reactome) by design grade, as Δ relative to the negative-control baseline.
## Usage: Rscript 13_c2_summary.R <human|mouse>
suppressMessages(library(data.table))
ORG <- commandArgs(trailingOnly=TRUE)[1]
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
g <- fread(sprintf("results/c2/gsea_summary_%s.tsv", ORG))
m <- fread(sprintf("results/main_%s_summary_v2.tsv", ORG))[, .(study, design_sex, design_age, sex_diff, age_smd, n, M1_jaccard, M0perm_jaccard, M1_n_deg_ref, M1_t_spearman, M2_jaccard, M2perm_jaccard)]
d <- merge(g, m, by="study")
lv <- c("balanced","mildly_imbalanced","partially_confounded","completely_confounded")
boot_ci <- function(x, R=2000) { x <- x[!is.na(x)]; if (length(x) < 3) return(c(NA,NA)); set.seed(1); b <- replicate(R, median(sample(x, length(x), TRUE))); round(quantile(b, c(.025,.975)), 3) }
mc <- function(x) { ci <- boot_ci(x); sprintf("%.3f [%.3f-%.3f]", median(x, na.rm=TRUE), ci[1], ci[2]) }
cat(sprintf("########## C2 pathway-level aggregation: %s (n=%d) ##########\n", ORG, nrow(d)))

for (col in c("H","R")) {
  nm <- if (col == "H") "Hallmark(50)" else "Reactome(1839)"
  cat(sprintf("\n== [%s] by sex grade: M0 vs M1 (real) / M0 vs M0perm (permuted) — Jaccard of significant pathways (padj<.05) ==\n", nm))
  j1 <- paste0("P1_", col, "_jaccard"); jp <- paste0("P0perm_", col, "_jaccard"); ns <- paste0("P1_", col, "_n_sig_ref"); na <- paste0("P1_", col, "_n_sig_alt")
  t1 <- paste0("P1_", col, "_top20_overlap"); tp <- paste0("P0perm_", col, "_top20_overlap"); s1 <- paste0("P1_", col, "_nes_spearman"); sp <- paste0("P0perm_", col, "_nes_spearman")
  f1 <- paste0("P1_", col, "_sign_flip_sig"); o1 <- paste0("P1_", col, "_top1_same"); op <- paste0("P0perm_", col, "_top1_same")
  x <- d[!is.na(get(j1)) & !is.na(get(jp))]
  print(x[, .(n = .N, sig_M0 = as.numeric(median(get(ns))), sig_M1 = as.numeric(median(get(na))),
              J_real = mc(get(j1)), J_perm = mc(get(jp)), dJ = round(median(get(j1) - get(jp)), 3),
              top20_real = round(median(get(t1)), 2), top20_perm = round(median(get(tp)), 2),
              NESsp_real = round(median(get(s1)), 3), NESsp_perm = round(median(get(sp)), 3),
              signflip = round(mean(get(f1) > 0), 2), top1_same_real = round(mean(get(o1), na.rm=TRUE), 2), top1_same_perm = round(mean(get(op), na.rm=TRUE), 2)),
          by=design_sex][order(match(design_sex, lv))])
}

cat("\n== Gene level vs pathway level (Hallmark): Δ comparison under partial confounding ==\n")
print(d[!is.na(M1_jaccard) & !is.na(P1_H_jaccard) & !is.na(P0perm_H_jaccard),
        .(n = .N, gene_dJ = round(median(M1_jaccard - M0perm_jaccard), 3), path_dJ = round(median(P1_H_jaccard - P0perm_H_jaccard), 3),
          gene_J_real = round(median(M1_jaccard), 3), path_J_real = round(median(P1_H_jaccard), 3)), by=design_sex][order(match(design_sex, lv))])

cat("\n== [Hallmark] by age grade: ref vs M2 (real) / ref vs M2perm (permuted) ==\n")
x <- d[!is.na(P2_H_jaccard) & !is.na(P2perm_H_jaccard)]
if (nrow(x)) print(x[, .(n = .N, J_real = mc(P2_H_jaccard), J_perm = mc(P2perm_H_jaccard), dJ = round(median(P2_H_jaccard - P2perm_H_jaccard), 3),
                         NESsp_real = round(median(P2_H_nes_spearman), 3), NESsp_perm = round(median(P2perm_H_nes_spearman), 3)), by=design_age][order(match(design_age, lv))])

cat("\n== C2-b: nature of DEGs lost from M0 to M1 (real vs permuted) ==\n")
print(d[, .(n = .N, n_deg0 = as.numeric(median(G1_n_deg0)), lost_real = as.numeric(median(G1_n_lost)), lost_perm = as.numeric(median(G0perm_n_lost, na.rm=TRUE)),
            lost_frac_real = round(median(G1_n_lost / pmax(G1_n_deg0, 1)), 3), lost_frac_perm = round(median(G0perm_n_lost / pmax(G0perm_n_deg0, 1), na.rm=TRUE), 3),
            boundary_real = round(median(G1_lost_frac_boundary, na.rm=TRUE), 2), boundary_perm = round(median(G0perm_lost_frac_boundary, na.rm=TRUE), 2),
            bigLFC_lost_real = round(median(G1_lost_frac_big_lfc, na.rm=TRUE), 2),
            lfc_lost = round(median(G1_lost_med_abs_lfc, na.rm=TRUE), 2), lfc_kept = round(median(G1_kept_med_abs_lfc, na.rm=TRUE), 2)),
        by=design_sex][order(match(design_sex, lv))])

cat("\n== Number of partially confounded studies whose 'pathway conclusion actually changed' (Hallmark: real J < permuted J − 0.2 & top1 changed) ==\n")
x <- d[design_sex == "partially_confounded" & !is.na(P1_H_jaccard) & !is.na(P0perm_H_jaccard)]
x[, changed := (P1_H_jaccard < P0perm_H_jaccard - 0.2)]
cat(sprintf("n=%d, changed by the J criterion %d (%.0f%%), of which the top1 pathway also changed %d\n", nrow(x), sum(x$changed), 100*mean(x$changed), sum(x$changed & !x$P1_H_top1_same, na.rm=TRUE)))
print(x[changed == TRUE][order(P1_H_jaccard)][1:min(10,.N), .(study, n, sex_diff, M1_jaccard, P1_H_n_sig_ref, P1_H_n_sig_alt, P1_H_jaccard, P0perm_H_jaccard, P1_H_nes_spearman, P1_H_top1_same)])
fwrite(d, sprintf("results/c2/c2_merged_%s.tsv", ORG), sep="\t")
