#!/usr/bin/env Rscript
# Semi-synthetic simulation: take a balanced large study as ground truth (full-data M2), induce sex/age imbalance by subsampling, then
# compare precision/recall/rank concordance of M0(~group), M1(+sex), M2(+sex+age).
# Usage: Rscript 06_semisynthetic_sim.R <SRP> [reps]
suppressMessages({library(recount3); library(edgeR); library(limma); library(data.table)})
args <- commandArgs(trailingOnly=TRUE); SRP <- args[1]; REPS <- if (length(args) >= 2) as.integer(args[2]) else 15
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
OUT <- file.path(PROJ, "results/simulation", SRP); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
panel <- fread(file.path(PROJ, "config/sex_linked_genes.tsv"))
set.seed(20260906)

## 1. Data: pilot pheno (finalised group/sex/age) + cached counts
ph <- fread(file.path(PROJ, "results/pilot", SRP, "pheno.tsv"))[!is.na(age) & sex %in% c("F","M")]
ph[, group := factor(group, levels=c("control","case"))]
bfc <- recount3_cache(file.path(PROJ, "data/raw/recount3_cache", SRP))
rse <- create_rse_manual(project=SRP, project_home="data_sources/sra", organism="human", annotation="gencode_v29", type="gene", bfc=bfc)
cts <- transform_counts(rse)[, ph$external_id]; sym <- rowData(rse)$gene_name
keep <- filterByExpr(DGEList(cts), model.matrix(~group, ph)); cts <- cts[keep, ]; sym <- sym[keep]
cat(SRP, ": samples", nrow(ph), " genes", nrow(cts), "\n"); print(table(ph$group, ph$sex))

fit_de <- function(idx, f) {
  d <- ph[idx]; X <- model.matrix(as.formula(paste0("~group", f)), d)
  if (qr(X)$rank < ncol(X)) return(NULL)   # fully confounded (e.g. all cases female, all controls male) -> cannot fit
  y <- calcNormFactors(DGEList(cts[, d$external_id])); v <- voom(y, X); fit <- eBayes(lmFit(v, X))
  tt <- topTable(fit, coef="groupcase", number=Inf, sort.by="none"); data.table(gene_id=rownames(tt), gene=sym, logFC=tt$logFC, P=tt$P.Value, FDR=tt$adj.P.Val)
}
## 2. Ground truth: full-data M2. Studies with thousands of DEGs are tightened to FDR<0.01 & |logFC|>0.5
full <- fit_de(seq_len(nrow(ph)), " + sex + age")
## Two truth sets: a lenient set for precision (full M2 FDR<0.05), a strict set for recall (FDR<0.01 & |logFC|>0.5)
truth_loose <- full[FDR < 0.05, gene_id]; truth <- full[FDR < 0.01 & abs(logFC) > 0.5, gene_id]; truth_top100 <- full[order(P)][1:100, gene_id]
full0 <- fit_de(seq_len(nrow(ph)), "")
cat("truth strict:", length(truth), " loose:", length(truth_loose), " (full M0 FDR<0.05:", sum(full0$FDR < 0.05), ")\n")
fwrite(full, file.path(OUT, "full_M2.tsv.gz"), sep="\t")

## 3. Imbalance-inducing subsampling
## sex_delta: difference in female proportion (case p=0.5+δ/2, control p=0.5-δ/2).
## age_shift: 0=random, 0.5=moderate (cases above the 25th percentile, controls below the 75th -> SMD≈1), 1=strong (cases from the upper half, controls from the lower half -> SMD≈2.5)
draw <- function(n_per, sex_delta, age_shift) {
  q <- quantile(ph$age, c(0.25, 0.5, 0.75))
  pick <- function(g, pf) {
    d <- ph[group == g]
    if (age_shift == 1) d <- if (g == "case") d[age >= q[2]] else d[age < q[2]]
    if (age_shift == 0.5) d <- if (g == "case") d[age >= q[1]] else d[age <= q[3]]
    nf <- round(n_per * pf); nm <- n_per - nf
    f <- d[sex == "F"]; m <- d[sex == "M"]
    if (nrow(f) < nf || nrow(m) < nm) return(NULL)
    c(sample(f$external_id, nf), sample(m$external_id, nm))
  }
  ids <- c(pick("case", 0.5 + sex_delta/2), pick("control", 0.5 - sex_delta/2))
  if (is.null(ids) || length(ids) < 2 * n_per) return(NULL); which(ph$external_id %in% ids)
}
metrics <- function(r, label) {
  if (is.null(r)) return(list(n_deg = NA_integer_, precision = NA_real_, precision_strict = NA_real_, recall = NA_real_, top100_vs_truth = NA_real_, sexgene_deg = NA_integer_))
  deg <- r[FDR < 0.05, gene_id]; top <- r[order(P)][1:100, gene_id]
  list(n_deg = length(deg), precision = if (length(deg)) mean(deg %in% truth_loose) else NA_real_, precision_strict = if (length(deg)) mean(deg %in% truth) else NA_real_, recall = mean(truth %in% deg),
       top100_vs_truth = mean(top %in% truth_top100), sexgene_deg = sum(r[gene_id %in% deg, gene] %in% panel$gene))
}
grid <- CJ(n_per = c(15L, 30L), sex_delta = c(0, 0.25, 0.5, 0.75, 1), age_shift = c(0, 0.5, 1), rep = seq_len(REPS))
res <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i]; idx <- draw(g$n_per, g$sex_delta, g$age_shift); if (is.null(idx)) next
  d <- ph[idx]; imb <- d[, .(pf = mean(sex == "F"), ma = mean(age), sa = sd(age)), by=group]
  obs_sex_diff <- abs(diff(imb$pf)); obs_smd <- abs(diff(imb$ma)) / sqrt(mean(imb$sa^2))
  r0 <- fit_de(idx, ""); r1 <- fit_de(idx, " + sex"); r2 <- fit_de(idx, " + sex + age")
  m0 <- metrics(r0); m1 <- metrics(r1); m2 <- metrics(r2)
  chg <- if (is.null(r2)) list(obs_jaccard_M0_M2 = NA_real_, obs_top100_M0_M2 = NA_real_, obs_tspearman_M0_M2 = NA_real_) else {
    d0 <- r0[FDR < 0.05, gene_id]; d2 <- r2[FDR < 0.05, gene_id]
    list(obs_jaccard_M0_M2 = length(intersect(d0, d2)) / max(1, length(union(d0, d2))),
         obs_top100_M0_M2 = length(intersect(r0[order(P)][1:100, gene_id], r2[order(P)][1:100, gene_id])) / 100,
         obs_tspearman_M0_M2 = cor(sign(r0$logFC) * -log10(r0$P), sign(r2$logFC) * -log10(r2$P), method="spearman")) }
  res[[i]] <- data.table(study = SRP, g, obs_sex_diff = round(obs_sex_diff, 2), obs_age_smd = round(obs_smd, 2),
                         as.data.table(setNames(c(m0, m1, m2), c(paste0("M0_", names(m0)), paste0("M1_", names(m1)), paste0("M2_", names(m2))))), as.data.table(chg))
  if (i %% 25 == 0) cat(i, "/", nrow(grid), "\n")
}
out <- rbindlist(res, fill=TRUE); out[, n_truth := length(truth)][, n_truth_loose := length(truth_loose)]
fwrite(out, file.path(OUT, "sim_results.tsv"), sep="\t")
cat("done:", nrow(out), "fits\n")
print(out[, .(reps = .N, M0_prec = round(mean(M0_precision, na.rm=TRUE), 2), M2_prec = round(mean(M2_precision, na.rm=TRUE), 2),
              M0_rec = round(mean(M0_recall), 3), M2_rec = round(mean(M2_recall), 3), M0_sexg = round(mean(M0_sexgene_deg), 1),
              obs_jacc = round(mean(obs_jaccard_M0_M2), 2)), by=.(n_per, sex_delta, age_shift)])
