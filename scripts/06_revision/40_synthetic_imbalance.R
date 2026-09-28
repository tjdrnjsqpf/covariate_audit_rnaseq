#!/usr/bin/env Rscript
## S (2026-09-18, review A priority 1 / C 4.2): within-study synthetic imbalance experiment.
## From one balanced-design study, draw subsamples (fixed size) that vary only the per-group sex composition, and at |Δf| ∈ {0,.2,.4,.6,.8}
## refit M0 / M1 / M0perm (full permutation) / M0permw (within-group permutation). Because study, tissue, effect size and n are fixed,
## the dose–response is identified within the study, and with the full-sample M1 fit as reference we measure "error" rather than "change".
## v2 (2026-09-22, review H/I): in addition to reference A (full-sample M1), add error vs reference C (M1 on the remaining samples excluding the subsample; out-of-sample), overlap vs a balanced draw of the same size (reference B),
##   and bias²/variance decomposition of the logFC error on reference A's DE genes (M0 vs M1). Usage: Rscript 40_synthetic_imbalance.R <human|mouse> <SRP> [reps=10] [outdir=results/revision/synthetic]
suppressMessages({library(recount3); library(edgeR); library(limma); library(data.table)})
args <- commandArgs(trailingOnly=TRUE); ORG <- args[1]; SRP <- args[2]; REPS <- if (length(args) >= 3) as.integer(args[3]) else 10L
OUTDIR <- if (length(args) >= 4) args[4] else "results/revision/synthetic"
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ); dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
IS_MOUSE <- ORG == "mouse"
panel <- if (IS_MOUSE) data.table(gene=c("Xist","Ddx3y","Eif2s3y","Kdm5d","Uty")) else fread("config/sex_linked_genes.tsv")
ph <- fread(sprintf("results/main_%s_v2/%s/pheno.tsv", ORG, SRP))[, .(external_id, group, sex)]
ph[, group := factor(group, levels=c("control","case"))]

## counts (same run summation as the DE script)
bfc <- recount3_cache(file.path(PROJ, "data/raw/recount3_cache", SRP))
rse <- create_rse_manual(project=SRP, project_home="data_sources/sra", organism=if (IS_MOUSE) "mouse" else "human", annotation=if (IS_MOUSE) "gencode_v23" else "gencode_v29", type="gene", bfc=bfc)
assay(rse, "counts") <- transform_counts(rse); cts <- assay(rse, "counts"); sym <- rowData(rse)$gene_name
sacc_col <- grep("^sra\\.sample_acc", names(colData(rse)), value=TRUE)[1]
sacc <- data.table(external_id = colData(rse)$external_id, sample_acc = as.character(colData(rse)[[sacc_col]])); sacc[is.na(sample_acc) | sample_acc == "", sample_acc := external_id]
if (anyDuplicated(sacc$sample_acc)) { sacc[, sample_acc := factor(sample_acc, levels=unique(sample_acc))]; M <- model.matrix(~ 0 + sample_acc, sacc); rownames(M) <- sacc$external_id
  rep_id <- sacc[, external_id[1], by=sample_acc]$V1; cts <- as.matrix(cts[, rownames(M)]) %*% M; colnames(cts) <- rep_id }
stopifnot(all(ph$external_id %in% colnames(cts))); cts <- cts[, ph$external_id]

## Cell sizes and subsample size m (identical across all steps)
cell <- ph[, .N, by=.(group, sex)]; gc <- function(g, s) { v <- cell[group == g & sex == s, N]; if (length(v)) v else 0L }
Fc <- gc("case","F"); Mc <- gc("case","M"); Fk <- gc("control","F"); Mk <- gc("control","M")
m <- min(Fc, Mc, Fk, Mk, 25L); if (m < 6) stop("min cell < 6: ", paste(Fc, Mc, Fk, Mk))
KEEP <- filterByExpr(DGEList(cts), model.matrix(~group, ph))   # common universe defined on the full study
fit <- function(d, design) { y <- calcNormFactors(DGEList(cts[KEEP, d$external_id])); v <- voom(y, design); f <- eBayes(lmFit(v, design))
  tt <- as.data.table(topTable(f, coef="groupcase", number=Inf, sort.by="none")); tt$gene_id <- rownames(cts)[KEEP]; tt$gene <- sym[KEEP]; tt }
## Reference: full sample, sex-adjusted (|Δf|≈0, maximum power)
REF <- fit(ph, model.matrix(~group + sex, ph)); ref_de <- REF[adj.P.Val < 0.05, gene_id]; ref_t <- setNames(REF$t, REF$gene_id)
jac <- function(a, b) { u <- length(union(a, b)); if (u == 0) NA_real_ else length(intersect(a, b)) / u }
cmp_ref <- function(tt, tag) { de <- tt[adj.P.Val < 0.05, gene_id]; top <- head(tt[order(adj.P.Val, -abs(t))], 100)$gene_id; reftop <- head(REF[order(adj.P.Val, -abs(t))], 100)$gene_id
  setNames(list(length(de), if (length(de)) length(intersect(de, ref_de)) / length(de) else NA_real_, if (length(ref_de)) length(intersect(de, ref_de)) / length(ref_de) else NA_real_,
                round(cor(tt$t, ref_t[tt$gene_id], method="spearman"), 4), length(intersect(top, reftop)) / 100, sum(tt[gene_id %in% de, gene] %in% panel$gene)),
           paste0(tag, c("_n_de","_precision","_recall","_t_sp_ref","_top100_ref","_sexgene_de"))) }
ref_share <- length(ref_de) / sum(KEEP)
out <- list(); set.seed(sum(utf8ToInt(SRP))); BAL <- list(); LFC <- list()
for (d in c(0, 0.2, 0.4, 0.6, 0.8)) { kc <- round(m * (0.5 + d/2)); kf <- round(m * (0.5 - d/2))
  for (r in seq_len(REPS)) {
    set.seed(sum(utf8ToInt(SRP)) + 1000L*which(c(0,.2,.4,.6,.8) == d) + r)
    pick <- function(g, s, k) if (k > 0) sample(ph[group == g & sex == s, external_id], k) else character(0)
    ids <- c(pick("case","F",kc), pick("case","M",m-kc), pick("control","F",kf), pick("control","M",m-kf))
    s <- ph[match(ids, external_id)]; s[, group := factor(group, levels=c("control","case"))]
    s[, sex_perm := sample(sex)]; s[, sex_permw := if (.N > 1) sample(sex) else sex, by=group]
    mm <- function(f) model.matrix(as.formula(paste0("~group", f)), s)
    M0 <- fit(s, mm("")); M1 <- fit(s, mm(" + sex"))
    Xp <- mm(" + sex_perm"); M0perm <- if (qr(Xp)$rank == ncol(Xp)) fit(s, Xp) else NULL
    Xw <- mm(" + sex_permw"); M0permw <- if (qr(Xw)$rank == ncol(Xw)) fit(s, Xw) else NULL
    de0 <- M0[adj.P.Val < 0.05, gene_id]; de1 <- M1[adj.P.Val < 0.05, gene_id]
    ## Reference C: M1 on the remaining samples excluding the subsample (out-of-sample reference; only when ≥4 per group and ≥8 in total)
    ## v3 (2026-09-23, review J B1): the complement set is imbalanced in the opposite direction to the draw, so fit M1 on a 'balanced complement subset' with equal sex counts within each group and use it as reference C (≥2 of each sex per group, ≥8 in total)
    comp <- ph[!external_id %in% ids]; comp[, group := factor(group, levels=c("control","case"))]; refC <- NULL
    df_comp <- if (nrow(comp) && all(table(comp$group) > 0)) abs(diff(comp[, mean(sex == "F"), by=group][order(group)]$V1)) else NA_real_
    set.seed(sum(utf8ToInt(SRP)) + 7000L*which(c(0,.2,.4,.6,.8) == d) + r); cb <- comp[, { k <- min(sum(sex == "F"), sum(sex == "M")); .SD[c(sample(which(sex == "F"), k), sample(which(sex == "M"), k))] }, by=group]
    if (nrow(cb) >= 8 && all(table(cb$group) >= 4) && uniqueN(cb$sex) == 2) { Xc <- model.matrix(~group + sex, cb); if (qr(Xc)$rank == ncol(Xc)) refC <- fit(cb, Xc) }
    cmpC <- function(tt, tag) { if (is.null(refC)) return(setNames(list(NA_real_, NA_real_, NA_real_, NA_integer_, NA_real_), paste0(tag, c("_precision_C","_recall_C","_t_sp_C","_n_refC_de","_top100_C"))))
      deC <- refC[adj.P.Val < 0.05, gene_id]; de <- tt[adj.P.Val < 0.05, gene_id]; rC <- setNames(refC$t, refC$gene_id)
      top <- function(d) d[order(adj.P.Val, -abs(t))][1:min(100, .N), gene_id]
      setNames(list(if (length(de)) length(intersect(de, deC)) / length(de) else NA_real_, if (length(deC)) length(intersect(de, deC)) / length(deC) else NA_real_, round(cor(tt$t, rC[tt$gene_id], method="spearman"), 4), length(deC), length(intersect(top(tt), top(refC))) / 100), paste0(tag, c("_precision_C","_recall_C","_t_sp_C","_n_refC_de","_top100_C"))) }
    ## Reference B: overlap with the DE set of a balanced draw of the same size (step 0, same r)
    if (d == 0) BAL[[r]] <- list(de0=de0, de1=de1)
    JB0 <- if (d == 0) NA_real_ else jac(de0, BAL[[r]]$de1); JB1 <- if (d == 0) NA_real_ else jac(de1, BAL[[r]]$de1)
    ## Save logFC (for reference A's DE genes; for the bias/variance decomposition)
    key <- as.character(d); if (is.null(LFC[[key]])) LFC[[key]] <- list(M0=list(), M1=list())
    LFC[[key]]$M0[[r]] <- M0$logFC[match(ref_de, M0$gene_id)]; LFC[[key]]$M1[[r]] <- M1$logFC[match(ref_de, M1$gene_id)]
    row <- c(list(study=SRP, organism=ORG, m_per_arm=m, target_df=d, realized_df=abs(kc - kf)/m, rep=r, n_ref_de=length(ref_de), ref_share=ref_share, n_comp=nrow(comp), df_comp=df_comp, n_compbal=nrow(cb),
                  J_M1=jac(de0, de1), J_M0perm=if (is.null(M0perm)) NA_real_ else jac(de0, M0perm[adj.P.Val < 0.05, gene_id]),
                  J_M0permw=if (is.null(M0permw)) NA_real_ else jac(de0, M0permw[adj.P.Val < 0.05, gene_id]),
                  t_sp_M0_M1=round(cor(M0$t, M1$t, method="spearman"), 4), J_M0_vs_balanced=JB0, J_M1_vs_balanced=JB1),
             cmp_ref(M0, "M0"), cmp_ref(M1, "M1"), cmpC(M0, "M0"), cmpC(M1, "M1"))
    out[[length(out)+1]] <- as.data.table(row) } }
## Bias²/variance decomposition: for each of reference A's DE genes, squared deviation of the across-draw mean from the reference (bias²) and the across-draw variance
lref <- REF$logFC[match(ref_de, REF$gene_id)]; BV <- list()
for (key in names(LFC)) for (mod in c("M0","M1")) { M <- do.call(cbind, LFC[[key]][[mod]]); mu <- rowMeans(M, na.rm=TRUE); v <- apply(M, 1, var, na.rm=TRUE)
  BV[[length(BV)+1]] <- data.table(study=SRP, organism=ORG, target_df=as.numeric(key), model=mod, n_ref_de=length(ref_de), bias2=mean((mu - lref)^2, na.rm=TRUE), variance=mean(v, na.rm=TRUE), mse=mean((mu - lref)^2, na.rm=TRUE) + mean(v, na.rm=TRUE), signed_bias=mean(mu - lref, na.rm=TRUE), rmse_per_draw=mean(sqrt(colMeans((M - lref)^2, na.rm=TRUE)))) }
fwrite(rbindlist(BV), file.path(OUTDIR, paste0(SRP, "_biasvar.tsv")), sep="\t")
res <- rbindlist(out); fwrite(res, file.path(OUTDIR, paste0(SRP, ".tsv")), sep="\t")
cat(SRP, ": m =", m, " ref DE =", length(ref_de), " rows =", nrow(res), "\n")
