#!/usr/bin/env Rscript
## Engine robustness check (advisor comment 4, 2026-10-10): refit M0 / M1 / M0perm / M0permw with limma-voom, edgeR quasi-likelihood and DESeq2
## on the stored phenotype (incl. the stored permuted labels sex_perm / sex_permw) and the shared gene universe (filterByExpr on ~group).
## Usage: Rscript 54_engine_check.R <human|mouse> <SRP>  -> results/response/engine/<SRP>.tsv
suppressMessages({library(recount3); library(edgeR); library(limma); library(DESeq2); library(data.table)})
a <- commandArgs(trailingOnly=TRUE); ORG <- a[1]; SRP <- a[2]; PROJ <- Sys.getenv("COV_AUDIT_PROJ", "/path/to/covariate_audit_rnaseq"); setwd(PROJ)
OUT <- "results/response/engine"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE); IS_MOUSE <- ORG == "mouse"
D <- sprintf("results/main_%s_v2/%s", ORG, SRP); ph <- fread(file.path(D, "pheno.tsv")); ph[, group := factor(group, levels=c("control","case"))]
stopifnot(all(c("sex","sex_perm","sex_permw") %in% names(ph)))
bfc <- recount3_cache(file.path(PROJ, "data/raw/recount3_cache", SRP))
rse <- create_rse_manual(project=SRP, project_home="data_sources/sra", organism=if (IS_MOUSE) "mouse" else "human", annotation=if (IS_MOUSE) "gencode_v23" else "gencode_v29", type="gene", bfc=bfc)
assay(rse, "counts") <- transform_counts(rse); cts <- assay(rse, "counts")
sacc_col <- grep("^sra\\.sample_acc", names(colData(rse)), value=TRUE)[1]; sacc <- data.table(external_id=colData(rse)$external_id, sample_acc=as.character(colData(rse)[[sacc_col]])); sacc[is.na(sample_acc) | sample_acc == "", sample_acc := external_id]
if (anyDuplicated(sacc$sample_acc)) { sacc[, sample_acc := factor(sample_acc, levels=unique(sample_acc))]; M <- model.matrix(~ 0 + sample_acc, sacc); rownames(M) <- sacc$external_id; rep_id <- sacc[, external_id[1], by=sample_acc]$V1; cts <- as.matrix(cts) %*% M; colnames(cts) <- rep_id }
stopifnot(all(ph$external_id %in% colnames(cts))); cts <- round(as.matrix(cts[, ph$external_id]))
KEEP <- filterByExpr(DGEList(cts), model.matrix(~group, ph)); cts <- cts[KEEP, ]
designs <- list(M0=~group, M1=~group + sex, M0perm=~group + sex_perm, M0permw=~group + sex_permw)
ok <- sapply(designs, function(f) { X <- model.matrix(f, ph); qr(X)$rank == ncol(X) })
fit_limma <- function(f) { X <- model.matrix(f, ph); y <- calcNormFactors(DGEList(cts)); v <- voom(y, X); tt <- topTable(eBayes(lmFit(v, X)), coef="groupcase", number=Inf, sort.by="none"); rownames(tt)[tt$adj.P.Val < 0.05] }
fit_edger <- function(f) { X <- model.matrix(f, ph); y <- estimateDisp(calcNormFactors(DGEList(cts)), X); q <- glmQLFit(y, X); tt <- topTags(glmQLFTest(q, coef="groupcase"), n=Inf, sort.by="none")$table; rownames(tt)[tt$FDR < 0.05] }
fit_deseq <- function(f) { cd <- as.data.frame(ph); for (k in c("sex","sex_perm","sex_permw")) cd[[k]] <- factor(cd[[k]]); dds <- DESeqDataSetFromMatrix(cts, cd, f); dds <- DESeq(dds, quiet=TRUE, fitType="parametric"); r <- results(dds, name="group_case_vs_control"); rownames(r)[!is.na(r$padj) & r$padj < 0.05] }
jac <- function(a, b) { u <- length(union(a, b)); if (u == 0) NA_real_ else length(intersect(a, b)) / u }
rows <- list()
for (eng in c("limma", "edgeR", "DESeq2")) { fn <- switch(eng, limma=fit_limma, edgeR=fit_edger, DESeq2=fit_deseq)
  sets <- lapply(names(designs), function(m) if (ok[[m]]) tryCatch(fn(designs[[m]]), error=function(e) { message(eng, " ", m, ": ", conditionMessage(e)); NULL }) else NULL); names(sets) <- names(designs)
  J <- function(m) if (is.null(sets[[m]]) || is.null(sets$M0)) NA_real_ else jac(sets$M0, sets[[m]])
  rows[[eng]] <- data.table(study=SRP, organism=ORG, engine=eng, n=nrow(ph), n_universe=nrow(cts), n_de_M0=length(sets$M0), n_de_M1=length(sets$M1), n_de_perm=length(sets$M0perm), n_de_permw=length(sets$M0permw),
    J_M1=J("M1"), J_perm=J("M0perm"), J_permw=J("M0permw")) }
res <- rbindlist(rows)[, `:=`(d_across=J_M1 - J_perm, d_within=J_M1 - J_permw)]
fwrite(res, file.path(OUT, paste0(SRP, ".tsv")), sep="\t"); print(res[, .(engine, n_de_M0, n_de_M1, J_M1, J_perm, J_permw, d_across, d_within)])
