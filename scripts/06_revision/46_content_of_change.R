#!/usr/bin/env Rscript
## Content of change + permutation distribution (review H/I v3 requests 1, 4, 5). For one study, newly fit B within-group and B across-group permutations, and compute
## (i) whether the genes dropped by the real labels are the same genes dropped by the permuted labels (set overlap, sex-chromosome fraction),
## (ii) Jaccard excluding sex chromosomes, (iii) t-rank perturbation and logFC shift (vs within-group permutation), (iv) pathway-level permutation distribution and per-study p-values.
## Usage: Rscript 46_content_of_change.R <human|mouse> <SRP> [B=20] [outdir=results/revision/content]
suppressMessages({library(recount3); library(edgeR); library(limma); library(data.table); library(fgsea); library(msigdbr)})
a <- commandArgs(trailingOnly=TRUE); ORG <- a[1]; SRP <- a[2]; B <- if (length(a) >= 3) as.integer(a[3]) else 20L; OUTDIR <- if (length(a) >= 4) a[4] else "results/revision/content"
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ); dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE); IS_MOUSE <- ORG == "mouse"
D <- sprintf("results/main_%s_v2/%s", ORG, SRP); ph <- fread(file.path(D, "pheno.tsv"))[, .(external_id, group, sex)]; ph[, group := factor(group, levels=c("control","case"))]
rd <- function(m) { x <- fread(file.path(D, paste0("de_", m, ".tsv.gz"))); x[, .(gene_id, gene, logFC, t, adj.P.Val)] }
de <- list(M0=rd("M0"), M1=rd("M1"), M0perm=rd("M0perm"), M0permw=rd("M0permw"))
## counts + chromosome
bfc <- recount3_cache(file.path(PROJ, "data/raw/recount3_cache", SRP))
rse <- create_rse_manual(project=SRP, project_home="data_sources/sra", organism=if (IS_MOUSE) "mouse" else "human", annotation=if (IS_MOUSE) "gencode_v23" else "gencode_v29", type="gene", bfc=bfc)
assay(rse, "counts") <- transform_counts(rse); cts <- assay(rse, "counts"); sym <- rowData(rse)$gene_name; chrom <- as.character(seqnames(rowRanges(rse)))
sacc_col <- grep("^sra\\.sample_acc", names(colData(rse)), value=TRUE)[1]; sacc <- data.table(external_id=colData(rse)$external_id, sample_acc=as.character(colData(rse)[[sacc_col]])); sacc[is.na(sample_acc) | sample_acc == "", sample_acc := external_id]
if (anyDuplicated(sacc$sample_acc)) { sacc[, sample_acc := factor(sample_acc, levels=unique(sample_acc))]; M <- model.matrix(~ 0 + sample_acc, sacc); rownames(M) <- sacc$external_id; rep_id <- sacc[, external_id[1], by=sample_acc]$V1; cts <- as.matrix(cts[, rownames(M)]) %*% M; colnames(cts) <- rep_id }
stopifnot(all(ph$external_id %in% colnames(cts))); cts <- cts[, ph$external_id]
KEEP <- filterByExpr(DGEList(cts), model.matrix(~group, ph)); gid <- rownames(cts)[KEEP]; stopifnot(setequal(gid, de$M0$gene_id))
XY <- data.table(gene_id=gid, gene=sym[KEEP], chr=chrom[KEEP])[, xy := chr %in% c("chrX","chrY","X","Y")]
fit <- function(design) { y <- calcNormFactors(DGEList(cts[KEEP, ])); v <- voom(y, design); f <- eBayes(lmFit(v, design)); tt <- as.data.table(topTable(f, coef="groupcase", number=Inf, sort.by="none")); tt$gene_id <- gid; tt$gene <- sym[KEEP]; tt[, .(gene_id, gene, logFC, t, adj.P.Val)] }
## Hallmark
GS <- { m <- as.data.table(msigdbr(species=if (IS_MOUSE) "Mus musculus" else "Homo sapiens", category="H")); split(m$gene_symbol, m$gs_name) }
keep_ids <- de$M0[!is.na(t) & !is.na(gene) & gene != ""][order(-abs(t))][!duplicated(gene), gene_id]
gsea <- function(d) { d <- d[gene_id %in% keep_ids & !is.na(t)]; rk <- setNames(d$t, d$gene); set.seed(1); r <- as.data.table(fgsea(GS, rk, minSize=15, maxSize=500, nproc=1))[, .(pathway, NES, padj)]; list(sig=r[padj < 0.05, pathway], top1=if (any(r$padj < 0.05)) r[order(padj, -abs(NES))][1, pathway] else NA_character_) }
jac <- function(a, b) { u <- length(union(a, b)); if (u == 0) NA_real_ else length(intersect(a, b)) / u }
deset <- function(d, noxy=FALSE) { s <- d[adj.P.Val < 0.05, gene_id]; if (noxy) s <- setdiff(s, XY[xy == TRUE, gene_id]); s }
## reference quantities
DE0 <- deset(de$M0); DE0n <- deset(de$M0, TRUE); t0 <- de$M0$t; l0 <- de$M0$logFC; ml0 <- median(abs(l0))
P0 <- gsea(de$M0)
summ_model <- function(d, tag) { s <- deset(d); L <- setdiff(DE0, s); G <- setdiff(s, DE0); p <- gsea(d)
  list(tag=tag, n_de=length(s), J=jac(DE0, s), J_noxy=jac(DE0n, deset(d, TRUE)), t_sp=cor(t0, d$t, method="spearman"), n_lost=length(L), n_gained=length(G),
       xy_lost=mean(XY[match(L, gene_id), xy]), xy_gained=mean(XY[match(G, gene_id), xy]), lfc_shift=median(abs(d$logFC - l0)) / ml0,
       xy_share_shift=sum((d$logFC - l0)[XY$xy]^2) / sum((d$logFC - l0)^2), PJ=jac(P0$sig, p$sig), top1_same=if (is.na(P0$top1)) NA else identical(P0$top1, p$top1), L=list(L), G=list(G), Psig=list(p$sig)) }
real <- summ_model(de$M1, "M1"); across1 <- summ_model(de$M0perm, "M0perm"); within1 <- summ_model(de$M0permw, "M0permw")
## B permutations of each kind
base <- sum(utf8ToInt(SRP)); mm <- function(v) model.matrix(~ group + v, data.frame(group=ph$group, v=v))
perm <- list()
for (b in seq_len(B)) { ph[, vw := { s <- sex; set.seed(base + 1000L + b + 17L * .GRP); if (.N > 1) sample(s) else s }, by=group]; vw <- ph$vw   # within-group permutation (per-group seed; := preserves row order — extracting $V1 misaligned labels in studies with interleaved groups, fixed 2026-09-23)
  Xw <- mm(vw); if (qr(Xw)$rank == ncol(Xw)) perm[[length(perm)+1]] <- c(summ_model(fit(Xw), "within"), b=b)
  set.seed(base + 2000L + b); va <- sample(ph$sex); Xa <- mm(va); if (qr(Xa)$rank == ncol(Xa)) perm[[length(perm)+1]] <- c(summ_model(fit(Xa), "across"), b=b) }
P <- rbindlist(lapply(perm, function(x) as.data.table(x[c("tag","b","n_de","J","J_noxy","t_sp","n_lost","n_gained","xy_lost","xy_gained","lfc_shift","xy_share_shift","PJ","top1_same")])), fill=TRUE)
P[, study := SRP]; P[, organism := ORG]
W <- perm[sapply(perm, `[[`, "tag") == "within"]; A <- perm[sapply(perm, `[[`, "tag") == "across"]
pair_mean <- function(lst, key) { if (length(lst) < 2) return(NA_real_); v <- c(); for (i in 1:(length(lst)-1)) for (j in (i+1):length(lst)) v <- c(v, jac(lst[[i]][[key]][[1]], lst[[j]][[key]][[1]])); mean(v, na.rm=TRUE) }
one <- function(x) if (length(x)) x else NA
row <- data.table(study=SRP, organism=ORG, n=nrow(ph), n_universe=length(gid), n_de_M0=length(DE0), n_xy_universe=sum(XY$xy), n_xy_de_M0=sum(XY[match(DE0, gene_id), xy]),
  J_M1=real$J, J_M1_noxy=real$J_noxy, J_across1=across1$J, J_within1=within1$J,
  J_within_mean=mean(sapply(W, `[[`, "J")), J_within_sd=sd(sapply(W, `[[`, "J")), J_across_mean=mean(sapply(A, `[[`, "J")), J_across_sd=sd(sapply(A, `[[`, "J")),
  J_within_noxy_mean=mean(sapply(W, `[[`, "J_noxy")), p_J_within=mean(sapply(W, `[[`, "J") <= real$J), p_J_across=mean(sapply(A, `[[`, "J") <= real$J),
  tsp_M1=real$t_sp, tsp_within_mean=mean(sapply(W, `[[`, "t_sp")), p_tsp_within=mean(sapply(W, `[[`, "t_sp") <= real$t_sp),
  n_lost_M1=real$n_lost, n_lost_within_mean=mean(sapply(W, `[[`, "n_lost")), n_gained_M1=real$n_gained, n_gained_within_mean=mean(sapply(W, `[[`, "n_gained")),
  J_lost_M1_vs_within=mean(sapply(W, function(w) jac(real$L[[1]], w$L[[1]]))), J_lost_within_vs_within=pair_mean(W, "L"),
  J_gained_M1_vs_within=mean(sapply(W, function(w) jac(real$G[[1]], w$G[[1]]))), J_gained_within_vs_within=pair_mean(W, "G"),
  xy_lost_M1=real$xy_lost, xy_lost_within_mean=mean(sapply(W, `[[`, "xy_lost"), na.rm=TRUE), xy_gained_M1=real$xy_gained, xy_gained_within_mean=mean(sapply(W, `[[`, "xy_gained"), na.rm=TRUE),
  lfc_shift_M1=real$lfc_shift, lfc_shift_within_mean=mean(sapply(W, `[[`, "lfc_shift")), xy_share_shift_M1=real$xy_share_shift, xy_share_shift_within_mean=mean(sapply(W, `[[`, "xy_share_shift")),
  PJ_M1=real$PJ, PJ_within_mean=mean(sapply(W, `[[`, "PJ"), na.rm=TRUE), PJ_across_mean=mean(sapply(A, `[[`, "PJ"), na.rm=TRUE),
  p_PJ_within=if (is.na(real$PJ)) NA_real_ else mean(sapply(W, `[[`, "PJ") <= real$PJ, na.rm=TRUE), p_PJ_across=if (is.na(real$PJ)) NA_real_ else mean(sapply(A, `[[`, "PJ") <= real$PJ, na.rm=TRUE),
  J_Psig_M1_vs_within=mean(sapply(W, function(w) jac(setdiff(P0$sig, real$Psig[[1]]), setdiff(P0$sig, w$Psig[[1]])))), 
  top1_changed_M1=if (is.na(real$top1_same)) NA else !real$top1_same, top1_changed_within_frac=mean(!sapply(W, `[[`, "top1_same"), na.rm=TRUE), top1_changed_across_frac=mean(!sapply(A, `[[`, "top1_same"), na.rm=TRUE), B_within=length(W), B_across=length(A))
fwrite(row, file.path(OUTDIR, paste0(SRP, ".tsv")), sep="\t"); fwrite(P, file.path(OUTDIR, paste0(SRP, "_perms.tsv")), sep="\t")
cat(SRP, ": J_M1", round(real$J,3), " within mean", round(row$J_within_mean,3), " p_within", row$p_J_within, " lost overlap real-vs-within", round(row$J_lost_M1_vs_within,3), " null", round(row$J_lost_within_vs_within,3), "\n")
