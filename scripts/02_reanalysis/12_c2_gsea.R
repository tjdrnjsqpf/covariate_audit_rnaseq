#!/usr/bin/env Rscript
## C2: does covariate adjustment change the "conclusion" (pathway level)?
## Run fgsea (Hallmark 50, Reactome) on the saved t statistics and compute the change in the set of significant pathways
## for M0 vs M1/M2 and the negative controls (M0perm/M2perm). Tests whether the gene-level metric (Jaccard) merely reflects movement of borderline genes.
## Also C2-b: are the DEGs lost from M0 to M1 borderline (adj.P 0.01–0.05) or substantive genes?
## Usage: Rscript 12_c2_gsea.R <human|mouse> [ncores]
suppressMessages({library(data.table); library(fgsea); library(msigdbr); library(parallel)})
args <- commandArgs(trailingOnly=TRUE); ORG <- args[1]; NC <- if (length(args) >= 2) as.integer(args[2]) else 16L
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ); dir.create("results/c2", showWarnings=FALSE)
SPECIES <- if (ORG == "human") "Homo sapiens" else "Mus musculus"
studies <- readLines(if (length(args) >= 3) args[3] else sprintf("config/c2_analysis_set_%s.txt", ORG))   # 3rd argument: study list file (for testing)
TAG <- if (length(args) >= 3) "_test" else ""
OUTD <- sprintf("results/main_%s_v2", ORG)

## Gene sets
gs_list <- function(cat, sub=NULL) { m <- as.data.table(msigdbr(species=SPECIES, category=cat, subcollection=sub)); split(m$gene_symbol, m$gs_name) }
GS <- list(H = gs_list("H"), R = gs_list("C2", "CP:REACTOME"))
cat("Gene sets: Hallmark", length(GS$H), " Reactome", length(GS$R), "\n")

## t statistics -> named rank vector (duplicate symbols: max |t|)
## 2026-09-18 (review D 3.2): the representative gene_id for duplicate symbols is chosen once by max |t| in M0 and applied identically to all models (model-independent)
keep_ids_of <- function(d0) { d0 <- d0[!is.na(t) & !is.na(gene) & gene != ""]; d0[order(-abs(t))][!duplicated(gene), gene_id] }
rank_of <- function(d, keep_ids) { d <- d[gene_id %in% keep_ids & !is.na(t)]; setNames(d$t, d$gene) }
run_gsea <- function(rk, sets) { set.seed(1); as.data.table(fgsea(sets, rk, minSize=15, maxSize=500, nproc=1))[, .(pathway, NES, padj)] }

## Pathway-level change metrics (reference vs alternative)
pmetric <- function(a, b, thr=0.05) {
  m <- merge(a, b, by="pathway", suffixes=c("_r","_a"))
  sr <- m[padj_r < thr, pathway]; sa <- m[padj_a < thr, pathway]; un <- length(union(sr, sa))
  ## 2026-09-18 (review D 3.4): padj ties are broken by -|NES|; top1_same is NA if the reference model has no significant pathway
  tr <- m[order(padj_r, -abs(NES_r))][1:min(20,.N), pathway]; ta <- m[order(padj_a, -abs(NES_a))][1:min(20,.N), pathway]
  list(n_sig_ref = length(sr), n_sig_alt = length(sa),
       jaccard = if (un == 0) NA_real_ else round(length(intersect(sr, sa)) / un, 3),
       top20_overlap = round(length(intersect(tr, ta)) / length(tr), 3),
       nes_spearman = round(cor(m$NES_r, m$NES_a, method="spearman", use="complete.obs"), 3),
       sign_flip_sig = sum(sign(m[pathway %in% sr, NES_r]) != sign(m[pathway %in% sr, NES_a])),
       top1_same = if (length(sr)) tr[1] == ta[1] else NA)
}
## C2-b: nature of the lost DEGs at gene level
gmetric <- function(a, b, fdr=0.05) {
  m <- merge(a[, .(gene_id, lfc0=logFC, p0=adj.P.Val, t0=t)], b[, .(gene_id, p1=adj.P.Val)], by="gene_id")
  lost <- m[p0 < fdr & p1 >= fdr]; kept <- m[p0 < fdr & p1 < fdr]; gained <- m[p0 >= fdr & p1 < fdr]
  list(n_deg0 = sum(m$p0 < fdr), n_lost = nrow(lost), n_gained = nrow(gained),
       lost_frac_boundary = if (nrow(lost)) round(mean(lost$p0 > 0.01), 3) else NA_real_,      # fraction of borderline genes (adj.P 0.01–0.05)
       lost_med_abs_lfc = if (nrow(lost)) round(median(abs(lost$lfc0)), 3) else NA_real_,
       kept_med_abs_lfc = if (nrow(kept)) round(median(abs(kept$lfc0)), 3) else NA_real_,
       lost_frac_big_lfc = if (nrow(lost)) round(mean(abs(lost$lfc0) > 1), 3) else NA_real_)      # fraction of lost DEGs with |logFC|>1
}

one <- function(s) {
  d <- file.path(OUTD, s); rd <- function(m) { f <- file.path(d, paste0("de_", m, ".tsv.gz")); if (file.exists(f)) fread(f) else NULL }
  de <- list(M0=rd("M0"), M1=rd("M1"), M0perm=rd("M0perm"), M0permw=rd("M0permw"), M0a=rd("M0a"), M2=rd("M2"), M2perm=rd("M2perm"))
  if (is.null(de$M0) || is.null(de$M1)) return(NULL)
  out <- list(study = s); long <- list()
  keep_ids <- keep_ids_of(de$M0)
  gs <- lapply(de[!sapply(de, is.null)], function(x) { rk <- rank_of(x, keep_ids); list(H = run_gsea(rk, GS$H), R = run_gsea(rk, GS$R)) })
  for (m in names(gs)) long[[m]] <- gs[[m]]$H[, .(study = s, model = m, pathway, NES, padj)]
  cmp <- list(P1 = c("M0","M1"), P0perm = c("M0","M0perm")); if (!is.null(de$M0permw)) cmp$P0permw <- c("M0","M0permw")
  if (!is.null(de$M2)) { ref2 <- if (is.null(de$M0a)) "M0" else "M0a"; cmp$P2 <- c(ref2, "M2"); if (!is.null(de$M2perm)) cmp$P2perm <- c(ref2, "M2perm") }
  for (k in names(cmp)) { r <- cmp[[k]][1]; a <- cmp[[k]][2]; if (is.null(gs[[r]]) || is.null(gs[[a]])) next
    for (col in c("H","R")) { mt <- pmetric(gs[[r]][[col]], gs[[a]][[col]]); for (n in names(mt)) out[[paste0(k, "_", col, "_", n)]] <- mt[[n]] } }
  ## C2-b
  gm <- gmetric(de$M0, de$M1); for (n in names(gm)) out[[paste0("G1_", n)]] <- gm[[n]]
  if (!is.null(de$M0perm)) { gm <- gmetric(de$M0, de$M0perm); for (n in names(gm)) out[[paste0("G0perm_", n)]] <- gm[[n]] }
  list(summ = as.data.table(out), long = rbindlist(long))
}

res <- mclapply(studies, function(s) tryCatch(one(s), error=function(e) { message(s, ": ", conditionMessage(e)); NULL }), mc.cores=NC)
ok <- res[!sapply(res, is.null)]
summ <- rbindlist(lapply(ok, `[[`, "summ"), fill=TRUE); long <- rbindlist(lapply(ok, `[[`, "long"))
fwrite(summ, sprintf("results/c2/gsea_summary_%s%s.tsv", ORG, TAG), sep="\t")
fwrite(long, sprintf("results/c2/hallmark_long_%s%s.tsv.gz", ORG, TAG), sep="\t")
cat("Done:", nrow(summ), "/", length(studies), "studies\n")
