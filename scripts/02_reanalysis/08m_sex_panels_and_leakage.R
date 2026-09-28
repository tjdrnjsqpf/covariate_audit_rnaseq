#!/usr/bin/env Rscript
# (1) Fit ~group+sex in each mouse pilot dataset -> save sex coefficients
# (2) Derive tissue-specific leave-one-out sex-biased autosomal panels
# (3) Compute extended-panel leakage in pilot M0 DEGs -> relation to sex_diff, detection sensitivity / false alarms
suppressMessages({library(recount3); library(edgeR); library(limma); library(data.table)})
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
PIL <- file.path(PROJ, "results/pilot_mouse"); OUT <- file.path(PROJ, "results/sex_panels_mouse"); dir.create(OUT, showWarnings=FALSE)
ps <- fread(file.path(PROJ, "results/pilot_mouse_summary.tsv"))[, .(study, tissue_category, sex_diff, design, models, n)]
source(file.path(PROJ, "scripts/02_reanalysis/sex_panel_utils.R"))
sexchr <- c("Xist","Ddx3y","Eif2s3y","Kdm5d","Uty")

## (1) Sex coefficients
for (st in ps$study) {
  f <- file.path(PIL, st, "sex_de.tsv.gz"); if (file.exists(f)) next
  ph <- fread(file.path(PIL, st, "pheno.tsv"))[sex %in% c("F","M")]; if (uniqueN(ph$sex) < 2 || min(table(ph$sex)) < 2) next
  bfc <- recount3_cache(file.path(PROJ, "data/raw/recount3_cache", st))
  rse <- tryCatch(create_rse_manual(project=st, project_home="data_sources/sra", organism="mouse", annotation="gencode_v23", type="gene", bfc=bfc), error=function(e) NULL); if (is.null(rse)) next
  cts <- transform_counts(rse)[, ph$external_id]; sym <- rowData(rse)$gene_name; chr <- as.character(seqnames(rowRanges(rse)))
  ph[, group := factor(group, levels=c("control","case"))][, sex := factor(sex, levels=c("F","M"))]
  X <- model.matrix(~group + sex, ph); if (qr(X)$rank < ncol(X)) next
  y <- DGEList(cts); keep <- filterByExpr(y, X); y <- calcNormFactors(y[keep, , keep.lib.sizes=FALSE]); v <- voom(y, X); fit <- eBayes(lmFit(v, X))
  tt <- topTable(fit, coef="sexM", number=Inf, sort.by="none")
  fwrite(data.table(gene_id=rownames(tt), gene=sym[keep], chr=chr[keep], logFC_M_vs_F=tt$logFC, P=tt$P.Value, FDR=tt$adj.P.Val), f, sep="\t")
  cat("sex DE:", st, nrow(ph), "samples\n")
}
sd_all <- rbindlist(lapply(ps$study, function(st) { f <- file.path(PIL, st, "sex_de.tsv.gz"); if (!file.exists(f)) return(NULL); fread(f)[, study := st] }))
sd_all <- merge(sd_all, ps[, .(study, tissue_category, n)], by="study")
fwrite(sd_all, file.path(OUT, "sex_de_all_pilot_mouse.tsv.gz"), sep="\t")

## (2) Tissue-specific signature (Stouffer meta-analysis, autosomal, top 300) — full and leave-one-out
tissues <- unique(ps$tissue_category)
cat("\n== tissue signatures (all datasets in tissue; autosomal, FDR<0.01, |mean logFC|>0.5, top300) ==\n")
sig_all <- lapply(setNames(tissues, tissues), function(tc) sex_signature(sd_all, ps[tissue_category == tc, study]))
for (tc in tissues) cat(sprintf("%-14s datasets=%2d signature=%3d  e.g. %s\n", tc, ps[tissue_category == tc, .N], nrow(sig_all[[tc]]), paste(head(sig_all[[tc]]$gene, 8), collapse=",")))
fwrite(rbindlist(lapply(tissues, function(tc) cbind(tissue=tc, sig_all[[tc]]))), file.path(OUT, "tissue_sex_signatures_mouse.tsv"), sep="\t")
sig_pan <- sex_signature(sd_all, ps$study); cat("pan-tissue signature:", nrow(sig_pan), "\n")

## (3) Pilot M0: (a) number of DEGs among the 5 sex-chromosome genes, (b) t-correlation with the LOO tissue signature, (c) correlation with the LOO pan-tissue signature
res <- rbindlist(lapply(ps$study, function(st) {
  f <- file.path(PIL, st, "de_M0.tsv.gz"); if (!file.exists(f)) return(NULL)
  m0 <- fread(f); deg <- m0[adj.P.Val < 0.05, gene]; tc <- ps[study == st, tissue_category]
  sig_t <- sex_signature(sd_all, setdiff(ps[tissue_category == tc, study], st)); sig_p <- sex_signature(sd_all, setdiff(ps$study, st))
  lt <- leakage_cor(m0, sig_t); lp <- leakage_cor(m0, sig_p)
  ## fraction of signature genes that are DEGs vs overall DEG fraction
  data.table(study=st, tissue=tc, n_deg=length(deg), deg_frac=length(deg)/nrow(m0), sexchr_in_deg=sum(deg %in% sexchr_genes),
             sig_t_n=nrow(sig_t), sig_t_cor=round(lt$cor, 3), sig_t_in_deg=sum(deg %in% sig_t$gene), sig_t_enrich=if (nrow(sig_t)) round(mean(sig_t$gene %in% deg) / max(1e-9, length(deg)/nrow(m0)), 2) else NA_real_,
             sig_p_n=nrow(sig_p), sig_p_cor=round(lp$cor, 3), sig_p_in_deg=sum(deg %in% sig_p$gene))
}))
res <- merge(res, ps[, .(study, sex_diff, design, models, n)], by="study")
## Direction correction: if the case group has more males, a positive correlation with signature z (M vs F) is expected, so also compute a sign-aligned "signed leakage"
ph_dir <- rbindlist(lapply(ps$study, function(st) { f <- file.path(PIL, st, "pheno.tsv"); if (!file.exists(f)) return(NULL); p <- fread(f)[sex %in% c("F","M")]; data.table(study=st, male_excess_case = p[group == "case", mean(sex == "M")] - p[group == "control", mean(sex == "M")]) }))
res <- merge(res, ph_dir, by="study", all.x=TRUE)
res[, sig_t_cor_signed := sig_t_cor * sign(male_excess_case)][, sig_p_cor_signed := sig_p_cor * sign(male_excess_case)]
fwrite(res, file.path(OUT, "pilot_leakage_signature.tsv"), sep="\t")
cat("\n== Spearman(sex_diff, metric), all with M0 result ==\n")
for (v in c("sexchr_in_deg","sig_t_cor_signed","sig_p_cor_signed","sig_t_in_deg","sig_t_enrich")) { ok <- !is.na(res[[v]]); cat(sprintf("  %-18s rho=%.2f (n=%d)\n", v, cor(res$sex_diff[ok], res[[v]][ok], method="spearman"), sum(ok))) }
cat("\n== detector by sex_diff bin: sexchr>=1 | |sig_t_cor|>=0.2 | |sig_p_cor|>=0.2 | signed sig_p_cor>=0.2 ==\n")
res[, sd_bin := cut(sex_diff, c(-0.01, 0.15, 0.3, 0.5, 0.75, 1.01), labels=c("0-0.15","0.15-0.3","0.3-0.5","0.5-0.75","0.75-1"))]
print(res[, .(n=.N, sexchr=round(mean(sexchr_in_deg >= 1), 2), sig_t=round(mean(abs(sig_t_cor) >= 0.2, na.rm=TRUE), 2), sig_p=round(mean(abs(sig_p_cor) >= 0.2, na.rm=TRUE), 2), sig_p_signed=round(mean(sig_p_cor_signed >= 0.2, na.rm=TRUE), 2), med_sig_p_signed=round(median(sig_p_cor_signed, na.rm=TRUE), 2)), by=sd_bin][order(sd_bin)])
print(res[order(-sex_diff), .(study, tissue, design, sex_diff, n_deg, sexchr_in_deg, sig_t_n, sig_t_cor_signed, sig_p_n, sig_p_cor_signed)])
