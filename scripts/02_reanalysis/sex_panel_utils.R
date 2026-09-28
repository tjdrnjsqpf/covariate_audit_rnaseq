# Tissue-specific sex-biased gene signature (Stouffer meta-analysis) + leakage statistic. Shared by 08m (pilot) and 06m (simulation).
sexchr_genes <- c("Xist","Ddx3y","Eif2s3y","Kdm5d","Uty")
## sd_all: study, gene, chr, logFC_M_vs_F, P, FDR, n (sample size). Returns: per-gene combined z (M vs F), mean logFC, n_datasets
sex_signature <- function(sd_all, studies, top_n=300, autosomal_only=TRUE) {
  d <- sd_all[study %in% studies & !chr %in% c("chrM")]
  if (autosomal_only) d <- d[!chr %in% c("chrX","chrY")]
  if (!length(unique(d$study))) return(data.table(gene=character(0), z=numeric(0), mlfc=numeric(0), nd=integer(0)))
  d[, z := sign(logFC_M_vs_F) * qnorm(pmax(P/2, 1e-300), lower.tail=FALSE)][, w := sqrt(n)]
  g <- d[, .(z = sum(w * z) / sqrt(sum(w^2)), mlfc = mean(logFC_M_vs_F), nd = .N, consist = max(mean(logFC_M_vs_F > 0), mean(logFC_M_vs_F < 0))), by=gene]
  g <- g[nd >= 2 & consist >= 0.7][order(-abs(z))]
  g[, fdr := p.adjust(2 * pnorm(-abs(z)), "BH")]
  head(g[fdr < 0.01 & abs(mlfc) > 0.5], top_n)
}
## Leakage statistic: correlation between the study's group t-statistics (M0) and signature z (restricted to signature genes). A large absolute value means sex signal is mixed into the group effect
leakage_cor <- function(de_table, sig) {   # de_table: gene, t
  m <- merge(de_table[, .(gene, t)], sig[, .(gene, z)], by="gene"); if (nrow(m) < 20) return(list(cor=NA_real_, n=nrow(m)))
  list(cor = cor(m$t, m$z, method="spearman"), n = nrow(m))
}
