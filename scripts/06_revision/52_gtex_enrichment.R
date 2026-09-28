#!/usr/bin/env Rscript
## Enrichment of tissue-matched GTEx v8 (Oliva 2020) sex-biased genes — review H item 3 / I item 1.
## For each study in the human analysis set: how strongly the gene sets lost/gained by adjustment (real labels vs within-group vs across-group permutation) and the logFC shifts
## are enriched for tissue-matched GTEx sex-biased genes (autosomal only; sex chromosomes separately), and how many sex-biased genes the unadjusted DE set contains.
## Usage (server): Rscript scripts/06_revision/52_gtex_enrichment.R  → results/revision/G_gtex_per_study.tsv
suppressMessages(library(data.table))
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"; setwd(PROJ)
GT <- fread("/var2/lsg/db/gtex_v8_sbgenes/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt")[, gid := sub("\\..*", "", gene)]
ANN <- fread("data/annot/gencode_v29_gene_chr.tsv")[, gid := sub("\\..*", "", gene_id)]; XY <- ANN[chr %in% c("chrX","chrY"), gene_id]
TIS <- list(blood="Whole_Blood", brain="^Brain_", gut="^(Colon_|Small_Intestine|Stomach|Esophagus_)", liver="^Liver$", lung="^Lung$", kidney="^Kidney_", heart="^Heart_",
            muscle="^Muscle_Skeletal", adipose="^Adipose_", skin="^Skin_")
ntis <- uniqueN(GT$tissue); shared <- GT[, .N, by=gid][N >= ceiling(ntis/4), gid]          # no tissue match → 'shared' set: sex-biased in at least 1/4 of the 44 tissues
setmatch <- function(tc) { if (!is.null(TIS[[tc]])) { s <- unique(GT[grepl(TIS[[tc]], tissue), gid]); list(set=s, kind=tc) } else list(set=shared, kind="shared") }
S <- fread("results/main_human_summary_v2.tsv"); A <- readLines("config/c2_analysis_set_human.txt"); S <- S[study %in% A]
rd <- function(D, m) fread(file.path(D, paste0("de_", m, ".tsv.gz")))[, .(gene_id, logFC, t, adj.P.Val)]
fisher_or <- function(inset, insb) { tb <- table(factor(inset, levels=c(FALSE,TRUE)), factor(insb, levels=c(FALSE,TRUE))); f <- fisher.test(tb); c(or=unname(f$estimate), p=f$p.value) }
rows <- list()
for (i in seq_len(nrow(S))) { st <- S$study[i]; D <- sprintf("results/main_human_v2/%s", st); if (!file.exists(file.path(D, "de_M0permw.tsv.gz"))) next
  de <- list(M0=rd(D,"M0"), M1=rd(D,"M1"), Mw=rd(D,"M0permw"), Ma=rd(D,"M0perm")); U <- de$M0$gene_id; ug <- sub("\\..*", "", U)
  sm <- setmatch(S$tissue_category[i]); sb <- ug %in% sm$set; xy <- U %in% XY; sba <- sb & !xy
  DE0 <- de$M0[adj.P.Val < 0.05, gene_id]; sets <- lapply(de[c("M1","Mw","Ma")], function(d) d[adj.P.Val < 0.05, gene_id])
  L <- lapply(sets, function(s) setdiff(DE0, s)); G <- lapply(sets, function(s) setdiff(s, DE0))
  fr <- function(g) if (length(g) == 0) NA_real_ else mean(sba[match(g, U)])
  orf <- function(g) if (length(g) < 5) c(or=NA_real_, p=NA_real_) else fisher_or(U %in% g, sba)
  l0 <- de$M0$logFC; sh <- function(d) { dd <- (d$logFC - l0)^2; c(share=sum(dd[sba]) / sum(dd), ratio_medabs=median(sqrt(dd[sba])) / median(sqrt(dd[!sba & !xy]))) }
  o0 <- orf(DE0); o1 <- orf(sets$M1); shr <- sh(de$M1); shw <- sh(de$Mw); sha <- sh(de$Ma)
  rows[[st]] <- data.table(study=st, tissue_category=S$tissue_category[i], set_kind=sm$kind, design_sex=S$design_sex[i], sex_diff=S$sex_diff[i], n=S$n[i],
    n_U=length(U), n_sb=sum(sb), n_sb_auto=sum(sba), frac_U_sb_auto=mean(sba), n_de0=length(DE0),
    frac_de0_sb=fr(DE0), or_de0_sb=o0["or"], p_de0_sb=o0["p"], frac_de1_sb=fr(sets$M1), or_de1_sb=o1["or"],
    n_lost_real=length(L$M1), n_lost_w=length(L$Mw), n_lost_a=length(L$Ma), frac_lost_real=fr(L$M1), frac_lost_w=fr(L$Mw), frac_lost_a=fr(L$Ma),
    or_lost_real=orf(L$M1)["or"], or_lost_w=orf(L$Mw)["or"], or_lost_a=orf(L$Ma)["or"],
    n_gain_real=length(G$M1), n_gain_w=length(G$Mw), frac_gain_real=fr(G$M1), frac_gain_w=fr(G$Mw), frac_gain_a=fr(G$Ma), or_gain_real=orf(G$M1)["or"], or_gain_w=orf(G$Mw)["or"],
    shift_share_real=shr["share"], shift_share_w=shw["share"], shift_share_a=sha["share"], shift_ratio_real=shr["ratio_medabs"], shift_ratio_w=shw["ratio_medabs"], shift_ratio_a=sha["ratio_medabs"],
    xy_frac_lost_real=if (length(L$M1)) mean(xy[match(L$M1, U)]) else NA_real_, xy_frac_lost_w=if (length(L$Mw)) mean(xy[match(L$Mw, U)]) else NA_real_)
  if (i %% 20 == 0) cat(i, "\n") }
R <- rbindlist(rows); fwrite(R, "results/revision/G_gtex_per_study.tsv", sep="\t"); cat("studies:", nrow(R), " set kinds:", paste(names(table(R$set_kind)), table(R$set_kind), collapse=", "), "\n")
