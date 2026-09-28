#!/usr/bin/env Rscript
## E9b + E10 (server part). Reads ONLY the frozen per-study DE tables; no refitting. One row per study.
## Usage (server):  nice -n 19 Rscript 39a_server_extract_mechanism_leakage.R <org> <frozen_dir> <summary.tsv> <out.tsv> [genelist.tsv] [cores=8]
##   frozen_dir e.g. results/main_human_v2_bak20260917 ; per study: de_M0.tsv.gz de_M1.tsv.gz de_M0perm.tsv.gz (gene_id, gene, logFC, t, adj.P.Val) + pheno.tsv (group, sex, sex_perm)
##   genelist.tsv (optional): columns gene [, tissue]; gene = symbol or Ensembl id (version ignored). With a tissue column the study's tissue_category is matched when present, else the union is used.
## Mechanism: SE = logFC / t (moderated SE). For OLS, SE(M1)/SE(M0) = sqrt(VIF) * sigma1/sigma0 and logFC(M0) − logFC(M1) = (difference in sex composition) × (gene's sex effect).
suppressMessages({library(data.table); library(parallel)})
a <- commandArgs(trailingOnly=TRUE); ORG <- a[1]; FRZ <- a[2]; SUMF <- a[3]; OUTF <- a[4]; GL <- if (length(a) >= 5 && nzchar(a[5]) && a[5] != "NA") a[5] else NA; CORES <- if (length(a) >= 6) as.integer(a[6]) else 8L
stopifnot(CORES <= 8); setDTthreads(1)
S <- fread(SUMF, select=c("study", "tissue_category", "design_sex", "valid", "single_sex"))
PANEL <- if (ORG == "human") c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "NLGN4Y", "TXLNGY", "USP9Y", "UTY", "ZFY") else c("Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
gl <- if (!is.na(GL)) { g <- fread(GL); if (!"tissue" %in% names(g)) g[, tissue := "all"]; g[, gene := sub("\\.[0-9]+$", "", gene)]; g[!gene %in% PANEL] } else NULL
studies <- S[valid == TRUE & single_sex == FALSE, study]; studies <- studies[file.exists(file.path(FRZ, studies, "de_M1.tsv.gz")) & file.exists(file.path(FRZ, studies, "de_M0.tsv.gz"))]
rd <- function(s, m) { f <- file.path(FRZ, s, sprintf("de_%s.tsv.gz", m)); if (!file.exists(f)) return(NULL); fread(f, select=c("gene_id", "gene", "logFC", "t", "adj.P.Val")) }
gini <- function(x) { x <- sort(x); n <- length(x); if (!n || sum(x) == 0) return(NA_real_); 2*sum(x*seq_len(n))/(n*sum(x)) - (n + 1)/n }
shift <- function(d, lab, ispanel) { ad <- abs(d); ss <- d^2; q99 <- quantile(ad, .99); top <- ad >= q99
  r <- list(med_abs=median(ad), p99_abs=unname(q99), sd=sd(d), top1pct_share_ss=sum(ss[top])/sum(ss), gini_abs=gini(ad), excess_kurtosis=mean((d - mean(d))^4)/var(d)^2 - 3,
            n_panel=sum(ispanel), panel_med_abs=if (any(ispanel)) median(ad[ispanel]) else NA_real_, panel_share_ss=sum(ss[ispanel])/sum(ss), panel_mean_pctile=if (any(ispanel)) mean(rank(ad)[ispanel])/length(ad) else NA_real_,
            panel_over_bg=if (any(ispanel)) median(ad[ispanel])/median(ad[!ispanel]) else NA_real_)
  setNames(r, paste0(lab, "_", names(r))) }
leak <- function(m, L, lab, B=200) { U <- nrow(m); inL <- m$inL_tmp; nL <- sum(inL); d0 <- m$p0 < 0.05; lost <- d0 & !(m$p1 < 0.05); kept <- d0 & m$p1 < 0.05; lostp <- if ("pp" %in% names(m)) d0 & !(m$pp < 0.05) else rep(NA, U); gained <- !d0 & m$p1 < 0.05
  set.seed(1); rnd <- if (sum(lost) > 0 && nL > 0) replicate(B, { z <- sample.int(U, nL); mean(lost[z])*nL/sum(lost) }) else NA_real_   # share of lost genes that fall in a random same-size gene set
  r <- list(n_list_in_universe=nL, bg_share=nL/U, n_lost=sum(lost), n_kept=sum(kept), n_gained=sum(gained), n_lostperm=sum(lostp), k_lost=sum(lost & inL), k_kept=sum(kept & inL), k_gained=sum(gained & inL), k_lostperm=sum(lostp & inL), k_M0de=sum(d0 & inL),
            share_lost=if (sum(lost)) mean(inL[lost]) else NA_real_, share_kept=if (sum(kept)) mean(inL[kept]) else NA_real_, share_gained=if (sum(gained)) mean(inL[gained]) else NA_real_, share_lostperm=if (sum(lostp, na.rm=TRUE)) mean(inL[which(lostp)]) else NA_real_,
            random_share_lost_mean=mean(rnd), random_share_lost_q975=if (all(is.na(rnd))) NA_real_ else unname(quantile(rnd, .975)))
  setNames(r, paste0(lab, "_", names(r))) }
one <- function(s) tryCatch({ d0 <- rd(s, "M0"); d1 <- rd(s, "M1"); dp <- rd(s, "M0perm")
  m <- merge(d0[, .(gene_id, gene, l0=logFC, t0=t, p0=adj.P.Val)], d1[, .(gene_id, l1=logFC, t1=t, p1=adj.P.Val)], by="gene_id"); if (!is.null(dp)) m <- merge(m, dp[, .(gene_id, lp=logFC, tp=t, pp=adj.P.Val)], by="gene_id", all.x=TRUE)
  m <- m[is.finite(t0) & is.finite(t1) & t0 != 0 & t1 != 0]; m[, `:=`(se0=l0/t0, se1=l1/t1)]; ispanel <- m$gene %in% PANEL
  ph <- fread(file.path(FRZ, s, "pheno.tsv")); g <- as.integer(ph$group == "case"); x <- as.integer(ph$sex == "F"); r <- suppressWarnings(cor(g, x)); rp <- if ("sex_perm" %in% names(ph)) suppressWarnings(cor(g, as.integer(ph$sex_perm == "F"))) else NA_real_
  out <- list(study=s, n=nrow(ph), n_genes=nrow(m), r_group_sex=r, vif=1/(1 - r^2), delta_f_signed=mean(x[g == 1]) - mean(x[g == 0]), r_group_sexperm=rp, vif_perm=1/(1 - rp^2),
              se_ratio_M1_median=median(m$se1/m$se0), se_ratio_M1_q1=unname(quantile(m$se1/m$se0, .25)), se_ratio_M1_q3=unname(quantile(m$se1/m$se0, .75)), se_ratio_M1_panel_median=if (any(ispanel)) median((m$se1/m$se0)[ispanel]) else NA_real_,
              med_abs_l0=median(abs(m$l0)), cor_l0_l1=cor(m$l0, m$l1), n_deg_M0=sum(m$p0 < 0.05), n_deg_M1=sum(m$p1 < 0.05))
  out <- c(out, shift(m$l1 - m$l0, "dM1", ispanel))
  if ("tp" %in% names(m)) { k <- is.finite(m$tp) & m$tp != 0; out$se_ratio_perm_median <- median((m$lp/m$tp)[k]/m$se0[k]); out$n_deg_perm <- sum(m$pp < 0.05, na.rm=TRUE); out <- c(out, shift((m$lp - m$l0)[k], "dPerm", ispanel[k])) }
  ## why are M0-DE genes lost under M1? log|t1/t0| = log|l1/l0| − log(se1/se0)
  lost <- m$p0 < 0.05 & !(m$p1 < 0.05); if (sum(lost) >= 5) { out$lost_med_log_lfc_ratio <- median(log(abs(m$l1/m$l0))[lost]); out$lost_med_log_se_ratio <- median(log(m$se1/m$se0)[lost]); out$lost_med_log_t_ratio <- median(log(abs(m$t1/m$t0))[lost])
    out$lost_frac_lfc_shrunk_gt_se_inflated <- mean((-log(abs(m$l1/m$l0)) > log(m$se1/m$se0))[lost]) }
  m[, inL_tmp := ispanel]; out <- c(out, leak(m, NULL, "xy"))
  if (!is.null(gl)) { tc <- S[study == s, tissue_category][1]; L <- if (!is.na(tc) && tc %in% gl$tissue) gl[tissue == tc, unique(gene)] else unique(gl$gene); out$list_tissue_matched <- !is.na(tc) && tc %in% gl$tissue
    m[, inL_tmp := (gene %in% L | sub("\\.[0-9]+$", "", gene_id) %in% L) & !ispanel]; out <- c(out, leak(m, NULL, "list")) }
  as.data.table(out) }, error=function(e) data.table(study=s, error=conditionMessage(e)))
res <- rbindlist(mclapply(studies, one, mc.cores=CORES, mc.preschedule=FALSE), fill=TRUE); res[, organism := ORG]
fwrite(res, paste0(OUTF, ".tmp"), sep="\t"); file.rename(paste0(OUTF, ".tmp"), OUTF); cat(ORG, nrow(res), "studies;", if ("error" %in% names(res)) sum(!is.na(res$error)) else 0, "errors\n")
