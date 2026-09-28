#!/usr/bin/env Rscript
## E6. Top-1 Hallmark pathway: recompute top1_same from hallmark_long under the OLD rule (rank by padj, ties -> file/alphabetical order, never NA)
##     and the NEW rule (rank by padj then −|NES|; NA when M0 has no pathway at padj < 0.05). Class-level "top pathway changed" rates old vs new.
## NOTE 2026-09-18: the local results/c2/hallmark_long_human.tsv.gz (09-08) lacks SRP114567; the server copy (09-10) has it. Reported numbers were produced with
##   HALLMARK_HUMAN=results/revision/inputs/hallmark_long_human_server_20260910.tsv.gz  (a verbatim copy of the server file).
##     "Conclusion changed" = P1_H_jaccard < P0perm_H_jaccard − margin at margins 0.1/0.15/0.2/0.3, also stratified by number of significant M0 pathways.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
HL  <- c(human=Sys.getenv("HALLMARK_HUMAN", "results/c2/hallmark_long_human.tsv.gz"), mouse=Sys.getenv("HALLMARK_MOUSE", "results/c2/hallmark_long_mouse.tsv.gz"))
GSS <- c(human=Sys.getenv("GSEA_HUMAN", "results/c2/gsea_summary_human.tsv"), mouse=Sys.getenv("GSEA_MOUSE", "results/c2/gsea_summary_mouse.tsv"))
D <- load_summary()
T1 <- rbindlist(lapply(ORGS, function(o) { h <- fread(HL[[o]])[model %in% c("M0", "M1", "M0perm")]; h[, idx := seq_len(.N)]   # idx keeps file order (the old tie-break)
  w <- dcast(h, study + pathway ~ model, value.var=c("NES", "padj")); w <- merge(w, h[model == "M0", .(study, pathway, idx)], by=c("study", "pathway")); setorder(w, study, idx)
  one <- function(x, alt) { if (!nrow(x) || all(is.na(x$padj_M0))) return(list(NA, NA, NA_integer_, NA_integer_))   # NA padj sorts last, as in 12_c2_gsea.R
    pa <- x[[paste0("padj_", alt)]]; na <- x[[paste0("NES_", alt)]]
    old <- x$pathway[order(x$padj_M0)[1]] == x$pathway[order(pa)[1]]
    new <- if (any(x$padj_M0 < 0.05, na.rm=TRUE)) x$pathway[order(x$padj_M0, -abs(x$NES_M0))[1]] == x$pathway[order(pa, -abs(na))[1]] else NA
    list(old, new, sum(x$padj_M0 == min(x$padj_M0, na.rm=TRUE), na.rm=TRUE), sum(x$padj_M0 < 0.05, na.rm=TRUE)) }
  r <- w[, { a <- one(.SD, "M1"); b <- one(.SD, "M0perm"); .(top1_same_M1_old=a[[1]], top1_same_M1_new=a[[2]], top1_same_perm_old=b[[1]], top1_same_perm_new=b[[2]], n_tied_at_min_padj_M0=a[[3]], n_sig_M0=a[[4]]) }, by=study]
  r[, organism := o]; r }))
G <- rbindlist(lapply(ORGS, function(o) { g <- fread(GSS[[o]]); g[, organism := o]; g[, .(organism, study, P1_H_n_sig_ref, P1_H_jaccard, P0perm_H_jaccard, P1_H_top1_same, P0perm_H_top1_same)] }))
X <- merge(D[, .(organism=as.character(organism), study, analysis, design_sex, n, sex_diff, M1_jaccard, M0perm_jaccard)], merge(G, T1, by=c("organism", "study"), all=TRUE), by=c("organism", "study"))
X <- X[analysis == TRUE & design_sex %in% CLS & !is.na(M1_jaccard) & !is.na(M0perm_jaccard) & !is.na(P1_H_jaccard) & !is.na(P0perm_H_jaccard)]   # Fig 2A set
X[, cls := factor(design_sex, levels=CLS)]; X[, organism := factor(organism, levels=ORGS)]
cat("agreement of recomputed OLD rule with stored gsea_summary top1_same (M1 / perm):", X[, mean(top1_same_M1_old == P1_H_top1_same, na.rm=TRUE)], X[, mean(top1_same_perm_old == P0perm_H_top1_same, na.rm=TRUE)], "\n")
cat("n_sig_M0 recomputed == stored:", X[, mean(n_sig_M0 == P1_H_n_sig_ref, na.rm=TRUE)], "\n")
rate <- function(v) { v <- v[!is.na(v)]; w <- wilson(sum(!v), length(v)); sprintf("%d/%d = %.1f%% [%.1f, %.1f]", w$k, w$n, w$pct, w$lo, w$hi) }
num <- function(v) { v <- v[!is.na(v)]; 100*mean(!v) }
R1 <- X[, .(n_studies=.N, n_M0_has_sig=sum(n_sig_M0 > 0, na.rm=TRUE), n_missing_from_hallmark_long=sum(is.na(n_sig_M0)), median_tied_at_min_padj=as.double(median(n_tied_at_min_padj_M0)),
            changed_M1_old=rate(top1_same_M1_old), changed_perm_old=rate(top1_same_perm_old), changed_M1_new=rate(top1_same_M1_new), changed_perm_new=rate(top1_same_perm_new),
            pct_M1_old=num(top1_same_M1_old), pct_perm_old=num(top1_same_perm_old), pct_M1_new=num(top1_same_M1_new), pct_perm_new=num(top1_same_perm_new)), keyby=.(organism, cls)]
R1[, `:=`(ratio_old=pct_M1_old/pct_perm_old, ratio_new=pct_M1_new/pct_perm_new, excess_pp_old=pct_M1_old - pct_perm_old, excess_pp_new=pct_M1_new - pct_perm_new)]; out(R1, "E6_top1_changed_old_vs_new.tsv")
## paired (McNemar) real vs permuted under the new rule
MC <- X[!is.na(top1_same_M1_new) & !is.na(top1_same_perm_new), { b <- sum(!top1_same_M1_new & top1_same_perm_new); c <- sum(top1_same_M1_new & !top1_same_perm_new)
  .(only_real_changed=b, only_perm_changed=c, mcnemar_exact_p=if (b + c) binom.test(b, b + c)$p.value else NA_real_) }, keyby=.(organism, cls)]; out(MC, "E6_top1_mcnemar_new_rule.tsv")
## conclusion changed at margins
MG <- c(0.1, 0.15, 0.2, 0.3)
cc <- function(d, by) rbindlist(lapply(MG, function(m) d[, { w <- wilson(sum(P1_H_jaccard < P0perm_H_jaccard - m), .N); wr <- wilson(sum(P0perm_H_jaccard < P1_H_jaccard - m), .N)
  .(margin=m, n_studies=.N, k_changed=w$k, pct_changed=w$pct, lo=w$lo, hi=w$hi, k_reverse=wr$k, pct_reverse=wr$pct) }, keyby=by]))
R2 <- cc(X, c("organism", "cls")); setkey(R2, organism, cls, margin); out(R2, "E6_conclusion_changed_margins.tsv")
X[, nsig_grp := cut(P1_H_n_sig_ref, c(-1, 0, 5, 15, 50), labels=c("0", "1-5", "6-15", ">15"))]
R3 <- cc(X, c("organism", "nsig_grp", "cls")); setkey(R3, organism, nsig_grp, cls, margin); out(R3, "E6_conclusion_changed_margins_by_nsig.tsv")
NSG <- X[, .(n_studies=.N, nsig_median=as.double(median(P1_H_n_sig_ref)), pct_nsig0=100*mean(P1_H_n_sig_ref == 0), pct_nsig_1to5=100*mean(P1_H_n_sig_ref %in% 1:5)), keyby=.(organism, cls)]; out(NSG, "E6_nsig_by_class.tsv")
out(X[, .(organism, study, design_sex, n_sig_M0, n_tied_at_min_padj_M0, top1_same_M1_old, top1_same_M1_new, top1_same_perm_old, top1_same_perm_new, P1_H_jaccard, P0perm_H_jaccard)], "E6_per_study_top1.tsv")
options(width=250); print(R1[, 1:10]); print(R1[, c(1, 2, 11:18)], digits=3); print(MC); print(dcast(R2, organism + cls + n_studies ~ margin, value.var="pct_changed"), digits=3); print(dcast(R2, organism + cls ~ margin, value.var="pct_reverse"), digits=3)
print(dcast(R3[margin == 0.2], organism + cls ~ nsig_grp, value.var=c("n_studies", "pct_changed")), digits=3); print(NSG, digits=3)
