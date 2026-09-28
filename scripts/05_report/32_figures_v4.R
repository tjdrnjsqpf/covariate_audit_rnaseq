#!/usr/bin/env Rscript
## v4 (2026-09-23) regenerate Figures 3 and 4 — incorporating chain 6 results.
##   Fig 3: (A) Δ vs two controls (B) false-DE fraction vs out-of-sample reference C (C) recall vs reference C (D) t correlation vs reference C (E) logFC bias²/variance (reference A) (F) leakage (fraction of draws)
##   Fig 4: (A) fraction with p ≤ 0.05 against the permutation distribution (B) pathway conclusion change (C) top-pathway change (D) SVA (E) paper audit
## Colours: red = real sex/adjusted, light grey = across-group permutation, black = within-group permutation, dark grey = descriptive/unadjusted. Species: left (Human) / right (Mouse).
## Run (Mac): LC_ALL=en_US.UTF-8 COV_AUDIT_PROJ="$(pwd)" Rscript scripts/05_report/32_figures_v4.R [content_per_study.tsv]
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
PROJ <- Sys.getenv("COV_AUDIT_PROJ", getwd()); setwd(PROJ)
CFILE <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(CFILE)) CFILE <- "results/revision/C_content_per_study.tsv"
FONT <- if (Sys.info()[["sysname"]] == "Darwin") "Helvetica Neue" else "sans"
PD <- "results/figures/draft_v4/panels"; SD <- "results/figures/draft_v4/source_data"; KD <- "results/figures/draft_v4/source_data_keynote"
for (d in c(PD, SD, KD)) dir.create(d, showWarnings=FALSE, recursive=TRUE)
RED <- "#B02418"; GRY <- "#C4C4C4"; BLK <- "#1A1A1A"; DGRY <- "#5F5F5F"; LRED <- "#E3A29B"
ORG <- c(human="Human", mouse="Mouse"); lv <- c("balanced","mildly_imbalanced","partially_confounded"); lab1 <- c(balanced="Balanced", mildly_imbalanced="Mild", partially_confounded="Confounded")
th <- theme_classic(base_size=9, base_family=FONT) + theme(axis.line=element_line(colour="black", linewidth=0.4), axis.ticks=element_line(colour="black", linewidth=0.4), axis.text=element_text(colour="black"),
  legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"), legend.key.size=unit(9,"pt"), legend.text=element_text(size=8),
  strip.background=element_blank(), strip.text=element_text(face="bold", size=9, hjust=0), panel.spacing.x=unit(14,"pt"), plot.margin=margin(4,8,4,4))
tg <- theme(plot.tag=element_text(face="bold", size=11, family=FONT), plot.tag.position=c(0, 1), plot.margin=margin(16, 14, 6, 8))
update_geom_defaults("text", list(family=FONT))
bci <- function(x, R=2000) { x <- x[is.finite(x)]; if (length(x) < 3) return(c(median(x), NA, NA)); set.seed(1); b <- replicate(R, median(sample(x, replace=TRUE))); c(median(x), unname(quantile(b,.025)), unname(quantile(b,.975))) }
sv <- function(p, f, w, h) { pdfdev <- if (Sys.info()[["sysname"]] == "Darwin") function(filename, width, height, ...) grDevices::quartz(file=filename, type="pdf", width=width, height=height) else cairo_pdf
  ggsave(file.path(PD, paste0(f, ".pdf")), p, width=w, height=h, device=pdfdev); tmp <- tempfile(fileext=".png")
  if (requireNamespace("ragg", quietly=TRUE)) ggsave(tmp, p, width=w, height=h, dpi=300, device=ragg::agg_png) else ggsave(tmp, p, width=w, height=h, dpi=300); invisible(file.copy(tmp, file.path(PD, paste0(f, ".png")), overwrite=TRUE)) }
src <- function(d, f) fwrite(d, file.path(SD, paste0(f, ".tsv")), sep="\t"); kn <- function(d, f) fwrite(d, file.path(KD, paste0(f, ".tsv")), sep="\t")
wide_err <- function(d, f, xcol, scol) { d <- copy(d)[, `:=`(err_minus=m-lo, err_plus=hi-m)]; kn(dcast(d, as.formula(sprintf("organism + %s ~ %s", xcol, scol)), value.var=c("m","err_minus","err_plus")), paste0(f, "_wide")) }
sx <- list(scale_x_continuous(breaks=c(0, 0.4, 0.8)), labs(x="Induced |Δf|"))
pl <- function(d, ylab, cols, dodge=0.04) ggplot(d, aes(target_df, m, colour=variable, group=variable)) + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.03, linewidth=0.4, position=position_dodge(dodge)) +
  geom_line(linewidth=0.5, position=position_dodge(dodge)) + geom_point(size=1.4, position=position_dodge(dodge)) + facet_wrap(~organism) + scale_colour_manual(values=cols) + labs(x="Induced |Δf|", y=ylab, colour=NULL) + th

## ================= Figure 3: synthetic v2 =================
S <- rbindlist(lapply(list.files("results/revision/synthetic", pattern="^[A-Z]+[0-9]+[.]tsv$", full.names=TRUE), fread), fill=TRUE); S[, organism := factor(ORG[organism], levels=ORG)]
BV <- rbindlist(lapply(list.files("results/revision/synthetic", pattern="_biasvar[.]tsv$", full.names=TRUE), fread), fill=TRUE); BV[, organism := factor(ORG[organism], levels=ORG)]
per <- S[, .(dA=mean(J_M1 - J_M0perm), dW=mean(J_M1 - J_M0permw), f0=mean(1 - M0_precision, na.rm=TRUE), f1=mean(1 - M1_precision, na.rm=TRUE), r0=mean(M0_recall_C, na.rm=TRUE), r1=mean(M1_recall_C, na.rm=TRUE),
             t0=mean(M0_t_sp_C, na.rm=TRUE), t1=mean(M1_t_sp_C, na.rm=TRUE)), by=.(organism, study, target_df)]
summ <- function(v, nm) per[, { z <- bci(get(v)); .(variable=nm, n=.N, m=z[1], lo=z[2], hi=z[3]) }, by=.(organism, target_df)]
## (A)
dA <- rbind(summ("dA", "vs. across-group control"), summ("dW", "vs. within-group control"))[, variable := factor(variable, levels=c("vs. across-group control", "vs. within-group control"))]
pA <- ggplot(dA, aes(target_df, m, colour=variable, shape=variable, group=variable)) + geom_hline(yintercept=0, linetype=2, linewidth=0.4, colour="grey45") + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.03, linewidth=0.4) + geom_line(linewidth=0.5) + geom_point(size=1.6, fill="white") +
  facet_wrap(~organism) + scale_colour_manual(values=c(GRY, BLK)) + scale_shape_manual(values=c(16, 21)) + labs(x="Induced |Δf|", y="Δ Jaccard, real label minus control", colour=NULL, shape=NULL, tag="(A)") + th + tg + guides(colour=guide_legend(nrow=2), shape=guide_legend(nrow=2))
sv(pA, "Fig3A_synthetic_delta_two_controls", 3.6, 3.0); src(dA, "Fig3A_synthetic_delta_two_controls"); wide_err(dA, "Fig3A_synthetic_delta_two_controls", "target_df", "variable")
## (B) fixed full-sample reference A (the out-of-sample balanced reference has too few samples at the tail steps); (C)(D) out-of-sample balanced reference C
mk <- function(v0, v1) rbind(summ(v0, "Unadjusted"), summ(v1, "Sex-adjusted"))[, variable := factor(variable, levels=c("Unadjusted", "Sex-adjusted"))]
dB <- mk("f0", "f1"); dC <- mk("r0", "r1"); dD <- mk("t0", "t1")
pB <- pl(dB, "False-DE fraction\n(absent from full-data reference)", c(DGRY, RED)) + labs(tag="(B)") + tg + guides(colour=guide_legend(nrow=2)) + sx
pC <- pl(dC, "Recall of out-of-sample\nbalanced reference DE set", c(DGRY, RED)) + labs(tag="(C)") + tg + theme(legend.position="none") + sx
pD <- pl(dD, "Spearman r of t-statistics with\nout-of-sample balanced reference", c(DGRY, RED)) + labs(tag="(D)") + tg + theme(legend.position="none") + sx
for (z in list(list(dB, "Fig3B_synthetic_false_de_fraction_refA"), list(dC, "Fig3C_synthetic_recall_refC"), list(dD, "Fig3D_synthetic_t_correlation_refC"))) { src(z[[1]], z[[2]]); wide_err(z[[1]], z[[2]], "target_df", "variable") }
sv(pB, "Fig3B_synthetic_false_de_fraction_refA", 3.6, 3.0); sv(pC, "Fig3C_synthetic_recall_refC", 3.6, 3.0); sv(pD, "Fig3D_synthetic_t_correlation_refC", 3.6, 3.0)
## (E) bias² and variance (reference A)
bv <- BV[, .(bias2=mean(bias2 - variance/10), variance=mean(variance)), by=.(organism, study, target_df, model)]   # correct the bias² floor (variance/10)
dE <- rbind(bv[, { z <- bci(bias2); .(component="Squared bias", n=.N, m=z[1], lo=z[2], hi=z[3]) }, by=.(organism, target_df, model)], bv[, { z <- bci(variance); .(component="Variance", n=.N, m=z[1], lo=z[2], hi=z[3]) }, by=.(organism, target_df, model)])
dE[, model := factor(fifelse(model == "M0", "Unadjusted", "Sex-adjusted"), levels=c("Unadjusted", "Sex-adjusted"))]; dE[, component := factor(component, levels=c("Squared bias", "Variance"))]
pE <- ggplot(dE, aes(target_df, m, colour=model, linetype=component, group=interaction(model, component))) + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.03, linewidth=0.4, show.legend=FALSE) + geom_line(linewidth=0.5) + geom_point(size=1.3, show.legend=FALSE) +
  facet_wrap(~organism, scales="free_y") + scale_colour_manual(values=c(DGRY, RED)) + scale_linetype_manual(values=c(2, 1)) + labs(x="Induced |Δf|", y="log-fold-change error vs reference\n(mean over reference DE genes)", colour=NULL, linetype=NULL, tag="(E)") + th + tg + sx +
  guides(colour="none", linetype=guide_legend(nrow=1, override.aes=list(colour="black", linewidth=0.6))) + theme(legend.key.width=unit(20, "pt"))
sv(pE, "Fig3E_synthetic_bias_variance", 3.6, 3.0); src(dE, "Fig3E_synthetic_bias_variance"); dEw <- copy(dE)[, variable := paste(model, component, sep=" / ")]; wide_err(dEw, "Fig3E_synthetic_bias_variance", "target_df", "variable")
## (F) leakage: share of draws pooled over studies
dF <- melt(S[, .(Unadjusted=mean(M0_sexgene_de > 0), `Sex-adjusted`=mean(M1_sexgene_de > 0), n_draws=.N), by=.(organism, target_df)], id.vars=c("organism","target_df","n_draws"))
pF <- ggplot(dF, aes(target_df, value, colour=variable, group=variable)) + geom_line(linewidth=0.5) + geom_point(size=1.4) + facet_wrap(~organism) + scale_colour_manual(values=c(DGRY, RED)) + scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
  labs(x="Induced |Δf|", y="Draws with sex-chromosome\ngenes among DE genes", colour=NULL, tag="(F)") + th + tg + theme(legend.position="none") + sx
sv(pF, "Fig3F_synthetic_leakage", 3.6, 3.0); src(dF, "Fig3F_synthetic_leakage"); kn(dcast(dF[, pct := round(100*value, 1)], organism + target_df ~ variable, value.var="pct"), "Fig3F_synthetic_leakage_wide")
F3 <- ((pA | pB | pC) + plot_layout(widths=c(1.35, 1, 1))) / ((pD | pE | pF) + plot_layout(widths=c(1, 1.35, 1))) + plot_layout(heights=c(1.08, 1))
sv(F3, "Fig3_composite", 7.1, 6.3)

## ================= Figure 4 =================
cls <- rbindlist(lapply(names(ORG), function(o) fread(sprintf("results/main_%s_summary_v2.tsv", o))[, .(study, design_sex, sv_sex_r2, valid, analysis=study %in% readLines(sprintf("config/c2_analysis_set_%s.txt", o)), organism=ORG[o])]))
C <- fread(CFILE); for (k in c("design_sex","sex_diff","tissue_category")) if (k %in% names(C)) C[, (k) := NULL]
C <- merge(C, cls[, .(study, design_sex)], by="study")[design_sex %in% lv]; C[, organism := factor(ORG[organism], levels=ORG)]; C[, cls := factor(lab1[design_sex], levels=lab1)]
## (A) permutation-distribution calls
a4 <- C[, .(n=.N, B=as.numeric(median(B_within)), `Below all within-group permutations`=mean(p_J_within <= 1/(B_within+1) + 1e-9, na.rm=TRUE), `Below all across-group permutations`=mean(p_J_across <= 1/(B_across+1) + 1e-9, na.rm=TRUE)), by=.(organism, cls)][order(organism, cls)]
a4l <- melt(a4, id.vars=c("organism","cls","n","B")); a4l[, variable := factor(variable, levels=c("Below all within-group permutations", "Below all across-group permutations"))]; a4l[, null_floor := 1/(B+1)]
thx4 <- theme(axis.text.x=element_text(angle=35, hjust=1, vjust=1))
q4a <- ggplot(a4l, aes(cls, value, fill=variable)) + geom_col(position=position_dodge(0.8), width=0.75) + geom_text(aes(y=pmax(value, null_floor + 0.01), label=sprintf("%.0f", 100*value)), position=position_dodge(0.8), vjust=-0.4, size=2.2) +   # lift labels above the null dashed line so they do not overlap it
  geom_errorbar(data=unique(a4l[, .(organism, cls, null_floor)]), mapping=aes(x=cls, ymin=null_floor, ymax=null_floor), width=0.8, linetype=2, linewidth=0.4, colour="grey45", inherit.aes=FALSE) +
  facet_wrap(~organism) + scale_fill_manual(values=c(BLK, GRY)) + scale_y_continuous(labels=scales::percent, expand=c(0,0)) + coord_cartesian(ylim=c(0, 0.7)) +
  labs(x=NULL, y="Studies whose real-label Jaccard falls\nbelow every permutation (p = 1/(B + 1))", fill=NULL, tag="(A)") + th + thx4 + tg + guides(fill=guide_legend(ncol=1))
sv(q4a, "Fig4A_permutation_calls", 3.6, 3.2); src(a4, "Fig4A_permutation_calls"); kn(copy(a4)[, `:=`(within_pct=round(100*`Below all within-group permutations`,1), across_pct=round(100*`Below all across-group permutations`,1), null_pct=round(100/(B+1),1))][, .(organism, cls, n, B, within_pct, across_pct, null_pct)], "Fig4A_permutation_calls_wide")
## (B)(C) from the v3 source data (single seeded permutation), (D)(E) unchanged
v3 <- "results/figures/draft_v3/source_data"
c3 <- fread(file.path(v3, "Fig4A_pathway_conclusion_changed.tsv")); c3l <- melt(c3, id.vars=c("organism","cls","n")); c3l[, variable := factor(variable, levels=unique(variable))]
c3l[, cls := factor(lab1[lv][match(cls, c("Balanced\n|Δf| < 0.15", "Mild\n0.15–0.3", "Confounded\n≥ 0.3"))], levels=lab1)]; c3l[, organism := factor(organism, levels=ORG)]
q4b <- ggplot(c3l, aes(cls, value, fill=variable)) + geom_col(position=position_dodge(0.8), width=0.75, colour=BLK, linewidth=0.3) + geom_text(aes(label=sprintf("%.0f", 100*value)), position=position_dodge(0.8), vjust=-0.4, size=2.2) +
  facet_wrap(~organism) + scale_fill_manual(values=c(GRY, BLK, "white"), labels=c("Real sex vs. across-group control", "Real sex vs. within-group control", "Within-group label vs. across-group control")) + scale_y_continuous(labels=scales::percent, expand=c(0,0)) + coord_cartesian(ylim=c(0,0.45)) +
  labs(x=NULL, y="Studies whose significant Hallmark\nset changed beyond the control", fill=NULL, tag="(B)") + th + thx4 + tg + guides(fill=guide_legend(ncol=1))
t3 <- fread(file.path(v3, "Fig4B_top_pathway_changed.tsv")); t3l <- melt(t3, id.vars=c("organism","cls","n")); t3l[, variable := factor(variable, levels=unique(variable))]
t3l[, cls := factor(lab1[lv][match(cls, c("Balanced\n|Δf| < 0.15", "Mild\n0.15–0.3", "Confounded\n≥ 0.3"))], levels=lab1)]; t3l[, organism := factor(organism, levels=ORG)]
q4c <- ggplot(t3l, aes(cls, value, fill=variable)) + geom_col(position=position_dodge(0.8), width=0.75) + geom_text(aes(label=sprintf("%.0f", 100*value)), position=position_dodge(0.8), vjust=-0.4, size=2.2) +
  facet_wrap(~organism) + scale_fill_manual(values=c(RED, BLK, GRY)) + scale_y_continuous(labels=scales::percent, expand=c(0,0)) + coord_cartesian(ylim=c(0,0.45)) + labs(x=NULL, y="Top-ranked Hallmark pathway changed", fill=NULL, tag="(C)") + th + thx4 + tg + guides(fill=guide_legend(ncol=1))
Sv <- fread(file.path(v3, "Fig4C_sva_capture.tsv"))[, cls := factor(lab1[design_sex], levels=lab1)][, organism := factor(organism, levels=ORG)]
q4d <- ggplot(Sv, aes(cls, frac)) + geom_col(width=0.62, fill=DGRY) + geom_text(aes(label=sprintf("%.0f%%", 100*frac)), vjust=-0.4, size=2.4) + facet_wrap(~organism) + scale_y_continuous(labels=scales::percent, limits=c(0, 0.62), expand=c(0,0)) +
  labs(x=NULL, y="Studies in which a surrogate\nvariable captures sex (R² > 0.5)", tag="(D)") + th + thx4 + tg
G <- fread(file.path(v3, "Fig4D_publication_audit.tsv")); IT <- c("Sex in DE model", "Sex in DE model — confounded designs only", "Age in DE model", "Batch / latent variable", "Group sex composition reported"); IT4 <- c("Sex in model", "Sex in model —\nconfounded designs", "Age in model", "Batch / latent variable", "Group sex reported")
G[, item := factor(IT4[match(item, IT)], levels=rev(IT4))]; G[, lab := sprintf("%.1f%% (%d / %d)", 100*frac, k, n)]; G[, organism := factor(organism, levels=ORG)]
q4e <- ggplot(G, aes(item, frac)) + geom_col(width=0.62, fill=DGRY) + geom_text(aes(label=lab), hjust=-0.08, size=2.2) + coord_flip() + facet_wrap(~organism) + scale_y_continuous(labels=scales::percent, breaks=seq(0,0.6,0.3), limits=c(0, 1.05), expand=c(0,0)) +
  labs(x=NULL, y="Source publications", tag="(E)") + th + tg + theme(panel.spacing.x=unit(8,"pt"))
for (z in list(list(c3, "Fig4B_pathway_conclusion_changed"), list(t3, "Fig4C_top_pathway_changed"), list(Sv, "Fig4D_sva_capture"), list(G[, .(organism, item=gsub("\n"," ",item), k, n, frac)], "Fig4E_publication_audit"))) src(z[[1]], z[[2]])
sv(q4b, "Fig4B_pathway_conclusion_changed", 3.6, 3.2); sv(q4c, "Fig4C_top_pathway_changed", 3.6, 3.2); sv(q4d, "Fig4D_sva_capture", 3.6, 3.0); sv(q4e, "Fig4E_publication_audit", 3.6, 3.0)
F4 <- wrap_elements(full=((q4a | q4b) / (q4c | q4d) + plot_layout(heights=c(1.1, 1)))) / wrap_elements(full=q4e) + plot_layout(heights=c(2.2, 0.9))
sv(F4, "Fig4_composite", 7.1, 8.6)
## Keynote copies of v3 tables for B–E
for (f in c("Fig4A_pathway_conclusion_changed_wide","Fig4B_top_pathway_changed_wide","Fig4C_sva_capture_wide","Fig4D_publication_audit_wide")) { p <- file.path("results/figures/draft_v3/source_data_keynote", paste0(f, ".tsv")); if (file.exists(p)) file.copy(p, file.path(KD, sub("^Fig4A", "Fig4B", sub("^Fig4B", "Fig4C", sub("^Fig4C", "Fig4D", sub("^Fig4D", "Fig4E", f))))), overwrite=TRUE) }
cat("content file:", CFILE, " studies:", nrow(C), "\n"); cat("panels:", length(list.files(PD, pattern="pdf$")), " source:", length(list.files(SD)), " keynote:", length(list.files(KD)), "\n")
