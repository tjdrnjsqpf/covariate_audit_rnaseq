#!/usr/bin/env Rscript
## Draft figures for the restructured manuscript (2026-09-22): per-panel PDF (vector, for Keynote) + PNG + per-panel source data TSV.
## Colour rules: red = real sex, light grey = full permutation (degrees-of-freedom cost), black = within-group permutation (df + collinearity cost), dark grey = descriptive statistics.
## Run (Mac): LC_ALL=en_US.UTF-8 COV_AUDIT_PROJ="$(pwd)" Rscript scripts/05_report/26_draft_figures_v3.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
PROJ <- Sys.getenv("COV_AUDIT_PROJ", "/var2/lsg/Claude_Code/covariate_audit_rnaseq"); setwd(PROJ)
FONT <- if (Sys.info()[["sysname"]] == "Darwin") "Helvetica Neue" else "sans"
PD <- "results/figures/draft_v3/panels"; SD <- "results/figures/draft_v3/source_data"; dir.create(PD, showWarnings=FALSE, recursive=TRUE); dir.create(SD, showWarnings=FALSE, recursive=TRUE)
RED <- "#B02418"; GRY <- "#C4C4C4"; BLK <- "#1A1A1A"; DGRY <- "#5F5F5F"; LRED <- "#E3A29B"
LAB <- c("Real sex", "Permuted across groups\n(df cost)", "Permuted within groups\n(df + collinearity cost)")
COL3 <- setNames(c(RED, GRY, BLK), LAB)
lv <- c("balanced","mildly_imbalanced","partially_confounded"); lab <- c(balanced="Balanced\n|Δf| < 0.15", mildly_imbalanced="Mild\n0.15–0.3", partially_confounded="Confounded\n≥ 0.3")
ORG <- c(human="Human", mouse="Mouse")
th <- theme_classic(base_size=9, base_family=FONT) + theme(axis.line=element_line(colour="black", linewidth=0.4), axis.ticks=element_line(colour="black", linewidth=0.4), axis.text=element_text(colour="black"),
  legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"), legend.key.size=unit(9,"pt"), legend.text=element_text(size=8),
  strip.background=element_blank(), strip.text=element_text(face="bold", size=9, hjust=0), panel.spacing.x=unit(14,"pt"), plot.margin=margin(4,8,4,4))
update_geom_defaults("text", list(family=FONT))
bci <- function(x, R=2000) { x <- x[!is.na(x)]; if (length(x) < 3) return(c(median(x), NA, NA)); set.seed(1); b <- replicate(R, median(sample(x, replace=TRUE))); c(median(x), unname(quantile(b,.025)), unname(quantile(b,.975))) }
sv <- function(p, f, w, h) { pdfdev <- if (Sys.info()[["sysname"]] == "Darwin") function(filename, width, height, ...) grDevices::quartz(file=filename, type="pdf", width=width, height=height) else cairo_pdf
  ggsave(file.path(PD, paste0(f, ".pdf")), p, width=w, height=h, device=pdfdev); tmp <- tempfile(fileext=".png")
  if (requireNamespace("ragg", quietly=TRUE)) ggsave(tmp, p, width=w, height=h, dpi=300, device=ragg::agg_png) else ggsave(tmp, p, width=w, height=h, dpi=300); invisible(file.copy(tmp, file.path(PD, paste0(f, ".png")), overwrite=TRUE)) }
src <- function(d, f) fwrite(d, file.path(SD, paste0(f, ".tsv")), sep="\t")

## ---------- data ----------
D <- rbindlist(lapply(names(ORG), function(o) { d <- fread(sprintf("results/main_%s_summary_v2.tsv", o)); A <- readLines(sprintf("config/c2_analysis_set_%s.txt", o))
  g <- fread(sprintf("results/c2/gsea_summary_%s.tsv", o))[, .(study, P1_H_jaccard, P0perm_H_jaccard, P0permw_H_jaccard, P1_H_top1_same, P0perm_H_top1_same, P0permw_H_top1_same)]
  d <- merge(d, g, by="study", all.x=TRUE); d[, organism := ORG[o]]; d[, analysis := study %in% A]; d }), fill=TRUE)
D[, organism := factor(organism, levels=ORG)]
A2 <- D[analysis == TRUE & design_sex %in% lv & !is.na(M1_jaccard) & !is.na(M0perm_jaccard) & !is.na(M0permw_jaccard)]
A2[, cls := factor(lab[design_sex], levels=lab[lv])]
S <- fread("results/revision/S_per_study_rung.tsv")[, organism := factor(ORG[organism], levels=ORG)]

## ================= Figure 1: prevalence, and the two costs of a term =================
## 1A |Δf| distribution
V <- D[valid == TRUE & single_sex == FALSE & !is.na(sex_diff)]; sh <- V[, .(n=.N, p=mean(design_sex == "partially_confounded")), by=organism]
## 1B: count the 20 bins [0,0.05) … [0.95,1.00] directly via integer indices (avoids floating-point boundary issues) and draw as bars — source table and figure share the same numbers
V[, bi := pmin(floor(sex_diff / 0.05 + 1e-9), 19)]
BINS <- data.table(bi=0:19, bin=sprintf("%.2f\u2013%.2f", seq(0, 0.95, 0.05), seq(0.05, 1, 0.05)), xmid=seq(0.025, 0.975, 0.05))
HB <- merge(CJ(bi=0:19, organism=levels(V$organism)), V[, .N, by=.(bi, organism)], by=c("bi","organism"), all.x=TRUE)[is.na(N), N := 0L]
HB <- merge(HB, BINS, by="bi")[, organism := factor(organism, levels=ORG)][order(organism, bi)]
p1a <- ggplot(HB, aes(xmid, N)) + geom_col(width=0.05, fill=DGRY, colour="white", linewidth=0.3) + geom_vline(xintercept=c(0.15, 0.3), linetype=2, linewidth=0.4, colour="grey45") +
  geom_text(data=sh, aes(x=1.0, y=Inf, label=sprintf("confounded: %.1f%% of %d", 100*p, n)), vjust=1.8, hjust=1, size=2.6, inherit.aes=FALSE) +
  facet_wrap(~organism, scales="free_y", ncol=1) + scale_x_continuous(breaks=seq(0, 1, 0.25)) + labs(x="|\u0394f| between compared groups", y="Studies") + th + theme(strip.text=element_text(face="bold", size=9, hjust=0))
sv(p1a, "Fig1B_imbalance_distribution", 3.4, 3.6)   # meant to sit vertically next to the Fig 1 schematic (Keynote)
src(V[, .(study, organism, sex_diff, design_sex, n, n_case, n_control, n_female, n_male)], "Fig1B_imbalance_distribution")
src(dcast(HB, bin ~ organism, value.var="N")[match(BINS$bin, bin)], "Fig1B_histogram_bins")   # 20-bin count table identical to the figure (for Keynote)

## 1B Jaccard by class: real / across / within
long <- rbind(A2[, .(study, organism, cls, type=LAB[1], J=M1_jaccard)], A2[, .(study, organism, cls, type=LAB[2], J=M0perm_jaccard)], A2[, .(study, organism, cls, type=LAB[3], J=M0permw_jaccard)])
long[, type := factor(type, levels=LAB)]; nn <- A2[, .N, by=.(organism, cls)][order(organism, cls)]
p1b <- ggplot(long, aes(cls, J, fill=type)) + geom_boxplot(outlier.size=0.35, outlier.alpha=0.6, width=0.75, position=position_dodge(0.85), linewidth=0.3, colour=BLK) +
  facet_wrap(~organism, labeller=as_labeller(setNames(nn[, sprintf("%s (n = %s)", organism[1], paste(N, collapse=" / ")), by=organism]$V1, levels(nn$organism)))) +
  scale_fill_manual(values=COL3) + coord_cartesian(ylim=c(0,1)) + labs(x="Sex imbalance class", y="Jaccard, unadjusted vs adjusted DE set", fill=NULL) + th
sv(p1b, "Fig2A_jaccard_three_controls", 7.1, 3.2); src(dcast(long, study + organism + cls ~ type, value.var="J"), "Fig2A_jaccard_three_controls")

## 1C dose–response: Δ vs across and vs within, by |Δf| bin
A3 <- copy(A2)[, bin := cut(sex_diff, c(-0.01,0.1,0.2,0.3,0.4,0.5,0.7,1), labels=c("0–0.1","0.1–0.2","0.2–0.3","0.3–0.4","0.4–0.5","0.5–0.7","0.7–1"))]
dr <- rbind(A3[, { z <- bci(M1_jaccard - M0perm_jaccard); .(control="vs. permuted across groups", n=.N, m=z[1], lo=z[2], hi=z[3]) }, by=.(organism, bin)],
            A3[, { z <- bci(M1_jaccard - M0permw_jaccard); .(control="vs. permuted within groups", n=.N, m=z[1], lo=z[2], hi=z[3]) }, by=.(organism, bin)])
nb <- dr[control == dr$control[1], .(organism, bin, n)]   # number of studies per bin (shared by both controls)
p1c <- ggplot(dr, aes(bin, m, colour=control, group=control, shape=control)) + geom_hline(yintercept=0, linetype=2, linewidth=0.4, colour="grey45") +
  geom_linerange(aes(ymin=lo, ymax=hi), linewidth=0.4, position=position_dodge(0.4)) + geom_line(data=dr[n >= 3], linewidth=0.5, position=position_dodge(0.4)) + geom_point(data=dr[n >= 3], size=1.7, position=position_dodge(0.4), fill="white") +
  geom_point(data=dr[n < 3], size=1.7, shape=1, position=position_dodge(0.4), show.legend=FALSE) +   # n < 3: hollow points, no line or interval
  geom_text(data=nb, aes(bin, -Inf, label=n), inherit.aes=FALSE, vjust=-0.4, size=2.2, colour="grey30") +
  facet_wrap(~organism) + scale_colour_manual(values=c(GRY, BLK)) + scale_shape_manual(values=c(16, 21)) + labs(x="|Δf| bin (number of studies below the axis)", y="Δ Jaccard (real sex − control)", colour=NULL, shape=NULL) + scale_y_continuous(expand=expansion(mult=c(0.18, 0.05))) + th   # bottom margin: keep the n labels from overlapping the points
sv(p1c, "Fig2B_dose_response_two_controls", 7.1, 2.8); src(dr[order(organism, control, bin)], "Fig2B_dose_response_two_controls")

## Figure 2 composite (to check the Keynote layout; font base_size 9 as in Fig 1B, width 7.1 in)
tg <- theme(plot.tag=element_text(face="bold", size=11, family=FONT), plot.tag.position=c(0, 1))
F2 <- (p1b + labs(tag="(A)") + tg) / (p1c + labs(tag="(B)") + tg) + plot_layout(heights=c(1.15, 1))
sv(F2, "Fig2_composite", 7.1, 6.0)

## ================= Figure 2: synthetic imbalance within studies =================
rung <- function(cols, labs_) { x <- melt(S[, c("organism","study","target_df", cols), with=FALSE], id.vars=c("organism","study","target_df")); x[, variable := factor(labs_[as.character(variable)], levels=labs_)]
  x[, { z <- bci(value); .(n=.N, m=z[1], lo=z[2], hi=z[3]) }, by=.(organism, target_df, variable)] }
pl <- function(dd, ylab, cols, shapes=c(16,21), hline=NULL) { p <- ggplot(dd, aes(target_df, m, colour=variable, group=variable, shape=variable)) + { if (!is.null(hline)) geom_hline(yintercept=hline, linetype=2, linewidth=0.4, colour="grey45") } +
  geom_linerange(aes(ymin=lo, ymax=hi), linewidth=0.4, position=position_dodge(0.04)) + geom_line(linewidth=0.5, position=position_dodge(0.04)) + geom_point(size=1.7, fill="white", position=position_dodge(0.04)) +
  facet_wrap(~organism) + scale_colour_manual(values=cols) + scale_shape_manual(values=shapes) + scale_x_continuous(breaks=c(0,.2,.4,.6,.8)) + labs(x="Induced |Δf| (fixed n per arm)", y=ylab, colour=NULL, shape=NULL) + th; p }
d2a <- rung(c("d_across","d_within"), c(d_across="vs. permuted across groups", d_within="vs. permuted within groups"))
sv(pl(d2a, "Δ Jaccard (real sex − control)", c(GRY, BLK), hline=0), "Fig3A_synthetic_delta_two_controls", 7.1, 2.6); src(d2a, "Fig3A_synthetic_delta_two_controls")
d2b <- rung(c("fd_M0","fd_M1"), c(fd_M0="Unadjusted", fd_M1="Sex-adjusted"))
sv(pl(d2b, "DE genes absent from balanced reference\n(false-DE fraction)", c(DGRY, RED)), "Fig3B_synthetic_false_de_fraction", 7.1, 2.6); src(d2b, "Fig3B_synthetic_false_de_fraction")
d2c <- rung(c("M0_recall","M1_recall"), c(M0_recall="Unadjusted", M1_recall="Sex-adjusted"))
sv(pl(d2c, "Recall of balanced-reference DE set", c(DGRY, RED)), "Fig3C_synthetic_recall", 7.1, 2.6); src(d2c, "Fig3C_synthetic_recall")
d2d <- rung(c("M0_t_sp_ref","M1_t_sp_ref"), c(M0_t_sp_ref="Unadjusted", M1_t_sp_ref="Sex-adjusted"))
sv(pl(d2d, "Spearman r of t-statistics\nwith balanced reference", c(DGRY, RED)), "Fig3D_synthetic_t_correlation", 7.1, 2.6); src(d2d, "Fig3D_synthetic_t_correlation")
lk <- S[, .(Unadjusted=mean(M0_sexgene_de > 0), `Sex-adjusted`=mean(M1_sexgene_de > 0), n=.N), by=.(organism, target_df)]; lkl <- melt(lk, id.vars=c("organism","target_df","n"), variable.name="variable", value.name="m")[, `:=`(lo=NA_real_, hi=NA_real_)]
sv(pl(lkl, "Studies with ≥ 1 sex-chromosome gene\namong DE genes", c(DGRY, RED)) + scale_y_continuous(labels=scales::percent), "Fig3E_synthetic_leakage", 7.1, 2.6); src(lk, "Fig3E_synthetic_leakage")
src(S, "Fig3_synthetic_per_study_rung_means")

## Figure 3 composite: A full width / B|C / D|E, legends only on A and B (C–E use the same series as B)
pA <- pl(d2a, "\u0394 Jaccard (real sex \u2212 control)", c(GRY, BLK), hline=0) + labs(tag="(A)") + tg
pB <- pl(d2b, "False-DE fraction\n(DE genes absent from reference)", c(DGRY, RED)) + labs(tag="(B)") + tg
pC <- pl(d2c, "Recall of reference DE set", c(DGRY, RED)) + labs(tag="(C)") + tg + theme(legend.position="none")
pD <- pl(d2d, "Spearman r of t-statistics\nwith reference", c(DGRY, RED)) + labs(tag="(D)") + tg + theme(legend.position="none")
pE <- pl(lkl, "Studies with sex-chromosome\ngenes among DE genes", c(DGRY, RED)) + scale_y_continuous(labels=scales::percent) + labs(tag="(E)") + tg + theme(legend.position="none")
## Panel margins: widen the gap between tag and panel and between panels
tg3 <- theme(plot.tag=element_text(face="bold", size=11, family=FONT), plot.tag.position=c(0, 1), plot.margin=margin(16, 14, 6, 8))
## Top row: A legend in two rows with short names, B legend also in two rows — so they do not intrude on the neighbouring panel tag
pA <- pA + tg3 + scale_colour_manual(values=c(GRY, BLK), labels=c("vs. across-group control", "vs. within-group control")) + scale_shape_manual(values=c(16, 21), labels=c("vs. across-group control", "vs. within-group control")) + guides(colour=guide_legend(nrow=2), shape=guide_legend(nrow=2))
pB <- pB + tg3 + guides(colour=guide_legend(nrow=2), shape=guide_legend(nrow=2))
## Bottom row: panels are narrow, so x ticks at 0 / 0.4 / 0.8 and short axis titles
sx <- list(scale_x_continuous(breaks=c(0, 0.4, 0.8)), labs(x="Induced |\u0394f|"))
pC <- pC + tg3 + sx; pD <- pD + tg3 + sx; pE <- pE + tg3 + sx
## Widths: top row A:B = 2:1, bottom row C:D:E = 1:1:1 → B = C = D = E, A twice as wide
pB <- pB + sx   # B also becomes narrow, so use the same axis settings as the bottom row
F3 <- ((pA | pB) + plot_layout(widths=c(2, 1))) / ((pC | pD | pE) + plot_layout(widths=c(1, 1, 1))) + plot_layout(heights=c(1.05, 1))
sv(F3, "Fig3_composite", 7.1, 6.0)

## ================= Figure 3: pathway level, SVA, publications =================
P <- A2[!is.na(P1_H_jaccard) & !is.na(P0perm_H_jaccard) & !is.na(P0permw_H_jaccard)]
c3 <- P[order(organism, cls)][, .(n=.N, `Real sex vs. across-group control`=mean(P1_H_jaccard < P0perm_H_jaccard - 0.2), `Real sex vs. within-group control`=mean(P1_H_jaccard < P0permw_H_jaccard - 0.2), `Within-group label vs. across-group control`=mean(P0permw_H_jaccard < P0perm_H_jaccard - 0.2)), by=.(organism, cls)]
c3l <- melt(c3, id.vars=c("organism","cls","n")); c3l[, variable := factor(variable, levels=unique(variable))]
p3a <- ggplot(c3l, aes(cls, value, fill=variable)) + geom_col(position=position_dodge(0.8), width=0.75) + geom_text(aes(label=sprintf("%.0f", 100*value)), position=position_dodge(0.8), vjust=-0.4, size=2.3) +
  facet_wrap(~organism, labeller=as_labeller(setNames(c3[, sprintf("%s (n = %s)", organism[1], paste(n, collapse=" / ")), by=organism]$V1, levels(c3$organism)))) +
  scale_fill_manual(values=c(GRY, BLK, LRED)) + scale_y_continuous(labels=scales::percent, expand=c(0,0)) + coord_cartesian(ylim=c(0,0.45)) + labs(x="Sex imbalance class", y="Studies whose significant Hallmark\nset changed beyond the control", fill=NULL) + th + guides(fill=guide_legend(ncol=1))
sv(p3a, "Fig4A_pathway_conclusion_changed", 7.1, 3.0); src(c3, "Fig4A_pathway_conclusion_changed")
t3 <- P[order(organism, cls)][, .(n=sum(!is.na(P1_H_top1_same)), `Real sex`=mean(!P1_H_top1_same, na.rm=TRUE), `Permuted within groups`=mean(!P0permw_H_top1_same, na.rm=TRUE), `Permuted across groups`=mean(!P0perm_H_top1_same, na.rm=TRUE)), by=.(organism, cls)]
t3l <- melt(t3, id.vars=c("organism","cls","n")); t3l[, variable := factor(variable, levels=unique(variable))]
p3b <- ggplot(t3l, aes(cls, value, fill=variable)) + geom_col(position=position_dodge(0.8), width=0.75) + geom_text(aes(label=sprintf("%.0f", 100*value)), position=position_dodge(0.8), vjust=-0.4, size=2.3) +
  facet_wrap(~organism) + scale_fill_manual(values=c(RED, BLK, GRY)) + scale_y_continuous(labels=scales::percent, expand=c(0,0)) + coord_cartesian(ylim=c(0,0.45)) + labs(x="Sex imbalance class", y="Top-ranked Hallmark pathway changed", fill=NULL) + th
sv(p3b, "Fig4B_top_pathway_changed", 7.1, 3.0); src(t3, "Fig4B_top_pathway_changed")
Sv <- D[analysis == TRUE & !is.na(sv_sex_r2) & design_sex %in% lv, .(frac=mean(sv_sex_r2 > 0.5), n=.N), by=.(organism, design_sex)][, cls := factor(lab[design_sex], levels=lab[lv])]
p3c <- ggplot(Sv, aes(cls, frac)) + geom_col(width=0.62, fill=DGRY) + geom_text(aes(label=sprintf("%.0f%%", 100*frac)), vjust=-0.4, size=2.5) + facet_wrap(~organism, labeller=as_labeller(setNames(Sv[order(cls), sprintf("%s (n = %s)", organism[1], paste(n, collapse=" / ")), by=organism]$V1, levels(Sv$organism)))) +
  scale_y_continuous(labels=scales::percent, limits=c(0, 0.62), expand=c(0,0)) + labs(x="Sex imbalance class", y="Studies in which a surrogate\nvariable captures sex (R² > 0.5)") + th
sv(p3c, "Fig4C_sva_capture", 7.1, 2.8); src(Sv[, .(organism, design_sex, n, frac)], "Fig4C_sva_capture")
IT <- c("Sex in DE model", "Sex in DE model —\nconfounded designs only", "Age in DE model", "Batch / latent variable", "Group sex composition reported")
G <- rbindlist(lapply(names(ORG), function(o) { g <- fread(sprintf("results/gate2/gate2_merged_%s.tsv", o))[, .(study, text_relevant, sex_adjusted, age_adjusted, batch_adj, sex_reported)]
  g <- merge(g, fread(sprintf("results/main_%s_summary_v2.tsv", o))[valid == TRUE, .(study, design_sex)], by="study", suffixes=c("", "_new"))[text_relevant == TRUE]; pc <- g[design_sex == "partially_confounded"]   # 2026-09-22: only papers in the current valid set (studies dropped by the filter fix are excluded)
  data.table(organism=ORG[o], item=IT, k=c(sum(g$sex_adjusted == "yes"), sum(pc$sex_adjusted == "yes"), sum(g$age_adjusted == "yes"), sum(g$batch_adj == "yes"), sum(g$sex_reported == "yes")), n=c(nrow(g), nrow(pc), nrow(g), nrow(g), nrow(g))) }))
G[, frac := k/n]; G[, lab := sprintf("%.1f%% (%d / %d)", 100*frac, k, n)]; G[, organism := factor(organism, levels=ORG)]; G[, item := factor(item, levels=rev(IT))]
p3d <- ggplot(G, aes(item, frac)) + geom_col(width=0.62, fill=DGRY) + geom_text(aes(label=lab), hjust=-0.08, size=2.4) + coord_flip() + facet_wrap(~organism, labeller=as_labeller(setNames(G[item == IT[1], sprintf("%s (n = %d publications)", organism, n)], levels(G$organism)))) +
  scale_y_continuous(labels=scales::percent, breaks=seq(0,0.6,0.2), limits=c(0, 0.82), expand=c(0,0)) + labs(x=NULL, y="Source publications") + th
sv(p3d, "Fig4D_publication_audit", 7.1, 2.8); src(G[, .(organism, item=gsub("\n"," ",item), k, n, frac)], "Fig4D_publication_audit")

## Figure 4 composite: (A)(B) / (C)(D). Half-width panels, so class labels on one line and slanted; n removed from the titles and given in the legend
tg4 <- theme(plot.tag=element_text(face="bold", size=11, family=FONT), plot.tag.position=c(0, 1), plot.margin=margin(16, 14, 6, 8))
lab1 <- c(balanced="Balanced", mildly_imbalanced="Mild", partially_confounded="Confounded")
thx4 <- theme(axis.text.x=element_text(angle=35, hjust=1, vjust=1))
c3s <- copy(c3l)[, cls := factor(lab1[as.character(factor(cls, levels=lab[lv], labels=lv))], levels=lab1)]
q4a <- ggplot(c3s, aes(cls, value, fill=variable)) + geom_col(position=position_dodge(0.8), width=0.75) + geom_text(aes(label=sprintf("%.0f", 100*value)), position=position_dodge(0.8), vjust=-0.4, size=2.2) +
  facet_wrap(~organism) + scale_fill_manual(values=c(GRY, BLK, LRED), labels=c("Real sex vs. across-group control", "Real sex vs. within-group control", "Within-group label vs. across-group control")) + scale_y_continuous(labels=scales::percent, expand=c(0,0)) + coord_cartesian(ylim=c(0,0.45)) +
  labs(x=NULL, y="Studies whose significant Hallmark\nset changed beyond the control", fill=NULL, tag="(A)") + th + thx4 + tg4 + guides(fill=guide_legend(ncol=1))
t3s <- copy(t3l)[, cls := factor(lab1[as.character(factor(cls, levels=lab[lv], labels=lv))], levels=lab1)]
q4b <- ggplot(t3s, aes(cls, value, fill=variable)) + geom_col(position=position_dodge(0.8), width=0.75) + geom_text(aes(label=sprintf("%.0f", 100*value)), position=position_dodge(0.8), vjust=-0.4, size=2.2) +
  facet_wrap(~organism) + scale_fill_manual(values=c(RED, BLK, GRY)) + scale_y_continuous(labels=scales::percent, expand=c(0,0)) + coord_cartesian(ylim=c(0,0.45)) +
  labs(x=NULL, y="Top-ranked Hallmark pathway changed", fill=NULL, tag="(B)") + th + thx4 + tg4 + guides(fill=guide_legend(ncol=1))
Sv2 <- copy(Sv)[, cls := factor(lab1[design_sex], levels=lab1)]
q4c <- ggplot(Sv2, aes(cls, frac)) + geom_col(width=0.62, fill=DGRY) + geom_text(aes(label=sprintf("%.0f%%", 100*frac)), vjust=-0.4, size=2.4) + facet_wrap(~organism) +
  scale_y_continuous(labels=scales::percent, limits=c(0, 0.62), expand=c(0,0)) + labs(x=NULL, y="Studies in which a surrogate\nvariable captures sex (R\u00b2 > 0.5)", tag="(C)") + th + thx4 + tg4
IT4 <- c("Sex in model", "Sex in model \u2014\nconfounded designs", "Age in model", "Batch / latent variable", "Group sex reported")
G4 <- copy(G)[, item := factor(IT4[match(as.character(item), IT)], levels=rev(IT4))]
q4d <- ggplot(G4, aes(item, frac)) + geom_col(width=0.62, fill=DGRY) + geom_text(aes(label=lab), hjust=-0.08, size=2.2) + coord_flip() + facet_wrap(~organism) +
  scale_y_continuous(labels=scales::percent, breaks=seq(0,0.6,0.3), limits=c(0, 1.05), expand=c(0,0)) + labs(x=NULL, y="Source publications", tag="(D)") + th + tg4 + theme(panel.spacing.x=unit(8,"pt"))
F4 <- (q4a | q4b) / ((q4c | q4d) + plot_layout(widths=c(1, 1.6))) + plot_layout(heights=c(1.15, 1))
sv(F4, "Fig4_composite", 7.1, 6.8)

## ================= Supplementary candidates =================
L <- A2[, .(n=.N, `Unadjusted (real label)`=mean(M1_sexgene_deg_ref > 0), `Permuted within groups`=mean(M0permw_sexgene_deg_alt > 0), `Sex-adjusted`=mean(M1_sexgene_deg_alt > 0)), by=.(organism, cls)]; Ll <- melt(L, id.vars=c("organism","cls","n"))
pS1 <- ggplot(Ll, aes(cls, value, fill=variable)) + geom_col(position=position_dodge(0.8), width=0.75) + geom_text(aes(label=sprintf("%.0f", 100*value)), position=position_dodge(0.8), vjust=-0.4, size=2.3) + facet_wrap(~organism) +
  scale_fill_manual(values=c(DGRY, BLK, RED)) + scale_y_continuous(labels=scales::percent, expand=c(0,0)) + coord_cartesian(ylim=c(0,0.7)) + labs(x="Sex imbalance class", y="Studies with ≥ 1 sex-chromosome gene\namong DE genes", fill=NULL) + th
sv(pS1, "FigS1_observed_leakage", 7.1, 2.8); src(L, "FigS1_observed_leakage")
E1 <- fread("results/revision/E1_single_sex_split.tsv"); src(E1, "FigS2_single_sex_split_source_E1")
E7 <- fread("results/revision/E7_sim_table.tsv"); src(E7, "FigS3_simulation_ground_truth_source_E7")
cat("panels:", length(list.files(PD, pattern="pdf$")), " source tables:", length(list.files(SD)), "\n")
