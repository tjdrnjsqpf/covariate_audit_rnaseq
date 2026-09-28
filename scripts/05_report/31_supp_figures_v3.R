#!/usr/bin/env Rscript
## Supplementary figures S1–S8 (restructured manuscript). Input: results/revision/E*.tsv (scripts/06_revision/30–38, recomputed on chain 5 data) + summary tables.
## Run (Mac): LC_ALL=en_US.UTF-8 COV_AUDIT_PROJ="$(pwd)" Rscript scripts/05_report/31_supp_figures_v3.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
PROJ <- Sys.getenv("COV_AUDIT_PROJ", "/var2/lsg/Claude_Code/covariate_audit_rnaseq"); setwd(PROJ)
FONT <- if (Sys.info()[["sysname"]] == "Darwin") "Helvetica Neue" else "sans"
PD <- "results/figures/draft_v3/panels"; SD <- "results/figures/draft_v3/source_data"; RV <- "results/revision"
RED <- "#B02418"; GRY <- "#C4C4C4"; BLK <- "#1A1A1A"; DGRY <- "#5F5F5F"; MGRY <- "#8C8C8C"; LRED <- "#E3A29B"
ORG <- c(human="Human", mouse="Mouse"); lv <- c("balanced","mildly_imbalanced","partially_confounded"); lab1 <- c(balanced="Balanced", mildly_imbalanced="Mild", partially_confounded="Confounded")
CLS3 <- setNames(c(GRY, MGRY, BLK), lab1)   # class series use three grey levels (red is reserved for real sex adjustment)
th <- theme_classic(base_size=9, base_family=FONT) + theme(axis.line=element_line(colour="black", linewidth=0.4), axis.ticks=element_line(colour="black", linewidth=0.4), axis.text=element_text(colour="black"),
  legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"), legend.key.size=unit(9,"pt"), legend.text=element_text(size=8),
  strip.background=element_blank(), strip.text=element_text(face="bold", size=9, hjust=0), panel.spacing.x=unit(14,"pt"), plot.margin=margin(6,10,4,4), plot.tag=element_text(face="bold", size=11), plot.tag.position=c(0,1))
update_geom_defaults("text", list(family=FONT))
sv <- function(p, f, w, h) { pdfdev <- if (Sys.info()[["sysname"]] == "Darwin") function(filename, width, height, ...) grDevices::quartz(file=filename, type="pdf", width=width, height=height) else cairo_pdf
  ggsave(file.path(PD, paste0(f, ".pdf")), p, width=w, height=h, device=pdfdev); tmp <- tempfile(fileext=".png"); ggsave(tmp, p, width=w, height=h, dpi=300, device=ragg::agg_png); invisible(file.copy(tmp, file.path(PD, paste0(f, ".png")), overwrite=TRUE)) }
src <- function(d, f) fwrite(d, file.path(SD, paste0(f, ".tsv")), sep="\t")
rd <- function(f) fread(file.path(RV, f))
orgf <- function(d) { d[, organism := factor(ORG[as.character(organism)], levels=ORG)]; d }

## ---------- S2: sex direction of single-sex studies ----------
e1 <- rd("E1_single_sex_split.tsv")
s2 <- rbind(e1[stratum == "all", .(organism, group="All single-sex studies", n_single_sex, pct_female_only, lo, hi)],
            e1[organism == "human" & stratum == "human_sex_specific_tissue(title+abstract)" & level == "excluding_sex_specific", .(organism, group="Excluding sex-specific tissues", n_single_sex, pct_female_only, lo, hi)],
            e1[organism == "mouse" & stratum == "design_type" & level %in% c("genotype_comparison","in_vivo_treatment_arms"), .(organism, group=fifelse(level == "genotype_comparison", "Genotype comparisons", "Treatment arms"), n_single_sex, pct_female_only, lo, hi)])
s2 <- orgf(s2); s2[, group := factor(group, levels=c("All single-sex studies","Excluding sex-specific tissues","Genotype comparisons","Treatment arms"))]; s2[, lbl := sprintf("%.0f%%  (n = %d)", pct_female_only, n_single_sex)]
pS2 <- ggplot(s2, aes(group, pct_female_only)) + geom_hline(yintercept=50, linetype=2, linewidth=0.4, colour="grey45") + geom_col(width=0.6, fill=DGRY) + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.18, linewidth=0.4) +
  geom_text(aes(y=hi + 3, label=lbl), size=2.4, hjust=0) + coord_flip(ylim=c(0, 100)) + facet_wrap(~organism, scales="free_y", ncol=1) + scale_x_discrete(limits=rev) +
  labs(x=NULL, y="Female-only studies among single-sex studies (%)") + th
sv(pS2, "FigS2_single_sex_direction", 7.1, 3.4); src(s2[, .(organism, group, n_single_sex, pct_female_only, lo, hi)], "FigS2_single_sex_direction")

## ---------- S3: semi-synthetic simulation (ground truth) ----------
e7 <- orgf(rd("E7_sim_table.tsv")[set == "complete_ladders"])
s3 <- melt(e7[, .(organism, sex_delta, `Unadjusted`=M0_false_deg, `Sex-adjusted`=M1_false_deg)], id.vars=c("organism","sex_delta"), variable.name="model", value.name="false_deg")
s3r <- melt(e7[, .(organism, sex_delta, `Unadjusted`=M0_recall, `Sex-adjusted`=M1_recall)], id.vars=c("organism","sex_delta"), variable.name="model", value.name="recall")
pl3 <- function(d, y, ylab) ggplot(d, aes(sex_delta, .data[[y]], colour=model, shape=model, group=model)) + geom_line(linewidth=0.5) + geom_point(size=1.8, fill="white") + facet_wrap(~organism, scales="free_y") +
  scale_colour_manual(values=c(DGRY, RED)) + scale_shape_manual(values=c(16, 21)) + scale_x_continuous(breaks=c(0,.25,.5,.75,1)) + labs(x="Induced sex imbalance, |Δf|", y=ylab, colour=NULL, shape=NULL) + th
pS3 <- (pl3(s3, "false_deg", "False DE genes (not in ground truth)") + labs(tag="(A)")) / (pl3(s3r, "recall", "Recall of ground-truth DE genes") + labs(tag="(B)") + theme(legend.position="none"))
sv(pS3, "FigS3_simulation_ground_truth", 7.1, 5.0); src(e7, "FigS3_simulation_ground_truth")

## ---------- S4: stratification by sample size ----------
e2 <- orgf(rd("E2_delta_by_class_nstrata.tsv")[nstr %in% c("<12","12-19","20-39",">=40")]); e2[, nstr := factor(nstr, levels=c("<12","12-19","20-39",">=40"))]; e2[, cls := factor(lab1[cls], levels=lab1)]
pS4 <- ggplot(e2, aes(nstr, delta_median, colour=cls, group=cls)) + geom_hline(yintercept=0, linetype=2, linewidth=0.4, colour="grey45") + geom_linerange(aes(ymin=lo, ymax=hi), linewidth=0.4, position=position_dodge(0.5)) +
  geom_line(linewidth=0.5, position=position_dodge(0.5)) + geom_point(size=1.8, position=position_dodge(0.5)) + geom_text(aes(y=0.72, label=n_studies), position=position_dodge(0.5), size=2.2, show.legend=FALSE) +
  facet_wrap(~organism) + scale_colour_manual(values=CLS3) + coord_cartesian(ylim=c(-0.85, 0.78)) + labs(x="Sample size stratum (n)", y="Δ Jaccard (real sex − across-group control)", colour="Class") + th
sv(pS4, "FigS4_sample_size_strata", 7.1, 3.2); src(e2[, .(organism, nstr, cls, n_studies, n_median, delta_median, lo, hi)], "FigS4_sample_size_strata")

## ---------- S5: dependence of the 'conclusion change' call on the margin ----------
e6 <- orgf(rd("E6_conclusion_changed_margins.tsv")); e6[, cls := factor(lab1[cls], levels=lab1)]
pS5 <- ggplot(e6, aes(factor(margin), pct_changed, colour=cls, group=cls)) + geom_linerange(aes(ymin=lo, ymax=hi), linewidth=0.4, position=position_dodge(0.4)) + geom_line(linewidth=0.5, position=position_dodge(0.4)) + geom_point(size=1.8, position=position_dodge(0.4)) +
  facet_wrap(~organism) + scale_colour_manual(values=CLS3) + labs(x="Margin below the across-group control (Jaccard)", y="Studies called 'pathway conclusion changed' (%)", colour="Class") + th
sv(pS5, "FigS5_conclusion_margin_sensitivity", 7.1, 3.0); src(e6, "FigS5_conclusion_margin_sensitivity")

## ---------- S6: SVA capture — by R² threshold ----------
e5 <- orgf(rd("E5_capture_by_class_thresholds.tsv")); e5[, cls := factor(lab1[cls], levels=lab1)]
s6 <- rbindlist(lapply(c("0.3","0.5","0.7"), function(t) e5[, .(organism, cls, n_studies, threshold=paste0("R² > ", t), pct=get(paste0("capture_r2_", t)), lo=get(paste0("lo_", t)), hi=get(paste0("hi_", t)))]))
pS6 <- ggplot(s6, aes(threshold, pct, colour=cls, group=cls)) + geom_linerange(aes(ymin=lo, ymax=hi), linewidth=0.4, position=position_dodge(0.4)) + geom_line(linewidth=0.5, position=position_dodge(0.4)) + geom_point(size=1.8, position=position_dodge(0.4)) +
  facet_wrap(~organism) + scale_colour_manual(values=CLS3) + labs(x="Threshold for 'a surrogate variable captures sex'", y="Studies (%)", colour="Class") + th
sv(pS6, "FigS6_sva_threshold_sensitivity", 7.1, 3.0); src(s6, "FigS6_sva_threshold_sensitivity")

## ---------- S7: chance expectation of imbalance ----------
e8 <- rd("E8_imbalance_null_exact.tsv")[source == "recomputed, tolerance 1e-9 both sides"]; e8 <- orgf(e8)
s7 <- melt(e8[, .(organism, threshold=paste0("|Δf| ≥ ", threshold), n_studies, Observed=observed, `Expected by chance`=expected, ratio, lo=ratio_exact_lo, hi=ratio_exact_hi, p=p_two_sided)], id.vars=c("organism","threshold","n_studies","ratio","lo","hi","p"), variable.name="kind", value.name="studies")
pS7 <- ggplot(s7, aes(threshold, studies, fill=kind)) + geom_col(position=position_dodge(0.7), width=0.62) + geom_text(aes(label=round(studies)), position=position_dodge(0.7), vjust=-0.4, size=2.4) +
  geom_text(data=unique(s7[, .(organism, threshold, ratio, lo, hi, studies=pmax(studies) )]), aes(threshold, y=Inf, label=sprintf("ratio %.2f [%.2f, %.2f]", ratio, lo, hi)), vjust=1.5, size=2.4, inherit.aes=FALSE) +
  facet_wrap(~organism, scales="free_y") + scale_fill_manual(values=c(DGRY, GRY)) + scale_y_continuous(expand=expansion(mult=c(0, 0.25))) + labs(x=NULL, y="Mixed-sex studies", fill=NULL) + th
sv(pS7, "FigS7_imbalance_chance_expectation", 7.1, 3.0); src(e8[, .(organism, threshold, n_studies, observed, expected, ratio, ratio_exact_lo, ratio_exact_hi, p_two_sided)], "FigS7_imbalance_chance_expectation")

## ---------- S8: mechanism — SE inflation = √VIF (per study) ----------
e10 <- rbindlist(lapply(names(ORG), function(o) { f <- file.path(RV, sprintf("E10_per_study_%s.tsv", o)); if (!file.exists(f)) return(NULL); d <- fread(f)[, .(study, vif, se_ratio_M1_median)]; d[, organism := ORG[o]]; d }))
if (nrow(e10)) { e10[, organism := factor(organism, levels=ORG)]; e10[, sqrt_vif := sqrt(vif)]
  pS8 <- ggplot(e10, aes(sqrt_vif, se_ratio_M1_median)) + geom_abline(slope=1, intercept=0, linetype=2, linewidth=0.4, colour="grey45") + geom_point(size=1.2, alpha=0.6, colour=DGRY) + facet_wrap(~organism) +
    labs(x="√VIF of the sex term (from the group × sex table)", y="Median SE ratio, sex-adjusted / unadjusted") + th
  sv(pS8, "FigS8_se_inflation_vs_vif", 7.1, 3.0); src(e10, "FigS8_se_inflation_vs_vif") }
cat("supplementary panels done\n")
