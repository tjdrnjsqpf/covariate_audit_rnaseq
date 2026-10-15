#!/usr/bin/env Rscript
## Figure S9 (Additional file 9): (A) share of projects with sex recorded by year, (B) single-sex share by year, (C) |df| >= 0.3 by era, (D) d_across / d_within by engine.
## Input: results/response/M_metadata_by_year.tsv (58), M_imbalance_by_era.tsv (60), E_engine_by_class.tsv (61). Usage (mac): LC_ALL=en_US.UTF-8 Rscript scripts/07_response/62_response_figures.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
PROJ <- Sys.getenv("COV_AUDIT_PROJ", getwd()); setwd(PROJ); OD <- "results/figures/response"; dir.create(OD, showWarnings=FALSE, recursive=TRUE)
FONT <- if (Sys.info()[["sysname"]] == "Darwin") "Helvetica Neue" else "sans"; RED <- "#B02418"; GRY <- "#C4C4C4"; BLK <- "#1A1A1A"; DGRY <- "#5F5F5F"
th <- theme_classic(base_size=9, base_family=FONT) + theme(axis.line=element_line(colour="black", linewidth=0.4), axis.ticks=element_line(colour="black", linewidth=0.4), axis.text=element_text(colour="black"),
  legend.position="top", legend.justification="left", legend.margin=margin(0,0,0,0), legend.box.spacing=unit(2,"pt"), legend.key.size=unit(9,"pt"), legend.text=element_text(size=8),
  strip.background=element_blank(), strip.text=element_text(face="bold", size=9, hjust=0), panel.spacing.x=unit(14,"pt"), plot.margin=margin(4,8,4,4))
tg <- theme(plot.tag=element_text(face="bold", size=11, family=FONT), plot.tag.position=c(0, 1), plot.margin=margin(14, 12, 6, 8))
ORG <- c(human="Human", mouse="Mouse")
sv <- function(p, f, w, h) { pdfdev <- if (Sys.info()[["sysname"]] == "Darwin") function(filename, width, height, ...) grDevices::quartz(file=filename, type="pdf", width=width, height=height) else cairo_pdf
  ggsave(file.path(OD, paste0(f, ".pdf")), p, width=w, height=h, device=pdfdev); tmp <- tempfile(fileext=".png")
  if (requireNamespace("ragg", quietly=TRUE)) ggsave(tmp, p, width=w, height=h, dpi=300, device=ragg::agg_png) else ggsave(tmp, p, width=w, height=h, dpi=300); invisible(file.copy(tmp, file.path(OD, paste0(f, ".png")), overwrite=TRUE)) }
wil <- function(k, n) { if (is.na(n) || n == 0) return(c(NA, NA)); z <- 1.96; p <- k/n; d <- 1 + z^2/n; c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d, (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d) }
## ---- A, B: by year — one SRA runinfo series (2009–2026) in red; the recount3-parser series (to 2019) as a gray dashed reference ----
Y <- fread("results/response/M_metadata_by_year.tsv")[!is.na(year) & year >= 2013 & year <= 2026]; Y[, organism := ORG[organism]]
Y[, src := fifelse(grepl("recount3", era), "recount3 attribute parser (to 2019)", "SRA runinfo, BioSample sex (2009–2026)")]
Y[, src := factor(src, c("SRA runinfo, BioSample sex (2009–2026)", "recount3 attribute parser (to 2019)"))]
Y[, c("lo_rec","hi_rec") := as.list(wil(round(frac_sex_any*n_studies), n_studies)), by=.(organism, era, year)]
Y[, c("lo_sin","hi_sin") := as.list(wil(n_single, n_eligible)), by=.(organism, era, year)]
yr_sc <- list(scale_colour_manual(values=c(RED, DGRY)), scale_fill_manual(values=c(RED, DGRY)), scale_linetype_manual(values=c("solid","22")), scale_x_continuous(breaks=seq(2013, 2026, 3)),
  guides(colour=guide_legend(nrow=2), fill="none", linetype="none"))
pA <- ggplot(Y, aes(year, frac_sex_any, colour=src, linetype=src)) + geom_ribbon(aes(ymin=lo_rec, ymax=hi_rec, fill=src), alpha=0.12, colour=NA) + geom_line(linewidth=0.5) + geom_point(size=1.2) +
  facet_wrap(~organism) + yr_sc + scale_y_continuous(labels=scales::percent, limits=c(0, NA)) + labs(x="Year of SRA release", y="Projects with sex recorded", colour=NULL) + th
pB <- ggplot(Y[n_eligible >= 20], aes(year, frac_single, colour=src, linetype=src)) + geom_ribbon(aes(ymin=lo_sin, ymax=hi_sin, fill=src), alpha=0.12, colour=NA) + geom_line(linewidth=0.5) + geom_point(size=1.2) +
  facet_wrap(~organism) + yr_sc + scale_y_continuous(labels=scales::percent, limits=c(0, 1)) + labs(x="Year of SRA release", y="Single-sex studies\n(sex recorded, 8–200 samples)", colour=NULL) + th
fwrite(Y[, .(organism, era, year, n_studies, frac_sex_recorded=frac_sex_any, n_eligible, n_single, frac_single)], file.path(OD, "source_A_B.tsv"), sep="\t")
## ---- C: |Δf| by era ----
pC <- NULL; fC <- "results/response/M_imbalance_by_era.tsv"
if (file.exists(fC)) { I <- fread(fC)[!is.na(n_with_group)]; if (!"year" %in% names(I)) I[, year := NA_integer_]; I[, organism := ORG[organism]]
  I[, period := fifelse(is.na(year), fifelse(grepl("recount3", era), "to 2019", "2020–26"), as.character(year))]
  I <- I[is.na(year)]; I[, period := factor(period, unique(c("to 2019", "2020–26", sort(unique(I[!is.na(year), period])))))]
  I[, c("lo","hi") := as.list(wil(round(frac_ge30*n_with_group), n_with_group)), by=.(organism, period)]; I[, kind := fifelse(is.na(year), "era", "year")]
  pC <- ggplot(I, aes(period, frac_ge30, colour=kind)) + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.25, linewidth=0.4) + geom_point(size=1.6) + facet_wrap(~organism) +
    scale_colour_manual(values=c(era=RED, year=DGRY), guide="none") + scale_y_continuous(labels=scales::percent, limits=c(0, NA)) + geom_text(aes(y=0, label=n_with_group), size=2.3, vjust=-0.5, colour=BLK, family=FONT) +
    labs(x="SRA release", y="Mixed-sex studies with |Δf| ≥ 0.3\n(metadata only)") + th + coord_cartesian(clip="off")
  fwrite(I[, .(organism, period, n_with_group, frac_ge30, lo, hi)], file.path(OD, "source_C.tsv"), sep="\t") }
## ---- D: engines ----
pD <- NULL; fD <- "results/response/E_engine_by_class.tsv"
if (file.exists(fD)) { E <- fread(fD); num <- function(s) { m <- regmatches(s, regexec("^(-?[0-9.]+)(?: \\[(-?[0-9.]+), (-?[0-9.]+)\\])?", s, perl=TRUE))[[1]]; as.numeric(m[2:4]) }
  D <- rbindlist(lapply(seq_len(nrow(E)), function(i) rbind(data.table(E[i, .(organism, design_sex, engine, n)], metric="Δ vs across-group control", t(num(E$d_across[i]))), data.table(E[i, .(organism, design_sex, engine, n)], metric="Δ vs within-group control", t(num(E$d_within[i]))))))
  setnames(D, c("V1","V2","V3"), c("m","lo","hi")); D <- D[n >= 10]; D[, `:=`(organism = ORG[organism], engine = factor(c(limma="limma-
voom", edgeR="edgeR
QL", DESeq2="DESeq2")[engine], c("limma-
voom","edgeR
QL","DESeq2")), metric = factor(metric, c("Δ vs across-group control","Δ vs within-group control")),
    cls = factor(fifelse(design_sex == "balanced", "Balanced", "Confounded"), c("Balanced","Confounded")))]
  D[, lab := factor(sprintf("%s, %s", organism, tolower(cls)))]; D[, lab := factor(lab, unique(lab[order(cls, organism)]))]
  pD <- ggplot(D, aes(engine, m, colour=metric)) + geom_hline(yintercept=0, linewidth=0.3, colour=GRY) + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.25, linewidth=0.4, position=position_dodge(0.5)) + geom_point(size=1.6, position=position_dodge(0.5)) +
    facet_wrap(~lab, nrow=2, scales="free_y") + scale_colour_manual(values=c(GRY, BLK)) + labs(x="DE engine (limma-voom, edgeR QL, DESeq2)", y="Δ Jaccard (median, 95% CI)", colour=NULL) + th + theme(panel.spacing.y=unit(8,"pt"), panel.spacing.x=unit(10,"pt"))
  fwrite(D[, .(organism, cls, engine, n, metric, m, lo, hi)], file.path(OD, "source_D.tsv"), sep="\t") }
top <- wrap_elements(pA + labs(tag="A") + tg) / wrap_elements(pB + labs(tag="B") + tg); sv(top, "response_AB", 6.3, 5)
if (!is.null(pC) && !is.null(pD)) { P <- (wrap_elements(pA + labs(tag="A") + tg) | wrap_elements(pB + labs(tag="B") + tg)) / (wrap_elements(pC + labs(tag="C") + tg) + wrap_elements(pD + labs(tag="D") + tg) + plot_layout(widths=c(1, 1.3))); sv(P, "response_figure", 7.2, 6.6) }
if (!is.null(pD)) sv(pD + labs(tag="D") + tg, "response_D", 4.5, 3.6)
cat("done\n")
