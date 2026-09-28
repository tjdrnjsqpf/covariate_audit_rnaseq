#!/usr/bin/env Rscript
# Mouse pilot: recount3 mouse counts -> sex (reported ∪ expression-recovered, from 04m) -> M0/M1(/M2) limma-voom -> change metrics.
# Usage: Rscript 03_pilot_de_one_project.R <SRP> <group_key> <case_level> <control_level> [exclude]
#   levels: exact lower-case match first, otherwise regex.
#   5th argument (optional): ";"-separated. "ex:key=value" (exclude) "in:key=value" (keep only matching samples) "cov:key" (categorical covariate added to every model).
#   "key=value" without a prefix is treated as ex:.
suppressMessages({library(recount3); library(edgeR); library(limma); library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
SRP <- args[1]; GKEY <- tolower(args[2]); CASE <- tolower(args[3]); CTRL <- tolower(args[4]); EXCL <- if (length(args) >= 5) args[5] else ""
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
OUTDIR <- if (length(args) >= 6) args[6] else "results/pilot_mouse"
OUT <- file.path(PROJ, OUTDIR, SRP); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
panel <- data.table(gene=c("Xist","Ddx3y","Eif2s3y","Kdm5d","Uty"), chromosome=c("X","Y","Y","Y","Y"))
bfc <- recount3_cache(file.path(PROJ, "data/raw/recount3_cache", SRP))   # per-project cache: avoids SQLite lock conflicts in parallel runs; kept under /var2

## 1. count
rse <- create_rse_manual(project=SRP, project_home="data_sources/sra", organism="mouse", annotation="gencode_v23", type="gene", bfc=bfc)
assay(rse, "counts") <- transform_counts(rse)
cts <- assay(rse, "counts"); sym <- rowData(rse)$gene_name
md <- as.data.table(colData(rse))[, .(external_id, sra.sample_attributes, sra.sample_title)]

## Pseudoreplication fix (2026-09-08): runs sharing the same sample_acc (biological sample) are technical replicates (re-sequencing / lane splits), so
## counts are summed to one column per sample. Column names keep the representative run (first external_id) so downstream code is unchanged.
sacc_col <- grep("^sra\\.sample_acc", names(colData(rse)), value=TRUE)[1]   # the .x suffix depends on whether recount_pred was merged
sacc <- data.table(external_id = colData(rse)$external_id, sample_acc = as.character(colData(rse)[[sacc_col]]))
sacc[is.na(sample_acc) | sample_acc == "", sample_acc := external_id]
n_runs_raw <- ncol(cts)
if (anyDuplicated(sacc$sample_acc)) {
  sacc[, sample_acc := factor(sample_acc, levels = unique(sample_acc))]
  M <- model.matrix(~ 0 + sample_acc, sacc); rownames(M) <- sacc$external_id
  rep_id <- sacc[, external_id[1], by = sample_acc]$V1
  cts <- as.matrix(cts[, rownames(M)]) %*% M; colnames(cts) <- rep_id   # keep double: coercing to integer makes deep-sequencing sums exceed 2^31 and become NA
  md <- md[external_id %in% rep_id]
}
n_after_run_collapse <- ncol(cts)

## 2. metadata
al <- md[, .(kv = as.character(unlist(strsplit(sra.sample_attributes, "|", fixed=TRUE)))), by=external_id]
al[, c("key","value") := tstrsplit(kv, ";;", fixed=TRUE, keep=1:2)][, key := tolower(trimws(key))][, value := tolower(trimws(value))][, kv := NULL]
match_level <- function(v, pat) { m <- v == pat; if (!any(m, na.rm=TRUE)) m <- grepl(pat, v, perl=TRUE); m & !is.na(m) }
grp <- al[key == GKEY, .(external_id, gval = value)]
grp[, group := fifelse(match_level(gval, CASE), "case", fifelse(match_level(gval, CTRL), "control", NA_character_))]
st <- fread(file.path(PROJ, "data/processed/mouse_sample_table_sexinferred.tsv.gz"), colClasses=list(character=c("study","external_id")), na.strings=c("","NA"))[study == SRP]
sx  <- st[, .(external_id, sex_meta = fifelse(sex %in% c("F","M"), sex, NA_character_), sex_final_tbl = sex_final)]   # sex = reported (replaced by expression sex in label-flipped projects)
ag  <- st[!is.na(age_wk), .(external_id, age = age_wk)]
grp <- unique(grp, by="external_id"); sx <- unique(sx, by="external_id"); ag <- unique(ag, by="external_id")   # 2026-09-18: prevent duplicate rows for the same sample
ph <- Reduce(function(a, b) merge(a, b, by="external_id", all.x=TRUE), list(data.table(external_id=colnames(cts)), grp, sx, ag))
COVAR <- NULL; SUBJ <- NULL; SEXFIX <- NULL; TISSUE <- NULL; IS_MOUSE <- TRUE
if (nzchar(EXCL)) for (e in strsplit(EXCL, ";", fixed=TRUE)[[1]]) {
  e <- trimws(e); if (!nzchar(e)) next
  typ <- if (grepl("^(ex|in|cov|subj|sexfix|tissue):", e)) sub(":.*$", "", e) else "ex"; body <- sub("^(ex|in|cov|subj|sexfix|tissue):", "", e)
  if (typ == "cov") { COVAR <- tolower(body); next }
  if (typ == "subj") { SUBJ <- tolower(body); next }
  if (typ == "sexfix") { SEXFIX <- tolower(body); next }
  if (typ == "tissue") { TISSUE <- tolower(body); next }   # selects the cell-composition panel   # sensitivity: expr=use expression sex on mismatch, drop=exclude mismatched samples
  opt_k <- tolower(strsplit(body, "=", fixed=TRUE)[[1]][1]); opt_v <- tolower(sub("^[^=]*=", "", body))   # variable names chosen not to clash with data.table column names
  hit <- al[key == opt_k & value == opt_v, external_id]
  cat("option", typ, opt_k, "=", opt_v, ": matched", length(hit), "samples\n")
  ph <- if (typ == "ex") ph[!external_id %in% hit] else ph[external_id %in% hit] }
if (!is.null(COVAR)) { cv <- al[key == COVAR, .(external_id, covar = value)]; ph <- merge(ph, cv, by="external_id", all.x=TRUE); ph <- ph[!is.na(covar)]; ph[, covar := factor(covar)] }
ph <- ph[!is.na(group)][, group := factor(group, levels=c("control","case"))]
n_before_subj <- nrow(ph)
if (!is.null(SUBJ)) { sj <- al[key == SUBJ, .(external_id, subj = value)]
  ph <- merge(ph, sj, by="external_id", all.x=TRUE)
  ph <- rbind(ph[is.na(subj)], ph[!is.na(subj)][, .SD[1], by=.(subj, group)], fill=TRUE) }
if (nrow(ph) < 8 || any(table(ph$group) < 4)) stop("too few samples before sex assignment: n=", nrow(ph))   # 2026-09-18: final check is below (after removing samples with undetermined sex)
cts <- cts[, ph$external_id]

## 3. Expression-based sex inference: XIST vs mean Y-linked logCPM
lc <- cpm(cts, log=TRUE, prior.count=1)
## QC: detect targeted sequencing (amplicon, TCR repertoire, etc.) rather than whole transcriptome -> number of genes with CPM>1 at the sample median
n_genes_detected <- median(colSums(cpm(cts) > 1))
if (n_genes_detected < 5000) warning(SRP, ": only ", n_genes_detected, " genes detected (median CPM>1) - likely targeted assay")
ph[, xist := lc[which(sym == "Xist")[1], external_id]]
ph[, y_score := colMeans(lc[sym %in% panel[chromosome == "Y", gene], external_id, drop=FALSE])]
## Call: "ambiguous" if the XIST-Y difference is below 2 logCPM (both low = low quality, both high = suspected pooling/contamination)
ph[, sex_expr := fifelse(xist - y_score >= 2, "F", fifelse(y_score - xist >= 2, "M", "ambiguous"))]
ph[, sex_mismatch := !is.na(sex_meta) & sex_expr != "ambiguous" & sex_meta != sex_expr]
ph[, sex := fifelse(is.na(sex_meta), sex_expr, sex_meta)]   # mouse: prefer the sex finalised in 04m (reported, or expression if flipped), otherwise expression
if (!is.null(SEXFIX)) { if (SEXFIX == "expr") ph[sex_expr != "ambiguous", sex := sex_expr] else if (SEXFIX == "drop") ph <- ph[sex_mismatch == FALSE] }
n_sex_ambiguous_excluded <- sum(ph$sex == "ambiguous")   # no metadata sex and undeterminable from expression -> number of samples excluded (2026-09-18: review D 1.4)
ph <- ph[sex != "ambiguous"]   # exclude samples with no metadata sex that cannot be called from expression
stopifnot(nrow(ph) >= 8, all(table(ph$group) >= 4))   # 2026-09-18 (review D 1.6): minimum size (≥8, ≥4 per group) is applied after removing samples with undetermined sex


## C4: RNA quality proxies (recount3 QC) — record between-group SMD and use as covariates in M4
cdf <- as.data.frame(colData(rse)); qcol <- function(pat) { cc <- grep(pat, names(cdf), value=TRUE)[1]; if (is.na(cc)) rep(NA_real_, nrow(cdf)) else suppressWarnings(as.numeric(cdf[[cc]])) }
qc <- data.table(external_id = cdf$external_id, qc_mito = qcol("aligned_reads\\.\\.chrm$"), qc_unique = qcol("uniquely_mapped_reads_\\._both$"),
                 qc_intron = qcol("intron_sum_\\.$"), qc_assigned = qcol("gene_fc\\.unique_\\.$"), qc_depth = log10(qcol("number_of_input_reads_both$")))
ph <- merge(ph, qc, by="external_id", all.x=TRUE, sort=FALSE)
smd <- function(x, g) { a <- x[g == "case"]; b <- x[g == "control"]; s <- sqrt((var(a, na.rm=TRUE) + var(b, na.rm=TRUE)) / 2); if (is.na(s) || s == 0) NA_real_ else (mean(a, na.rm=TRUE) - mean(b, na.rm=TRUE)) / s }
QC_SMD <- sapply(c("qc_mito","qc_unique","qc_intron","qc_assigned","qc_depth"), function(k) round(smd(ph[[k]], ph$group), 3))


## C3: tissue-specific cell-composition proxy scores (mean z of marker genes). Panel chosen by option tissue:<category>, otherwise generic (immune/stromal).
PANELS <- list(
  blood   = list(a = c("FCGR3B","CSF3R","CXCR2","S100A8","S100A9","MNDA","FPR1"), b = c("CD3D","CD3E","CD2","IL7R","CD247","CD19","MS4A1","CD79A","CD79B","NKG7","GNLY","KLRD1")),
  brain   = list(a = c("RBFOX3","SNAP25","SYT1","GAD1","SLC17A7","SYP","STMN2"), b = c("GFAP","AQP4","MBP","PLP1","MOG","AIF1","CX3CR1","OLIG2")),
  liver   = list(a = c("ALB","APOA1","TTR","SERPINA1","HP","APOB","FGB"), b = c("PTPRC","CD68","PECAM1","COL1A1","ACTA2","LYZ")),
  adipose = list(a = c("ADIPOQ","LEP","PLIN1","FABP4","CIDEA","LPL"), b = c("PTPRC","CD68","COL1A1","PECAM1","LYZ","CD3E")),
  muscle  = list(a = c("ACTA1","MYH1","MYH2","MYH7","TNNT1","CKM","MB"), b = c("PTPRC","COL1A1","PECAM1","PDGFRA","CD68")),
  lung    = list(a = c("SFTPC","SFTPB","AGER","SCGB1A1","NKX2-1","SFTPA1"), b = c("PTPRC","CD68","MARCO","CD3E","LYZ","COL1A1")),
  kidney  = list(a = c("LRP2","UMOD","SLC34A1","SLC12A1","AQP2","NPHS2"), b = c("PTPRC","CD68","COL1A1","PECAM1","LYZ")),
  gut     = list(a = c("EPCAM","VIL1","KRT20","CDH17","LGR5","KRT8"), b = c("PTPRC","CD3E","CD79A","COL1A1","LYZ","PECAM1")),
  heart   = list(a = c("MYH6","MYH7","TNNT2","ACTC1","NPPA","MYL2"), b = c("PTPRC","COL1A1","PECAM1","PDGFRA","CD68")),
  generic = list(a = c("PTPRC","CD3E","CD68","CD79A","LYZ","CD14"), b = c("COL1A1","COL1A2","DCN","ACTA2","PECAM1","LUM")))
MOUSE_BLOOD <- list(a = c("S100a8","S100a9","Ly6g","Csf3r","Cxcr2","Mmp9","Retnlg"), b = c("Cd3d","Cd3e","Cd2","Il7r","Cd19","Ms4a1","Cd79a","Cd79b","Nkg7","Gzma","Klrd1"))
to_mouse <- function(g) paste0(substr(g, 1, 1), tolower(substring(g, 2)))
tissue_key <- function(t) { t <- tolower(if (is.null(t)) "" else t)
  if (grepl("blood|immune|pbmc", t)) "blood" else if (grepl("brain|cns|neur", t)) "brain" else if (grepl("liver", t)) "liver" else if (grepl("adipos|fat", t)) "adipose"
  else if (grepl("muscle|musculo", t)) "muscle" else if (grepl("lung", t)) "lung" else if (grepl("kidney|renal", t)) "kidney" else if (grepl("gut|intest|colon", t)) "gut" else if (grepl("heart|cardiac", t)) "heart" else "generic" }
CELL_PANEL <- tissue_key(TISSUE); pn <- PANELS[[CELL_PANEL]]
if (IS_MOUSE) pn <- if (CELL_PANEL == "blood") MOUSE_BLOOD else lapply(pn, to_mouse)
cell_score <- function(genes) { idx <- which(sym %in% genes); if (length(idx) < 3) return(rep(NA_real_, nrow(ph))); z <- t(scale(t(lc[idx, ph$external_id, drop=FALSE]))); colMeans(z, na.rm=TRUE) }
ph[, cell_a := cell_score(pn$a)]; ph[, cell_b := cell_score(pn$b)]
CELL_SMD <- c(cell_a_smd = round(smd(ph$cell_a, ph$group), 3), cell_b_smd = round(smd(ph$cell_b, ph$group), 3))
## C7: blood globin fraction (proxy for globin depletion and differences in its efficiency)
GLOB <- if (IS_MOUSE) c("Hbb-bs","Hbb-bt","Hbb-b1","Hbb-b2","Hba-a1","Hba-a2") else c("HBB","HBA1","HBA2","HBD")
ph[, globin_frac := as.numeric(colSums(cts[sym %in% GLOB, external_id, drop=FALSE]) / colSums(cts[, external_id, drop=FALSE]))]
GLOBIN <- c(globin_smd = round(smd(ph$globin_frac, ph$group), 3), globin_max = round(max(ph$globin_frac, na.rm=TRUE), 4), globin_range = round(diff(range(ph$globin_frac, na.rm=TRUE)), 4))
## C6: technical batch (SRA metadata). Cramér's V between group and batch variables.
bcols <- c(platform = "^sra\\.platform_model$", layout = "^sra\\.library_layout$", center = "^sra\\.run_center_name$", pubdate = "^sra\\.run_published$", readlen = "average_input_read_length$")
bt <- data.table(external_id = cdf$external_id)
for (nm in names(bcols)) { cc <- grep(bcols[[nm]], names(cdf), value=TRUE)[1]; bt[[nm]] <- if (is.na(cc)) NA_character_ else as.character(cdf[[cc]]) }
bt[, pubmonth := substr(pubdate, 1, 7)]; bt[, readlen := as.character(round(suppressWarnings(as.numeric(readlen)) / 10) * 10)]
ph <- merge(ph, bt[, .(external_id, platform, layout, center, pubmonth, readlen)], by="external_id", all.x=TRUE, sort=FALSE)
cramer <- function(x, g) { ok <- !is.na(x); if (sum(ok) < 4) return(NA_real_); tab <- table(x[ok], g[ok]); if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct=FALSE)$statistic); round(as.numeric(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1)))), 3) }
BATCH <- sapply(c("platform","layout","center","pubmonth","readlen"), function(k) cramer(ph[[k]], ph$group))
batch_var <- if (all(is.na(BATCH))) NA_character_ else names(BATCH)[which.max(BATCH)]; batch_v_max <- if (is.na(batch_var)) NA_real_ else BATCH[[batch_var]]
intron_range <- round(diff(range(ph$qc_intron, na.rm=TRUE)), 3)

## 4. Imbalance metrics
imb <- ph[, .(n = .N, pf = mean(sex == "F"), age_m = mean(age, na.rm=TRUE), age_sd = sd(age, na.rm=TRUE)), by=group]
sex_diff <- abs(diff(imb$pf)); age_smd <- abs(diff(imb$age_m)) / sqrt(mean(imb$age_sd^2, na.rm=TRUE))

## 5. DE
## B4 fix: compute filterByExpr once with the M0 design so every model uses the same gene universe
KEEP <- filterByExpr(DGEList(cts[, ph$external_id]), model.matrix(~group, ph))
run_de <- function(d, design) {
  y <- calcNormFactors(DGEList(cts[KEEP, d$external_id]))
  v <- voom(y, design); fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef="groupcase", number=Inf, sort.by="none"); tt$gene <- sym[KEEP]; tt$gene_id <- rownames(tt); as.data.table(tt)
}
cvt <- if (is.null(COVAR)) "" else " + covar"   # categorical covariates the original authors already had in the design (cell subtype etc.) are included in every model
mm <- function(f, d) model.matrix(as.formula(paste0("~group", f, cvt)), d)
set.seed(sum(utf8ToInt(SRP)))   # reproducible per-study seed
ph[, sex_perm := sample(sex)]                 # negative control: shuffle sex independently of group
res <- list(M0 = run_de(ph, mm("", ph)))
full_rank <- function(d, f) { X <- mm(f, d); qr(X)$rank == ncol(X) && all(colSums(X != 0) > 0) }
sex_ok <- uniqueN(ph$sex) == 2 && full_rank(ph, " + sex")   # estimable even if one group is single-sex, as long as the other group has both sexes
if (sex_ok) {
  res$M1 <- run_de(ph, mm(" + sex", ph))
  Xp <- mm(" + sex_perm", ph)
  if (qr(Xp)$rank == ncol(Xp)) res$M0perm <- run_de(ph, Xp)   # cost of adding a covariate per se (baseline)
}
pha <- ph[!is.na(age)]
age_ok <- sex_ok && nrow(pha) >= 6 && all(table(pha$group) >= 3) && sd(pha$age) > 0 && full_rank(pha, " + sex + age")
if (age_ok) {
  if (nrow(pha) < nrow(ph)) res$M0a <- run_de(pha, mm("", pha))   # M0 on the same subset when age is missing
  res$M2 <- run_de(pha, mm(" + sex + age", pha))
  pha[, age_perm := sample(age)]                 # negative control: shuffle age independently of group (cost of the extra age df per se)
  Xp <- mm(" + sex + age_perm", pha)
  if (qr(Xp)$rank == ncol(Xp)) res$M2perm <- run_de(pha, Xp)
}


## C5: SVA (M3) — do data-driven latent variables absorb the sex effect?
n_sv <- 0L; sv_sex_r2 <- NA_real_
if (sex_ok && nrow(ph) >= 10) tryCatch({
  suppressMessages(library(sva))
  v0 <- voom(calcNormFactors(DGEList(cts[KEEP, ph$external_id])), mm("", ph))
  nsv <- tryCatch(num.sv(v0$E, mm("", ph), method="be"), error=function(e) 0L)
  nsv <- min(nsv, 5L, floor((nrow(ph) - 4) / 2))
  if (nsv >= 1) {
    svo <- sva(v0$E, mm("", ph), model.matrix(~1, ph), n.sv=nsv)
    SV <- svo$sv; colnames(SV) <- paste0("sv", seq_len(ncol(SV))); n_sv <- ncol(SV)
    sv_sex_r2 <- round(max(apply(SV, 2, function(s) summary(lm(s ~ ph$sex))$r.squared)), 3)
    res$M3 <- run_de(ph, cbind(mm("", ph), SV))
  }
}, error=function(e) message("SVA failed: ", conditionMessage(e)))
## C4: M4 = M1 + quality covariates (mitochondrial fraction, unique alignment rate)
if (sex_ok && all(!is.na(ph$qc_mito)) && all(!is.na(ph$qc_unique)) && sd(ph$qc_mito) > 0 && sd(ph$qc_unique) > 0) {
  X4 <- cbind(mm(" + sex", ph), qc_mito = as.numeric(scale(ph$qc_mito)), qc_unique = as.numeric(scale(ph$qc_unique)))
  if (qr(X4)$rank == ncol(X4)) res$M4 <- tryCatch(run_de(ph, X4), error=function(e) { message("M4 failed: ", conditionMessage(e)); NULL })
}
## C3: M5 = M1 + cell-composition scores (tissue-specific panel)
if (sex_ok && all(!is.na(ph$cell_a)) && all(!is.na(ph$cell_b))) {
  X5 <- cbind(mm(" + sex", ph), cell_a = ph$cell_a, cell_b = ph$cell_b)
  if (qr(X5)$rank == ncol(X5)) res$M5 <- tryCatch(run_de(ph, X5), error=function(e) { message("M5 failed: ", conditionMessage(e)); NULL })
}
## C6: M6 = M1 + batch factor (batch variable most strongly entangled with group; not possible if fully confounded (V=1))
if (sex_ok && !is.na(batch_var) && batch_v_max < 0.999) {
  bf <- factor(ph[[batch_var]]); if (nlevels(bf) >= 2) { X6 <- cbind(mm(" + sex", ph), model.matrix(~bf)[, -1, drop=FALSE])
    if (qr(X6)$rank == ncol(X6) && nrow(ph) - ncol(X6) >= 4) res$M6 <- tryCatch(run_de(ph, X6), error=function(e) { message("M6 failed: ", conditionMessage(e)); NULL }) }
}
## C7: M7 = M1 + logit(globin fraction) — only in studies where globin is detected (max >= 1%)
if (sex_ok && GLOBIN[["globin_max"]] >= 0.01 && sd(ph$globin_frac) > 0) {
  X7 <- cbind(mm(" + sex", ph), globin = qlogis(pmin(pmax(ph$globin_frac, 1e-5), 1 - 1e-5)))
  if (qr(X7)$rank == ncol(X7)) res$M7 <- tryCatch(run_de(ph, X7), error=function(e) { message("M7 failed: ", conditionMessage(e)); NULL })
}


## Negative controls for the alternative axes: the same covariates randomly permuted across samples (baseline for the cost of adding a term; same logic as M0perm for sex)
pidx <- sample(nrow(ph))
if (!is.null(res$M4)) { Xq <- cbind(mm(" + sex", ph), qc_mito = as.numeric(scale(ph$qc_mito))[pidx], qc_unique = as.numeric(scale(ph$qc_unique))[pidx]); if (qr(Xq)$rank == ncol(Xq)) res$M4perm <- tryCatch(run_de(ph, Xq), error=function(e) { message("M4perm failed: ", conditionMessage(e)); NULL }) }
if (!is.null(res$M5)) { Xc <- cbind(mm(" + sex", ph), cell_a = ph$cell_a[pidx], cell_b = ph$cell_b[pidx]); if (qr(Xc)$rank == ncol(Xc)) res$M5perm <- tryCatch(run_de(ph, Xc), error=function(e) { message("M5perm failed: ", conditionMessage(e)); NULL }) }
if (!is.null(res$M6)) { bfp <- factor(ph[[batch_var]])[pidx]; Xb <- cbind(mm(" + sex", ph), model.matrix(~bfp)[, -1, drop=FALSE]); if (qr(Xb)$rank == ncol(Xb) && nrow(ph) - ncol(Xb) >= 4) res$M6perm <- tryCatch(run_de(ph, Xb), error=function(e) { message("M6perm failed: ", conditionMessage(e)); NULL }) }
if (!is.null(res$M7)) { Xg <- cbind(mm(" + sex", ph), globin = qlogis(pmin(pmax(ph$globin_frac, 1e-5), 1 - 1e-5))[pidx]); if (qr(Xg)$rank == ncol(Xg)) res$M7perm <- tryCatch(run_de(ph, Xg), error=function(e) { message("M7perm failed: ", conditionMessage(e)); NULL }) }

## W (2026-09-18, review A/C): within-group permutation control — shuffle sex only within each group, keeping |Δf| (group-sex collinearity) intact while breaking only the sex–expression link.
## Full permutation (M0perm) puts only the df cost into the baseline; within-group permutation (M0permw) includes df + collinearity cost. A separate seed keeps the existing random stream unchanged.
set.seed(sum(utf8ToInt(SRP)) + 1L)
ph[, sex_permw := if (.N > 1) sample(sex) else sex, by=group]
if (sex_ok) { Xw <- mm(" + sex_permw", ph); if (qr(Xw)$rank == ncol(Xw)) res$M0permw <- tryCatch(run_de(ph, Xw), error=function(e) { message("M0permw failed: ", conditionMessage(e)); NULL }) }

## 6. Change metrics
metric <- function(a, b, fdr=0.05) {
  m <- merge(a[, .(gene_id, gene, t0 = t, lfc0 = logFC, p0 = adj.P.Val)], b[, .(gene_id, t1 = t, lfc1 = logFC, p1 = adj.P.Val)], by="gene_id")
  d0 <- m[p0 < fdr, gene_id]; d1 <- m[p1 < fdr, gene_id]; un <- length(union(d0, d1))
  t0 <- head(m[order(p0, -abs(t0))], 100)$gene_id; t1 <- head(m[order(p1, -abs(t1))], 100)$gene_id   # 2026-09-18: ties broken by |t|; guards against a universe smaller than 100
  list(n_deg_ref = length(d0), n_deg_alt = length(d1), jaccard = if (un == 0) NA_real_ else round(length(intersect(d0, d1)) / un, 3),
       top100_overlap = length(intersect(t0, t1)) / 100, t_spearman = round(cor(m$t0, m$t1, method="spearman"), 3),
       sign_flip_deg = sum(sign(m[gene_id %in% d0, lfc0]) != sign(m[gene_id %in% d0, lfc1])),
       sexgene_deg_ref = sum(m[gene_id %in% d0, gene] %in% panel$gene), sexgene_deg_alt = sum(m[gene_id %in% d1, gene] %in% panel$gene))
}
summ <- data.table(study = SRP, group_key = GKEY, options = EXCL, n_runs_raw = n_runs_raw, n_after_run_collapse = n_after_run_collapse, n_before_subject_collapse = n_before_subj, n = nrow(ph), n_case = sum(ph$group == "case"), n_control = sum(ph$group == "control"),
                   frac_female = round(mean(ph$sex == "F"), 2), sex_diff = round(sex_diff, 3), age_smd = round(age_smd, 3),
                   n_genes_detected = n_genes_detected, n_female = sum(ph$sex == "F"), n_male = sum(ph$sex == "M"),
                   n_female_case = sum(ph$sex == "F" & ph$group == "case"), n_female_control = sum(ph$sex == "F" & ph$group == "control"), n_sex_meta_missing = sum(is.na(ph$sex_meta)), n_sex_mismatch = sum(ph$sex_mismatch), n_sex_ambiguous = sum(ph$sex_expr == "ambiguous"), n_sex_ambiguous_excluded = n_sex_ambiguous_excluded, n_age_missing = sum(is.na(ph$age)),
                   cell_panel = CELL_PANEL, cell_a_smd = CELL_SMD[["cell_a_smd"]], cell_b_smd = CELL_SMD[["cell_b_smd"]], globin_smd = GLOBIN[["globin_smd"]], globin_max = GLOBIN[["globin_max"]], globin_range = GLOBIN[["globin_range"]], batch_var = batch_var, batch_v_max = batch_v_max, batch_v_platform = BATCH[["platform"]], batch_v_layout = BATCH[["layout"]], batch_v_center = BATCH[["center"]], batch_v_pubmonth = BATCH[["pubmonth"]], batch_v_readlen = BATCH[["readlen"]], intron_range = intron_range, n_sv = n_sv, sv_sex_r2 = sv_sex_r2, qc_mito_smd = QC_SMD[["qc_mito"]], qc_unique_smd = QC_SMD[["qc_unique"]], qc_intron_smd = QC_SMD[["qc_intron"]], qc_assigned_smd = QC_SMD[["qc_assigned"]], qc_depth_smd = QC_SMD[["qc_depth"]], models = paste(names(res), collapse=","))
if (!is.null(res$M1)) { mt <- metric(res$M0, res$M1); for (k in names(mt)) summ[[paste0("M1_", k)]] <- mt[[k]] }
if (!is.null(res$M0perm)) { mt <- metric(res$M0, res$M0perm); for (k in names(mt)) summ[[paste0("M0perm_", k)]] <- mt[[k]] }
if (!is.null(res$M0permw)) { mt <- metric(res$M0, res$M0permw); for (k in names(mt)) summ[[paste0("M0permw_", k)]] <- mt[[k]] }
if (!is.null(res$M2)) { ref <- if (is.null(res$M0a)) res$M0 else res$M0a; mt <- metric(ref, res$M2); for (k in names(mt)) summ[[paste0("M2_", k)]] <- mt[[k]] }
if (!is.null(res$M2perm)) { ref <- if (is.null(res$M0a)) res$M0 else res$M0a; mt <- metric(ref, res$M2perm); for (k in names(mt)) summ[[paste0("M2perm_", k)]] <- mt[[k]] }
for (mx in c("M3","M4","M5","M6","M7","M4perm","M5perm","M6perm","M7perm")) if (!is.null(res[[mx]])) {
  mt <- metric(res$M0, res[[mx]]); for (k in names(mt)) summ[[paste0(mx, "_", k)]] <- mt[[k]]
  if (!is.null(res$M1)) { mt <- metric(res$M1, res[[mx]]); for (k in names(mt)) summ[[paste0(mx, "v1_", k)]] <- mt[[k]] }   # vs M1: does the extra covariate change more on top of sex adjustment?
}

fwrite(ph, file.path(OUT, "pheno.tsv"), sep="\t")
for (m in intersect(names(res), c("M0","M1","M0perm","M0permw","M0a","M2","M2perm"))) fwrite(res[[m]], file.path(OUT, paste0("de_", m, ".tsv.gz")), sep="\t")   # also save the perm tables (C2 pathway-level negative control)
fwrite(summ, file.path(OUT, "summary.tsv"), sep="\t")
print(t(summ))
