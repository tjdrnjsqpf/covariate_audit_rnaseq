#!/usr/bin/env Rscript
# Parse recount3 human SRA metadata and select pilot candidates that satisfy the blood / sex / age / case-control criteria.
# Usage: /var2/lsg/miniforge3/envs/recount3/bin/Rscript scripts/00_pilot/02_parse_metadata_blood_candidates.R
suppressMessages({library(data.table)})
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
RAW  <- file.path(PROJ, "data/raw/recount3_metadata")
PROC <- file.path(PROJ, "data/processed"); dir.create(PROC, showWarnings=FALSE, recursive=TRUE)
RES  <- file.path(PROJ, "results"); dir.create(RES, showWarnings=FALSE)

## ---- 1. Read all metadata (required columns only) ----
keep <- c("external_id","study","sra.study_title","sra.study_abstract","sra.experiment_title",
          "sra.library_strategy","sra.library_source","sra.library_selection","sra.library_layout",
          "sra.library_construction_protocol","sra.platform_model","sra.sample_attributes",
          "sra.sample_name","sra.sample_title","sra.run_published")
files <- list.files(RAW, pattern="^sra\\.sra\\..*\\.MD\\.gz$", full.names=TRUE)
cat("metadata files:", length(files), "\n")
read1 <- function(f) {
  d <- tryCatch(fread(f, sep="\t", quote="", header=TRUE, showProgress=FALSE, colClasses="character"),
                error=function(e) NULL)
  if (is.null(d)) return(NULL)
  # raw file columns lack the "sra." prefix (recount3::read_metadata adds it) -> add it here
  pre <- setdiff(names(d), c("rail_id","external_id","study")); setnames(d, pre, paste0("sra.", pre))
  miss <- setdiff(keep, names(d)); for (m in miss) d[[m]] <- NA_character_
  d[, ..keep]
}
md <- rbindlist(lapply(files, read1), use.names=TRUE, fill=TRUE)
cat("samples:", nrow(md), " projects:", uniqueN(md$study), "\n")
fwrite(md, file.path(PROC, "recount3_human_sra_metadata_min.tsv.gz"), sep="\t")

## ---- 2. sample_attributes → long ----
attr_long <- md[!is.na(sra.sample_attributes) & sra.sample_attributes != "",
                .(kv = as.character(unlist(strsplit(sra.sample_attributes, "|", fixed=TRUE)))), by=.(study, external_id)]
attr_long[, c("key","value") := tstrsplit(kv, ";;", fixed=TRUE, keep=1:2)]
attr_long[, key := tolower(trimws(key))][, value := trimws(value)][, kv := NULL]
attr_long <- attr_long[!is.na(key) & !is.na(value)]
fwrite(attr_long, file.path(PROC, "sample_attributes_long.tsv.gz"), sep="\t")

## ---- 3. Extract sex / age / tissue ----
sex_keys <- c("sex","gender","sex/gender","gender/sex","patient sex","donor sex","subject sex","biological sex","donor gender")
sx <- attr_long[key %in% sex_keys]
sx[, sex := fifelse(grepl("^(f|female|woman|w)$", tolower(value)), "F",
             fifelse(grepl("^(m|male|man)$", tolower(value)), "M", NA_character_))]
sx <- sx[!is.na(sex), .(sex = sex[1]), by=.(study, external_id)]

ag <- attr_long[grepl("^age", key) & !grepl("stage|group|category|range|bin|onset|diagnosis|death|at death", key)]
parse_age <- function(v) {
  v <- tolower(v)
  num <- suppressWarnings(as.numeric(sub("^.*?(-?[0-9]+\\.?[0-9]*).*$", "\\1", v)))
  # range "45-50" -> midpoint
  rng <- regmatches(v, regexpr("([0-9]+)\\s*[-–to]+\\s*([0-9]+)", v))
  if (length(rng)) { p <- as.numeric(unlist(strsplit(rng, "[^0-9]+"))); p <- p[!is.na(p)]; if (length(p)==2) num <- mean(p) }
  if (grepl("month", v)) num <- num/12
  if (grepl("week", v)) num <- num/52
  if (grepl("day", v)) num <- num/365
  if (grepl(">\\s*89|89\\+|90\\+", v)) num <- 90
  num
}
ag[, age := sapply(value, parse_age)]
ag <- ag[!is.na(age) & age >= 0 & age <= 110, .(age = age[1]), by=.(study, external_id)]

tissue_keys <- c("tissue","source_name","cell type","sample type","tissue type","cell_type","source name","tissue/cell type","biomaterial","sample_type","tissue source","material")
ts <- attr_long[key %in% tissue_keys, .(tissue_txt = paste(unique(value), collapse=" ; ")), by=.(study, external_id)]
blood_re <- "blood|pbmc|peripheral blood|leukocyte|leucocyte|monocyte|neutrophil|lymphocyte|granulocyte|buffy|cd4|cd8|cd14|cd19|t[- ]cell|b[- ]cell|nk[- ]cell|t lymph|b lymph|paxgene|tempus"
ts[, is_blood := grepl(blood_re, tolower(tissue_txt))]

cl <- attr_long[key %in% c("cell line","cell_line","cell line name")]
cl <- cl[!tolower(value) %in% c("","none","na","n/a","not applicable","no","-"), .(has_cell_line = TRUE), by=.(study, external_id)]

## ---- 4. Merge sample table ----
s <- md[, .(study, external_id, gsm = sra.sample_name, title = sra.sample_title,
            strategy = sra.library_strategy, source = sra.library_source, layout = sra.library_layout,
            protocol = sra.library_construction_protocol, platform = sra.platform_model,
            study_title = sra.study_title, abstract = sra.study_abstract, exp_title = sra.experiment_title,
            published = sra.run_published)]
s <- merge(s, sx, by=c("study","external_id"), all.x=TRUE)
s <- merge(s, ag, by=c("study","external_id"), all.x=TRUE)
s <- merge(s, ts, by=c("study","external_id"), all.x=TRUE)
s <- merge(s, cl, by=c("study","external_id"), all.x=TRUE)
s[is.na(has_cell_line), has_cell_line := FALSE][is.na(is_blood), is_blood := FALSE]
sc_re <- "single[- ]cell|single[- ]nucle|scrna|snrna|10x genomics|chromium|smart-seq|drop-seq|inDrop|cel-seq|mars-seq|seq-well"
s[, is_sc := grepl(sc_re, tolower(paste(protocol, study_title, exp_title, abstract)))]
small_re <- "mirna|micro rna|microrna|small rna|small-rna|piRNA|cage|rip-seq|clip|ribo-seq|ribosome profiling|circrna|cfrna|cell-free|exosom"
s[, is_small := grepl(small_re, tolower(paste(protocol, study_title, exp_title)))]
fwrite(s, file.path(PROC, "sample_table.tsv.gz"), sep="\t")

## ---- 5. Search candidate group keys per project ----
skip_keys <- c(sex_keys, tissue_keys, "age","source_name","sample name","sample_name","title","id","sample id","patient id","subject id","individual","donor","donor id","subject","replicate","batch","lane","library","barcode","run","biosample","description","sample description","molecule","organism","strain","isolate","ethnicity","race","bmi","cell type","tissue")
grp_re <- "disease|condition|diagnos|group|status|phenotype|case|control|treatment|state|cohort|outcome|severity|infect|response|type|class|subtype"
attr_long[, nchar_v := nchar(value)]
gk <- attr_long[!key %in% skip_keys & !grepl("^age", key) & nchar_v <= 60,
                .(n_lev = uniqueN(value), n = .N, min_lev = min(table(value)), levs = paste(names(sort(table(value), dec=TRUE))[1:min(4, uniqueN(value))], collapse=" / "),
                  cnts = paste(sort(table(value), dec=TRUE)[1:min(4, uniqueN(value))], collapse="/")), by=.(study, key)]
gk <- gk[n_lev >= 2 & n_lev <= 6 & min_lev >= 3]
gk[, pri := grepl(grp_re, key) * 10 + (n_lev == 2) * 5 + pmin(min_lev, 20)/20]
setorder(gk, study, -pri)
best <- gk[, .SD[1], by=study][, .(study, group_key = key, group_levels = levs, group_counts = cnts, group_nlev = n_lev, group_min_n = min_lev)]

## ---- 6. Project summary ----
p <- s[, .(n = .N, n_sex = sum(!is.na(sex)), n_age = sum(!is.na(age)), n_both = sum(!is.na(sex) & !is.na(age)),
           frac_female = mean(sex == "F", na.rm=TRUE), age_min = suppressWarnings(min(age, na.rm=TRUE)), age_max = suppressWarnings(max(age, na.rm=TRUE)),
           frac_blood = mean(is_blood), any_cell_line = any(has_cell_line), any_sc = any(is_sc), any_small = any(is_small),
           strategy = names(sort(table(strategy), dec=TRUE))[1], source = names(sort(table(source), dec=TRUE))[1],
           layout = names(sort(table(layout), dec=TRUE))[1], platform = names(sort(table(platform), dec=TRUE))[1],
           gse_hint = sub(":.*$", "", exp_title[1]), study_title = study_title[1], abstract = substr(abstract[1], 1, 400),
           first_published = substr(min(published, na.rm=TRUE), 1, 10)), by=study]
p <- merge(p, best, by="study", all.x=TRUE)
fwrite(p, file.path(RES, "project_summary_all.tsv"), sep="\t")
cat("projects summarized:", nrow(p), "\n")
cat("with sex>=90%:", sum(p$n_sex/p$n >= .9), " with age>=90%:", sum(p$n_age/p$n >= .9), " both:", sum(p$n_both/p$n >= .9), "\n")

## ---- 7. Blood pilot candidates ----
cand <- p[frac_blood >= 0.8 & n >= 12 & n <= 400 & n_both/n >= 0.9 & !any_cell_line & !any_sc & !any_small &
          strategy == "RNA-Seq" & source == "TRANSCRIPTOMIC" & !is.na(group_key) & group_min_n >= 5]
## Compute sex/age imbalance between the two largest groups
imb <- function(st, key) {
  a <- attr_long[study == st & key == key]
  x <- merge(a[, .(external_id, grp = value)], s[study == st, .(external_id, sex, age)], by="external_id")
  top2 <- names(sort(table(x$grp), dec=TRUE))[1:2]; x <- x[grp %in% top2]
  pf <- x[, .(pf = mean(sex == "F", na.rm=TRUE), ma = mean(age, na.rm=TRUE), sa = sd(age, na.rm=TRUE), n = .N), by=grp]
  sex_dp <- abs(diff(pf$pf)); smd <- abs(diff(pf$ma)) / sqrt(mean(pf$sa^2, na.rm=TRUE))
  fp <- tryCatch(fisher.test(table(x$grp, x$sex))$p.value, error=function(e) NA_real_)
  data.table(study = st, n_top2 = sum(pf$n), sex_diff = round(sex_dp, 2), sex_fisher_p = signif(fp, 2), age_smd = round(smd, 2))
}
if (nrow(cand)) {
  im <- rbindlist(Map(imb, cand$study, cand$group_key))
  cand <- merge(cand, im, by="study")
  cand[, design := fifelse(sex_diff >= 0.99, "completely_confounded", fifelse(sex_diff >= 0.3 | age_smd >= 0.8, "partially_confounded", "balanced"))]
  setorder(cand, -n_top2)
}
out_cols <- c("study","gse_hint","n","n_top2","n_both","frac_female","age_min","age_max","group_key","group_levels","group_counts",
              "sex_diff","sex_fisher_p","age_smd","design","layout","platform","first_published","study_title","abstract")
fwrite(cand[, ..out_cols], file.path(RES, "pilot_candidates_blood.tsv"), sep="\t")
cat("blood pilot candidates:", nrow(cand), "\n"); print(table(cand$design))
## Save the relaxed criterion (sex only, no age) separately
cand2 <- p[frac_blood >= 0.8 & n >= 12 & n <= 400 & n_sex/n >= 0.9 & n_age/n < 0.9 & !any_cell_line & !any_sc & !any_small &
           strategy == "RNA-Seq" & source == "TRANSCRIPTOMIC" & !is.na(group_key) & group_min_n >= 5]
fwrite(cand2, file.path(RES, "pilot_candidates_blood_sexonly.tsv"), sep="\t")
cat("blood sex-only candidates:", nrow(cand2), "\n")
