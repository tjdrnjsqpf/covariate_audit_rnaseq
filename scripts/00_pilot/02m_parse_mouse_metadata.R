#!/usr/bin/env Rscript
# Parse recount3 mouse SRA metadata + prevalence survey (sex/age/strain reporting rates, fraction of single-sex studies, between-group imbalance in mixed-sex studies)
suppressMessages(library(data.table))
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
RAW  <- file.path(PROJ, "data/raw/recount3_metadata_mouse"); PROC <- file.path(PROJ, "data/processed"); RES <- file.path(PROJ, "results")
keep <- c("external_id","study","sra.study_title","sra.study_abstract","sra.experiment_title","sra.library_strategy","sra.library_source",
          "sra.library_layout","sra.library_construction_protocol","sra.platform_model","sra.sample_attributes","sra.sample_name","sra.sample_title","sra.run_published")
read1 <- function(f) { d <- tryCatch(fread(f, sep="\t", quote="", header=TRUE, showProgress=FALSE, colClasses="character"), error=function(e) NULL); if (is.null(d)) return(NULL)
  pre <- setdiff(names(d), c("rail_id","external_id","study")); setnames(d, pre, paste0("sra.", pre)); for (m in setdiff(keep, names(d))) d[[m]] <- NA_character_; d[, ..keep] }
files <- list.files(RAW, pattern="^sra\\.sra\\..*\\.MD\\.gz$", full.names=TRUE)
md <- rbindlist(lapply(files, read1), use.names=TRUE, fill=TRUE); cat("mouse files:", length(files), " samples:", nrow(md), " projects:", uniqueN(md$study), "\n")
al <- md[!is.na(sra.sample_attributes) & sra.sample_attributes != "", .(kv = as.character(unlist(strsplit(sra.sample_attributes, "|", fixed=TRUE)))), by=.(study, external_id)]
al[, c("key","value") := tstrsplit(kv, ";;", fixed=TRUE, keep=1:2)][, key := tolower(trimws(key))][, value := tolower(trimws(value))][, kv := NULL]; al <- al[!is.na(key) & !is.na(value)]
fwrite(al, file.path(PROC, "mouse_sample_attributes_long.tsv.gz"), sep="\t")

## Sex
sx <- al[key %in% c("sex","gender","sex/gender","animal sex","mouse sex")]
sx[, sex := fifelse(grepl("^(f|female|females|w)$", value), "F", fifelse(grepl("^(m|male|males)$", value), "M", fifelse(grepl("mixed|both|pool", value), "mixed", NA_character_)))]
sx <- sx[!is.na(sex), .(sex = sex[1]), by=.(study, external_id)]
## Age -> weeks. E##=embryonic (stage flag instead of negative values), P##=postnatal days
ag <- al[grepl("^age|^developmental stage|^dev_stage|^stage$", key) & !grepl("group|category|bin|onset|death", key)]
parse_wk <- function(v) {
  n <- suppressWarnings(as.numeric(sub("^.*?([0-9]+\\.?[0-9]*).*$", "\\1", v)))
  rng <- regmatches(v, regexpr("([0-9]+\\.?[0-9]*)\\s*(-|–|to)\\s*([0-9]+\\.?[0-9]*)", v)); if (length(rng)) { p <- as.numeric(unlist(regmatches(rng, gregexpr("[0-9]+\\.?[0-9]*", rng)))); if (length(p) == 2) n <- mean(p) }
  if (grepl("^e[0-9]|embryo|gestation|^gd|dpc", v)) return(c(NA_real_, "embryonic"))
  if (grepl("^p[0-9]|postnatal day|pnd", v)) return(c(n/7, "postnatal_day"))
  if (grepl("day", v)) return(c(n/7, "day"))
  if (grepl("month|mo\\b|mos", v)) return(c(n*4.33, "month"))
  if (grepl("year|yr", v)) return(c(n*52, "year"))
  if (grepl("week|wk", v)) return(c(n, "week"))
  if (grepl("adult", v) & is.na(n)) return(c(NA_real_, "adult_unspecified"))
  return(c(n, "unitless"))
}
pw <- t(sapply(ag$value, parse_wk)); ag[, age_wk := suppressWarnings(as.numeric(pw[,1]))][, age_unit := pw[,2]]
ag <- ag[, .(age_wk = age_wk[1], age_unit = age_unit[1]), by=.(study, external_id)]
## strain / genotype / treatment / tissue
st <- al[key %in% c("strain","strain background","genetic background","background strain","mouse strain"), .(strain = value[1]), by=.(study, external_id)]
ts <- al[key %in% c("tissue","source_name","cell type","sample type","tissue type","cell_type","organ","source name"), .(tissue_txt = paste(unique(value), collapse=" ; ")), by=.(study, external_id)]
tissue_cat <- function(x) { x <- tolower(x)
  fifelse(grepl("blood|pbmc|leukocyte|monocyte|neutrophil|lymphocyte|splenocyte|bone marrow|bmdm|macrophage|t cell|b cell", x), "blood_immune",
  fifelse(grepl("brain|cortex|hippocamp|cerebell|striatum|hypothalam|spinal|retina|neuron|glia", x), "brain_cns",
  fifelse(grepl("liver|hepat", x), "liver", fifelse(grepl("colon|intestin|ileum|gut|stomach|jejun|duoden|cecum", x), "gut",
  fifelse(grepl("lung|airway|trachea", x), "lung", fifelse(grepl("heart|cardiac|ventric|myocard", x), "heart",
  fifelse(grepl("kidney|renal", x), "kidney", fifelse(grepl("muscle|gastrocnemius|quadricep|tibialis|soleus", x), "muscle",
  fifelse(grepl("adipose|fat|wat|bat", x), "adipose", fifelse(grepl("skin|epiderm|dermis", x), "skin",
  fifelse(grepl("tumor|tumour|carcinoma|cancer|melanoma|glioma|lymphoma|leukemia", x), "tumor",
  fifelse(grepl("embryo|esc|es cell|ipsc|blastocyst|organoid|mef|fibroblast|cell line|3t3|hek|raw264|c2c12|nih", x), "cell_model_embryo",
  fifelse(grepl("testis|ovary|uterus|placenta|prostate|mammary|spermat|oocyte", x), "reproductive", fifelse(is.na(x) | x == "", "unknown", "other")))))))))))))) }
ts[, tissue_category := tissue_cat(tissue_txt)]
s <- md[, .(study, external_id, gsm = sra.sample_name, strategy = sra.library_strategy, source = sra.library_source, protocol = sra.library_construction_protocol,
            study_title = sra.study_title, exp_title = sra.experiment_title, abstract = sra.study_abstract)]
for (x in list(sx, ag, st, ts)) s <- merge(s, x, by=c("study","external_id"), all.x=TRUE)
s[, is_sc := grepl("single[- ]cell|single[- ]nucle|scrna|snrna|10x genomics|chromium|smart-seq|drop-seq|indrop|cel-seq|mars-seq", tolower(paste(protocol, study_title, exp_title, abstract)))]
s[, is_small := grepl("mirna|micro rna|microrna|small rna|small-rna|pirna|cage|rip-seq|clip|ribo-seq|ribosome profiling", tolower(paste(protocol, study_title, exp_title)))]
fwrite(s, file.path(PROC, "mouse_sample_table.tsv.gz"), sep="\t")

## group key (same scoring as human)
skip <- c("sex","gender","age","strain","tissue","source_name","cell type","sample type","tissue type","cell_type","organ","title","id","sample id","replicate","batch","lane","library","barcode","run","biosample","description","organism","isolate","source name","developmental stage")
good_re <- "genotype|treatment|condition|disease|diet|group|infect|status|phenotype|knockout|ko|wt|mutant|exposure|model|challenge|stimul"
bad_re  <- "date|update|flowcell|lane|run|batch|barcode|library|cohort|time|day|week|hour|replicate|^id$|_id$| id$|name|accession|biosample|insdc|checklist|protocol|kit|platform|instrument|read|strand|layout|rin|passage|well|plate|age|sex|gender|strain"
al[, nv := nchar(value)]
gk <- al[!key %in% skip & nv <= 60, .(n_lev = uniqueN(value), n = .N, min_lev = min(table(value)), levs = paste(names(sort(table(value), dec=TRUE))[1:min(4, uniqueN(value))], collapse=" / ")), by=.(study, key)][n_lev >= 2 & n_lev <= 6 & min_lev >= 2]
gk[, score := 10*grepl(good_re, key) - 10*grepl(bad_re, key) + 5*(n_lev == 2) + pmin(min_lev, 10)/10]; setorder(gk, study, -score)
best <- gk[, .(group_key = key[1], group_levels = levs[1], group_nlev = n_lev[1], group_min_n = min_lev[1], group_score = score[1]), by=study]

## Project summary
mode1 <- function(x) { t <- table(x); if (!length(t)) NA_character_ else names(sort(t, decreasing=TRUE))[1] }
p <- s[, .(n = .N, n_sex = sum(!is.na(sex)), n_F = sum(sex == "F", na.rm=TRUE), n_M = sum(sex == "M", na.rm=TRUE), n_mixedpool = sum(sex == "mixed", na.rm=TRUE),
           n_age = sum(!is.na(age_wk)), n_age_any = sum(!is.na(age_unit)), n_embryo = sum(age_unit == "embryonic", na.rm=TRUE), age_med_wk = suppressWarnings(median(age_wk, na.rm=TRUE)),
           n_strain = sum(!is.na(strain)), strain_top = mode1(strain), tissue_category = mode1(tissue_category), tissue_example = mode1(tissue_txt),
           strategy = mode1(strategy), source = mode1(source), any_sc = any(is_sc), any_small = any(is_small), study_title = study_title[1]), by=study]
p <- merge(p, best, by="study", all.x=TRUE)
p[, sex_design := fifelse(n_sex/n < 0.9, "unreported", fifelse(n_F > 0 & n_M == 0, "female_only", fifelse(n_M > 0 & n_F == 0, "male_only", fifelse(n_F > 0 & n_M > 0, "mixed", "pooled/other"))))]
## Mixed-sex studies: difference in sex proportion between the two largest groups of the best group key
imb <- function(st, gkey) { a <- al[study == st & key == gkey, .(external_id, grp = value)]; x <- merge(a, s[study == st, .(external_id, sex)], by="external_id")[sex %in% c("F","M")]
  top2 <- names(sort(table(x$grp), dec=TRUE))[1:2]; x <- x[grp %in% top2]; if (uniqueN(x$grp) < 2) return(NA_real_); pf <- x[, mean(sex == "F"), by=grp]$V1; abs(diff(pf)) }
mixed <- p[sex_design == "mixed" & !is.na(group_key) & group_score > 0]
mixed[, sex_diff := mapply(imb, study, group_key)]
p <- merge(p, mixed[, .(study, sex_diff)], by="study", all.x=TRUE)
p[, design := fifelse(is.na(sex_diff), NA_character_, fifelse(sex_diff >= 0.99, "completely_confounded", fifelse(sex_diff >= 0.3, "partially_confounded", fifelse(sex_diff >= 0.15, "mildly_imbalanced", "balanced"))))]
fwrite(p, file.path(RES, "mouse_project_summary_all.tsv"), sep="\t")

## Prevalence figures
base <- p[strategy == "RNA-Seq" & source == "TRANSCRIPTOMIC" & !any_sc & !any_small & n >= 6 & !tissue_category %in% c("cell_model_embryo")]
cat("\n== MOUSE Aim1 (bulk RNA-seq, tissue/primary, n>=6):", nrow(base), "projects,", sum(base$n), "samples ==\n")
cat(sprintf("sex >=90%%: %d (%.1f%%)\n", sum(base$n_sex/base$n >= .9), 100*mean(base$n_sex/base$n >= .9)))
cat(sprintf("age(any) >=90%%: %d (%.1f%%)   age numeric >=90%%: %d (%.1f%%)\n", sum(base$n_age_any/base$n >= .9), 100*mean(base$n_age_any/base$n >= .9), sum(base$n_age/base$n >= .9), 100*mean(base$n_age/base$n >= .9)))
cat(sprintf("strain >=90%%: %d (%.1f%%)\n", sum(base$n_strain/base$n >= .9), 100*mean(base$n_strain/base$n >= .9)))
cat(sprintf("sex+age+strain all >=90%%: %d (%.1f%%)\n", sum(base$n_sex/base$n >= .9 & base$n_age_any/base$n >= .9 & base$n_strain/base$n >= .9), 100*mean(base$n_sex/base$n >= .9 & base$n_age_any/base$n >= .9 & base$n_strain/base$n >= .9)))
cat("\n-- sex design (all base projects) --\n"); print(round(100*prop.table(table(base$sex_design)), 1))
cat("-- sex design among sex-reported --\n"); print(round(100*prop.table(table(base[sex_design != "unreported", sex_design])), 1))
cat("\n-- mixed-sex projects with a group key: imbalance class --\n"); print(table(base$design, useNA="ifany"))
cat("\n-- by tissue: projects / sex reported% / female_only% / male_only% / mixed% --\n")
bt <- base[, .(projects=.N, sex_rep=round(100*mean(n_sex/n >= .9)), female_only=round(100*mean(sex_design=="female_only")), male_only=round(100*mean(sex_design=="male_only")), mixed=round(100*mean(sex_design=="mixed")), age_rep=round(100*mean(n_age_any/n >= .9)), strain_rep=round(100*mean(n_strain/n >= .9))), by=tissue_category][order(-projects)]
print(bt); fwrite(bt, file.path(RES, "mouse_aim1_by_tissue.tsv"), sep="\t")
cat("\n-- age unit distribution (samples with age) --\n"); print(sort(table(s$age_unit), decreasing=TRUE))
cat("-- top strains --\n"); print(head(sort(table(s$strain), decreasing=TRUE), 10))
