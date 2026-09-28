#!/usr/bin/env Rscript
# All-tissue version of 02b: select candidates without restricting to blood and assign tissue_category.
# Usage: /var2/lsg/miniforge3/envs/recount3/bin/Rscript scripts/00_pilot/02b_select_blood_candidates.R
suppressMessages(library(data.table))
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
PROC <- file.path(PROJ, "data/processed"); RES <- file.path(PROJ, "results")
s  <- fread(file.path(PROC, "sample_table.tsv.gz"), colClasses=list(character=c("study","external_id","sex","gsm")), na.strings=c("","NA"))
al <- fread(file.path(PROC, "sample_attributes_long.tsv.gz"), colClasses="character", na.strings=c("","NA"))[!is.na(key) & !is.na(value)]

## ---- 1. Score candidate group variables ----
sex_keys <- c("sex","gender","sex/gender","gender/sex","patient sex","donor sex","subject sex","biological sex","donor gender")
tissue_keys <- c("tissue","source_name","cell type","sample type","tissue type","cell_type","source name","tissue/cell type","biomaterial","sample_type","tissue source","material")
skip_keys <- c(sex_keys, tissue_keys, "title","id","sample id","patient id","subject id","individual","donor","donor id","subject","replicate","batch","lane","library","barcode","run","biosample","description","sample description","molecule","organism","strain","isolate","ethnicity","race","bmi")
good_re <- "disease|diagnos|condition|status|phenotype|health|group|case|control|patient|clinical|severity|infect|outcome|state|karyotype|genotype|mutation|responder|response"
bad_re  <- "date|update|flowcell|lane|run|batch|barcode|library|country|location|site|center|centre|cohort|time|day|week|hour|visit|vaccin|stimul|dose|replicate|^id$|_id$| id$|name|accession|biosample|insdc|ena |sra|checklist|protocol|kit|platform|instrument|read|strand|layout|rin|concentr|volume|passage|well|plate|tube|extraction|collection|age|sex|gender|ethnic|race|bmi|weight|height|smok"
al[, value := tolower(trimws(value))]
al[, nv := nchar(value)]
gk <- al[!key %in% skip_keys & nv <= 60,
         .(n_lev = uniqueN(value), n = .N, min_lev = min(table(value)),
           levs = paste(names(sort(table(value), dec=TRUE))[1:min(4, uniqueN(value))], collapse=" / "),
           cnts = paste(sort(table(value), dec=TRUE)[1:min(4, uniqueN(value))], collapse="/")), by=.(study, key)]
gk <- gk[n_lev >= 2 & n_lev <= 6 & min_lev >= 3]
gk[, score := 10*grepl(good_re, key) - 10*grepl(bad_re, key) + 3*(key %in% c("treatment","treatment group","treatmentgroup")) + 5*(n_lev == 2) + pmin(min_lev, 20)/20]
## Penalise values that look like dates/IDs
gk[, score := score - 10*grepl("^[0-9]{4}-[0-9]{2}|^[A-Z0-9]{8,}$|^\\[?[0-9]+\\]?$", sub(" /.*$", "", levs))]
setorder(gk, study, -score)
best <- gk[, .(group_key = key[1], group_levels = levs[1], group_counts = cnts[1], group_nlev = n_lev[1], group_min_n = min_lev[1], group_score = score[1],
               alt_keys = paste(head(key[-1], 3), collapse=" | ")), by=study]

## ---- 2. Project summary ----
mode1 <- function(x) { t <- table(x); if (!length(t)) NA_character_ else names(sort(t, decreasing=TRUE))[1] }
p <- s[, .(n = .N, n_sex = sum(!is.na(sex)), n_age = sum(!is.na(age)), n_both = sum(!is.na(sex) & !is.na(age)),
           frac_female = round(mean(sex == "F", na.rm=TRUE), 2), age_min = suppressWarnings(round(min(age, na.rm=TRUE), 1)), age_max = suppressWarnings(round(max(age, na.rm=TRUE), 1)),
           frac_blood = mean(is_blood), tissue_example = mode1(tissue_txt),
           any_cell_line = any(has_cell_line), any_sc = any(is_sc), any_small = any(is_small),
           strategy = mode1(strategy), source = mode1(source),
           layout = mode1(layout), platform = mode1(platform),
           gsm_first = gsm[1], study_title = study_title[1], abstract = substr(abstract[1], 1, 500),
           first_published = substr(min(published, na.rm=TRUE), 1, 10)), by=study]
p <- merge(p, best, by="study", all.x=TRUE)
p[, pediatric := is.finite(age_max) & age_max < 18]
## Design flags: in-vitro treatment / repeated measures on the same subject / cell-subtype comparison -> recommended for exclusion from the case-control pilot
invitro_re <- "stimulat|in vitro|in-vitro|cultur|vehicle|dmso|xenograft|transduc|expanded|conditioning|treated with|incubat|infected with|exposed to"
repeat_re  <- "pre or post|pre/post|pre- and post|timepoint|time point|time course|timecourse|visit|longitudinal|before and after|follow-up"
subtype_re <- "cell_subtype|cell subtype|subset|sorting|sorted|cell type|cell_type"
kv <- al[, .(txt = paste(unique(key), collapse=" ; ")), by=study]
p <- merge(p, kv, by="study", all.x=TRUE)
p[, flag_invitro := grepl(invitro_re, tolower(paste(study_title, group_key, group_levels, tissue_example)))]
p[, flag_repeat  := grepl(repeat_re,  tolower(paste(group_key, group_levels, txt)))]
p[, flag_subtype := grepl(subtype_re, tolower(paste(group_key, group_levels)))]
p[, txt := NULL]
tissue_cat <- function(x) { x <- tolower(x)
  fifelse(grepl("blood|pbmc|leukocyte|leucocyte|monocyte|neutrophil|lymphocyte|granulocyte|buffy|cd4|cd8|cd14|cd19|t[- ]cell|b[- ]cell|nk[- ]cell|paxgene|platelet", x), "blood",
  fifelse(grepl("brain|cortex|hippocamp|cerebell|striatum|amygdala|thalam|hypothalam|substantia|putamen|caudate|nucleus accumbens|white matter|grey matter|gray matter|spinal|dorsolateral", x), "brain",
  fifelse(grepl("colon|rectum|rectal|ileum|ileal|intestin|duoden|jejun|gut|bowel|gastric|stomach|esophag|oesophag", x), "gut",
  fifelse(grepl("liver|hepat", x), "liver",
  fifelse(grepl("lung|bronch|airway|nasal|trachea|sputum|alveol", x), "lung",
  fifelse(grepl("skin|epiderm|dermis|keratinocyte|psoria|lesion", x), "skin",
  fifelse(grepl("muscle|vastus|biceps|quadricep|myo", x), "muscle",
  fifelse(grepl("adipose|fat|subcutaneous|omental", x), "adipose",
  fifelse(grepl("kidney|renal|glomerul", x), "kidney",
  fifelse(grepl("heart|cardiac|ventric|atri|myocard", x), "heart",
  fifelse(grepl("tumor|tumour|carcinoma|cancer|adenoma|glioma|melanoma|lymphoma|leukemia|leukaemia|sarcoma|neoplasm|metasta", x), "tumor",
  fifelse(grepl("synovi|joint|cartilage|bone|tendon", x), "musculoskeletal",
  fifelse(grepl("placenta|endometri|ovar|uter|cervi|prostate|testis|testic|breast|mammary", x), "reproductive",
  fifelse(grepl("fibroblast|ipsc|stem cell|organoid|cell line|hek|hela", x), "cell_model",
  fifelse(is.na(x) | x == "", "unknown", "other"))))))))))))))) }
p[, tissue_category := tissue_cat(tissue_example)]
fwrite(p, file.path(RES, "project_summary_all.tsv"), sep="\t")

## ---- 2b. For review: top 5 candidate group keys per project ----
gk_top <- gk[, head(.SD, 5), by=study][, .(study, key, n_lev, min_lev, score = round(score, 1), levs, cnts)]
fwrite(gk_top, file.path(RES, "group_key_candidates_top5.tsv"), sep="\t")  # shared across all projects

## ---- 3. Aim 1 preliminary prevalence figures ----
base <- p[strategy == "RNA-Seq" & source == "TRANSCRIPTOMIC" & !any_sc & !any_small & n >= 12]
cat("== Aim1 preliminary (human bulk RNA-seq, n>=12) ==\n")
cat("projects:", nrow(base), "\n")
cat("sex >=90%:", sum(base$n_sex/base$n >= .9), sprintf("(%.1f%%)", 100*mean(base$n_sex/base$n >= .9)), "\n")
cat("age >=90%:", sum(base$n_age/base$n >= .9), sprintf("(%.1f%%)", 100*mean(base$n_age/base$n >= .9)), "\n")
cat("both>=90%:", sum(base$n_both/base$n >= .9), sprintf("(%.1f%%)", 100*mean(base$n_both/base$n >= .9)), "\n")
cat("blood (>=50%) projects:", sum(base$frac_blood >= .5), " of which both>=90%:", sum(base$frac_blood >= .5 & base$n_both/base$n >= .9), "\n")
cat("\n== by tissue_category (n>=12, no cell line): projects / sex>=90% / age>=90% / both>=90% ==\n")
bt <- base[any_cell_line == FALSE, .(projects=.N, sex=sum(n_sex/n>=.9), age=sum(n_age/n>=.9), both=sum(n_both/n>=.9)), by=tissue_category][order(-projects)]
print(bt); fwrite(bt, file.path(RES, "aim1_metadata_by_tissue.tsv"), sep="\t")

## ---- 4. Imbalance computation (two largest groups) ----
imb <- function(st, gkey) {
  a <- al[study == st & key == gkey, .(external_id, grp = value)]
  x <- merge(a, s[study == st, .(external_id, sex, age)], by="external_id")
  top2 <- names(sort(table(x$grp), dec=TRUE))[1:2]; x <- x[grp %in% top2]
  pf <- x[, .(pf = mean(sex == "F", na.rm=TRUE), ma = mean(age, na.rm=TRUE), sa = sd(age, na.rm=TRUE), n = .N, nf = sum(sex == "F", na.rm=TRUE), nm = sum(sex == "M", na.rm=TRUE)), by=grp]
  sex_dp <- abs(diff(pf$pf)); smd <- abs(diff(pf$ma)) / sqrt(mean(pf$sa^2, na.rm=TRUE))
  fp <- tryCatch(fisher.test(table(x$grp, x$sex))$p.value, error=function(e) NA_real_)
  data.table(study = st, top2_levels = paste(top2, collapse=" vs "), n_top2 = sum(pf$n),
             sex_table = paste0(pf$grp, ":F", pf$nf, "/M", pf$nm, collapse="; "),
             age_table = paste0(pf$grp, ":", round(pf$ma, 1), "±", round(pf$sa, 1), collapse="; "),
             sex_diff = round(sex_dp, 2), sex_fisher_p = signif(fp, 2), age_smd = round(smd, 2))
}
classify <- function(d) d[, design := fifelse(sex_diff >= 0.99, "completely_confounded",
                                       fifelse(sex_diff >= 0.3 | (!is.na(age_smd) & age_smd >= 0.8), "partially_confounded",
                                       fifelse(sex_diff >= 0.15 | (!is.na(age_smd) & age_smd >= 0.4), "mildly_imbalanced", "balanced")))]

## ---- 5. Blood candidates: strict / relaxed / sex-only ----
common <- quote(!tissue_category %in% c("cell_model") & n >= 12 & n <= 500 & !any_cell_line & !any_sc & !any_small &
                strategy == "RNA-Seq" & source == "TRANSCRIPTOMIC" & !is.na(group_key) & group_score > 0 & group_min_n >= 5)
strict  <- p[eval(common) & n_both/n >= 0.9]
relaxed <- p[eval(common) & n_both/n >= 0.6 & n_both/n < 0.9]
sexonly <- p[eval(common) & n_sex/n >= 0.9 & n_age/n < 0.6]
addimb <- function(d, tier) { if (!nrow(d)) return(d[, tier := character(0)]); d <- merge(d, rbindlist(Map(imb, d$study, d$group_key)), by="study"); classify(d); d[, tier := tier]; d }
strict <- addimb(strict, "strict"); relaxed <- addimb(relaxed, "relaxed_age"); sexonly <- addimb(sexonly, "sex_only")
cand <- rbindlist(list(strict, relaxed, sexonly), fill=TRUE)
setorder(cand, tier, -n_top2)
out_cols <- c("study","tier","design","flag_invitro","flag_repeat","flag_subtype","gsm_first","n","n_top2","n_sex","n_age","frac_female","age_min","age_max","pediatric",
              "group_key","group_score","top2_levels","group_counts","sex_table","age_table","sex_diff","sex_fisher_p","age_smd",
              "alt_keys","tissue_example","layout","platform","first_published","study_title","abstract")
out_cols <- c(out_cols[1:3], "tissue_category", out_cols[-(1:3)])
fwrite(cand[, ..out_cols], file.path(RES, "candidates_all_tissues_v1.tsv"), sep="\t")
cat("\n== candidates by tissue x tier ==\n"); print(table(cand$tissue_category, cand$tier))
cat("\n== strict candidates by tissue x design ==\n"); print(table(cand[tier=="strict", tissue_category], cand[tier=="strict", design]))
