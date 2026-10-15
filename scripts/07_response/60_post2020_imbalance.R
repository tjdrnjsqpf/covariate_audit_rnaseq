#!/usr/bin/env Rscript
## Advisor comment Q5, step 2: distribution of the between-group sex-ratio difference |df| in recount3 (to 2019) versus SRA 2020-2026 under one metadata-only rule.
## Both eras use recorded sex only (no expression-based inference); the grouping key follows the rules of 02_parse (human) and 02m_parse (mouse).
## Input: data/processed/{,mouse_}sample_attributes_long.tsv.gz, results/response/biosample_attrs_{sp}_2020plus.tsv.gz (59), M_projects_{sp}_2020plus.tsv (58)
## Output: results/response/M_imbalance_by_era.tsv, M_imbalance_projects.tsv. Usage: Rscript scripts/07_response/60_post2020_imbalance.R [recount3|post2020|both]
suppressMessages(library(data.table)); PROJ <- Sys.getenv("COV_AUDIT_PROJ", "/path/to/covariate_audit_rnaseq"); setwd(PROJ)
mode <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(mode)) mode <- "both"
sex_keys <- c("sex","gender","sex/gender","gender/sex","patient sex","donor sex","subject sex","biological sex","donor gender")
H <- list(skip = c(sex_keys, "tissue","tissue type","organ","body site","cell type","cell line","age","source name","sample name","title","id","sample id","patient id","subject id","individual","donor","donor id","subject","replicate","batch","lane","library","barcode","run","biosample","description","sample description","molecule","organism","strain","isolate","ethnicity","race","bmi","isolation source","geo loc name","collection date","env broad scale","env local scale","env medium","lat lon","biomaterial provider","sample type","dev stage","host","organism part"),
          grp_re = "disease|condition|diagnos|group|status|phenotype|case|control|treatment|state|cohort|outcome|severity|infect|response|type|class|subtype", min_lev = 3)
M <- list(skip = c("sex","gender","age","strain","tissue","source name","cell type","sample type","tissue type","organ","title","id","sample id","replicate","batch","lane","library","barcode","run","biosample","description","organism","isolate","developmental stage","dev stage","isolation source","geo loc name","collection date","env broad scale","env local scale","env medium","lat lon","biomaterial provider","host","organism part","breed","cell line"),
          good_re = "genotype|treatment|condition|disease|diet|group|infect|status|phenotype|knockout|ko|wt|mutant|exposure|model|challenge|stimul",
          bad_re = "date|update|flowcell|lane|run|batch|barcode|library|cohort|time|day|week|hour|replicate|^id$|_id$| id$|name|accession|biosample|insdc|checklist|protocol|kit|platform|instrument|read|strand|layout|rin|passage|well|plate|age|sex|gender|strain", min_lev = 2)
norm_key <- function(k) gsub("_", " ", trimws(tolower(k)))
## al: study, sid, key, value (key before normalization); returns per-study sex design, grouping key, |Δf|
classify <- function(al, sp) {
  al[, key := norm_key(key)]; al <- al[nchar(value) <= 60 & value != ""]
  sx <- al[key %in% sex_keys]; sx[, sex := fifelse(grepl("^(f|female|woman|w)$", tolower(value)), "F", fifelse(grepl("^(m|male|man)$", tolower(value)), "M", NA_character_))]
  sx <- sx[!is.na(sex), .(sex = sex[1]), by=.(study, sid)]
  ns <- al[, .(n = uniqueN(sid)), by=study]; p <- merge(ns, sx[, .(n_sex=.N, n_F=sum(sex=="F"), n_M=sum(sex=="M")), by=study], by="study", all.x=TRUE)
  p[is.na(n_sex), `:=`(n_sex=0L, n_F=0L, n_M=0L)]
  p[, sex_design := fifelse(n_sex/n < 0.9, "unreported", fifelse(n_F > 0 & n_M == 0, "female_only", fifelse(n_M > 0 & n_F == 0, "male_only", "mixed")))]
  R <- if (sp == "human") H else M
  gk <- al[!key %in% R$skip & !grepl("^age", key), .(n_lev = uniqueN(value), min_lev = min(table(value))), by=.(study, key)][n_lev >= 2 & n_lev <= 6 & min_lev >= R$min_lev]
  if (sp == "human") gk[, score := grepl(R$grp_re, key) * 10 + (n_lev == 2) * 5 + pmin(min_lev, 20)/20] else gk[, score := 10*grepl(R$good_re, key) - 10*grepl(R$bad_re, key) + 5*(n_lev == 2) + pmin(min_lev, 10)/10]
  setorder(gk, study, -score); best <- gk[, .SD[1], by=study][score > 0, .(study, group_key = key, group_nlev = n_lev, group_min_n = min_lev)]
  p <- merge(p, best, by="study", all.x=TRUE)
  imb <- function(st, gkey) { a <- al[study == st & key == gkey, .(sid, grp = value)]; x <- merge(a, sx[study == st, .(sid, sex)], by="sid")
    top2 <- names(sort(table(x$grp), dec=TRUE))[1:2]; x <- x[grp %in% top2]; if (uniqueN(x$grp) < 2) return(NA_real_); pf <- x[, mean(sex == "F"), by=grp]$V1; abs(diff(pf)) }
  mx <- p[sex_design == "mixed" & !is.na(group_key) & n >= 8 & n <= 200]; mx[, sex_diff := mapply(imb, study, group_key)]
  p <- merge(p, mx[, .(study, sex_diff)], by="study", all.x=TRUE)
  p[, design := fifelse(is.na(sex_diff), NA_character_, fifelse(sex_diff >= 0.99, "completely_confounded", fifelse(sex_diff >= 0.3, "partially_confounded", fifelse(sex_diff >= 0.15, "mildly_imbalanced", "balanced"))))]
  p }
out <- list()
if (mode %in% c("recount3","both")) for (sp in c("human","mouse")) {
  al <- fread(sprintf("data/processed/%ssample_attributes_long.tsv.gz", if (sp == "mouse") "mouse_" else ""), showProgress=FALSE); setnames(al, "external_id", "sid")
  st <- fread(sprintf("results/%sproject_summary_all.tsv", if (sp == "mouse") "mouse_" else ""), select=c("study","strategy","source"))[strategy == "RNA-Seq" & source == "TRANSCRIPTOMIC", study]
  p <- classify(al[study %in% st], sp); p[, `:=`(organism = sp, era = "recount3 (-2019)", year = NA_integer_)]; out[[paste0(sp, "_r3")]] <- p; cat(sp, "recount3 done:", nrow(p), "\n") }
if (mode %in% c("post2020","both")) for (sp in c("human","mouse")) {
  f <- sprintf("results/response/biosample_attrs_%s_2020plus.tsv.gz", sp); if (!file.exists(f)) next
  al <- fread(f, quote="", showProgress=FALSE); setnames(al, "biosample", "sid"); yr <- fread(sprintf("results/response/M_projects_%s_2020plus.tsv", sp), select=c("study","year"))
  p <- classify(al, sp); p <- merge(p, yr, by="study", all.x=TRUE); p[, `:=`(organism = sp, era = "2020-2026")]; out[[paste0(sp, "_post")]] <- p; cat(sp, "post2020 done:", nrow(p), "\n") }
P <- rbindlist(out, use.names=TRUE, fill=TRUE); fwrite(P, "results/response/M_imbalance_projects.tsv", sep="\t")
summ <- function(p) { m <- p[!is.na(sex_diff)]; if (nrow(m) == 0) return(data.table(n_sexrec=NA_integer_)); data.table(n_sexrec=sum(p$n_sex/p$n >= 0.9 & p$n >= 8 & p$n <= 200), n_mixed=sum(p$sex_design == "mixed" & p$n >= 8 & p$n <= 200), n_with_group=nrow(m),
  frac_ge30=mean(m$sex_diff >= 0.3), frac_ge15=mean(m$sex_diff >= 0.15), frac_ge99=mean(m$sex_diff >= 0.99), median_df=median(m$sex_diff),
  ci_ge30 = { b <- binom.test(sum(m$sex_diff >= 0.3), nrow(m))$conf.int; sprintf("%.3f-%.3f", b[1], b[2]) }) }
S <- P[, summ(.SD), by=.(organism, era)]; if (any(P$era == "2020-2026")) S <- rbind(S, P[era == "2020-2026", summ(.SD), by=.(organism, era, year)][order(organism, year)], fill=TRUE)
fwrite(S, "results/response/M_imbalance_by_era.tsv", sep="\t"); print(S)
