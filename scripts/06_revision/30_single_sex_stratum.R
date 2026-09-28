#!/usr/bin/env Rscript
## E1. Single-sex stratum: female-only vs male-only by species, tissue, (human) excluding sex-specific tissue, (mouse) design_type. Wilson CIs.
## Run: LC_ALL=en_US.UTF-8 Rscript scripts/06_revision/30_single_sex_stratum.R
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
D <- load_summary(); S <- D[valid == TRUE & single_sex == TRUE]
S[, which := fifelse(n_female > 0 & n_male == 0, "female_only", fifelse(n_male > 0 & n_female == 0, "male_only", NA_character_))]; stopifnot(!anyNA(S$which))
## title/abstract text to flag sex-specific tissue or condition (human)
ps <- fread("results/project_summary_all.tsv", select=c("study", "study_title", "abstract", "tissue_example"))
S <- merge(S, ps, by="study", all.x=TRUE)
RX <- "breast|mammar|ovar|uter[iu]|endometri|myometri|cervi|placent|pregnan|gestation|preeclampsia|trophoblast|decidua|fallopian|vagin|menstr|prostat|\\btest(is|es|ic)|sperm|seminal"
S[, txt := tolower(paste(study_title, abstract, tissue_example, note, group_key))]
S[, sex_specific := organism == "human" & (tissue_category == "reproductive" | grepl(RX, txt))]
S[, sex_specific_strict := organism == "human" & (tissue_category == "reproductive" | grepl(RX, tolower(paste(study_title, tissue_example, group_key))))]
S[, sex_source := fifelse(n_sex_meta_missing >= n, "inferred_from_expression_only", "metadata_available")]
row <- function(d, stratum, level) { if (!nrow(d)) return(NULL); k <- sum(d$which == "female_only"); n <- nrow(d); w <- wilson(k, n)
  data.table(organism=d$organism[1], stratum=stratum, level=level, n_single_sex=n, female_only=k, male_only=n - k, pct_female_only=w$pct, lo=w$lo, hi=w$hi,
             binom_p_vs_half=binom.test(k, n)$p.value) }
R <- rbindlist(c(
  lapply(ORGS, function(o) row(S[organism == o], "all", "all")),
  lapply(ORGS, function(o) rbindlist(lapply(sort(unique(S[organism == o, tissue_category])), function(t) row(S[organism == o & tissue_category == t], "tissue_category", ifelse(t == "", "(blank)", t))))),
  list(row(S[organism == "human" & sex_specific == TRUE],  "human_sex_specific_tissue(title+abstract)", "sex_specific"),
       row(S[organism == "human" & sex_specific == FALSE], "human_sex_specific_tissue(title+abstract)", "excluding_sex_specific"),
       row(S[organism == "human" & sex_specific_strict == TRUE],  "human_sex_specific_tissue(title+tissue only)", "sex_specific"),
       row(S[organism == "human" & sex_specific_strict == FALSE], "human_sex_specific_tissue(title+tissue only)", "excluding_sex_specific"),
       row(S[organism == "human" & sex_specific == FALSE & tissue_category != "tumor" & !grepl("tumor_tissue", scope_flags)], "human_excl_sex_specific_and_tumor", "excluding")),
  lapply(ORGS, function(o) rbindlist(lapply(sort(unique(S[organism == o, design_type])), function(t) row(S[organism == o & design_type == t], "design_type", t)))),
  lapply(ORGS, function(o) rbindlist(lapply(sort(unique(S[organism == o, sex_source])), function(t) row(S[organism == o & sex_source == t], "sex_source", t))))))
out(R, "E1_single_sex_split.tsv")
## single-sex prevalence (denominator = valid studies) by species, with Wilson CI; and mouse design_type prevalence
V <- D[valid == TRUE]
P <- rbind(V[, c(.(stratum="all", level="all"), wilson(sum(single_sex), .N)), by=organism],
           V[, c(.(stratum="design_type"), wilson(sum(single_sex), .N)), by=.(organism, level=design_type)][, .(organism, stratum, level, k, n, pct, lo, hi)],
           V[, c(.(stratum="tissue_category"), wilson(sum(single_sex), .N)), by=.(organism, level=tissue_category)][, .(organism, stratum, level, k, n, pct, lo, hi)])
out(P, "E1_single_sex_prevalence.tsv")
out(S[organism == "human", .(study, which, n, tissue_category, design_type, sex_specific, sex_specific_strict, sex_source, study_title)], "E1_human_single_sex_listing.tsv")
print(R[stratum %in% c("all") | grepl("human_", stratum) | stratum == "design_type" | stratum == "sex_source"], digits=3)
## mouse: genotype vs treatment difference
m <- S[organism == "mouse" & design_type %in% c("genotype_comparison", "in_vivo_treatment_arms")]; print(table(m$design_type, m$which)); print(fisher.test(table(m$design_type, m$which)))
