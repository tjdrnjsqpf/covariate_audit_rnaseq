#!/usr/bin/env Rscript
## Pseudoreplication audit
## Murad 2026 (bioRxiv 10.64898/2026.08.05.742891) reported that 15.1% of studies pass the replicate threshold when replicates
## are counted by accession but only 4.1% when counted by individual (3.73x). Here we count directly at the
## run / sample_acc / donor levels to check whether our pipeline is exposed to the same inflation.
## Usage: Rscript 21_pseudoreplication_audit.R <human|mouse> [summary.tsv]
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
ORG  <- if (length(args) >= 1) args[1] else "human"
SUMF <- if (length(args) >= 2) args[2] else sprintf("results/main_%s_summary.tsv", ORG)
PROJ <- "/var2/lsg/Claude_Code/covariate_audit_rnaseq"
setwd(PROJ); dir.create("results/qc", showWarnings=FALSE)

acc <- fread(sprintf("data/processed/accession_levels_%s.tsv.gz", ORG), na.strings=c("","NA"))
att <- fread(if (ORG == "human") "data/processed/sample_attributes_long.tsv.gz"
             else "data/processed/mouse_sample_attributes_long.tsv.gz", na.strings=c("","NA"))
att <- att[study %in% unique(acc$study)]
## Read each study's group variable (group_key) and exclude it from the individual-identifier candidates
cfgf <- sprintf("config/main_%s_runs.tsv", ORG)
GKEY <- if (file.exists(cfgf)) { g <- fread(cfgf, colClasses="character"); paste(g$study, tolower(trimws(g$group_key))) } else character()
att[, key := tolower(trimws(key))]

## ---- Level 1: run -> sample_acc (technical replicates) ----
lv <- acc[, .(n_run = .N,
              n_sample_acc = uniqueN(sample_acc[!is.na(sample_acc)]),
              n_experiment = uniqueN(experiment_acc[!is.na(experiment_acc)])), by=study]
lv[n_sample_acc == 0, n_sample_acc := NA_integer_]

## Additional duplicates caught by explicit replicate suffixes in titles (_rep1, technical replicate 2, etc.)
acc[, stem := tolower(fifelse(is.na(sample_title), "", sample_title))]
acc[, stem := gsub("[ _.:-]*(technical|tech|biological|bio)?[ _.:-]*(replicate|replicates|rep|repeat|run|lane)[ _.:-]*[0-9]+[ _.:-]*$", "", stem)]
st <- acc[stem != "", .(n_title_stem = uniqueN(stem), n_titled = .N), by=study]
lv <- merge(lv, st, by="study", all.x=TRUE)

## ---- Level 2: sample_acc -> donor (individual identifier) ----
## Accept only keys that point to an individual. Group labels such as 'patient group' cause over-merging and are excluded.
ID_STRICT <- "^(submitted_)?(subject|donor|individual|patient|participant|volunteer|animal|mouse|rat)([ _.:-]?(id|ids|no|num|number|name|code|label))?$"
ID_LOOSE  <- "(subject|donor|individual|patient|participant|volunteer|animal)[ _.:-]?(id|no|num|number|code|label)"
ID_BAD    <- "group|status|type|class|categor|condition|disease|diagnos|treat|sex|gender|age|strain|genotype|tissue|cell|time|visit|batch|stage|grade|response|outcome|cohort|arm"
PRIO <- c("submitted_subject_id","donor","subject","individual","patient","participant","volunteer","animal","mouse","rat")

cand <- att[!is.na(value) & value != "" &
            (grepl(ID_STRICT, key) | grepl(ID_LOOSE, key)) & !grepl(ID_BAD, key)]
cand[, prio := {
  p <- rep(99L, .N)
  for (i in seq_along(PRIO)) { idx <- which(p == 99L)
    if (length(idx)) { hit <- grepl(PRIO[i], key[idx], fixed=FALSE); p[idx[hit]] <- i } }
  p }]
## Also check whether the values themselves look like individual identifiers (fixed 2026-09-10)
## Even if the key name is patient/individual, a value such as "Kerataconus patient" is a group/tissue label, not an individual identifier.
cand[, val_len := nchar(value)][, val_sp := lengths(regmatches(value, gregexpr(" ", value)))]
ck <- cand[, .(cov = uniqueN(external_id), nd = uniqueN(value),
               med_len = as.numeric(median(val_len)), frac_wordy = as.numeric(mean(val_sp >= 2))), by=.(study, key, prio)]
ck <- merge(ck, lv[, .(study, n_run)], by="study")
## A column used as the group variable in the design cannot be an individual identifier
if (exists("GKEY") && length(GKEY)) ck <- ck[!paste(study, key) %in% GKEY]
## Acceptance criteria: covers >=80% of samples, has >=4 distinct values (a study with only 2–3 "individuals" is more likely a group label),
##            and values are not descriptive text (median length ≤ 20 chars, fewer than half of the values contain >=2 spaces)
ck <- ck[cov >= 0.8 * n_run & nd >= 4 & med_len <= 20 & frac_wordy < 0.5]
setorder(ck, study, prio, -cov, nd)
best <- ck[, .SD[1], by=study][, .(study, donor_key = key, n_donor = nd)]
## False-positive flag: if the values are short serial numbers like 1..k or single characters, they are likely renumbered within each group (model, litter)
## (e.g. SRP161570 'mouse number' 1..21 repeated in each of 6 infection models). Such keys underestimate the number of individuals.
vv <- merge(cand, best, by.x=c("study","key"), by.y=c("study","donor_key"))[, .(study, value)]
seqflag <- vv[, {
  u <- unique(value); num <- suppressWarnings(as.integer(u))
  .(donor_key_suspect = (all(!is.na(num)) && max(num) <= 2 * length(num) && min(num) >= 0) ||
                        all(nchar(gsub("^(mouse|animal|rat|subject|donor|patient)[ _-]*", "", u)) <= 2)) }, by=study]
best <- merge(best, seqflag, by="study", all.x=TRUE)
lv <- merge(lv, best, by="study", all.x=TRUE)
lv[is.na(donor_key_suspect), donor_key_suspect := FALSE]
lv[, donor_resolved := !is.na(n_donor) & donor_key_suspect == FALSE]   # suspect keys are treated as unresolved

## ---- Inflation ratios ----
lv[, r_run_per_sample  := round(n_run / n_sample_acc, 3)]
lv[, r_sample_per_donor := round(n_sample_acc / n_donor, 3)]
lv[, r_run_per_donor   := round(n_run / n_donor, 3)]

## ---- Compare with the n actually used by our pipeline ----
if (file.exists(SUMF)) {
  sm <- fread(SUMF)[, .(study, our_n = n, our_n_before = n_before_subject_collapse,
                        our_options = options, our_ncase = n_case, our_nctrl = n_control)]
  lv <- merge(lv, sm, by="study", all.x=TRUE)
  lv[, used_subj_collapse := grepl("subj:", fifelse(is.na(our_options), "", as.character(our_options)))]
  lv[, overcount_vs_donor := !is.na(our_n) & !is.na(n_donor) & our_n > n_donor]
  lv[, overcount_vs_sample := !is.na(our_n) & !is.na(n_sample_acc) & our_n > n_sample_acc]
}
fwrite(lv, sprintf("results/qc/pseudorep_%s.tsv", ORG), sep="\t")

## ================= Report =================
N <- nrow(lv)
cat(sprintf("\n########## Pseudoreplication audit: %s (%d studies) ##########\n", ORG, N))

cat("\n=== Level 1: run -> sample_acc (technical replicates) ===\n")
infl1 <- lv[!is.na(n_sample_acc) & n_run > n_sample_acc]
cat(sprintf("Studies with more runs than sample_acc: %d / %d (%.1f%%)\n", nrow(infl1), N, 100*nrow(infl1)/N))
cat(sprintf("Total runs %d vs total sample_acc %d — population-level inflation %.2fx\n",
            sum(lv$n_run), sum(lv$n_sample_acc, na.rm=TRUE), sum(lv$n_run)/sum(lv$n_sample_acc, na.rm=TRUE)))
if (nrow(infl1)) print(infl1[order(-r_run_per_sample)][1:min(10,.N), .(study, n_run, n_sample_acc, r_run_per_sample)])

cat("\n=== Level 2: sample_acc -> donor (multiple samples from the same individual) ===\n")
cat(sprintf("Studies with a resolved individual identifier: %d / %d (%.1f%%)  [excluded as suspect serial-number keys: %d]\n",
            sum(lv$donor_resolved), N, 100*mean(lv$donor_resolved), sum(lv$donor_key_suspect)))
res <- lv[donor_resolved == TRUE]
if (nrow(res)) {
  infl2 <- res[n_sample_acc > n_donor]
  cat(sprintf("  of which studies with samples > individuals: %d (%.1f%% of resolved)\n", nrow(infl2), 100*nrow(infl2)/nrow(res)))
  cat(sprintf("  run/donor ratio among resolved studies: median %.2f, 90th percentile %.2f, max %.2f\n",
      median(res$r_run_per_donor), quantile(res$r_run_per_donor, .9), max(res$r_run_per_donor)))
  cat("  Most extreme cases:\n")
  print(res[order(-r_run_per_donor)][1:min(10,.N), .(study, donor_key, n_run, n_sample_acc, n_donor, r_run_per_donor)])
}
## Manski-style bounds: for unresolved studies the number of individuals lies within [1, n_sample_acc]
lo <- nrow(lv[donor_resolved == TRUE & n_sample_acc > n_donor])
hi <- lo + sum(!lv$donor_resolved)
cat(sprintf("\nBounds on the population-wide fraction of studies with 'samples > individuals': %.1f%% ~ %.1f%%\n", 100*lo/N, 100*hi/N))
cat("  (lower = only resolved and actually inflated studies; upper = assumes every unresolved study is inflated)\n")

cat("\n=== Explicit replicate suffixes in titles ===\n")
tt <- lv[!is.na(n_title_stem) & n_title_stem < n_titled]
cat(sprintf("Studies merged via rep/replicate numbers in titles: %d (%.1f%%)\n", nrow(tt), 100*nrow(tt)/N))

if ("our_n" %in% names(lv)) {
  cat("\n=== Actual exposure of our analysis ===\n")
  a <- lv[!is.na(our_n)]
  cat(sprintf("Studies with results: %d,  collapsed via the subj: option %d (%.1f%%)\n",
              nrow(a), sum(a$used_subj_collapse), 100*mean(a$used_subj_collapse)))
  cat(sprintf("Our n > number of sample_acc (technical replicates counted as independent samples): %d (%.1f%%)\n",
              sum(a$overcount_vs_sample), 100*mean(a$overcount_vs_sample)))
  ad <- a[donor_resolved == TRUE]
  cat(sprintf("Of %d studies with resolved individuals, our n > number of individuals: %d (%.1f%%)\n",
              nrow(ad), sum(ad$overcount_vs_donor), 100*mean(ad$overcount_vs_donor)))
  MINN <- 8L
  cat(sprintf("Studies falling below the minimum sample size (n>=%d) when recounted per individual: %d / %d (%.1f%%)\n",
              MINN, sum(ad$n_donor < MINN), nrow(ad), 100*mean(ad$n_donor < MINN)))
  if (sum(ad$overcount_vs_donor)) {
    cat("\n  Candidates to exclude or collapse in the reanalysis (our n exceeds the number of individuals):\n")
    print(ad[overcount_vs_donor == TRUE][order(-our_n/n_donor)][1:min(15,.N),
             .(study, donor_key, our_n, n_donor, n_sample_acc, used_subj_collapse, ratio = round(our_n/n_donor,2))])
  }
}
cat("\nTable saved: results/qc/pseudorep_", ORG, ".tsv\n", sep="")
