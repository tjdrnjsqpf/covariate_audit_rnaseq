#!/usr/bin/env Rscript
## Advisor comment Q5, step 1: sex recording and single-sex shares of SRA human/mouse bulk RNA-seq studies by year, from the NCBI runinfo table (2009-2019 and 2020-2026) and from the recount3 attribute parser (to 2019).
## Input: results/response/sra_runinfo_{sp}_{2020plus,to2019}.slim.tsv (55, 57), results/project_summary_all.tsv, results/mouse_project_summary_all.tsv, data/raw/recount3_metadata_mouse/sra.sra.*.MD.gz
## Output: results/response/M_metadata_by_year.tsv, M_metadata_by_era.tsv. Usage (server): Rscript scripts/07_response/58_post2020_metadata.R
suppressMessages(library(data.table)); PROJ <- Sys.getenv("COV_AUDIT_PROJ", "/path/to/covariate_audit_rnaseq"); setwd(PROJ)
summ <- function(p) p[, .(n_studies=.N, frac_sex_any=mean(n_sex > 0), frac_sex_90=mean(n_sex/n >= 0.9),
                          n_eligible=sum(elig), frac_single=mean(single[elig]), frac_female_only=sum(fonly[elig])/sum(single[elig]), n_single=sum(single[elig]))]
out <- list()
## ---- SRA runinfo series ('Sex' = the BioSample sex attribute): 2020–2026 and, when fetched, 2009–2019 (same source, so the two eras join) ----
for (run in list(list(suf="2020plus", era="SRA 2020-2026", y0=2020, y1=2026), list(suf="to2019", era="SRA 2009-2019 (runinfo)", y0=2009, y1=2019))) for (sp in c("human","mouse")) {
  f <- sprintf("results/response/sra_runinfo_%s_%s.slim.tsv", sp, run$suf); if (!file.exists(f)) next
  d <- fread(f, quote="", showProgress=FALSE)
  d[, sx := fifelse(tolower(Sex) %in% c("female","f"), "F", fifelse(tolower(Sex) %in% c("male","m"), "M", NA_character_))]
  b <- d[, .(sx = sx[!is.na(sx)][1], year = as.integer(min(substr(ReleaseDate, 1, 4)))), by=.(SRAStudy, BioSample)]
  p <- b[, .(n=.N, n_sex=sum(!is.na(sx)), n_F=sum(sx == "F", na.rm=TRUE), n_M=sum(sx == "M", na.rm=TRUE), year=min(year)), by=.(study=SRAStudy)]
  p[, `:=`(elig = n >= 8 & n <= 200 & n_sex/n >= 0.9, single = (n_F > 0 & n_M == 0) | (n_M > 0 & n_F == 0), fonly = n_F > 0 & n_M == 0)]
  p <- p[year >= run$y0 & year <= run$y1]
  out[[paste(sp, run$suf)]] <- rbind(cbind(organism=sp, era=run$era, p[, summ(.SD), by=year][order(year)]), cbind(organism=sp, era=run$era, year=NA_integer_, summ(p)))
  fwrite(p, sprintf("results/response/M_projects_%s_%s.tsv", sp, run$suf), sep="\t")
}
## ---- recount3 (–2019) ----
h <- fread("results/project_summary_all.tsv", select=c("study","n","n_sex","frac_female","first_published","strategy","source"))
h <- h[strategy == "RNA-Seq" & source == "TRANSCRIPTOMIC"]
h[, `:=`(year = as.integer(substr(first_published, 1, 4)), n_F = round(frac_female * n_sex), n_M = n_sex - round(frac_female * n_sex))]
h[is.na(n_F), `:=`(n_F = 0L, n_M = 0L)]
h[, `:=`(elig = n >= 8 & n <= 200 & n_sex/n >= 0.9, single = (n_F > 0 & n_M == 0) | (n_M > 0 & n_F == 0), fonly = n_F > 0 & n_M == 0)]
m <- fread("results/mouse_project_summary_all.tsv", select=c("study","n","n_sex","n_F","n_M","strategy","source"))[strategy == "RNA-Seq" & source == "TRANSCRIPTOMIC"]
fl <- list.files("data/raw/recount3_metadata_mouse", "^sra\\.sra\\..*\\.MD\\.gz$", full.names=TRUE)
dt <- rbindlist(lapply(fl, function(f) { x <- fread(f, select="run_published", showProgress=FALSE, quote=""); data.table(study=sub(".*sra\\.sra\\.(.*)\\.MD\\.gz", "\\1", f), first_published=min(x$run_published, na.rm=TRUE)) }))
m <- merge(m, dt, by="study", all.x=TRUE); m[, year := as.integer(substr(first_published, 1, 4))]
m[, `:=`(elig = n >= 8 & n <= 200 & n_sex/n >= 0.9, single = (n_F > 0 & n_M == 0) | (n_M > 0 & n_F == 0), fonly = n_F > 0 & n_M == 0)]
for (sp in c("human","mouse")) { p <- if (sp == "human") h else m; p <- p[!is.na(year) & year <= 2019]
  out[[paste0(sp, "_r3")]] <- rbind(cbind(organism=sp, era="recount3 parser (-2019)", p[, summ(.SD), by=year][order(year)]), cbind(organism=sp, era="recount3 parser (-2019)", year=NA_integer_, summ(p))) }
R <- rbindlist(out); fwrite(R, "results/response/M_metadata_by_year.tsv", sep="\t"); fwrite(R[is.na(year)], "results/response/M_metadata_by_era.tsv", sep="\t")
print(R[is.na(year)]); print(R[!is.na(year) & year >= 2012, .(organism, era, year, n_studies, frac_sex_any=round(frac_sex_any,3), frac_sex_90=round(frac_sex_90,3), n_eligible, frac_single=round(frac_single,3), frac_female_only=round(frac_female_only,3))])
