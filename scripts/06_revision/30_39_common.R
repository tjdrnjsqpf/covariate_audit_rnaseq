## Shared helpers for revision analyses 30-39 (Task E). No refitting; reads existing summary tables.
## Paths are parameterised so the scripts can be re-run after the refit:
##   COV_AUDIT_PROJ   project root (default: local Dropbox path if it exists, else server path)
##   SUMMARY_HUMAN / SUMMARY_MOUSE   summary tables (default results/main_{org}_summary_v2.tsv)
##   ASET_HUMAN / ASET_MOUSE         analysis-set lists (default config/c2_analysis_set_{org}.txt)
##   OUTDIR                          output dir (default results/revision)
suppressMessages(library(data.table))
.loc <- Sys.getenv("COV_AUDIT_PROJ", unset = getwd())   # local checkout (set COV_AUDIT_PROJ); falls back to the server path below
PROJ <- Sys.getenv("COV_AUDIT_PROJ", if (dir.exists(.loc)) .loc else "/var2/lsg/Claude_Code/covariate_audit_rnaseq"); setwd(PROJ)
ORGS <- c("human", "mouse")
SUMMARY <- c(human=Sys.getenv("SUMMARY_HUMAN", "results/main_human_summary_v2.tsv"), mouse=Sys.getenv("SUMMARY_MOUSE", "results/main_mouse_summary_v2.tsv"))
ASET    <- c(human=Sys.getenv("ASET_HUMAN", "config/c2_analysis_set_human.txt"),    mouse=Sys.getenv("ASET_MOUSE", "config/c2_analysis_set_mouse.txt"))
OUTDIR  <- Sys.getenv("OUTDIR", "results/revision"); dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
CLS <- c("balanced", "mildly_imbalanced", "partially_confounded")
load_summary <- function() { D <- rbindlist(lapply(ORGS, function(o) { d <- fread(SUMMARY[[o]]); d[, organism := o]; d[, analysis := study %in% readLines(ASET[[o]])]; d }), fill=TRUE)
  stopifnot(!anyDuplicated(D[, paste(organism, study)])); D[, organism := factor(organism, levels=ORGS)]; D }
## Wilson score interval
wilson <- function(k, n, conf=0.95) { z <- qnorm(1 - (1 - conf)/2); p <- k/n; den <- 1 + z^2/n; ctr <- (p + z^2/(2*n))/den; hw <- z*sqrt(p*(1 - p)/n + z^2/(4*n^2))/den
  data.table(k=k, n=n, pct=100*p, lo=100*pmax(0, ctr - hw), hi=100*pmin(1, ctr + hw)) }
## bootstrap CI of a median (same recipe as 24_letter_figures.R: 1000 resamples, seed 1)
bci <- function(x, B=1000) { x <- x[!is.na(x)]; if (length(x) < 3) return(list(med=if (length(x)) median(x) else NA_real_, lo=NA_real_, hi=NA_real_)); set.seed(1)
  b <- replicate(B, median(sample(x, replace=TRUE))); list(med=median(x), lo=unname(quantile(b, .025)), hi=unname(quantile(b, .975))) }
nstrata <- function(n) cut(n, c(-Inf, 11, 19, 39, Inf), labels=c("<12", "12-19", "20-39", ">=40"))
out <- function(x, f) { fwrite(x, file.path(OUTDIR, f), sep="\t"); cat("wrote", file.path(OUTDIR, f), "\n") }
