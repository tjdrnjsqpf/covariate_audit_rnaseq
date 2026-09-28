#!/usr/bin/env Rscript
## E7. Semi-synthetic simulation, compact table for the main text: false-DEG fraction / precision / recall of M0 vs M1 by induced sex imbalance,
##     and the exact definition + value of the reported correlation between the observable change metric and actual error.
## NOTE 2026-09-18: the local results/simulation/sim_results_all.tsv holds only the first 3 human datasets (670 fits); the server copy has all 6 (1,050 fits, matches sim_summary.tsv).
##   Reported numbers: SIM_HUMAN=results/revision/inputs/sim_results_all_human_server_20260906_1050fits.tsv (verbatim copy of the server file).
## Inputs: results/simulation{,_mouse}/sim_results_all.tsv (one row per subsample fit). Truth = full-data adjusted model (human M2 ~group+sex+age; mouse M1 ~group+sex);
##   precision is against the loose truth (FDR<0.05), recall against the strict truth (FDR<0.01 & |logFC|>0.5). Only age_shift == 0 rows are used (pure sex imbalance).
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
SIM <- c(human=Sys.getenv("SIM_HUMAN", "results/simulation/sim_results_all.tsv"), mouse=Sys.getenv("SIM_MOUSE", "results/simulation_mouse/sim_results_all.tsv"))
d <- rbindlist(lapply(ORGS, function(o) { x <- fread(SIM[[o]]); x[, organism := o]; x }), fill=TRUE); d[, organism := factor(organism, levels=ORGS)]
s <- d[age_shift == 0]
for (m in c("M0", "M1")) { s[, (paste0(m, "_false")) := get(paste0(m, "_n_deg"))*(1 - get(paste0(m, "_precision")))]; s[get(paste0(m, "_n_deg")) == 0, (paste0(m, "_false")) := 0] }
## cell = dataset × n_per × sex_delta; pooled over replicates; then macro-averaged over cells within species
cell <- s[, .(reps=.N, M0_deg=mean(M0_n_deg), M1_deg=mean(M1_n_deg), M0_false=mean(M0_false), M1_false=mean(M1_false), M0_fdp=sum(M0_false)/sum(M0_n_deg), M1_fdp=sum(M1_false, na.rm=TRUE)/sum(M1_n_deg, na.rm=TRUE),
              M0_recall=mean(M0_recall), M1_recall=mean(M1_recall), M0_sexgene=mean(M0_sexgene_deg), M1_sexgene=mean(M1_sexgene_deg)), keyby=.(organism, study, n_per, sex_delta)]
out(cell, "E7_sim_cells.tsv")
## to keep rungs comparable, use only dataset×n_per combinations present at every delta 0..0.75 (M1) — complete ladders
lad <- cell[sex_delta <= 0.75 & !is.na(M1_deg), .N, by=.(organism, study, n_per)][N == 4]
tab <- function(cc, lab) cc[, .(set=lab, n_cells=.N, M0_deg=mean(M0_deg), M1_deg=mean(M1_deg), M0_false_deg=mean(M0_false), M1_false_deg=mean(M1_false), M0_false_frac=mean(M0_fdp, na.rm=TRUE), M1_false_frac=mean(M1_fdp, na.rm=TRUE),
                           M0_precision=1 - mean(M0_fdp, na.rm=TRUE), M1_precision=1 - mean(M1_fdp, na.rm=TRUE), M0_recall=mean(M0_recall), M1_recall=mean(M1_recall), M0_sexgene=mean(M0_sexgene), M1_sexgene=mean(M1_sexgene)), keyby=.(organism, sex_delta)]
T1 <- rbind(tab(merge(cell, lad[, .(organism, study, n_per)], by=c("organism", "study", "n_per")), "complete_ladders"), tab(cell, "all_cells")); out(T1, "E7_sim_table.tsv")
## fold increase of false DEGs in M0, delta 0.75 (and 1) vs 0, per dataset×n_per
FI <- dcast(cell, organism + study + n_per ~ sex_delta, value.var=c("M0_false", "M1_false")); FI[, `:=`(M0_fold_0.75=M0_false_0.75/M0_false_0, M0_fold_1=M0_false_1/M0_false_0, M1_fold_0.75=M1_false_0.75/M1_false_0)]; out(FI, "E7_sim_false_deg_fold.tsv")
## the reported correlation: Spearman, per dataset, over ALL subsample fits (all n_per, sex_delta, age_shift, replicates; complete obs), between
##   x = obs_jaccard_M0_M2 = Jaccard of the M0 and adjusted-model DE sets (raw Jaccard, NOT Δ = real − permuted; no permutation was run in the simulation), and
##   y = M0_precision = share of M0 DE genes that are in the loose truth set.  Positive = lower Jaccard goes with lower precision (more error).
cr <- function(x, lab) x[, { ok <- !is.na(obs_jaccard_M0_M2) & !is.na(M0_precision); f <- M0_n_deg*(1 - M0_precision)
  .(subset=lab, n_fits=sum(ok), rho_jaccard_vs_M0_precision=if (sum(ok) > 3) cor(obs_jaccard_M0_M2[ok], M0_precision[ok], method="spearman") else NA_real_,
    rho_jaccard_vs_M0_false_deg=if (sum(ok) > 3) cor(obs_jaccard_M0_M2[ok], f[ok], method="spearman") else NA_real_,
    rho_jaccard_vs_induced_delta=if (sum(ok) > 3) cor(obs_jaccard_M0_M2[ok], sex_delta[ok], method="spearman") else NA_real_) }, keyby=.(organism, study)]
CR <- rbind(cr(d, "all_fits(as reported)"), cr(d[age_shift == 0], "age_shift=0 only"), cr(d[age_shift == 0 & M0_n_deg >= 10], "age_shift=0, M0 DE>=10")); out(CR, "E7_sim_obs_vs_error_correlation.tsv")
options(width=250); print(T1[set == "complete_ladders"], digits=3); print(lad[, .N, by=organism]); print(FI[, .(organism, study, n_per, M0_false_0, M0_false_0.75, M0_fold_0.75, M1_false_0, M1_false_0.75)], digits=3); print(CR, digits=2)
