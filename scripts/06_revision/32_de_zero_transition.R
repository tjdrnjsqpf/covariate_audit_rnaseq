#!/usr/bin/env Rscript
## E3. (M0 DE = 0) × (M1 DE = 0) transition table by |Δf| class and species, over all valid mixed-sex studies with M1 estimable.
## M1_n_deg_ref = DE count under M0 (reference), M1_n_deg_alt = DE count under M1; same for M0perm_ (permuted-sex control).
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
D <- load_summary(); V <- D[valid == TRUE & single_sex == FALSE & design_sex %in% CLS & !is.na(M1_n_deg_ref) & !is.na(M1_n_deg_alt)]
V[, cls := factor(design_sex, levels=CLS)]
tr <- function(ref, alt) fifelse(ref == 0 & alt == 0, "0->0", fifelse(ref == 0 & alt > 0, "0->pos", fifelse(ref > 0 & alt == 0, "pos->0", "pos->pos")))
V[, t_real := tr(M1_n_deg_ref, M1_n_deg_alt)]; V[, t_perm := tr(M0perm_n_deg_ref, M0perm_n_deg_alt)]
tab <- function(col, lab) { x <- V[!is.na(get(col)), .N, by=.(organism, cls, transition=get(col))]; x <- dcast(x, organism + cls ~ transition, value.var="N", fill=0)
  for (k in c("0->0", "0->pos", "pos->0", "pos->pos")) if (!k %in% names(x)) x[, (k) := 0L]
  x[, n := `0->0` + `0->pos` + `pos->0` + `pos->pos`]; x[, model := lab]
  g <- x[, wilson(`0->pos`, `0->0` + `0->pos`)]; l <- x[, wilson(`pos->0`, `pos->0` + `pos->pos`)]
  x[, `:=`(M0_zero=`0->0` + `0->pos`, pct_gain_given_M0zero=g$pct, gain_lo=g$lo, gain_hi=g$hi, M0_pos=`pos->0` + `pos->pos`, pct_loss_given_M0pos=l$pct, loss_lo=l$lo, loss_hi=l$hi)]
  x[, mcnemar_p := mapply(function(b, c) if (b + c == 0) NA_real_ else binom.test(b, b + c)$p.value, `0->pos`, `pos->0`)]; x[] }
R <- rbind(tab("t_real", "M1_real_sex"), tab("t_perm", "M0perm_permuted_sex")); setcolorder(R, c("organism", "cls", "model", "n")); setkey(R, organism, cls, model)
out(R, "E3_de_zero_transition.tsv")
## sizes of the gained sets (0 -> pos): how many DE genes appear
G <- V[t_real == "0->pos", .(n_studies=.N, M1_deg_median=as.double(median(M1_n_deg_alt)), M1_deg_max=max(M1_n_deg_alt), n_with_ge10=sum(M1_n_deg_alt >= 10)), keyby=.(organism, cls)]; out(G, "E3_gained_set_sizes.tsv")
options(width=250); print(R, digits=3); print(G)
