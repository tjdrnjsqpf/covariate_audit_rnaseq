#!/usr/bin/env Rscript
## E4. Wilson CIs for every Fig 2C percentage (recomputed from gate2_merged_*.tsv, so it updates with the data) and
##     the ceiling on the true sex-adjustment rate given 0 of 14 resolved 'unclear' records had adjusted.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "30_39_common.R"))
GATE2 <- c(human=Sys.getenv("GATE2_HUMAN", "results/gate2/gate2_merged_human.tsv"), mouse=Sys.getenv("GATE2_MOUSE", "results/gate2/gate2_merged_mouse.tsv"))
K_UNCLEAR_CHECKED <- as.integer(Sys.getenv("UNCLEAR_CHECKED", "14")); K_UNCLEAR_ADJ <- as.integer(Sys.getenv("UNCLEAR_ADJUSTED", "0"))   # hand verification: 0 of 14 resolved unclear-sex records had adjusted
W <- list(); C <- list()
for (o in ORGS) { g <- fread(GATE2[[o]])[text_relevant == TRUE]; pc <- g[design_sex == "partially_confounded"]; pcr <- pc[sex_reported == "yes"]
  it <- list(c("Sex in DE model", sum(g$sex_adjusted == "yes"), nrow(g)), c("Age in DE model", sum(g$age_adjusted == "yes"), nrow(g)), c("Batch / latent variable", sum(g$batch_adj == "yes"), nrow(g)),
             c("Group sex composition reported", sum(g$sex_reported == "yes"), nrow(g)), c("Sex in DE model - confounded designs only", sum(pc$sex_adjusted == "yes"), nrow(pc)),
             c("Sex composition reported - confounded designs only", nrow(pcr), nrow(pc)), c("Sex in DE model - confounded AND composition reported", sum(pcr$sex_adjusted == "yes"), nrow(pcr)),
             c("Sex model 'unclear'", sum(g$sex_adjusted == "unclear"), nrow(g)))
  for (x in it) { k <- as.integer(x[2]); n <- as.integer(x[3]); w <- wilson(k, n); cp <- binom.test(k, n)$conf.int
    W[[length(W) + 1]] <- data.table(organism=o, item=x[1], k=k, n=n, pct=w$pct, wilson_lo=w$lo, wilson_hi=w$hi, exact_lo=100*cp[1], exact_hi=100*cp[2]) }
  ## ceiling: true rate = P(yes) + P(unclear) * theta, theta = adjustment rate among 'unclear' papers; 0/14 observed.
  n <- nrow(g); ky <- sum(g$sex_adjusted == "yes"); ku <- sum(g$sex_adjusted == "unclear"); m <- K_UNCLEAR_CHECKED; stopifnot(K_UNCLEAR_ADJ == 0)
  th <- c(rule_of_three=3/m, exact_one_sided_95=1 - 0.05^(1/m), exact_two_sided_95=1 - 0.025^(1/m), wilson_two_sided_95=wilson(0, m)$hi/100, all_unclear_adjusted=1)
  C[[o]] <- data.table(organism=o, n_pubs=n, k_yes=ky, k_unclear=ku, pct_yes=100*ky/n, pct_unclear=100*ku/n, bound=names(th), theta_upper_pct=100*th, ceiling_pct=100*(ky + ku*th)/n) }
W <- rbindlist(W); C <- rbindlist(C); out(W, "E4_fig2c_wilson.tsv"); out(C, "E4_adjustment_rate_ceiling.tsv")
## 2/32 vs overall 19/217 (and mouse 0/18 vs 5/183): Fisher exact, confounded vs not confounded
FT <- rbindlist(lapply(ORGS, function(o) { g <- fread(GATE2[[o]])[text_relevant == TRUE]; t <- table(pc=g$design_sex == "partially_confounded", adj=g$sex_adjusted == "yes"); f <- fisher.test(t)
  data.table(organism=o, adj_pc=t["TRUE", "TRUE"], n_pc=sum(t["TRUE", ]), adj_other=t["FALSE", "TRUE"], n_other=sum(t["FALSE", ]), odds_ratio=f$estimate, or_lo=f$conf.int[1], or_hi=f$conf.int[2], fisher_p=f$p.value) })); out(FT, "E4_confounded_vs_other_fisher.tsv")
options(width=200); print(W, digits=3); print(C, digits=3); print(FT, digits=3)
