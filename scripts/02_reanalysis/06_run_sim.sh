#!/usr/bin/env bash
PROJ=/var2/lsg/Claude_Code/covariate_audit_rnaseq; R=/var2/lsg/miniforge3/envs/cov_audit/bin/Rscript
mkdir -p $PROJ/logs/sim
for s in "$@"; do
  ( $R $PROJ/scripts/02_reanalysis/06_semisynthetic_sim.R $s 10 > $PROJ/logs/sim/$s.log 2>&1 && echo "done $s" || echo "FAIL $s" ) &
done
wait; echo ALL_FINISHED
