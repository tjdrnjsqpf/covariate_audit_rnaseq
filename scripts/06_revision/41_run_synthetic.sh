#!/usr/bin/env bash
## S runner: run 40_synthetic_imbalance.R in parallel over the eligible-study list. Usage: bash 41_run_synthetic.sh <human|mouse> <list.txt> [parallel=8]
set -uo pipefail; PROJ=/var2/lsg/Claude_Code/covariate_audit_rnaseq; R=/var2/lsg/miniforge3/envs/cov_audit/bin/Rscript
ORG=$1; LIST=$2; P=${3:-8}; LOGD=$PROJ/logs/synthetic_$ORG; mkdir -p "$LOGD" "$PROJ/results/revision/synthetic"
cat "$LIST" | xargs -P "$P" -I{} bash -c '[ -f '"$PROJ"'/results/revision/synthetic/{}.tsv ] && { echo "skip {}"; exit 0; }; nice -n 10 '"$R"' '"$PROJ"'/scripts/06_revision/40_synthetic_imbalance.R '"$ORG"' {} 10 > '"$LOGD"'/{}.log 2>&1 && echo "done {}" || echo "FAIL {}"'
echo "SYNTH_FINISHED ok=$(ls $PROJ/results/revision/synthetic/*.tsv 2>/dev/null | wc -l)"
