#!/usr/bin/env bash
## Usage: bash 56_run_engine.sh <list.tsv: study<TAB>organism> [parallel=8]
set -uo pipefail; PROJ=/path/to/covariate_audit_rnaseq; R=/var2/lsg/miniforge3/envs/cov_audit/bin/Rscript; LIST=$1; P=${2:-8}; mkdir -p $PROJ/logs/engine $PROJ/results/response/engine
cat "$LIST" | xargs -P "$P" -L 1 bash -c 's=$0; o=$1; [ -f '"$PROJ"'/results/response/engine/$s.tsv ] && { echo "skip $s"; exit 0; }; nice -n 10 '"$R"' '"$PROJ"'/scripts/07_response/54_engine_check.R $o $s > '"$PROJ"'/logs/engine/$s.log 2>&1 && echo "done $s" || echo "FAIL $s"'
echo "ENGINE_FINISHED ok=$(ls $PROJ/results/response/engine/*.tsv 2>/dev/null | wc -l)"
