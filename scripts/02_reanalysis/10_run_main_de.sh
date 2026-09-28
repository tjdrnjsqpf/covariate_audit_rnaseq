#!/usr/bin/env bash
# Main reanalysis runner: run 03_pilot_de_one_project.R for each row of the config TSV.
# Usage: bash 10_run_main_de.sh <config.tsv> <outdir> [parallel] [limit]
set -uo pipefail
PROJ=/var2/lsg/Claude_Code/covariate_audit_rnaseq
R=/var2/lsg/miniforge3/envs/cov_audit/bin/Rscript
CFG=$1; OUTDIR=$2; P=${3:-6}; LIMIT=${4:-0}
LOGD=$PROJ/logs/$(basename "$OUTDIR"); mkdir -p "$LOGD" "$PROJ/$OUTDIR"
rows=$(tail -n +2 "$PROJ/$CFG")
[ "$LIMIT" -gt 0 ] && rows=$(echo "$rows" | head -n "$LIMIT")
echo "$rows" | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5}' | \
xargs -P "$P" -d '\n' -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r s k c t e <<< "{}"
  [ -f '"$PROJ/$OUTDIR"'/$s/summary.tsv ] && { echo "skip $s"; exit 0; }
  '"$R"' '"$PROJ"'/scripts/00_pilot/03_pilot_de_one_project.R "$s" "$k" "$c" "$t" "$e" "'"$OUTDIR"'" > '"$LOGD"'/$s.log 2>&1 \
    && echo "done $s" || echo "FAIL $s"'
echo "ALL_FINISHED  ok=$(ls $PROJ/$OUTDIR/*/summary.tsv 2>/dev/null | wc -l)"
