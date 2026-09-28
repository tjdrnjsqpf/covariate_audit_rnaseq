#!/usr/bin/env bash
# Run the 03 script in parallel for each row of config/pilot_core_runs.tsv (default 4 concurrent jobs).
set -uo pipefail
PROJ=/var2/lsg/Claude_Code/covariate_audit_rnaseq
R=/var2/lsg/miniforge3/envs/cov_audit/bin/Rscript
P=${1:-4}
mkdir -p "$PROJ/logs/pilot"
tail -n +2 "${CFG:-$PROJ/config/pilot_core_runs.tsv}" | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5}' | \
xargs -P "$P" -d '\n' -I{} bash -c '
  IFS=$'"'"'\t'"'"' read -r s k c t e <<< "{}"
  echo "[$(date +%H:%M:%S)] start $s"
  '"$R"' '"$PROJ"'/scripts/00_pilot/03_pilot_de_one_project.R "$s" "$k" "$c" "$t" "$e" > '"$PROJ"'/logs/pilot/$s.log 2>&1 \
    && echo "[$(date +%H:%M:%S)] done $s" || echo "[$(date +%H:%M:%S)] FAIL $s"'
echo "ALL_FINISHED"
