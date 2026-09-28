#!/usr/bin/env bash
## Usage: bash 47_run_content.sh <human|mouse> <list.txt> [parallel=8] [B=20]
set -uo pipefail; PROJ=/var2/lsg/Claude_Code/covariate_audit_rnaseq; R=/var2/lsg/miniforge3/envs/cov_audit/bin/Rscript
ORG=$1; LIST=$2; P=${3:-8}; B=${4:-20}; LOGD=$PROJ/logs/content_$ORG; mkdir -p "$LOGD" "$PROJ/results/revision/content"
cat "$LIST" | xargs -P "$P" -I{} bash -c '[ -f '"$PROJ"'/results/revision/content/{}.tsv ] && { echo "skip {}"; exit 0; }; nice -n 10 '"$R"' '"$PROJ"'/scripts/06_revision/46_content_of_change.R '"$ORG"' {} '"$B"' > '"$LOGD"'/{}.log 2>&1 && echo "done {}" || echo "FAIL {}"'
echo "CONTENT_FINISHED ok=$(ls $PROJ/results/revision/content/*.tsv 2>/dev/null | grep -vc _perms)"
