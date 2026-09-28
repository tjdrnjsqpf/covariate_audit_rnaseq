#!/usr/bin/env bash
# Download sample metadata (sra.sra.*.MD.gz) and recount_pred for all recount3 mouse SRA projects.
# Usage: bash scripts/00_pilot/01_download_recount3_metadata.sh
set -euo pipefail
PROJ=/var2/lsg/Claude_Code/covariate_audit_rnaseq
OUT=$PROJ/data/raw/recount3_metadata_mouse_mouse
LIST=$PROJ/data/raw/recount3_mouse_sra_projects.tsv
BASE=http://duffel.rail.bio/recount3/mouse/data_sources/sra/metadata
mkdir -p "$OUT" "$PROJ/logs"

# 1) Project list
/var2/lsg/miniforge3/envs/recount3/bin/Rscript -e '
suppressMessages(library(recount3))
hp <- available_projects(organism="mouse")
sra <- subset(hp, file_source=="sra" & project_type=="data_sources")
write.table(sra[, c("project","n_samples")], "'"$LIST"'", sep="\t", quote=FALSE, row.names=FALSE)
cat("projects:", nrow(sra), "\n")' 2>&1 | grep -v -E "caching|adding rname|^$"

# 2) Build URL list (sra.sra + recount_pred); skip files already downloaded
tail -n +2 "$LIST" | cut -f1 | while read -r p; do
  sub=${p: -2}
  for t in sra recount_pred; do
    f="$OUT/sra.$t.$p.MD.gz"
    [ -s "$f" ] || echo "$BASE/$sub/$p/sra.$t.$p.MD.gz $f"
  done
done > "$OUT/../urls_todo_mouse.txt"
echo "to download: $(wc -l < "$OUT/../urls_todo_mouse.txt")"

# 3) Parallel download
cat "$OUT/../urls_todo_mouse.txt" | xargs -P 24 -n 2 sh -c 'curl -sSL --retry 3 --retry-delay 2 -o "$1" "$0" || echo "FAIL $0"' \
  2>&1 | tee "$PROJ/logs/download_metadata_mouse.log" | grep -c FAIL || true
echo "done: $(ls "$OUT" | wc -l) files"
