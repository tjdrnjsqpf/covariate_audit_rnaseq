#!/usr/bin/env bash
# Extract the three-level run / experiment / sample accession structure of each project.
# recount3 gene_sums columns are per run (external_id), so multiple runs that share the same
# sample_acc are technical replicates, not independent samples.
# Usage: bash 20_extract_accession_levels.sh <human|mouse>
set -uo pipefail
PROJ=/var2/lsg/Claude_Code/covariate_audit_rnaseq
ORG=${1:-human}
if [ "$ORG" = human ]; then MD=$PROJ/data/raw/recount3_metadata;       CFG=$PROJ/config/main_human_runs.tsv
else                       MD=$PROJ/data/raw/recount3_metadata_mouse; CFG=$PROJ/config/main_mouse_runs.tsv; fi
OUT=$PROJ/data/processed/accession_levels_$ORG.tsv
printf 'study\texternal_id\tsample_acc\texperiment_acc\tsample_title\trun_alias\n' > "$OUT"
n=0
tail -n +2 "$CFG" | cut -f1 | sort -u | while read -r s; do
  f=$MD/sra.sra.$s.MD.gz
  [ -f "$f" ] || { echo "MISSING $s" >&2; continue; }
  zcat "$f" | awk -F'\t' -v s="$s" '
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    { for (k in h) if (h[k]>NF) ;
      t=$h["sample_title"]; a=$h["run_alias"];
      gsub(/[\t\r\n]/, " ", t); gsub(/[\t\r\n]/, " ", a);
      printf "%s\t%s\t%s\t%s\t%s\t%s\n", s, $h["external_id"], $h["sample_acc"], $h["experiment_acc"], t, a }'
  n=$((n+1))
done >> "$OUT"
gzip -f "$OUT"
echo "wrote $OUT.gz  rows=$(zcat $OUT.gz | wc -l)  studies=$(zcat $OUT.gz | tail -n +2 | cut -f1 | sort -u | wc -l)"
