#!/usr/bin/env bash
# Stream recount3 mouse gene_sums and extract only the counts of Xist + 4 Y-linked genes and the per-sample totals (no files saved).
# Usage: bash 03m_extract_mouse_sex_genes.sh [parallel]
set -uo pipefail
PROJ=/var2/lsg/Claude_Code/covariate_audit_rnaseq
OUT=$PROJ/data/processed/mouse_sexgenes; mkdir -p $OUT $PROJ/logs
LIST=$PROJ/data/processed/mouse_base_projects.txt
P=${1:-16}
extract() {
  p=$1; f=$OUT/$p.tsv; [ -s "$f" ] && return 0
  u="http://duffel.rail.bio/recount3/mouse/data_sources/sra/gene_sums/${p: -2}/$p/sra.gene_sums.$p.M023.gz"
  curl -sSL --retry 3 --retry-delay 3 "$u" | zcat 2>/dev/null | awk -F'\t' -v proj="$p" '
    /^##/ {next}
    $1=="gene_id" {n=NF; for(i=2;i<=NF;i++) id[i]=$i; next}
    { g=$1; sub(/\..*$/,"",g)
      for(i=2;i<=NF;i++) tot[i]+=$i
      if(g=="ENSMUSG00000086503"){for(i=2;i<=NF;i++) xist[i]=$i}
      else if(g=="ENSMUSG00000069045"){for(i=2;i<=NF;i++) ddx3y[i]=$i}
      else if(g=="ENSMUSG00000069049"){for(i=2;i<=NF;i++) eif2s3y[i]=$i}
      else if(g=="ENSMUSG00000056673"){for(i=2;i<=NF;i++) kdm5d[i]=$i}
      else if(g=="ENSMUSG00000068457"){for(i=2;i<=NF;i++) uty[i]=$i} }
    END { if(n>0) for(i=2;i<=n;i++) printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", proj, id[i], tot[i]+0, xist[i]+0, ddx3y[i]+0, eif2s3y[i]+0, kdm5d[i]+0, uty[i]+0 }' > "$f.tmp"
  if [ -s "$f.tmp" ]; then mv "$f.tmp" "$f"; else rm -f "$f.tmp"; echo "FAIL $p"; fi
}
export -f extract; export OUT
cat "$LIST" | xargs -P "$P" -I{} bash -c 'extract {}' 2>&1 | tee -a $PROJ/logs/03m_extract.log | grep -c FAIL || true
echo "done: $(ls $OUT/*.tsv 2>/dev/null | wc -l) / $(wc -l < $LIST)"
echo ALL_FINISHED
