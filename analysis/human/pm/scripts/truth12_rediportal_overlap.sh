#!/bin/bash
#$ -N t12_redi
#$ -cwd
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -l mem=8G
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
# REDIportal membership (position only, no strand) of every spectrum12 truth-gate site.
# Usage: qsub truth12_rediportal_overlap.sh [eval|eval_depth]   (default eval)
# Output: human/pm/gt_eval_v6/spectrum12/<eval dir>/truth12_sites_rediportal.tsv.gz (chr, pos, in_rediportal)
set -euo pipefail
V6=${SCRNA_DIR:-.}/human/pm/gt_eval_v6
EV=${1:-eval}
IN=$V6/spectrum12/$EV/truth12_sites.tsv.gz
OUT=$V6/spectrum12/$EV/truth12_sites_rediportal.tsv.gz
awk -F'\t' 'NR == FNR { if (FNR > 1) q[$2 FS $3] = 1; next } ($1 FS $3) in q { print $1 "\t" $3 }' \
  <(zcat "$IN") "$V6"/rediportal/bed/chr*.bed | sort -u > "${OUT%.gz}.hits"
zcat "$IN" | awk -F'\t' 'NR == FNR { h[$1 FS $2] = 1; next }
  FNR == 1 { print "chr\tpos\tin_rediportal"; next }
  { k = $2 FS $3; if (!(k in s)) { s[k] = 1; print $2 "\t" $3 "\t" ((k in h) ? "TRUE" : "FALSE") } }' "${OUT%.gz}.hits" - |
  gzip > "$OUT"
rm "${OUT%.gz}.hits"
zcat "$OUT" | awk 'NR > 1 { n++; t += ($3 == "TRUE") } END { print n " positions, " t " in REDIportal" }'
