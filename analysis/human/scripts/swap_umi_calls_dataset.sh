#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=4G
#$ -o ${SCRNA_DIR:-.}/tmp/s5/swap_umi_calls_dataset.$JOB_ID.log
#$ -j y
# swap_umi_calls_dataset.sh   (v10)
# Non-ground-truth UMI datasets (jzlxv3, jzlcel): replace the unseeded umi_tools outputs in human/<ds>/results_v8 with
# the seeded run (run_pipeline_dataset_v8.sh UMI_ONLY=1 -> results_v8_umi).
#   1. archive (mv) every umi_tools output of results_v8 to human/<ds>_unseeded/results_v8/<same path>, logged in
#      HANDOFF/MOVED.tsv: <unit>/0{3_variants,4_filtered}/<caller>/*UmiDup, <sample>/02_barcode_filtered/*UmiDup,
#      <sample>/02_dedup/*_umiDup.bam{,.bai} and <sample>/02_dedup/subsample (umi_tools subsample of the 3' bias test)
#   2. copy (cp -a) the same paths from results_v8_umi into results_v8 (the Cellsnp-lite hybrid calls sit inside
#      04_filtered/cellsnp/*UmiDup in both trees)
# Refuses unless every seeded path has an unseeded counterpart (jzlxv3 cellsnp hybrid excepted: new) and no archive
# destination exists. DRY_RUN=1: checks and a summary only. Usage: qsub -q slice -v DATASETS="jzlxv3 jzlcel" human/scripts/swap_umi_calls_dataset.sh
set -euo pipefail
BASEDIR=${SCRNA_DIR:-.}
cd "$BASEDIR"
DATASETS=${DATASETS:-"jzlxv3 jzlcel"}

umi_paths() {   # relative umi_tools paths under results tree $1
  (cd "$1" && {
    find . -mindepth 4 -maxdepth 4 -type d -name '*UmiDup' -path './*/0[34]_*/*'
    find . -mindepth 3 -maxdepth 3 -type d -name '*UmiDup' -path './*/02_barcode_filtered/*'
    find . -mindepth 3 -maxdepth 3 -type f -path './*/02_dedup/*' \( -name '*_umiDup.bam' -o -name '*_umiDup.bam.bai' \)
    find . -mindepth 3 -maxdepth 3 -type d -path './*/02_dedup/subsample'
  } | sed 's|^\./||' | sort)
}

for ds in $DATASETS; do
  PUB=human/$ds/results_v8; NEW=human/$ds/results_v8_umi; ARC=human/${ds}_unseeded/results_v8
  [ -d "$PUB" ] && [ -d "$NEW" ] || { echo "ERROR: $PUB or $NEW missing" >&2; exit 1; }
  old=$(mktemp -p tmp); new=$(mktemp -p tmp)
  umi_paths "$PUB" > "$old"; umi_paths "$NEW" | grep -v '/02_dedup/subsample$' > "$new"
  missing=$(comm -13 "$old" "$new"); extra=$(comm -23 "$old" "$new" | grep -v '/02_dedup/subsample$' || true)
  echo "== $ds: unseeded $(wc -l < "$old"), seeded $(wc -l < "$new"), seeded-only $(printf '%s' "$missing" | grep -c . || true), unseeded-only $(printf '%s' "$extra" | grep -c . || true)"
  [ -z "$missing" ] || { echo "ERROR: seeded paths without an unseeded counterpart:" >&2; echo "$missing" >&2; exit 1; }
  [ -z "$extra" ] || { echo "ERROR: unseeded paths the seeded run did not make:" >&2; echo "$extra" >&2; exit 1; }
  while read -r r; do [ -e "$ARC/$r" ] && { echo "ERROR: archive exists: $ARC/$r" >&2; exit 1; }; done < "$old"

  # cellsnp hybrid: the seeded run makes it for every *UmiDup BAM type; record which ones the unseeded tree lacked
  (cd "$NEW" && find . -mindepth 5 -maxdepth 5 -type d -name hybrid -path './*/04_filtered/cellsnp/*UmiDup/*') | sed 's|^\./||' | sort |
    while read -r h; do [ -d "$PUB/$h" ] || echo "new (no unseeded counterpart): $PUB/$h"; done

  if [ "${DRY_RUN:-0}" = 1 ]; then cut -d/ -f2- "$old" | sed 's|/[^/]*$||' | sort | uniq -c; rm -f "$old" "$new"; continue; fi
  now=$(date '+%Y-%m-%d %H:%M %Z')
  while read -r r; do
    mkdir -p "$(dirname "$ARC/$r")"; mv "$PUB/$r" "$ARC/$r"
    printf "%s\t%s\t%s\t%s\n" "$PUB/$r" "$ARC/$r" "$now" "v10 non-GT swap (unseeded umi_tools): $ds" >> HANDOFF/MOVED.tsv
  done < "$old"
  echo "$ds: archived $(wc -l < "$old") paths to $ARC"
  while read -r r; do cp -a "$NEW/$r" "$PUB/$r"; done < "$new"
  echo "$ds: placed $(wc -l < "$new") paths from $NEW"
  n_after=$(umi_paths "$PUB" | wc -l)
  [ "$n_after" -eq "$(wc -l < "$new")" ] || { echo "ERROR: $ds has $n_after umi_tools paths after placing" >&2; exit 1; }
  rm -f "$old" "$new"
done
echo "done $(date)"
