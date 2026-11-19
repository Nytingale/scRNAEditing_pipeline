#!/bin/bash
# TSO-clipped teratoma rerun, cell-type tier. Same validator settings and split as the v6 cell-type BAMs
# (validate_samplesheet.R --dataset pm_gt --min-cells 0; split per teratoma, CB/CR prefixed, merged across T1-T4):
#   raw   samplesheet_pm_gt_persample.csv  (02_gt_bams/T*)           -> _celltype_frags,           02_gt_bams_celltype
#   rmdup _celltype_frags_rmdup_v6 sheet   (02_gt_bams_rmdup_v6/T*)  -> _celltype_frags_rmdup_v6,  02_gt_bams_celltype_rmdup_v6
#   umi   _celltype_frags_umi_v6 sheet     (02_gt_bams_umi_v6/T*)    -> _celltype_frags_umi_v6,    02_gt_bams_celltype_umi_v6
#   (rmdup / umi sheets = samplesheet_pm_gt_persample.csv with the BAM folder replaced)
# Barcode tables go to <plan>/barcodes (gt_barcodes_celltype/ is a frozen input and is not rewritten).
# STAGE=launch (held on the three splits): checks every sentinel and BAM, then runs the v8 cell-type tier on the
# 23 NTC cell types (the scored subset, >= 20 truth sites, depends on the rebuilt truth) -> samplesheet_pm_gt_celltype_v8_NTC.csv.
#   bash human/pm/scripts/tso_rerun_celltype.sh            (STAGE=submit, from the login node)
#$ -S /bin/bash
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -j y
#$ -o ${SCRNA_DIR:-.}/human/pm/logs/
#$ -l h_vmem=4G
#$ -l h_rt=2:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
P=$B/human/pm
NC=$B/human/pm_noclip
SELF=$P/scripts/tso_rerun_celltype.sh
STAGE=${STAGE:-submit}
declare -A PLAN=([raw]=_celltype_frags [rmdup]=_celltype_frags_rmdup_v6 [umi]=_celltype_frags_umi_v6)
declare -A OUTB=([raw]=02_gt_bams_celltype [rmdup]=02_gt_bams_celltype_rmdup_v6 [umi]=02_gt_bams_celltype_umi_v6)
declare -A BAMD=([rmdup]=02_gt_bams_rmdup_v6 [umi]=02_gt_bams_umi_v6)
# PLANS=umi LAUNCH=0: the umi_tools split only, no cell-type pipeline launch (v10 rebuild 2026-11-12)
PLANS=${PLANS:-raw rmdup umi}
LAUNCH=${LAUNCH:-1}

if [ "$STAGE" = "submit" ]; then
  set +u; source ~/.bashrc >/dev/null 2>&1 || true; mamba activate scRNAEditing; set -u
  ids=""
  for k in $PLANS; do
    D=$P/${PLAN[$k]}; mkdir -p "$D"
    if [ "$k" = raw ]; then S=$P/samplesheet_pm_gt_persample.csv
    else S=$D/split_sheet_${k}_v6.csv; sed "s#/human/pm/02_gt_bams/#/human/pm/${BAMD[$k]}/#" $P/samplesheet_pm_gt_persample.csv > "$S"; fi
    Rscript $B/nextflow/bin/validate_samplesheet.R "$S" --dataset pm_gt --min-cells 0 \
      --tsv-out "$D/CELLTYPE_SCOPE_${k}.tsv" --emit-plan "$D" > "$D/validate.log" 2>&1
    n=$(wc -l < "$D/_expected.txt")
    rm -f "$D"/_ok_* "$D"/_frag_*.csv
    j=$(qsub -terse -q slice -t 1-$n -N "split_ct_$k" \
        -v PLANDIR="$D",OUTBAM="$P/${OUTB[$k]}",OUTBC="$D/barcodes",FRAGDIR="$D" $P/scripts/split_gt_bams_celltype.sh | cut -d. -f1)
    echo "$k: $n units, split job $j"; ids="${ids:+$ids,}$j"
  done
  [ "$LAUNCH" = 1 ] && qsub -terse -N tso_ct_launch -hold_jid "$ids" -v STAGE=launch "$SELF"
  exit 0
fi

if [ "$STAGE" = "launch" ]; then
  for k in raw rmdup umi; do
    D=$P/${PLAN[$k]}
    while read -r u; do [ -e "$D/_ok_$u" ] || { echo "ERROR: $k split unit $u has no sentinel" >&2; exit 1; }; done < "$D/_expected.txt"
  done
  SRC=$P/samplesheet_pm_gt_celltype_as_sample_dedupv6.csv
  OUT=$P/samplesheet_pm_gt_celltype_v8_NTC.csv
  awk -F, 'FNR == 1 || $1 ~ /^NTC_ctrl__/' "$SRC" > "$OUT"
  n=$(($(wc -l < "$OUT") - 1)); [ "$n" -eq 23 ] || { echo "ERROR: $n units in $OUT (expected 23)" >&2; exit 1; }
  miss=$(awk -F, 'NR > 1 { for (i = 1; i <= NF; i++) if ($i ~ /\.bam$/ && system("test -s " $i " && test -s " $i ".bai")) print $i }' "$OUT")
  [ -z "$miss" ] || { echo "ERROR: missing cell-type BAMs:" >&2; echo "$miss" >&2; exit 1; }
  for k in raw rmdup umi; do
    for u in $(cut -d, -f1 "$OUT" | tail -n +2); do
      cmp -s <(sort "$P/${PLAN[$k]}/barcodes/$u.tsv") <(sort "$P/gt_barcodes_celltype/$u.tsv") || echo "WARNING: $k barcodes differ for $u"
    done
  done
  qsub -terse -v TIER=celltype,SHEET_OVERRIDE="$OUT" $P/scripts/run_pipeline_pm_gt_v8.sh
  exit 0
fi
echo "ERROR: unknown STAGE $STAGE" >&2; exit 1
