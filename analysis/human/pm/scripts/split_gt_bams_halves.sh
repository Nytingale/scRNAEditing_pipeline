#!/bin/bash
# split_gt_bams_halves.sh   (2026-10-13)
# Split-by-cells test: cut the NTC BAMs into cell halves A / B (gt_split_cells.R barcodes) by CB tag.
#   half B (called): unprocessed, coordinate dedup, UMI dedup BAMs of T1-T4 NTC and the 11 cell-type-tier units
#   half A (truth) : UMI dedup BAMs of T1-T4 NTC
# Per-teratoma halves are then merged into merged_NTC_ctrl_<H> (dedup per GEM well, then merge, as
# merge_gt_umi_bams_v6.sh); the cell-type BAMs carry teratoma-prefixed CB tags (gt_barcodes_celltype format).
# Output: gt_eval_v8/split_cells/bams/{raw,rmdup,umi}/<unit>_<H>.bam(.bai)
#   bash split_gt_bams_halves.sh setup                 # task list
#   qsub -t 1-$(wc -l < tasks.tsv) split_gt_bams_halves.sh
#   qsub -hold_jid <that job> -v STAGE=merge split_gt_bams_halves.sh
#$ -S /bin/bash
#$ -N split_halves
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 4
#$ -l h_vmem=4G
#$ -l h_rt=24:00:00
set -euo pipefail
PM=${SCRNA_DIR:-.}/human/pm
S=$PM/gt_eval_v8/split_cells
BC=$S/barcodes
OUT=$S/bams
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"
declare -A SRC_T=([raw]=02_gt_bams [rmdup]=02_gt_bams_rmdup_v6 [umi]=02_gt_bams_umi_v6)
declare -A SRC_C=([raw]=02_gt_bams_celltype [rmdup]=02_gt_bams_celltype_rmdup_v6 [umi]=02_gt_bams_celltype_umi_v6)
STAGE=${STAGE:-${1:-run}}

if [ "$STAGE" = "setup" ]; then
  mkdir -p "$OUT"/{raw,rmdup,umi} "$BC/cb"
  : > "$S/tasks.tsv"
  for n in 1 2 3 4; do
    for t in raw rmdup umi; do printf 'T%s_NTC_ctrl\tB\t%s\t%s\n' $n $t "$PM/${SRC_T[$t]}/T${n}_NTC_ctrl.bam" >> "$S/tasks.tsv"; done
    printf 'T%s_NTC_ctrl\tA\tumi\t%s\n' $n "$PM/${SRC_T[umi]}/T${n}_NTC_ctrl.bam" >> "$S/tasks.tsv"
  done
  while read -r u; do
    for t in raw rmdup umi; do printf '%s\tB\t%s\t%s\n' "$u" $t "$PM/${SRC_C[$t]}/$u.bam" >> "$S/tasks.tsv"; done
  done < "$PM/gt_eval_v8/unit_cov/celltype_keep_min20.txt"
  for f in "$BC"/*_[AB].tsv; do cut -f1 "$f" > "$BC/cb/$(basename "$f" .tsv).txt"; done
  awk -F'\t' '{ if (system("test -s " $4)) print "MISSING " $4 }' "$S/tasks.tsv"
  echo "tasks: $(wc -l < "$S/tasks.tsv")"
  exit 0
fi

if [ "$STAGE" = "merge" ]; then
  for h in A B; do for t in raw rmdup umi; do
    ins=("$OUT/$t"/T[1-4]_NTC_ctrl_$h.bam)
    [ -s "${ins[0]}" ] || continue
    [ ${#ins[@]} -eq 4 ] || { echo "ERROR: $t $h has ${#ins[@]} teratoma BAMs" >&2; exit 1; }
    o="$OUT/$t/merged_NTC_ctrl_$h.bam"
    samtools merge -f -@ 4 -o "$o.tmp.bam" "${ins[@]}" && mv "$o.tmp.bam" "$o" && samtools index -@ 4 "$o"
    sum=0; for b in "${ins[@]}"; do sum=$(( sum + $(samtools idxstats "$b" | awk '{s += $3 + $4} END {print s}') )); done
    m=$(samtools idxstats "$o" | awk '{s += $3 + $4} END {print s}')
    [ "$m" -eq "$sum" ] || { echo "ERROR: $o has $m reads, teratomas sum to $sum" >&2; exit 1; }
    echo "merged $t $h: $m reads"
  done; done
  exit 0
fi

read -r U H T IN < <(sed -n "${SGE_TASK_ID}p" "$S/tasks.tsv")
o="$OUT/$T/${U}_$H.bam"
samtools view -@ 4 -b -D "CB:$BC/cb/${U}_$H.txt" -o "$o.tmp.bam" "$IN"
mv "$o.tmp.bam" "$o"
samtools index -@ 4 "$o"
n_in=$(samtools idxstats "$IN" | awk '{s += $3 + $4} END {print s}')
n_out=$(samtools idxstats "$o" | awk '{s += $3 + $4} END {print s}')
echo "$U $H $T: $n_out of $n_in mapped+unmapped reads"
