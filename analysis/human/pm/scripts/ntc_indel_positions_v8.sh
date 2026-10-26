#!/bin/bash
# ntc_indel_positions_v8.sh   (2026-10-23; Reviewer 2, major 3: false positives near microindels)
# Indel positions in the reads of the merged NTC arm (UMI-deduplicated BAM, as the ground-truth pileups), one task per
# chromosome. Reads: primary, mapped, not QC-fail or duplicate (-F 0xF04), MAPQ >= 20. From each CIGAR, every I and D
# op is recorded at its reference position: D = first deleted reference base; I = reference base after the insertion.
# Output: human/pm/gt_eval_v8/indel_sj/indels/<chr>.tsv.gz  (chr, pos, type I/D, len, n_reads)
#   qsub -t 1-23 human/pm/scripts/ntc_indel_positions_v8.sh
#$ -S /bin/bash
#$ -N ntc_indels
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 1
#$ -l h_vmem=3G
#$ -l h_rt=4:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
BAM=$B/human/pm/02_gt_bams_umi/merged_NTC_ctrl.bam
OUT=$B/human/pm/gt_eval_v8/indel_sj/indels
mkdir -p "$OUT"
CHRS=(chr{1..22} chrX)
chr=${CHRS[$(( ${SGE_TASK_ID:-1} - 1 ))]}
ST=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/samtools
$ST view -F 0xF04 -q 20 "$BAM" "$chr" | awk -F'\t' -v OFS='\t' '
  $6 ~ /[ID]/ {
    p = $4; c = $6
    while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
      n = substr(c, 1, RLENGTH - 1) + 0; op = substr(c, RLENGTH, 1); c = substr(c, RLENGTH + 1)
      if (op == "D") { k[$3 "\t" p "\tD\t" n]++ }
      else if (op == "I") { k[$3 "\t" p "\tI\t" n]++ }
      if (op ~ /[MDN=X]/) p += n
    }
  }
  END { for (x in k) print x, k[x] }' | sort -k2,2n | gzip > "$OUT/$chr.tsv.gz.tmp"
mv "$OUT/$chr.tsv.gz.tmp" "$OUT/$chr.tsv.gz"
echo "$chr done: $(zcat "$OUT/$chr.tsv.gz" | wc -l) indel position x length rows"
