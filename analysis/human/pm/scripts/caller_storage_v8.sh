#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/
#$ -l h_vmem=4G
# ============================================================================
# caller_storage_v8.sh   (2026-10-16; review R1 minor 9: storage)
# Disk space of each completed variant-calling task of the v8 teratoma runs: bytes of the regular files in its
# Nextflow work directory, excluding Nextflow's own .command.* / .exitcode files. Inputs are symlinks and
# are not counted (find -type f does not follow them).
# Array task 1-3 = arm, sample, celltype.
# Output: human/pm/gt_eval_v8/caller_storage/storage_<tier>.tsv (tier, process, tag, hash, status, bytes)
#   bytes = NA when the work directory no longer exists.
# Usage : qsub -t 1-3 -N caller_storage human/pm/scripts/caller_storage_v8.sh   (from the repo root)
# ============================================================================
set -euo pipefail
PM=human/pm
TIERS=(arm sample celltype)
RES=(results_gt_v8 results_gt_persample_v8 results_gt_celltype_v8)
WORK=(work_gt_v8 work_gt_persample_v8 work_gt_celltype_v8)
i=$(( ${SGE_TASK_ID:-1} - 1 ))
tier=${TIERS[$i]}
trace=$PM/${RES[$i]}/pipeline_trace.txt
work=$PM/${WORK[$i]}
out_dir=$PM/gt_eval_v8/caller_storage
mkdir -p "$out_dir"
out=$out_dir/storage_${tier}.tsv
tmp=$out.tmp

printf "tier\tprocess\ttag\thash\tstatus\tbytes\n" > "$tmp"
awk -F'\t' 'NR == 1 { for (k = 1; k <= NF; k++) h[$k] = k; next }
  $h["process"] ~ /^(BCFTOOLS_MPILEUP|CELLSNP_LITE|GATK_HAPLOTYPECALLER|JACUSA2_CALL1|REDITOOLS3|SCALLELE_CHROM|SCALLELE_GATHER)$/ &&
  $h["status"] == "COMPLETED" { print $h["process"] "\t" $h["tag"] "\t" $h["hash"] "\t" $h["status"] }' "$trace" |
while IFS=$'\t' read -r process tag hash status; do
  dir=$(compgen -G "$work/$hash*" | head -1 || true)
  if [[ -n "$dir" && -d "$dir" ]]; then
    bytes=$(find "$dir" -maxdepth 1 -type f ! -name '.command.*' ! -name '.exitcode' -printf '%s\n' | awk '{ s += $1 } END { printf "%.0f", s }')
  else
    bytes=NA
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tier" "$process" "$tag" "$hash" "$status" "$bytes" >> "$tmp"
done
mv "$tmp" "$out"
echo "done $tier: $(($(wc -l < "$out") - 1)) tasks, $(awk -F'\t' 'NR > 1 && $6 == "NA"' "$out" | wc -l) missing work dirs"
