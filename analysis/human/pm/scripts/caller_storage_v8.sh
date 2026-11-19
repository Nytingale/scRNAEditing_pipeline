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
# Array task 1-3 = arm, sample, celltype; 4-6 = the same tiers of the UMI-aware MarkDuplicates runs (2026-11-16).
# Output: human/pm/gt_eval_v8/caller_storage/storage_<tier>.tsv (tier, process, tag, hash, status, bytes)
#   bytes = NA when the work directory no longer exists.
# Array task 7 (2026-11-16): results_gt_v8/pipeline_trace.txt was overwritten on 2026-11-08 by a fully cached
#   rerun, so the arm-tier run time, CPU and peak memory are read back from each task's .command.trace
#   (tasks and hashes of storage_arm.tsv) -> caller_storage/trace_arm_from_work.tsv
#   (process, tag, hash, status, realtime "<ms>ms", %cpu, peak_rss "<KB> KB").
# Array task 8 (2026-11-16): the Cellsnp-lite hybrid-stranding second pass on the coordinate and umi_tools BAMs was
#   replayed outside Nextflow (cellsnp_hybrid_v8.sh fwd, SGE job 2041193); its elapsed and CPU time and memory are read
#   from the job logs (logs/csphyb_fwd_teratoma.o2041193.<id>), storage = the rewritten call file
#   -> caller_storage/cellsnp_hybrid_fwd.tsv (id, tier, rel_dir, elapsed_s, cpu_s, mem_gb, bytes).
# Array tasks 9-14 (2026-11-17): the Cellsnp-lite-only reruns with the fixed module (run_pipeline_pm_gt_v8.sh CSNP=1),
#   arm, sample, celltype, then their UMI-aware MarkDuplicates runs -> storage_<tier>_csnp.tsv.
# Usage : qsub -t 1-8 -N caller_storage human/pm/scripts/caller_storage_v8.sh   (from the repo root; 7 after 1)
#         Slurm: sbatch --array=9-14 --wrap 'SGE_TASK_ID=$SLURM_ARRAY_TASK_ID bash human/pm/scripts/caller_storage_v8.sh'
# ============================================================================
set -euo pipefail
PM=human/pm
TIERS=(arm sample celltype arm_umimd sample_umimd celltype_umimd)
RES=(results_gt_v8 results_gt_persample_v8 results_gt_celltype_v8
     results_gt_v8_umimd results_gt_persample_v8_umimd results_gt_celltype_v8_umimd)
WORK=(work_gt_v8 work_gt_persample_v8 work_gt_celltype_v8 work_gt_v8_umimd work_gt_persample_v8_umimd work_gt_celltype_v8_umimd)
for k in 0 1 2 3 4 5; do
  TIERS+=("${TIERS[$k]}_csnp"); RES+=("${RES[$k]}_csnp"); WORK+=("${WORK[$k]}_csnp")
done
i=$(( ${SGE_TASK_ID:-1} - 1 ))
if (( i == 7 )); then
  C=$PM/gt_eval_v8/cellsnp_hybrid
  out=$PM/gt_eval_v8/caller_storage/cellsnp_hybrid_fwd.tsv
  printf "id\ttier\trel_dir\telapsed_s\tcpu_s\tmem_gb\tbytes\n" > "$out.tmp"
  awk -F'\t' 'NR > 1 { print $1 "\t" $2 "\t" $5 }' "$C/manifest.tsv" |
  while IFS=$'\t' read -r id tier rel; do
    log=logs/csphyb_fwd_teratoma.o2041193.$id
    [[ "$(cut -f2 "$C/status/fwd/$id.tsv")" == 0 ]] || { echo "fwd $id not completed" >&2; exit 1; }
    bytes=$(find "$C/calls/$id" -maxdepth 1 -type f -printf '%s\n' | awk '{ s += $1 } END { printf "%.0f", s }')
    awk -v id="$id" -v t="$tier" -v r="$rel" -v b="$bytes" '$2 == "ELAPSED_TIME_S" { e = $3 } $2 == "CPU_USED_S" { c = $3 }
      $2 == "MEM_USED_G" { m = $3 } END { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", id, t, r, e, c, m, b }' "$log" >> "$out.tmp"
  done
  mv "$out.tmp" "$out"
  echo "done cellsnp hybrid fwd: $(($(wc -l < "$out") - 1)) tasks"
  exit 0
fi
if (( i == 6 )); then
  out=$PM/gt_eval_v8/caller_storage/trace_arm_from_work.tsv
  printf "process\ttag\thash\tstatus\trealtime\t%%cpu\tpeak_rss\n" > "$out.tmp"
  awk -F'\t' 'NR > 1 { print $2 "\t" $3 "\t" $4 "\t" $5 }' $PM/gt_eval_v8/caller_storage/storage_arm.tsv |
  while IFS=$'\t' read -r process tag hash status; do
    dir=$(compgen -G "$PM/work_gt_v8/$hash*" | head -1)
    awk -F= -v p="$process" -v t="$tag" -v h="$hash" -v s="$status" '{ v[$1] = $2 }
      END { printf "%s\t%s\t%s\t%s\t%sms\t%.1f%%\t%s KB\n", p, t, h, s, v["realtime"], v["%cpu"] / 10, v["peak_rss"] }' \
      "$dir/.command.trace" >> "$out.tmp"
  done
  mv "$out.tmp" "$out"
  echo "done arm trace: $(($(wc -l < "$out") - 1)) tasks"
  exit 0
fi
(( i >= 8 )) && i=$(( i - 2 ))
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
