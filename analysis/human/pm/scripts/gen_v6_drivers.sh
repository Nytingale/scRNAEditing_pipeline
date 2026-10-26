#!/bin/bash
# Generate the sample-tier and cell-type-tier v6 pipeline drivers from the
# merged-arm v6 driver: only the header note, job name, samplesheet, outdir,
# work-dir and run-dir change; all pipeline parameters are identical.
set -eo pipefail
S=${SCRNA_DIR:-.}/human/pm/scripts
SRC=$S/run_pipeline_pm_gt_merged_full_v6.sh

gen () {  # out tag job sheet note
  local out=$1 tag=$2 job=$3 sheet=$4 note=$5
  awk -v note="$note" -v out="$(basename "$out")" '
    NR == 3 { print "# " out "   (2026-09-28; copy of run_pipeline_pm_gt_merged_full_v6.sh)"; print note; print "# (merged-arm header retained below for provenance)"; next }
    { print }' "$SRC" \
  | sed -e "s|nf_run_gt_v6|nf_run_gt_${tag}_v6|g" \
        -e "s|results_gt_v6|results_gt_${tag}_v6|g" \
        -e "s|work_gt_v6|work_gt_${tag}_v6|g" \
        -e "s|samplesheet_pm_gt_merged_dedupv6.csv|${sheet}|g" \
        -e "s|^#\\\$ -N nf_pm_gt_v6|#\$ -N ${job}|" \
        -e "s|(MERGED ARMS, FULL GRID)|(${tag} tier v6, FULL GRID)|" > "$out"
  chmod +x "$out"
  echo "== $out"; grep -n "_${tag}_v6\|${sheet}\|^#\\\$ -N" "$out"
  if grep -q "results_gt_v6\|work_gt_v6\|merged_dedupv6" "$out"; then echo "STALE merged path in $out" >&2; exit 1; fi
}

gen $S/run_pipeline_pm_gt_persample_v6.sh persample nf_pm_ps_v6 samplesheet_pm_gt_persample_dedupv6.csv \
"# SAMPLE TIER (teratoma x arm, 20 units; each unit is one GEM well). sampleRmDup /
# sampleUmiDup are PREBUILT from the per-teratoma v6 dedup BAMs
# (02_gt_bams_{rmdup,umi}_v6/T*_*.bam), the same inputs the arm tier merged."

gen $S/run_pipeline_pm_gt_celltype_v6.sh celltype nf_pm_ct_v6 samplesheet_pm_gt_celltype_as_sample_dedupv6.csv \
"# CELL-TYPE TIER (arm x cell type, pooled over T1-T4). sampleRmDup / sampleUmiDup
# are PREBUILT: the per-teratoma v6 dedup BAMs are split by cell type and merged
# across teratomas (02_gt_bams_celltype_{rmdup,umi}_v6/), i.e. dedup literally per
# GEM well, then merge -- never dedup on prefixed composite barcodes."
