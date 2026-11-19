#!/bin/bash
# tso_rerun_score.sh  [submit|keep|rank]   (2026-10-30)
# TSO-clipped teratoma: scoring inputs and the headline v8 scoring (scoring_pwcov / rediportal_pwcov /
# mask_compare_pwcov), with the unchanged scripts. Data-independent inputs are linked from human/pm_noclip.
#   submit (login node, after the collect jobs given as $1, comma list)
#     v6 scan     run_gt_truth_rebuild_v6.sh scan (90) -> cands -> truthpile merged_NTC_ctrl BQ30
#                 (gt_eval_v6/pileups/bq30/pileup_merged_NTC_ctrl.tsv, candidates/candidates_annot.tsv: REDIportal universe)
#     REDIportal  rediportal_pileup_v6.sh setup + 28 units x 24 chromosomes -> gt_eval_v6/rediportal/pileup/
#     unit_cov    gt_v8_unit_cov.sh setup + 28 units + collect -> gt_eval_v8/unit_cov/ko_truth_unit_cov.tsv.gz
#     all 23      run_gt_v8_rescore.sh tasks 6 (KO mismatch celltype) and 14 (REDIportal mismatch), GT_UNIT_KEEP empty
#                 -> gt_eval_v8/_ct_all23/
#     keep        gt_v8_celltype_covered_positives.R -> gt_v8_celltype_keep_pw.R gate_cov 30 -> unit_cov/celltype_keep_pwcov.txt
#     headline    run_gt_v8_rescore.sh KO 1-6, 10-12 (none, mismatch, both), REDIportal 13-14, 16, mask 17-19,
#                 GT_UNIT_RULE=gate_cov; the clustered state (7-9, 15) needs the nomm_rerun_v8.sh replay first
#     rank        gt_rank_v7.R per state + mask_compare_pwcov/mismatch, gt_v8_compare.R (SCN=scoring_pwcov)
#   inputs <truth job> (login node; v10 rebuild 2026-11-12): the v6 scan, REDIportal, unit_cov, all 23 and keep steps
#                 of submit without the headline rescore, queued after the v8 truth copy (run_gt_truth_tso.sh copy),
#                 which gt_v8_unit_cov.sh setup reads. The scored call sets do not change these inputs (per-unit
#                 coverage and truth counts only), so they need no caller rerun.
#$ -S /bin/bash
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -j y
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=4:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
NC=$B/human/pm_noclip
S=$PM/scripts
V6=$PM/gt_eval_v6
V8=$PM/gt_eval_v8
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
SELF=$S/tso_rerun_score.sh
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"   # bedtools 2.31.0
cd "$B"
MODE=${1:-}; shift || true

scoring_inputs() {  # $1 = hold for the all-23 scoring; sets jk (keep job)
  local HOLD=$1
  mkdir -p "$V6/rediportal" "$V8/unit_cov" "$V8/_ct_all23"
  if [ -e "$NC/gt_eval_v6/rediportal/rediportal_dbsnp_common_sites.txt" ]; then
    ln -sfn "$NC/gt_eval_v6/rediportal/rediportal_dbsnp_common_sites.txt" "$V6/rediportal/rediportal_dbsnp_common_sites.txt"
  fi
  js=$(qsub -terse -q slice -N tsos_v6scan -t 1-"$(wc -l < "$V6/scan/regions.tsv")" $S/run_gt_truth_rebuild_v6.sh scan | cut -d. -f1)
  jc=$(qsub -terse -q slice -N tsos_v6cands -l h_vmem=32G -hold_jid "$js" $S/run_gt_truth_rebuild_v6.sh cands)
  k=$(ls "$PM"/02_gt_bams_umi_v6/*.bam | sed 's#.*/##' | sort | grep -n '^merged_NTC_ctrl\.bam$' | cut -d: -f1)
  jb=$(qsub -terse -q slice -N tsos_v6pile -l h_vmem=32G -t "$k" -hold_jid "$jc" $S/run_gt_truth_rebuild_v6.sh truthpile | cut -d. -f1)
  bash $S/rediportal_pileup_v6.sh setup
  # REDIportal positions in dbSNP156 common, as after_rediportal_pileup_v6.sh (when not linked above)
  [ -s "$V6/rediportal/rediportal_dbsnp_common_sites.txt" ] || cat "$V6"/rediportal/bed/chr*.bed | sort -k1,1 -k2,2n -S 8G \
    | bedtools intersect -u -a - -b "$V6/snp/dbSNP156_common.bed" \
    | awk -v OFS='\t' '{print $1 ":" $3}' > "$V6/rediportal/rediportal_dbsnp_common_sites.txt"
  jr=$(qsub -terse -N tsos_redi -t 1-"$(wc -l < "$V6/rediportal/tasks.tsv")" -tc 150 $S/rediportal_pileup_v6.sh | cut -d. -f1)
  bash $S/gt_v8_unit_cov.sh setup
  ju=$(qsub -terse -N tsos_ucov -t 1-"$(wc -l < "$V8/unit_cov/tasks.tsv")" $S/gt_v8_unit_cov.sh | cut -d. -f1)
  juc=$(qsub -terse -q slice -N tsos_ucov_collect -l h_vmem=8G -hold_jid "$ju" -b n $S/gt_v8_unit_cov.sh collect)
  H="$HOLD,$jb,$jr,$juc"
  A23="GT_UNIT_KEEP=,SC=$V8/_ct_all23/scoring,RD=$V8/_ct_all23/rediportal"
  ja=$(qsub -terse -N tsos_all23 -t 6-6 -hold_jid "$H" -v "$A23" $S/run_gt_v8_rescore.sh | cut -d. -f1)
  jar=$(qsub -terse -N tsos_all23r -t 14-14 -hold_jid "$H" -v "$A23" $S/run_gt_v8_rescore.sh | cut -d. -f1)
  jk=$(qsub -terse -N tsos_keep -hold_jid "$ja,$jar" "$SELF" keep)
  echo "v6 scan $js | cands $jc | truthpile $jb (task $k) | redi $jr | ucov $ju | ucov collect $juc | all23 $ja,$jar | keep $jk"
}

case "$MODE" in
inputs)
  qsub -terse -N tsos_inputs -l h_vmem=8G -hold_jid "${1:?truth copy job}" "$SELF" inputs_run
  ;;
inputs_run)
  export PATH="/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
  scoring_inputs "$JOB_ID"
  ;;
submit)
  scoring_inputs "${1:?collect job ids}"
  HL="GT_UNIT_RULE=gate_cov,SC=$V8/scoring_pwcov,RD=$V8/rediportal_pwcov,MC=$V8/mask_compare_pwcov,GT_UNIT_KEEP=$V8/unit_cov/celltype_keep_pwcov.txt"
  ids=""
  for t in 1-6 10-12 13-14 16-16 17-19; do
    ids="$ids,$(qsub -terse -N tsos_rescore -t $t -hold_jid "$jk" -v "$HL" $S/run_gt_v8_rescore.sh | cut -d. -f1)"
  done
  jrk=$(qsub -terse -N tsos_rank -hold_jid "${ids#,}" "$SELF" rank)
  echo "rescore ${ids#,} | rank $jrk"
  ;;
keep)
  $RS $B/human/scripts/gt_v8_celltype_covered_positives.R "$V8/_ct_all23/scoring/mismatch" "$V8/_ct_all23/rediportal/scoring_mismatch"
  $RS $B/human/scripts/gt_v8_celltype_keep_pw.R gate_cov 30
  echo "kept: $(wc -l < "$V8/unit_cov/celltype_keep_pwcov.txt")"; cat "$V8/unit_cov/celltype_keep_pwcov.txt"
  ;;
rank)
  # qsub -v splits on commas, so a submitted list comes colon-separated
  ST=${RANK_STATES:-none,mismatch,both}; ST=${ST//:/,}
  for st in ${ST//,/ }; do $RS $B/human/scripts/gt_rank_v7.R "$V8/scoring_pwcov/$st"; done
  $RS $B/human/scripts/gt_rank_v7.R "$V8/mask_compare_pwcov/mismatch"
  SCN=scoring_pwcov MC=mask_compare_pwcov GT_V8_STATES=$ST $RS $B/human/scripts/gt_v8_compare.R
  ;;
clustered)
  # held on the nomm_rerun_v8.sh array ($1): clustered-state scoring, then ranks over all four states
  HL="GT_UNIT_RULE=gate_cov,SC=$V8/scoring_pwcov,RD=$V8/rediportal_pwcov,MC=$V8/mask_compare_pwcov,GT_UNIT_KEEP=$V8/unit_cov/celltype_keep_pwcov.txt"
  H="${1:?nomm array job},$(qstat -u "$USER" | awk '$3 == "tsos_rank" { print $1 }' | paste -sd,)"
  a=$(qsub -terse -N tsos_rescore_cl -t 7-9 -hold_jid "${H%,}" -v "$HL" $S/run_gt_v8_rescore.sh | cut -d. -f1)
  b=$(qsub -terse -N tsos_rescore_cl -t 15-15 -hold_jid "${H%,}" -v "$HL" $S/run_gt_v8_rescore.sh | cut -d. -f1)
  r=$(qsub -terse -N tsos_rank4 -hold_jid "$a,$b" -v RANK_STATES=none:mismatch:clustered:both "$SELF" rank)
  echo "clustered rescore $a,$b | rank (4 states) $r"
  ;;
windows)
  # window-filter comparisons of the response (Response Figures 8-9, 17) need the clustered pass on:
  # held on nomm_rerun_v8.sh REPLAY=clustered ($1) and REPLAY=both ($2); rescore both states, rank all four
  HL="GT_UNIT_RULE=gate_cov,SC=$V8/scoring_pwcov,RD=$V8/rediportal_pwcov,MC=$V8/mask_compare_pwcov,GT_UNIT_KEEP=$V8/unit_cov/celltype_keep_pwcov.txt"
  a=$(qsub -terse -N tsos_rescore_cl -t 7-9 -hold_jid "${1:?clustered replay job}" -v "$HL" $S/run_gt_v8_rescore.sh | cut -d. -f1)
  b=$(qsub -terse -N tsos_rescore_cl -t 15-15 -hold_jid "$1" -v "$HL" $S/run_gt_v8_rescore.sh | cut -d. -f1)
  c=$(qsub -terse -N tsos_rescore_both -t 10-12 -hold_jid "${2:?both replay job}" -v "$HL" $S/run_gt_v8_rescore.sh | cut -d. -f1)
  d=$(qsub -terse -N tsos_rescore_both -t 16-16 -hold_jid "$2" -v "$HL" $S/run_gt_v8_rescore.sh | cut -d. -f1)
  r=$(qsub -terse -N tsos_rank4 -hold_jid "$a,$b,$c,$d" -v RANK_STATES=none:mismatch:clustered:both "$SELF" rank)
  echo "rescore clustered $a,$b | both $c,$d | rank (4 states) $r"
  ;;
*) echo "ERROR: usage tso_rerun_score.sh submit <collect jobs> | keep | rank" >&2; exit 1 ;;
esac
