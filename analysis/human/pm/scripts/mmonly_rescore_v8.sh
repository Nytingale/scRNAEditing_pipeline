#!/bin/bash
# mmonly_rescore_v8.sh <rerun job ids, comma-separated>   (2026-10-15)
# Submits the swap (mmonly_swap_v8.sh) held on the mmonly_rerun_v8.sh arrays, then every analysis that reads the
# main call state (*_step5_nomm), in parallel where independent:
#   KO mismatch x 3 tiers, REDIportal mismatch, mask-compare mismatch (run_gt_v8_rescore.sh 4-6, 14, 17-19)
#   -> ranks of all states + compare (run_gt_v8_rank_compare.sh) -> gt_v8_state_caller_summary.R
#   split cells (run_gt_split_score.sh 1-6) -> gt_split_compare.R
#   strong-evidence truths (run_gt_v8_strong_truth.sh 1-6) -> gt_rank_v7.R x 2
#   sweep (run_gt_v8_sweep_rank.sh score 1-45, rank) -> gt_v8_sweep_rank_compare.R
#   germline (run_gt_v8_germline.sh score 1-3)
#   replicates: k-filter, then k2 / k3 / k4 / jacusa2rep scoring in parallel, summary, reachable check
#   gt_v8_ag_percent.R, gt_vartype_counts_v8.R (8 threads), gt_v8_dedup_inputs.R
# mmonly_swap_v8.sh exits 100 if any rerun row failed, which leaves every dependent job held.
set -euo pipefail
HOLD=$1
B=${SCRNA_DIR:-.}
S=$B/human/pm/scripts
V8=$B/human/pm/gt_eval_v8
RJ=$S/rjob_v8.sh
cd "$B"
q() { qsub -terse "$@" | cut -d. -f1; }

jsw=$(q -N mm_swap -l h_vmem=4G -hold_jid "$HOLD" "$S/mmonly_swap_v8.sh")
ko=$(q -N mm_ko -t 4-6 -hold_jid "$jsw" "$S/run_gt_v8_rescore.sh")
redi=$(q -N mm_redi -t 14 -hold_jid "$jsw" "$S/run_gt_v8_rescore.sh")
mask=$(q -N mm_mask -t 17-19 -hold_jid "$jsw" "$S/run_gt_v8_rescore.sh")
rk=$(q -N mm_rank -l h_vmem=16G -hold_jid "$ko,$mask" "$RJ" "$S/run_gt_v8_rank_compare.sh")
sc=$(q -N mm_statesum -l h_vmem=16G -hold_jid "$rk,$redi" "$RJ" human/scripts/gt_v8_state_caller_summary.R)
sp=$(q -N mm_split -t 1-6 -hold_jid "$jsw" "$S/run_gt_split_score.sh")
spc=$(q -N mm_splitcmp -l h_vmem=16G -hold_jid "$sp" "$RJ" human/scripts/gt_split_compare.R)
st=$(q -N mm_strong -t 1-6 -hold_jid "$jsw" "$S/run_gt_v8_strong_truth.sh")
str=$(q -N mm_strongrk -l h_vmem=16G -hold_jid "$st" "$RJ" "$S/rank_strong_truth_v8.sh")
sw=$(q -N mm_sweep -t 1-45 -hold_jid "$jsw" "$S/run_gt_v8_sweep_rank.sh" score)
swr=$(q -N mm_sweeprk -l h_vmem=8G -hold_jid "$sw" "$S/run_gt_v8_sweep_rank.sh" rank)
swc=$(q -N mm_sweepcmp -l h_vmem=16G -hold_jid "$swr" "$RJ" human/scripts/gt_v8_sweep_rank_compare.R)
gl=$(q -N mm_germ -t 1-3 -l h_vmem=48G -hold_jid "$jsw" -v STAGE=score "$S/run_gt_v8_germline.sh")
kf=$(q -N mm_kfilter -hold_jid "$jsw" "$RJ" human/scripts/replicates_kfilter_v8.R)
kd=$(q -N mm_kdone -l h_vmem=1G -hold_jid "$kf" -b y touch "$V8/replicates/results_k2_arm/.done")
rs=""
for s in k2 k3 k4 jacusa2rep; do
  j=$(q -N mm_rep_$s -hold_jid "$kd" -v STRATS=$s "$S/replicates_score_v8.sh"); rs=${rs:+$rs,}$j
done
rsum=$(q -N mm_repsum -l h_vmem=16G -hold_jid "$rs,$ko,$redi" "$RJ" human/scripts/replicates_summary_v8.R)
rch=$(q -N mm_reach -hold_jid "$rsum" "$S/replicates_reachable_v8.sh")
rch2=$(q -N mm_reach2 -v REACH_MIN_TER=2 -hold_jid "$rsum" "$S/replicates_reachable_v8.sh")
ark=$(q -N mm_armrank -l h_vmem=8G -hold_jid "$rch,$rch2" "$RJ" "$S/replicates_arm_rank_v8.sh")
ag=$(q -N mm_ag -l h_vmem=16G -hold_jid "$jsw" "$RJ" human/scripts/gt_v8_ag_percent.R)
vt=$(q -N mm_vartype -pe threaded 8 -l h_vmem=4G -hold_jid "$jsw" "$RJ" human/scripts/gt_vartype_counts_v8.R 8)
dd=$(q -N mm_dedup -hold_jid "$jsw" "$RJ" human/scripts/gt_v8_dedup_inputs.R)
cat <<EOF
swap $jsw
KO $ko | REDIportal $redi | mask $mask -> rank $rk -> state summary $sc
split $sp -> $spc | strong $st -> $str | sweep $sw -> $swr -> $swc | germline $gl
replicates kfilter $kf -> $rs -> summary $rsum -> reachable $rch, >=2 teratomas $rch2 -> dataset-tier rank $ark
ag% $ag | vartype $vt | dedup $dd
EOF
