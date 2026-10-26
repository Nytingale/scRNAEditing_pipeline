#!/bin/bash
# rank_strong_truth_v8.sh   (2026-10-15): gt_rank_v7.R on both strong-NTC-evidence truths (after run_gt_v8_strong_truth.sh)
set -euo pipefail
V8=${SCRNA_DIR:-.}/human/pm/gt_eval_v8
for a in 10 20; do Rscript human/scripts/gt_rank_v7.R "$V8/sensitivity/strong_alt$a/scoring/mismatch"; done
