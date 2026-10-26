#!/bin/bash
# replicates_arm_rank_v8.sh   (2026-10-15)
# Dataset-tier placement of replicate-aware JACUSA2 among the 58 ranked configurations, then its supplementary figure.
#   qsub -hold_jid <replicates summary + reachable> -N mm_armrank human/pm/scripts/rjob_v8.sh human/pm/scripts/replicates_arm_rank_v8.sh
set -euo pipefail
Rscript human/scripts/replicates_arm_rank_v8.R
Rscript human/scripts/replicates_ranksum_v8.R
Rscript publication/teratoma_gt/scripts/TERATOMA_v8_replicates_arm_rank_plot.R
REPL_SET=reach Rscript publication/teratoma_gt/scripts/TERATOMA_v8_replicates_arm_rank_plot.R
convert -density 150 -background white -flatten publication/teratoma_gt/figures/pm_gt/TERATOMA_v8_replicates_arm_rank.svg \
  review/response_figures_2026-10-15/TERATOMA_v8_replicates_arm_rank.png
