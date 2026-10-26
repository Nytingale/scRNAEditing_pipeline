#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_compare.R   (2026-10-13)
# Collates the v8 rescoring (run_gt_v8_rescore.sh; ranks from gt_rank_v7.R per state) into
# gt_eval_v8/compare_config_rank.tsv, in the format of the earlier copy from gt_eval_v7
# (state, source_variant). Run gt_rank_v7.R on gt_eval_v8/scoring/<state> and gt_eval_v8/mask_compare/mismatch first.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
V8 = file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8")
# GT_V8_STATES: comma-separated subset (interim collation before every state is scored)
SCN = Sys.getenv("SCN", "scoring"); SFX = sub("^scoring", "", SCN)   # SCN=scoring_pwcov: headline denominator; MC mask_compare_pwcov
MC = Sys.getenv("MC", "mask_compare")
STATES = strsplit(Sys.getenv("GT_V8_STATES", "none,mismatch,clustered,both"), ",", fixed = TRUE)[[1]]
cf = rbindlist(lapply(STATES, function(s)
  fread(file.path(V8, SCN, s, "config_rank.tsv"))[, `:=`(state = s, source_variant = "v7_nomm_pgp1")]))
stopifnot(cf[, .N, by = .(state, gate, t)][, all(N == N[1])])
fwrite(cf, file.path(V8, paste0("compare_config_rank", SFX, ".tsv")), sep = "\t")
cat("compare_config_rank.tsv:", nrow(cf), "rows\n")
mk = file.path(V8, MC, "mismatch", "config_rank.tsv")
if (file.exists(mk)) {
  fwrite(fread(mk)[, `:=`(state = "mismatch", variant = "v7_mask_nomm_pgp1")],
         file.path(V8, MC, "compare_config_rank.tsv"), sep = "\t")
  cat("mask_compare/compare_config_rank.tsv written\n")
}
