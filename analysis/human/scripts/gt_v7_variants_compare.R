#!/usr/bin/env Rscript
# =============================================================================
# gt_v7_variants_compare.R   (2026-10-10)
# Collates the v7 truth-variant scoring (run_gt_v7_variants_score.sh / gt_rank_v7.R) across the
# four call-set window states, and checks that REDItools3 calls are in the RNA frame in every
# state (caller strand mode: share of calls at positives that are RNA AG on the A-to-I strand,
# next to JACUSA2; the pre-fix parse put minus-strand A-to-I at RNA TC, ~50%).
# Input : gt_eval_v7/variants/scoring/<state>/<variant>/{caller_rank,config_rank}.tsv,
#         .../<state>/<variant>/arm/gt_v7_calls_NTC.tsv.gz, gt_eval_v7/variants/truth/truth_variant_counts.tsv
# Output: gt_eval_v7/variants/
#   compare_caller_rank.tsv   state x variant x gate x t x caller (best configuration)
#   compare_config_rank.tsv   every reported configuration
#   reditools3_frame_check.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
V = file.path(B, "human/pm/gt_eval_v7/variants")
STATES = c("none", "mismatch", "clustered", "both")
VARS = c("v7", "v7_mask", "v7_nomm", "v7_mask_nomm", "v7_pgp1", "v7_mask_pgp1", "v7_nomm_pgp1", "v7_mask_nomm_pgp1")
rd = function(f) rbindlist(lapply(STATES, function(s) rbindlist(lapply(VARS, function(v)
  fread(file.path(V, "scoring", s, v, f))[, `:=`(state = s, variant = v)]))))
cr = rd("caller_rank.tsv"); cf = rd("config_rank.tsv")
cnt = fread(file.path(V, "truth", "truth_variant_counts.tsv"))
cr = merge(cr, cnt[, .(variant, gate, t, n_pos = positive, n_neg = informative_negative)], by = c("variant", "gate", "t"))
fwrite(cr, file.path(V, "compare_caller_rank.tsv"), sep = "\t")
fwrite(cf, file.path(V, "compare_config_rank.tsv"), sep = "\t")
print(dcast(cr[gate == "expected" & t == "t50" & variant %in% VARS],
            caller ~ state + variant, value.var = "caller_rank"))

fc = rbindlist(lapply(STATES, function(s) {
  x = fread(file.path(V, "scoring", s, "v7", "arm", "gt_v7_calls_NTC.tsv.gz"))
  tr = fread(file.path(V, "truth", "v7", "gt_v7_truth_sites.tsv.gz"), select = c("gate", "site", "class_t50"))
  x = merge(x, tr, by = c("gate", "site"))[class_t50 == "positive" & gate == "expected"]
  x[, c("caller", "bam_type", "strand_mode", "snp") := tstrsplit(config, "__", fixed = TRUE)]
  x[caller %in% c("reditools3", "jacusa2") & strand_mode == "caller",
    .(state = s, n_calls = .N, rna_ag = mean(varType == "AG"), rna_tc = mean(varType == "TC"),
      ag_on_minus = sum(varType == "AG" & strand == "-"), tc_on_minus = sum(varType == "TC" & strand == "-")),
    by = .(caller, bam_type)]
}))
fwrite(fc, file.path(V, "reditools3_frame_check.tsv"), sep = "\t")
print(fc)
r3 = fc[caller == "reditools3"]
stopifnot(nrow(r3) > 0, all(r3$rna_ag > 0.8), all(r3$ag_on_minus > 5 * r3$tc_on_minus))
message("REDItools3 RNA-frame check passed in all states")
