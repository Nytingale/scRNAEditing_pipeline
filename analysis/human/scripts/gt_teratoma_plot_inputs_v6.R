#!/usr/bin/env Rscript
# =============================================================================
# gt_teratoma_plot_inputs_v6.R   (2026-09-29)
#
# Writes the v6 sweep scoring (GROUND_TRUTH_DEFINITION.md: positives at
# t = 25/50/75%, fixed informative negatives) in the table layout that the
# publication/teratoma_gt plot scripts read (the v5 gt_tier_scoring schema), so
# the result figures can be regenerated against the current truth.
#
#   tier names    arm -> "dataset" (the plot scripts' name for the arm tier), sample, celltype
#   tier value    mean over units (gt_score_v6.R convention); n_* are unit means
#   error         sample: SE across teratomas; celltype: SD across cell types;
#                 arm: site-bootstrap SD over resolved sites (arm_bootstrap_se.tsv)
#   ag_percent    100 x A>G sites / all sites (unique site x varType) in the SAME
#                 *_pre100.tsv.gz call files that are scored; median over units
#   reported      configurations with bam_type != "sample" only (as scored)
#
# Inputs : gt_eval_v6/scoring/sweep_infneg/{unit_scores_real.tsv.gz, arm_bootstrap_se.tsv, class_counts.tsv}
#          human/pm/results_gt{,_persample,_celltype}_v6/<NTC unit>/04_filtered/**/*_pre100.tsv.gz
# Outputs: gt_eval_v6/scoring/sweep_infneg/teratoma_inputs/
#          gt_tier_scoring.tsv, gt_tier_scoring_with_se.tsv, gt_tier_replicate_scores.tsv   (t = 0.50)
#          invariance/d{0.25,0.50,0.75}/{gt_tier_scoring.tsv, gt_tier_scoring_with_se.tsv, gt_tier_replicate_scores.tsv}
#          invariance/gt_tier_replicate_scores_all_depletions.tsv
#          ag_percent_units.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
# Window-filter call sets (HANDOFF 2026-09-29 UPDATE-22) reuse this script through three overrides:
#   GT_SCORING_DIR  scoring root (default gt_eval_v6/scoring)
#   GT_RES_TEMPLATE results dir per tier relative to human/pm, %s = arm|sample|celltype
#   GT_CALL_RE      call-file regex (default _pre100\.tsv\.gz$)
#   GT_SWEEP_SUBDIR sweep folder under the scoring root (default sweep_infneg; the depth-gated
#                   truth uses sweep_infneg_depth)
SW = file.path(Sys.getenv("GT_SCORING_DIR", file.path(PM, "gt_eval_v6", "scoring")),
               Sys.getenv("GT_SWEEP_SUBDIR", "sweep_infneg"))
OUT = file.path(SW, "teratoma_inputs")
COMBO = c("caller", "bam_type", "strand_mode", "snp_subtract")
RES = c(arm = "results_gt_v6", sample = "results_gt_persample_v6", celltype = "results_gt_celltype_v6")
if (nzchar(Sys.getenv("GT_RES_TEMPLATE"))) RES[] = sprintf(Sys.getenv("GT_RES_TEMPLATE"), names(RES))
CALL_RE = Sys.getenv("GT_CALL_RE", "_pre100\\.tsv\\.gz$")
TIER_OUT = c(arm = "dataset", sample = "sample", celltype = "celltype")
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

u = fread(file.path(SW, "unit_scores_real.tsv.gz"))
bt = fread(file.path(SW, "arm_bootstrap_se.tsv"))
cc = fread(file.path(SW, "class_counts.tsv"))
u[, (COMBO) := tstrsplit(config, "__", fixed = TRUE)]
u[, dep := as.numeric(sub("^t", "", t)) / 100]

# ---- A-to-G percentage of each scored call file ---------------------------------------
ag = rbindlist(lapply(names(RES), function(tr) {
  uu = unique(u[tier == tr, .(unit, config, caller, bam_type, strand_mode, snp_subtract)])
  rbindlist(lapply(seq_len(nrow(uu)), function(i) {
    x = uu[i]
    dir = file.path(PM, RES[[tr]], x$unit, "04_filtered", x$caller, x$bam_type, x$strand_mode, paste0("snp_", x$snp_subtract))
    f = list.files(dir, CALL_RE, full.names = TRUE, recursive = TRUE)
    stopifnot(length(f) == 1)
    if (file.size(f) == 0) return(data.table(tier = tr, unit = x$unit, config = x$config, n_sites = 0L, n_ag = 0L))
    d = unique(fread(f, select = c("site", "varType"), showProgress = FALSE))
    data.table(tier = tr, unit = x$unit, config = x$config, n_sites = nrow(d), n_ag = sum(d$varType == "AG"))
  }))
}))
ag[, ag_percent := fifelse(n_sites > 0, 100 * n_ag / n_sites, NA_real_)]
dir.create(file.path(OUT, "invariance"), recursive = TRUE, showWarnings = FALSE)
fwrite(ag, file.path(OUT, "ag_percent_units.tsv"), sep = "\t")
msg("A-to-G percentage: ", nrow(ag), " call files; median by tier: ",
    paste(ag[, .(m = round(median(ag_percent, na.rm = TRUE), 1)), by = tier][, paste(tier, m)], collapse = ", "))
u = merge(u, ag[, .(tier, unit, config, ag_percent)], by = c("tier", "unit", "config"), all.x = TRUE)

# ---- tables per depletion ------------------------------------------------------------
se = function(x) { x = x[is.finite(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }
sdv = function(x) { x = x[is.finite(x)]; if (length(x) > 1) sd(x) else NA_real_ }
write_d = function(tn) {
  x = u[t == tn]; d0 = x$dep[1]
  n_pos = cc[abs(t - d0) < 1e-9, positive]
  x[, tier_out := TIER_OUT[tier]]
  x[, `:=`(teratoma = fifelse(tier == "sample", sub("_NTC_ctrl$", "", unit), ""),
           celltype = fifelse(tier == "celltype", sub("^NTC_ctrl__", "", unit), ""))]
  s = x[, .(n_discovered = mean(n_calls), n_true = mean(tp), n_fp = mean(fp),
            precision = mean(precision, na.rm = TRUE), recall_vs_union = mean(recall), f1 = mean(f1, na.rm = TRUE),
            ag_percent = median(ag_percent, na.rm = TRUE), n_reps = .N,
            precision_sd_units = if (.BY$tier == "celltype") sdv(precision) else se(precision),
            recall_sd_units = if (.BY$tier == "celltype") sdv(recall) else se(recall),
            f1_sd_units = if (.BY$tier == "celltype") sdv(f1) else se(f1)),
        by = c("tier", "tier_out", COMBO)]
  s[, n_true_total := n_pos]
  b = bt[t == tn][, (COMBO) := tstrsplit(config, "__", fixed = TRUE)]
  s = merge(s, b[, c(COMBO, "precision_err", "recall_err", "f1_err"), with = FALSE], by = COMBO, all.x = TRUE)
  s[, `:=`(precision_se = fifelse(tier == "arm", precision_err, precision_sd_units),
           recall_se = fifelse(tier == "arm", recall_err, recall_sd_units),
           f1_se = fifelse(tier == "arm", f1_err, f1_sd_units),
           se_method = fcase(tier == "arm", "site_bootstrap_resolved", tier == "sample", "teratoma_replicate_se",
                             tier == "celltype", "celltype_sd"))]
  s[, `:=`(precision_mean = precision, recall_mean = recall_vs_union, f1_mean = f1)]
  s[, tier := tier_out]
  base = c("tier", COMBO, "n_discovered", "n_true", "n_fp", "precision", "recall_vs_union", "f1",
           "n_true_total", "ag_percent")
  setorderv(s, c("tier", "f1"), c(1, -1))
  rep = x[tier != "arm", .(caller, bam_type, strand_mode, snp_subtract, n_discovered = n_calls, n_true = tp, n_fp = fp,
                           precision, recall_vs_union = recall, f1, tier = tier_out, teratoma, celltype, depletion = d0)]
  dd = file.path(OUT, "invariance", sprintf("d%.2f", d0)); dir.create(dd, showWarnings = FALSE)
  fwrite(s[, ..base], file.path(dd, "gt_tier_scoring.tsv"), sep = "\t")
  fwrite(s[, c(base, "n_reps", "precision_mean", "precision_se", "recall_mean", "recall_se",
               "f1_mean", "f1_se", "se_method"), with = FALSE], file.path(dd, "gt_tier_scoring_with_se.tsv"), sep = "\t")
  fwrite(rep, file.path(dd, "gt_tier_replicate_scores.tsv"), sep = "\t")
  msg(sprintf("d = %.2f: %d tier x config rows, %d replicate rows, %d positives", d0, nrow(s), nrow(rep), n_pos))
  rep
}
reps = rbindlist(lapply(sort(unique(u$t)), write_d))
fwrite(reps, file.path(OUT, "invariance", "gt_tier_replicate_scores_all_depletions.tsv"), sep = "\t")
for (f in c("gt_tier_scoring.tsv", "gt_tier_scoring_with_se.tsv", "gt_tier_replicate_scores.tsv"))
  file.copy(file.path(OUT, "invariance", "d0.50", f), file.path(OUT, f), overwrite = TRUE)
msg("done -> ", OUT)
