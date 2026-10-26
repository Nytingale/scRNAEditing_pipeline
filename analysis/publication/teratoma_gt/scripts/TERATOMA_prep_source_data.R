#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_prep_source_data.R      (stage 1 of 2 -- ANALYSIS ONLY, writes TSVs)
#
# Assembles the per-combination benchmark table: one row per unique
# caller x bam_type x strand_mode x snp_subtract per tier, with precision,
# recall, F1 and the A-to-G percentage.
#
# Reads
#   SCORING_DIR/gt_tier_scoring.tsv (+ _with_se)   (TERATOMA_common.R; v6 sweep scoring,
#       human/scripts/gt_teratoma_plot_inputs_v6.R, with ag_percent included)
#   human/pm/results_gt_persample/summary/collated/
#       all_results_summary.tsv   <- A-to-G only when the scoring table lacks ag_percent (v5)
#
# Usage: Rscript TERATOMA_prep_source_data.R [DATASET]
# ============================================================================

suppressPackageStartupMessages(library(data.table))
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "TERATOMA_common.R"))

argv <- commandArgs(trailingOnly = TRUE)
teratoma_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
msg <- function(...) cat("  ", ..., "\n", sep = "")
wr <- function(dt, stem) {
  f <- file.path(REPORT_DIR, paste0(stem, ".tsv")); fwrite(dt, f, sep = "\t")
  msg(stem, ".tsv  (", nrow(dt), " rows)"); invisible(f)
}

if (!file.exists(SCORING_TSV))
  stop("missing ", SCORING_TSV, "\nRun gt_eval_classify_v3.R with --discovery first.")

sc <- fread(SCORING_TSV)
numc <- intersect(c("precision", "recall_vs_union", "evaluability_rate",
                    "n_discovered", "n_true", "n_fp", "n_evaluable"), names(sc))
sc[, (numc) := lapply(.SD, as.numeric), .SDcols = numc]
sc[, f1 := fifelse(precision + recall_vs_union > 0,
                   2 * precision * recall_vs_union / (precision + recall_vs_union),
                   NA_real_)]

# Attach SE when gt_eval_tier_errorbars.R has been run.
if (file.exists(SCORING_SE_TSV)) {
  se <- fread(SCORING_SE_TSV)
  se_cols <- intersect(c("precision_se", "recall_se", "f1_se",
                         "precision_mean", "recall_mean", "f1_mean",
                         "se_method", "n_reps"), names(se))
  keyc <- intersect(c("tier", COMBO_DIMS), names(se))
  sc <- merge(sc, se[, c(keyc, se_cols), with = FALSE], by = keyc, all.x = TRUE)
  msg("SE columns joined from ", basename(SCORING_SE_TSV), " (",
      sum(is.finite(sc$f1_se)), " / ", nrow(sc), " rows with f1_se)")
} else {
  msg("WARNING: ", SCORING_SE_TSV, " not found; bars will have no error bars. ",
      "Run human/pm/scripts/gt_eval_tier_errorbars.R first.")
}

# ── A-to-G percentage, joined from the pipeline's collated summary ──────────
# The GT discovery table is filtered to varType == "AG" and so cannot supply
# this; it has to come from the pipeline's own per-run metrics.
#
# >>> UNITS. `ag_pct_step7` in the collated summary is a FRACTION (0-1),
# despite the name. It is converted to a PERCENTAGE here, because that is how
# A-to-G enrichment is reported throughout the manuscript, and stored under an
# unambiguous column name. <<<
if ("ag_percent" %in% names(sc)) {
  msg("A-to-G percentage already in the scoring table (v6: from the scored call files)")
} else if (file.exists(COLLATED)) {
  cl <- fread(COLLATED)
  keyc <- intersect(COMBO_DIMS, names(cl))
  if ("ag_pct_step7" %in% names(cl) && length(keyc)) {
    ag <- cl[, .(ag_percent = 100 * as.numeric(median(as.numeric(ag_pct_step7),
                                                      na.rm = TRUE))), by = keyc]
    sc <- merge(sc, ag, by = intersect(keyc, names(sc)), all.x = TRUE)
    msg("A-to-G percentage joined for ", sum(!is.na(sc$ag_percent)), " of ",
        nrow(sc), " rows; range ",
        round(min(sc$ag_percent, na.rm = TRUE), 1), "-",
        round(max(sc$ag_percent, na.rm = TRUE), 1), "%")
  }
} else msg("WARNING: collated summary not found; A-to-G percentage unavailable")

setorder(sc, tier, -f1)
wr(sc, "TERATOMA_method_performance")

for (tr in unique(sc$tier))
  msg("best in tier '", tr, "': ",
      sc[tier == tr][1, paste(caller, bam_type, strand_mode, snp_subtract,
                              "| F1 =", round(f1, 3))])

# Per-caller best combination, the form that survives "which reditools3?"
best <- sc[order(tier, -f1), .SD[1], by = .(tier, caller)]
wr(best[, c("tier", COMBO_DIMS, "n_discovered", "precision", "recall_vs_union",
            "f1", "ag_percent"), with = FALSE], "TERATOMA_best_per_caller")

cat("\nprep complete -> ", REPORT_DIR, "\n", sep = "")
