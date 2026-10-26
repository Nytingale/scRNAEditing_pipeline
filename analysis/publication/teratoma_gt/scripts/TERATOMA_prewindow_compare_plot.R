#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_prewindow_compare_plot.R
# Every method configuration scored twice against the v6 sweep truth: on its final
# call set (after the mismatch-window filters, the benchmark call set) and on the
# pre-window call set (step 4: after the coverage / editing filters, before the
# variable-region and clustering windows). HANDOFF 2026-09-29 UPDATE-21.
#
# Inputs
#   human/pm/gt_eval_v6/scoring/sweep_infneg/{unit_scores_real.tsv.gz,caller_best.tsv}
#   human/pm/gt_eval_v6/prewindow/scoring/sweep_infneg/{unit_scores_real.tsv.gz,caller_best.tsv}
#   call files for the A-to-G percentage: results_gt*_v6/<NTC unit>/04_filtered/**/_pre100
#   and gt_eval_v6/prewindow/results_<tier>/<NTC unit>/04_filtered/**/_step4_cov_edit_filtered
# Outputs (figures/pm_gt/)
#   TERATOMA_scatter_prewindow_vs_final_t{25,50,75}.svg   per-configuration metrics
# The best-F1 slope chart was replaced by TERATOMA_window_state_ranksum_plot.R (four window-filter states).
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_prewindow_compare_plot.R"

PM   <- file.path(REPO_DIR, "human", "pm")
V6   <- file.path(PM, "gt_eval_v6")
SW   <- c(final = file.path(V6, "scoring", "sweep_infneg"),
          prewindow = file.path(V6, "prewindow", "scoring", "sweep_infneg"))
RES  <- list(final = c(arm = "results_gt_v6", sample = "results_gt_persample_v6", celltype = "results_gt_celltype_v6"),
             prewindow = c(arm = "prewindow/results_arm", sample = "prewindow/results_sample",
                           celltype = "prewindow/results_celltype"))
CALL_RE <- c(final = "_pre100\\.tsv\\.gz$", prewindow = "_step4_cov_edit_filtered\\.tsv\\.gz$")
UNIT_RE <- c(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")
TIER_LAB <- c(arm = "Arm", sample = "Sample", celltype = "Cell type")
UNREPORTED_BAM <- "sample"
UNREPORTED_CALLER_BAM <- list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))

# ---- scores ----------------------------------------------------------------
sc <- rbindlist(lapply(names(SW), function(s) {
  u <- fread(file.path(SW[[s]], "unit_scores_real.tsv.gz"))
  u[, .(precision = mean(precision, na.rm = TRUE), recall = mean(recall), f1 = mean(f1, na.rm = TRUE),
        n_calls = mean(n_calls)), by = .(tier, t, config, caller)][, set := s]
}))

# ---- A-to-G percentage (median over units, as gt_teratoma_plot_inputs_v6.R) ---
ag <- rbindlist(lapply(names(RES), function(s) rbindlist(lapply(names(RES[[s]]), function(tr) {
  root <- file.path(PM, if (s == "final") RES[[s]][[tr]] else file.path("gt_eval_v6", RES[[s]][[tr]]))
  units <- grep(UNIT_RE[[tr]], list.dirs(root, recursive = FALSE, full.names = FALSE), value = TRUE)
  rbindlist(lapply(units, function(u) {
    base <- file.path(root, u, "04_filtered")
    f <- list.files(base, CALL_RE[[s]], recursive = TRUE, full.names = TRUE)
    p <- tstrsplit(sub(paste0("^", base, "/"), "", f), "/", fixed = TRUE)
    data.table(set = s, tier = tr, unit = u,
               config = paste(p[[1]], p[[2]], p[[3]], sub("^snp_", "", p[[4]]), sep = "__"),
               ag_percent = vapply(f, function(x) {
                 if (file.size(x) == 0) return(NA_real_)
                 v <- fread(x, select = "varType", showProgress = FALSE)$varType
                 100 * mean(v == "AG")
               }, numeric(1)))
  }))
}))))
ag <- ag[, .(ag_percent = median(ag_percent, na.rm = TRUE)), by = .(set, tier, config)]
sc <- merge(sc, ag, by = c("set", "tier", "config"), all.x = TRUE)
sc[, c("caller_", "bam_type", "strand_mode", "snp_subtract") := tstrsplit(config, "__", fixed = TRUE)]
sc[, reported := !(bam_type %in% UNREPORTED_BAM |
                   mapply(function(c, b) b %in% UNREPORTED_CALLER_BAM[[c]], caller, bam_type) &
                   caller %in% names(UNREPORTED_CALLER_BAM))]
sc <- sc[reported == TRUE][, caller_ := NULL]
sc[, `:=`(method = std_caller(caller), Tier = factor(TIER_LAB[tier], TIER_LAB))]
fwrite(sc, file.path(SRC_DIR, "TERATOMA_prewindow_vs_final_scores.tsv"), sep = "\t")

# ---- (a) per-configuration scatter, one figure per threshold -------------------
MET <- c(recall = "Recall", precision = "Precision", f1 = "F1", ag_percent = "A-to-G %")
w <- dcast(melt(sc, id.vars = c("set", "tier", "Tier", "t", "config", "method", "strand_mode"),
                measure.vars = names(MET), variable.name = "metric"),
           tier + Tier + t + config + method + strand_mode + metric ~ set, value.var = "value")
w[, Metric := factor(MET[as.character(metric)], MET)]
for (tt in c("t25", "t50", "t75")) {
  d <- w[t == tt]
  STEM <- paste0("TERATOMA_scatter_prewindow_vs_final_", tt)
  p <- ggplot(d, aes(prewindow, final, colour = method, shape = strand_mode)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 1.8, alpha = 0.8) +
    facet_wrap(~ Metric + Tier, scales = "free", ncol = 3,
               labeller = labeller(.multi_line = FALSE)) +
    scale_colour_manual_cvd(values = cmap_method, name = NULL) +
    scale_shape_manual(values = c(annotation = 16, hybrid = 17, caller = 15),
                       labels = label_strandMode, name = NULL) +
    guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2)) +
    labs(x = "Before window filters", y = "Final call set",
         title = sprintf("Before vs after the mismatch-window filters (positives >= %s%% depletion)",
                         sub("^t", "", tt))) +
    theme_manuscript() + theme(legend.position = "bottom", legend.box = "vertical")
  save_figure(p, STEM, width = 9.5, height = 12.5,
              sources = c(file.path(SW, "unit_scores_real.tsv.gz"), file.path(V6, "prewindow", "manifest.tsv")),
              data = d, grid = "both", script = SCRIPT)
}
