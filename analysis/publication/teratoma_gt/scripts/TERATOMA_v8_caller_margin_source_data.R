#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_caller_margin_source_data.R   (2026-11-02; user: how much BCFtools wins by vs how much it loses by)
# Source data for the caller margin figure: per v8 truth set (6 KO + 3 REDIportal, main call state), each caller's
# best configuration by mean F1 across tiers, and the margin of the winning configuration over the best configuration
# of a different caller, absolute and as % of the winner's mean F1.
# Input : human/pm/gt_eval_v8/state_summary_pwcov/config_ranks.tsv (human/scripts/gt_v8_state_caller_summary.R)
# Output: source_data/pm_gt/TERATOMA_v8_caller_margin.{best,margin}.tsv  (drawn by TERATOMA_v8_caller_margin_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_caller_margin"

SRC <- file.path(V8_DIR, paste0("state_summary", V8_SFX), "config_ranks.tsv")
r <- fread(SRC)[state == V8_STATE]
stopifnot(setequal(r$truth, V8_TRUTH_IDS), r[, .N, by = truth][, all(N == as.integer(Sys.getenv("GT_N_CONFIG", "58")))])

best <- r[order(-mean_f1, overall_rank)][, .SD[1], by = .(truth, caller)][, .(truth, caller, config, mean_f1, overall_rank)]
# winner = best configuration overall; runner-up = best configuration of any other caller (2026-11-19)
m <- best[order(-mean_f1, overall_rank)][, {
  w <- .SD[1]; o <- .SD[caller != w$caller][1]
  .(winner = w$caller, winner_config = w$config, winner_f1 = w$mean_f1,
    other = o$caller, other_config = o$config, other_f1 = o$mean_f1,
    margin = w$mean_f1 - o$mean_f1, margin_pct = 100 * (w$mean_f1 - o$mean_f1) / w$mean_f1)
}, by = truth]
stopifnot(m[, all(margin >= 0)])

SHORT <- setNames(c(sprintf("Power-gated, \u2265 %d%%", c(25, 50, 75)), sprintf("Depth-gated, \u2265 %d%%", c(25, 50, 75)),
                    "\u2265 1 variant read", "\u2265 3 variant reads", "\u2265 3 variant reads, \u2265 5% editing"), V8_TRUTH_IDS)
m[, `:=`(Truth = factor(SHORT[truth], rev(SHORT)),
         Family = factor(fifelse(grepl("^redi_", truth), "REDIportal", "Editing reduction in KO"),
                         c("Editing reduction in KO", "REDIportal")),
         Winner = caller_factor(winner), Other = caller_factor(other))]
best[, method := caller_factor(caller)]
print(m[, .(truth, winner, other, winner_f1 = round(winner_f1, 4), other_f1 = round(other_f1, 4),
            margin = round(margin, 4), margin_pct = round(margin_pct, 2))])

write_sd(best, SRC_DIR, STEM, "best")
write_sd(m, SRC_DIR, STEM, "margin")
