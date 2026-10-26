#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_window_state_ranksum_plot.R   (2026-10-05)
# Caller rank by cross-tier F1 rank sum under the four window-filter states, one panel per
# truth set (the nine truth sets of NINE_TRUTHS). Replaces the two-state slope chart of
# TERATOMA_prewindow_compare_plot.R section b.
#   none      = step 4, no window filter
#   mismatch  = step 5, "Remove if near other MM types" only (pass 1)
#   clustered = "Clustered" only (nomm reruns, pipeline default from 2026-10-01)
#   both      = step 6, pass 1 then pass 2 (the published v6 benchmark call sets)
# Rank sum as TERATOMA_cross_tier_grid_plot.R: F1 from the tier precision and recall_vs_union,
# configurations ranked within each truth x state x tier (ties = min), ranks summed over the
# three tiers, overall rank by rank sum with ties broken by the higher mean F1. A caller is placed
# by its best-ranked configuration; the figure shows the callers' order (1-5). The
# gt_expected_d0.50 / both table must equal the grid source data (checked below).
#
# Inputs : human/pm/gt_eval_v6/{prewindow/scoring,prewindow/scoring_step5,nomm/scoring,scoring}/
#            sweep_infneg{,_depth}/teratoma_inputs/invariance/d0.{25,50,75}/gt_tier_scoring.tsv
#          human/pm/gt_eval_v6/rediportal/scoring{_none,_mismatch,_clustered,}/gt_tier_scoring_<set>.tsv
# Outputs: figures/pm_gt/TERATOMA_slope_window_state_ranksum.svg
#          source_data/pm_gt/TERATOMA_window_state_config_ranksum.tsv   every configuration
#          source_data/pm_gt/TERATOMA_window_state_rank_cor.tsv        state-pair rank correlations
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_window_state_ranksum_plot.R"

V6 <- file.path(REPO_DIR, "human", "pm", "gt_eval_v6")
STATE_SCORING <- list(none = file.path(V6, "prewindow", "scoring"), mismatch = file.path(V6, "prewindow", "scoring_step5"),
                      clustered = file.path(V6, "nomm", "scoring"), both = file.path(V6, "scoring"))
STATE_REDI <- c(none = "scoring_none", mismatch = "scoring_mismatch", clustered = "scoring_clustered", both = "scoring")
UNREPORTED_BAM <- "sample"
UNREPORTED_CALLER_BAM <- list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))
STATE_LAB <- c(none = "No local\ncontext filters", mismatch = "+ Remove if near\nother MM types", clustered = "+ Clustered",
               both = "+ Remove if near\nother MM types\n+ Clustered")
cmap_state <- setNames(c("grey60", "#FDAE6B", "#6BAED6", "#54278F"), names(STATE_LAB))
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
is_unreported <- function(caller, bam_type)
  bam_type %in% UNREPORTED_BAM |
    mapply(function(c, b) b %in% UNREPORTED_CALLER_BAM[[c]], caller, bam_type) & caller %in% names(UNREPORTED_CALLER_BAM)
grid_f1 <- function(p, r) fifelse(p + r > 0, 2 * p * r / (p + r), NA_real_)

# ---- truth sets (order = panel order) --------------------------------------------------------
DEPL <- c("0.25", "0.50", "0.75")
REDI_ID <- c("cov10_alt1", "alt3", "alt3_e05")
TRUTH <- data.table(
  truth = c(paste0("gt_expected_d", DEPL), paste0("gt_depth_d", DEPL), paste0("redi_", REDI_ID)),
  label = c(sprintf("Editing reduction in KO \u2265 %d%%\nexpected KO alt reads \u2265 3", c(25, 50, 75)),
            sprintf("Editing reduction in KO \u2265 %d%%\nNTC and KO\ncoverage \u2265 10", c(25, 50, 75)),
            "REDIportal\n\u2265 10 reads, \u2265 1 G", "REDIportal\n\u2265 10 reads, \u2265 3 G",
            "REDIportal\n\u2265 10 reads, \u2265 3 G, \u2265 5% editing level"))
truth_file <- function(id, s) {
  if (startsWith(id, "redi_"))
    return(file.path(V6, "rediportal", STATE_REDI[[s]], paste0("gt_tier_scoring_", sub("^redi_", "", id), ".tsv")))
  g <- if (grepl("^gt_expected_", id)) "sweep_infneg" else "sweep_infneg_depth"
  file.path(STATE_SCORING[[s]], g, "teratoma_inputs", "invariance", sub("^gt_[a-z]+_", "", id), "gt_tier_scoring.tsv")
}
SRC <- unlist(lapply(TRUTH$truth, function(id) vapply(names(STATE_SCORING), function(s) truth_file(id, s), "")))
stopifnot(all(file.exists(SRC)))

# ---- per-configuration rank sums -------------------------------------------------------------
sc <- rbindlist(lapply(TRUTH$truth, function(id) rbindlist(lapply(names(STATE_SCORING), function(s)
  fread(truth_file(id, s))[, c("tier", COMBO, "precision", "recall_vs_union", "f1"), with = FALSE][
    , `:=`(truth = id, state = s)]))))
sc[tier == "dataset", tier := "arm"]
stopifnot(!any(is_unreported(sc$caller, sc$bam_type)), sc[, .N, by = .(truth, state, tier)][, uniqueN(N) == 1])
N_CONFIG <- sc[, .N, by = .(truth, state, tier)]$N[1]
sc[, f1_grid := grid_f1(precision, recall_vs_union)]
sc[, f1_rank := frank(-f1_grid, ties.method = "min"), by = .(truth, state, tier)]
rs <- sc[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1_grid, na.rm = TRUE)), by = c("truth", "state", COMBO)]
rs[, overall_rank := frank(list(rank_sum, -mean_f1), ties.method = "min"), by = .(truth, state)]
setorder(rs, truth, state, overall_rank)
fwrite(rs, file.path(SRC_DIR, "TERATOMA_window_state_config_ranksum.tsv"), sep = "\t")

GRID_SRC <- file.path(SRC_DIR, "TERATOMA_cross_tier_grid.source_data.tsv")
chk <- merge(rs[truth == "gt_expected_d0.50" & state == "both"], fread(GRID_SRC), by = COMBO, suffixes = c("", ".grid"))
stopifnot(nrow(chk) == N_CONFIG, chk[, all(rank_sum == rank_sum.grid & overall_rank == overall_rank.grid)])

# ---- caller level: best-ranked configuration -------------------------------------------------
cb <- rs[order(overall_rank), .SD[1], by = .(truth, state, caller)]
cb[, caller_rank := frank(overall_rank, ties.method = "min"), by = .(truth, state)]
stopifnot(cb[, !anyDuplicated(caller_rank), by = .(truth, state)]$V1)
cb[, `:=`(method = std_caller(caller), State = factor(STATE_LAB[state], STATE_LAB),
          Truth = factor(TRUTH$label[match(truth, TRUTH$truth)], TRUTH$label))]
print(dcast(cb, truth + method ~ state, value.var = "caller_rank")[, c("truth", "method", names(STATE_LAB)), with = FALSE])

# ---- rank correlation of per-configuration F1 between states --------------------------------
f1 <- copy(sc)[, config := do.call(paste, c(.SD, sep = "__")), .SDcols = COMBO]
pairs <- combn(names(STATE_SCORING), 2, simplify = FALSE)
spear <- function(x, y) cor(x, y, method = "spearman", use = "complete.obs")
w <- dcast(f1, truth + tier + config ~ state, value.var = list("f1", "f1_grid"))
rc <- rbindlist(lapply(pairs, function(p) w[, .(state_a = p[1], state_b = p[2], n_config = .N,
  rho_f1 = spear(get(paste0("f1_", p[1])), get(paste0("f1_", p[2]))),
  rho_f1_grid = spear(get(paste0("f1_grid_", p[1])), get(paste0("f1_grid_", p[2])))), by = .(truth, tier)]))
cw <- dcast(cb, truth + caller ~ state, value.var = "caller_rank")
rc_caller <- rbindlist(lapply(pairs, function(p) cw[, .(state_a = p[1], state_b = p[2],
  rho_caller_rank = spear(get(p[1]), get(p[2])),
  same_top_caller = caller[which.min(get(p[1]))] == caller[which.min(get(p[2]))]), by = truth]))
rc <- merge(rc, rc_caller, by = c("truth", "state_a", "state_b"), all.x = TRUE)
rc[, tier := factor(tier, c("arm", "sample", "celltype"))]
setorder(rc, truth, tier, state_a, state_b)
fwrite(rc, file.path(SRC_DIR, "TERATOMA_window_state_rank_cor.tsv"), sep = "\t")
print(dcast(rc[truth == "gt_expected_d0.50"], state_a + state_b ~ tier, value.var = c("rho_f1", "rho_f1_grid")), digits = 3)

# ---- figure: one panel per truth set ---------------------------------------------------------
bands <- data.table(State = factor(STATE_LAB, STATE_LAB), state = names(STATE_LAB))
p <- ggplot(cb, aes(State, caller_rank, colour = method, group = method)) +
  geom_rect(data = bands, aes(xmin = as.integer(State) - 0.46, xmax = as.integer(State) + 0.46, fill = state),
            ymin = -Inf, ymax = Inf, alpha = 0.12, inherit.aes = FALSE) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 6) +
  geom_text(aes(label = caller_rank), colour = "white", size = 3.4, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ Truth, ncol = 3) +
  scale_x_discrete(limits = levels(cb$State)) +
  scale_y_reverse(breaks = 1:5, expand = expansion(add = 0.4)) +
  scale_fill_manual(values = cmap_state, guide = "none") +
  scale_colour_manual_cvd(values = cmap_method, name = NULL) +
  labs(x = NULL, y = "Caller rank by cross-tier F1 rank sum") +
  theme_manuscript(base_size = 16, facets = TRUE) + grid_y +
  theme(axis.text.x = element_text(size = 12, lineheight = 0.9, angle = 35, hjust = 1, vjust = 1),
        strip.text = element_text(size = 13, lineheight = 0.95), legend.position = "bottom",
        legend.text = element_text(size = 14), panel.spacing = unit(6, "pt")) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
save_figure(p, "TERATOMA_slope_window_state_ranksum", width = 11, height = 12, sources = c(SRC, GRID_SRC),
            data = cb, grid = "y", script = SCRIPT)
