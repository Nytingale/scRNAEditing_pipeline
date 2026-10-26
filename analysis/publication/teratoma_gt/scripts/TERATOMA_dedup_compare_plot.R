#!/usr/bin/env Rscript
# =============================================================================
# TERATOMA_dedup_compare_plot.R
#
# UMI vs coordinate de-duplication on the v6 data (per-teratoma dedup BAMs,
# v6 ground truth, NTC-only scoring). Unprocessed BAMs are benchmarked and
# kept in the source tables but are not a reported configuration; they appear
# only as the denominator of "reads removed".
#
#   (a) Reads removed by de-duplication, per arm (points = teratomas)
#   (b) Ground-truth sites with >= 10 reads (BQ/MQ 30), per arm
#   (c) F1 per configuration, UMI vs coordinate (arm and sample tiers)
#   (d) Ground-truth sites called in NTC by either BAM (any configuration)
#   (e) Mismatch types among sites called in NTC (all 12 types, pre100)
#
# Inputs: human/pm/gt_eval_v6/dedup_compare/ (run_dedup_compare_v6.sh),
#         SCORING_DIR/gt_tier_scoring_with_se.tsv (v6 sweep, t = 50%, informative negatives;
#         gt_teratoma_plot_inputs_v6.R) and human/pm/gt_eval_v6/scoring/{arm,sample}/
#         gt_v6_calls_NTC.tsv.gz (panel d; TRUE = positive at t = 50%)
# Usage:  Rscript TERATOMA_dedup_compare_plot.R [DATASET]
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "TERATOMA_common.R"))

argv <- commandArgs(trailingOnly = TRUE)
teratoma_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT <- caller_script()

V6    <- file.path(REPO_DIR, "human", "pm", "gt_eval_v6")
DD    <- file.path(V6, "dedup_compare")
READS <- file.path(DD, "dedup_read_counts.tsv")
DEPTH <- file.path(DD, "dedup_truth_site_depth.tsv.gz")
# Scores: v6 sweep at t = 50% with informative negatives (GROUND_TRUTH_DEFINITION.md);
# error = site bootstrap (arm) / SE across teratomas (sample)
SCORE <- file.path(SCORING_DIR, "gt_tier_scoring_with_se.tsv")
CALLS <- c(arm    = file.path(V6, "scoring", "arm",    "gt_v6_calls_NTC.tsv.gz"),
           sample = file.path(V6, "scoring", "sample", "gt_v6_calls_NTC.tsv.gz"))
DISC  <- c(arm    = file.path(DD, "discovery", "arm",    "gt_eval_discovery.tsv.gz"),
           sample = file.path(DD, "discovery", "sample", "gt_eval_discovery.tsv.gz"))
N_TRUTH <- 1692L
MIN_COV <- 10L
BS <- 13

DEDUP_LAB  <- c(coord_dedup = "Coordinate", umi_dedup = "UMI")
BT_TO_DD   <- c(sampleRmDup = "coord_dedup", sampleUmiDup = "umi_dedup")
cmap_dd    <- c(Coordinate = unname(cmap_proc["De-duplicated (coordinate)"]),
                UMI        = unname(cmap_proc["De-duplicated (UMI)"]))
tera_arms  <- c("NTC_ctrl", "AAVS1_ctrl", "ADAR1_KO", "ADARB1_KO", "ADARB2_KO")
arm_lab    <- c(NTC_ctrl = "NTC", AAVS1_ctrl = "AAVS1", ADAR1_KO = "ADAR1\nKO",
                ADARB1_KO = "ADARB1\nKO", ADARB2_KO = "ADARB2\nKO")
tier_lab   <- c(arm = "Arm (T1-T4 merged)", sample = "Sample (mean of T1-T4)")
zcat <- function(f, ...) fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t",
                               showProgress = FALSE, ...)
split_unit <- function(d) {
  d[, `:=`(teratoma = sub("_.*$", "", unit), arm = sub("^T[0-9]+_", "", unit))]
  d[, arm := factor(arm, levels = tera_arms)]
}
dd_theme <- function(facets = FALSE) theme_manuscript(base_size = BS, facets = facets)

# ── (a) Reads removed ───────────────────────────────────────────────────────
rc <- dcast(fread(READS), unit ~ bam_type, value.var = "n_reads")
rc <- melt(rc, id.vars = c("unit", "unprocessed"),
           measure.vars = names(DEDUP_LAB), variable.name = "bam_type",
           value.name = "n_dedup")
rc[, `:=`(pct_removed = 100 * (1 - n_dedup / unprocessed),
          dedup = factor(DEDUP_LAB[as.character(bam_type)], levels = DEDUP_LAB))]
split_unit(rc)
rc_arm <- rc[, .(pct_removed = mean(pct_removed)), by = .(arm, dedup)]

pd <- position_dodge(width = 0.75)
p.a <- ggplot(rc, aes(arm, pct_removed, fill = dedup)) +
  geom_col(data = rc_arm, position = pd, width = 0.7, colour = "grey20", linewidth = 0.2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75,
                                             seed = manuscript_seed),
             size = 1.6, alpha = 0.8, show.legend = FALSE) +
  scale_x_discrete(labels = arm_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = cmap_dd, name = "De-duplication") +
  labs(x = NULL, y = "Reads removed (%)", title = "Reads removed") +
  dd_theme() + grid_y

# ── (b) Ground-truth sites with >= MIN_COV reads ────────────────────────────
dep <- zcat(DEPTH)
dep <- dep[bam_type %in% names(DEDUP_LAB),
           .(n_cov = sum(depth >= MIN_COV), mean_depth = mean(depth)),
           by = .(unit, bam_type)]
dep[, dedup := factor(DEDUP_LAB[bam_type], levels = DEDUP_LAB)]
split_unit(dep)
dep_arm <- dep[, .(n_cov = mean(n_cov)), by = .(arm, dedup)]

p.b <- ggplot(dep, aes(arm, n_cov, fill = dedup)) +
  geom_col(data = dep_arm, position = pd, width = 0.7, colour = "grey20", linewidth = 0.2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75,
                                             seed = manuscript_seed),
             size = 1.6, alpha = 0.8, show.legend = FALSE) +
  scale_x_discrete(labels = arm_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), labels = shorten_axis) +
  scale_fill_manual(values = cmap_dd, name = "De-duplication") +
  labs(x = NULL, y = paste0("Ground-truth sites with \u2265 ", MIN_COV, " reads"),
       title = paste0("Coverage at the ", format(N_TRUTH, big.mark = ","),
                      " ground-truth sites")) +
  dd_theme() + grid_y

# ── (c) F1 per configuration, UMI vs coordinate ─────────────────────────────
sc <- fread(SCORE)[tier %in% c("dataset", "sample")]
sc[, tier := fifelse(tier == "dataset", "arm", tier)]
sc[, `:=`(f1 = fifelse(precision + recall_vs_union > 0,
                       2 * precision * recall_vs_union / (precision + recall_vs_union), NA_real_),
          recall = recall_vs_union, precision_err = precision_se, recall_err = recall_se, f1_err = f1_se)]
sc <- sc[bam_type %in% names(BT_TO_DD)]
sc[, dd := BT_TO_DD[bam_type]]
key <- c("tier", "caller", "strand_mode", "snp_subtract")
metrics <- c("precision", "recall", "f1")
sc_w <- dcast(sc, as.formula(paste(paste(key, collapse = " + "), "~ dd")),
              value.var = c(metrics, paste0(metrics, "_err")))
sc_w[, `:=`(tier_lab = factor(tier_lab[tier], levels = tier_lab),
            caller_lab = label_caller[caller],
            strand_lab = factor(label_strandMode[strand_mode], levels = label_strandMode))]
lim_f1 <- range(sc_w[, c(f1_coord_dedup - f1_err_coord_dedup, f1_umi_dedup + f1_err_umi_dedup,
                         f1_coord_dedup + f1_err_coord_dedup, f1_umi_dedup - f1_err_umi_dedup)],
                na.rm = TRUE)
lim_f1 <- c(floor(lim_f1[1] * 20) / 20, ceiling(lim_f1[2] * 20) / 20)

p.c <- ggplot(sc_w, aes(f1_coord_dedup, f1_umi_dedup, colour = caller_lab, shape = strand_lab)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = f1_umi_dedup - f1_err_umi_dedup,
                    ymax = f1_umi_dedup + f1_err_umi_dedup), width = 0, alpha = 0.5) +
  geom_errorbar(aes(xmin = f1_coord_dedup - f1_err_coord_dedup,
                     xmax = f1_coord_dedup + f1_err_coord_dedup), width = 0, alpha = 0.5, orientation = "y") +
  geom_point(size = 2.4) +
  facet_wrap(~ tier_lab, nrow = 1) +
  coord_equal(xlim = lim_f1, ylim = lim_f1) +
  scale_colour_manual(values = cmap_method, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(x = "F1, coordinate de-duplication", y = "F1, UMI de-duplication",
       title = "F1 per configuration (\u00b1 SE)") +
  dd_theme(facets = TRUE) + grid_both +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))

# Precision and recall for the same configurations (supplementary panel)
pr <- rbindlist(lapply(c("precision", "recall"), function(m) {
  sc_w[, .(tier_lab, caller_lab, strand_lab, metric = m,
           coord = get(paste0(m, "_coord_dedup")), umi = get(paste0(m, "_umi_dedup")),
           coord_err = get(paste0(m, "_err_coord_dedup")), umi_err = get(paste0(m, "_err_umi_dedup")))]
}))
pr[, metric := factor(metric, levels = c("precision", "recall"), labels = c("Precision", "Recall"))]
p.c2 <- ggplot(pr, aes(coord, umi, colour = caller_lab, shape = strand_lab)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = umi - umi_err, ymax = umi + umi_err), width = 0, alpha = 0.5) +
  geom_errorbar(aes(xmin = coord - coord_err, xmax = coord + coord_err), width = 0, alpha = 0.5, orientation = "y") +
  geom_point(size = 2.2) +
  facet_grid(metric ~ tier_lab, scales = "free") +
  scale_colour_manual(values = cmap_method, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(x = "Coordinate de-duplication", y = "UMI de-duplication",
       title = "Precision and recall per configuration (\u00b1 SE)") +
  dd_theme(facets = TRUE) + grid_both +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))

# ── (d) Ground-truth sites called in NTC, any configuration ────────────────
calls <- rbindlist(lapply(names(CALLS), function(t) {
  d <- zcat(CALLS[[t]])[class == "TRUE"]
  d[, bam_type := sub("^[^_]+__([^_]+)__.*$", "\\1", config)]
  d <- unique(d[bam_type %in% names(BT_TO_DD), .(unit, site, bam_type)])
  d[, tier := t]
}))
ov <- dcast(calls, tier + unit + site ~ bam_type, fun.aggregate = length)
ov[, set := fifelse(sampleRmDup > 0 & sampleUmiDup > 0, "Both",
            fifelse(sampleRmDup > 0, "Coordinate only", "UMI only"))]
ov <- ov[, .N, by = .(tier, unit, set)]
ov[, unit_lab := fifelse(tier == "arm", "Arm", sub("_.*$", "", unit))]
ov[, `:=`(unit_lab = factor(unit_lab, levels = c("Arm", "T1", "T2", "T3", "T4")),
          set = factor(set, levels = c("UMI only", "Coordinate only", "Both")))]
cmap_ov <- c(Both = "grey60", "Coordinate only" = cmap_dd[["Coordinate"]],
             "UMI only" = cmap_dd[["UMI"]])

p.d <- ggplot(ov, aes(unit_lab, N, fill = set)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = shorten_axis) +
  scale_fill_manual(values = cmap_ov, name = NULL, breaks = rev(levels(ov$set))) +
  labs(x = NULL, y = "Ground-truth sites called",
       title = "Ground-truth sites called in NTC") +
  dd_theme() + grid_y

# ── (e) Mismatch spectrum of sites called in NTC ────────────────────────────
all12 <- all_vartypes()
spec <- rbindlist(lapply(names(DISC), function(t) {
  d <- zcat(DISC[[t]], select = c("bam_type", "site", "varType"))
  rbindlist(lapply(c("sample", names(BT_TO_DD)), function(bt) {
    vt <- d[bam_type == bt & varType %in% all12, .N, by = .(site, varType)]
    setorder(vt, site, -N, varType)
    vt <- vt[, .SD[1], by = site]
    s <- vt[, .(n = .N), by = varType]
    s <- merge(data.table(varType = all12), s, by = "varType", all.x = TRUE)
    s[is.na(n), n := 0L]
    s[, `:=`(tier = t, bam_type = bt, n_sites = nrow(vt), pct = 100 * n / nrow(vt))]
  }))
}))
spec_plot <- spec[bam_type %in% names(BT_TO_DD)]
spec_plot[, `:=`(dedup = factor(DEDUP_LAB[BT_TO_DD[bam_type]], levels = DEDUP_LAB),
                 varType = factor(varType, levels = all12),
                 tier_lab = factor(tier_lab[tier], levels = tier_lab))]

p.e <- ggplot(spec_plot, aes(varType, pct, fill = dedup)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75,
           colour = "grey20", linewidth = 0.2) +
  facet_wrap(~ tier_lab, ncol = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = cmap_dd, name = "De-duplication") +
  labs(x = "Mismatch type (transcript strand)", y = "Sites called (%)",
       title = "Mismatch types among sites called in NTC") +
  dd_theme(facets = TRUE) + grid_y

# ── Report tables (unprocessed kept here, not plotted) ─────────────────────
fwrite(rc, file.path(REPORT_DIR, "DEDUP_v6_reads_removed.tsv"), sep = "\t")
fwrite(zcat(DEPTH)[, .(n_cov = sum(depth >= MIN_COV), mean_depth = mean(depth),
                       n_sites = .N), by = .(unit, bam_type)],
       file.path(REPORT_DIR, "DEDUP_v6_truth_site_coverage.tsv"), sep = "\t")
fwrite(spec, file.path(REPORT_DIR, "DEDUP_v6_mismatch_spectrum_NTC.tsv"), sep = "\t")
fwrite(sc_w[, !c("tier_lab", "caller_lab", "strand_lab")],
       file.path(REPORT_DIR, "DEDUP_v6_scores_umi_vs_coord.tsv"), sep = "\t")

# ── Compose ─────────────────────────────────────────────────────────────────
leg_dd <- cowplot::get_plot_component(p.a + theme(legend.position = "bottom"),
                                      "guide-box-bottom", return_all = TRUE)
nol <- theme(legend.position = "none")
row1 <- cowplot::plot_grid(p.a + nol, p.b + nol, p.d, labels = c("a", "b", "d"),
                           nrow = 1, label_size = 16, align = "h", axis = "tb",
                           rel_widths = c(1, 1, 0.8))
row2 <- cowplot::plot_grid(p.c, p.e + nol, labels = c("c", "e"), nrow = 1,
                           label_size = 16, rel_widths = c(1.25, 1))
p_final <- cowplot::plot_grid(row1, leg_dd, row2, ncol = 1, rel_heights = c(1, 0.07, 1.25))

STEM <- "TERATOMA_dedup_compare"
srcs <- c(READS, DEPTH, SCORE, CALLS, DISC)
save_figure(p_final, STEM, width = 16, height = 12, sources = srcs, script = SCRIPT)
save_figure(p.a,  paste0(STEM, "_a_reads_removed"), 6.5, 4.8, sources = READS, data = rc, grid = "y", script = SCRIPT)
save_figure(p.b,  paste0(STEM, "_b_truth_coverage"), 6.5, 4.8, sources = DEPTH, data = dep, grid = "y", script = SCRIPT)
save_figure(p.c,  paste0(STEM, "_c_f1"), 9, 6, sources = SCORE, data = sc_w, grid = "both", script = SCRIPT)
save_figure(p.c2, paste0(STEM, "_c2_precision_recall"), 9, 8.5, sources = SCORE, data = pr, grid = "both", script = SCRIPT)
save_figure(p.d,  paste0(STEM, "_d_truth_overlap"), 5, 4.8, sources = CALLS, data = ov, grid = "y", script = SCRIPT)
save_figure(p.e,  paste0(STEM, "_e_mismatch_spectrum"), 8, 6, sources = DISC, data = spec_plot, grid = "y", script = SCRIPT)

message("Done -> figures/", DATASET, "/", STEM, "*.svg")
