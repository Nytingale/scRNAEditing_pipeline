#!/usr/bin/env Rscript
# =============================================================================
# TERATOMA_v8_dedup_compare_plot.R   (2026-10-13)
#
# UMI vs coordinate de-duplication vs unprocessed BAMs, ground truth v8 (gt_eval_v8/README.md).
# Layout of TERATOMA_dedup_compare_plot.R (v6); panel b also carries the unprocessed BAM, which
# replaces TERATOMA_dedup_audit_plot.R (depth at the v3 TRUE BED).
#   (a) Reads removed by de-duplication, per arm (points = teratomas)
#   (b) v8 truth sites (expected gate, t = 50%) with >= 10 reads, per arm and BAM type
#   (c) F1 per configuration, UMI vs coordinate, v8 t = 50% (expected and depth gates), per tier
#   (d) v8 truth sites with a correct call in NTC by either BAM (any reported configuration)
#   (e) Mismatch types among sites called in NTC (all 12 types; v8 call state)
# Unprocessed BAMs are not a reported configuration (scored, not ranked).
#
# Inputs: human/pm/gt_eval_v8/dedup_compare/ (gt_v8_dedup_inputs.R)
#         human/pm/gt_eval_v8/scoring/mismatch/<tier>/gt_v7_scores_reported.tsv, gt_v7_calls_NTC.tsv.gz
# Usage:  Rscript TERATOMA_v8_dedup_compare_plot.R [DATASET]
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "TERATOMA_common.R"))
argv <- commandArgs(trailingOnly = TRUE)
teratoma_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT <- caller_script()

DD    <- file.path(V8_DIR, "dedup_compare")
READS <- file.path(DD, "dedup_read_counts.tsv")
DEPTH <- file.path(DD, "dedup_truth_site_depth.tsv.gz")
TIERS <- c("arm", "sample", "celltype")
SCORE <- setNames(file.path(v8_scoring_dir(), TIERS, "gt_v7_scores_reported.tsv"), TIERS)
CALLS <- setNames(file.path(v8_scoring_dir(), c("arm", "sample"), "gt_v7_calls_NTC.tsv.gz"), c("arm", "sample"))
DISC  <- setNames(file.path(DD, "discovery", c("arm", "sample"), "gt_eval_discovery.tsv.gz"), c("arm", "sample"))
stopifnot(all(file.exists(c(READS, DEPTH, SCORE, CALLS, DISC, V8_TRUTH_SITES))))
truth <- fread(V8_TRUTH_SITES)
POS50 <- truth[gate == "expected" & class_t50 == "positive", site]
N_TRUTH <- length(POS50)
stopifnot(N_TRUTH == 1462L)
MIN_COV <- 10L
BS <- 13

DEDUP_LAB <- c(unprocessed = "Unprocessed", coord_dedup = "Coordinate", umi_dedup = "UMI")
BT_TO_DD  <- c(sampleRmDup = "coord_dedup", sampleUmiDup = "umi_dedup")
cmap_dd   <- c(Unprocessed = unname(cmap_proc["Unprocessed"]),
               Coordinate  = unname(cmap_proc["De-duplicated (coordinate)"]),
               UMI         = unname(cmap_proc["De-duplicated (UMI)"]))
tera_arms <- c("NTC_ctrl", "AAVS1_ctrl", "ADAR1_KO", "ADARB1_KO", "ADARB2_KO")
arm_lab   <- c(NTC_ctrl = "NTC", AAVS1_ctrl = "AAVS1", ADAR1_KO = "ADAR1\nKO",
               ADARB1_KO = "ADARB1\nKO", ADARB2_KO = "ADARB2\nKO")
tier_lab  <- c(arm = "Arm (T1-T4 merged)", sample = "Sample (mean of T1-T4)", celltype = sprintf("Cell type (mean of %d)", V8_N_CELLTYPES))
gate_lab  <- V8_GATE_LAB
zcat <- function(f, ...) fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t", showProgress = FALSE, ...)
split_unit <- function(d) {
  d[, `:=`(teratoma = sub("_.*$", "", unit), arm = sub("^T[0-9]+_", "", unit))]
  d[, arm := factor(arm, levels = tera_arms)]
}
dd_theme <- function(facets = FALSE) theme_manuscript(base_size = BS, facets = facets)
pd <- position_dodge(width = 0.8)
jd <- position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = manuscript_seed)

# ── (a) Reads removed ───────────────────────────────────────────────────────
rc <- dcast(fread(READS), unit ~ bam_type, value.var = "n_reads")
rc <- melt(rc, id.vars = c("unit", "unprocessed"), measure.vars = c("coord_dedup", "umi_dedup"),
           variable.name = "bam_type", value.name = "n_dedup")
rc[, `:=`(pct_removed = 100 * (1 - n_dedup / unprocessed),
          dedup = factor(DEDUP_LAB[as.character(bam_type)], levels = DEDUP_LAB))]
split_unit(rc)
rc_arm <- rc[, .(pct_removed = mean(pct_removed)), by = .(arm, dedup)]
p.a <- ggplot(rc, aes(arm, pct_removed, fill = dedup)) +
  geom_col(data = rc_arm, position = pd, width = 0.75, colour = "grey20", linewidth = 0.2) +
  geom_point(position = jd, size = 1.6, alpha = 0.8, show.legend = FALSE) +
  scale_x_discrete(labels = arm_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = cmap_dd, name = "BAM", drop = FALSE) +
  labs(x = NULL, y = "Reads removed (%)", title = "Reads removed") +
  dd_theme() + grid_y

# ── (b) v8 truth sites with >= MIN_COV reads, three BAM types ───────────────
dep <- zcat(DEPTH)[, .(n_cov = sum(depth >= MIN_COV), mean_depth = mean(depth), n_sites = .N), by = .(unit, bam_type)]
dep[, dedup := factor(DEDUP_LAB[bam_type], levels = DEDUP_LAB)]
split_unit(dep)
dep_arm <- dep[, .(n_cov = mean(n_cov), mean_depth = mean(mean_depth)), by = .(arm, dedup)]
p.b <- ggplot(dep, aes(arm, n_cov, fill = dedup)) +
  geom_col(data = dep_arm, position = pd, width = 0.75, colour = "grey20", linewidth = 0.2) +
  geom_point(position = jd, size = 1.4, alpha = 0.8, show.legend = FALSE) +
  scale_x_discrete(labels = arm_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), labels = shorten_axis) +
  scale_fill_manual(values = cmap_dd, name = "BAM") +
  labs(x = NULL, y = paste0("Ground-truth sites with \u2265 ", MIN_COV, " reads (per teratoma)"),
       title = paste0("Coverage at the ", format(N_TRUTH, big.mark = ","), " ground-truth sites")) +
  dd_theme() + grid_y

# ── (c) F1 per configuration, UMI vs coordinate ─────────────────────────────
sc <- rbindlist(lapply(TIERS, function(ti) fread(SCORE[[ti]])[, tier := ti]))[t == "t50" & bam_type %in% names(BT_TO_DD)]
sc[, dd := BT_TO_DD[bam_type]]
key <- c("tier", "gate", "caller", "strand_mode", "snp_subtract")
sc_w <- dcast(sc, as.formula(paste(paste(key, collapse = " + "), "~ dd")),
              value.var = c("precision", "recall", "f1", "f1_err"))
sc_w[, `:=`(tier_lab = factor(tier_lab[tier], levels = tier_lab), gate_lab = factor(gate_lab[gate], levels = gate_lab),
            caller_lab = caller_factor(caller),
            strand_lab = factor(label_strandMode[strand_mode], levels = label_strandMode))]
stopifnot(!anyNA(sc_w$f1_coord_dedup), !anyNA(sc_w$f1_umi_dedup))
sc_w[, `:=`(err_c = fifelse(tier == "sample", f1_err_coord_dedup, NA_real_),
            err_u = fifelse(tier == "sample", f1_err_umi_dedup, NA_real_))]
p.c <- ggplot(sc_w, aes(f1_coord_dedup, f1_umi_dedup, colour = caller_lab, shape = strand_lab)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = f1_umi_dedup - err_u, ymax = f1_umi_dedup + err_u), width = 0, alpha = 0.5, na.rm = TRUE) +
  geom_errorbar(aes(xmin = f1_coord_dedup - err_c, xmax = f1_coord_dedup + err_c), width = 0, alpha = 0.5,
                orientation = "y", na.rm = TRUE) +
  geom_point(size = 2.2) +
  facet_grid(gate_lab ~ tier_lab, scales = "free", labeller = labeller(gate_lab = label_wrap_gen(14))) +
  scale_colour_manual(values = cmap_method, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(x = "F1, coordinate de-duplication", y = "F1, UMI de-duplication",
       title = "F1 per configuration, ground truth \u2265 50% editing reduction in KO (sample: \u00b1 SE across teratomas)") +
  dd_theme(facets = TRUE) + grid_both +
  theme(aspect.ratio = 1, legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.position = "bottom")

# ── (d) v8 truth sites with a correct call in NTC, any reported configuration ─
rep_cfg <- unique(sc$config)
calls <- rbindlist(lapply(names(CALLS), function(ti) {
  d <- zcat(CALLS[[ti]])[gate == "expected" & correct == TRUE & site %chin% POS50]
  d[, bam_type := sub("^[^_]+__([^_]+)__.*$", "\\1", config)]
  unique(d[config %in% rep_cfg, .(unit, site, bam_type)])[, tier := ti]
}))
ov <- dcast(calls, tier + unit + site ~ bam_type, fun.aggregate = length)
ov[, set := fifelse(sampleRmDup > 0 & sampleUmiDup > 0, "Both", fifelse(sampleRmDup > 0, "Coordinate only", "UMI only"))]
ov <- ov[, .N, by = .(tier, unit, set)]
ov[, unit_lab := factor(fifelse(tier == "arm", "Arm", sub("_.*$", "", unit)), levels = c("Arm", "T1", "T2", "T3", "T4"))]
ov[, set := factor(set, levels = c("UMI only", "Coordinate only", "Both"))]
cmap_ov <- c(Both = "grey60", "Coordinate only" = cmap_dd[["Coordinate"]], "UMI only" = cmap_dd[["UMI"]])
p.d <- ggplot(ov, aes(unit_lab, N, fill = set)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = shorten_axis) +
  scale_fill_manual(values = cmap_ov, name = NULL, breaks = rev(levels(ov$set))) +
  labs(x = NULL, y = "Ground-truth sites called correctly", title = "Ground-truth sites called in NTC") +
  dd_theme() + grid_y + theme(legend.position = "bottom", legend.direction = "vertical")

# ── (e) Mismatch spectrum of sites called in NTC (v8 call state) ────────────
all12 <- all_vartypes()
spec <- rbindlist(lapply(names(DISC), function(ti) {
  d <- zcat(DISC[[ti]], select = c("bam_type", "site", "varType"))
  rbindlist(lapply(c("sample", names(BT_TO_DD)), function(bt) {
    vt <- d[bam_type == bt & varType %in% all12, .N, by = .(site, varType)]
    setorder(vt, site, -N, varType)
    vt <- vt[, .SD[1], by = site]
    s <- merge(data.table(varType = all12), vt[, .(n = .N), by = varType], by = "varType", all.x = TRUE)
    s[is.na(n), n := 0L]
    s[, `:=`(tier = ti, bam_type = bt, n_sites = nrow(vt), pct = 100 * n / nrow(vt))]
  }))
}))
spec_plot <- spec[bam_type %in% names(BT_TO_DD)]
spec_plot[, `:=`(dedup = factor(DEDUP_LAB[BT_TO_DD[bam_type]], levels = DEDUP_LAB),
                 varType = factor(varType, levels = all12), tier_lab = factor(tier_lab[tier], levels = tier_lab))]
p.e <- ggplot(spec_plot, aes(varType, pct, fill = dedup)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "grey20", linewidth = 0.2) +
  facet_wrap(~ tier_lab, ncol = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = cmap_dd, name = "BAM") +
  labs(x = "Mismatch type (RNA strand)", y = "Sites called (%)",
       title = "Mismatch types among sites called in NTC") +
  dd_theme(facets = TRUE) + grid_y

# ── Report tables ───────────────────────────────────────────────────────────
fwrite(rc, file.path(REPORT_DIR, "DEDUP_v8_reads_removed.tsv"), sep = "\t")
fwrite(dep, file.path(REPORT_DIR, "DEDUP_v8_truth_site_coverage.tsv"), sep = "\t")
fwrite(spec, file.path(REPORT_DIR, "DEDUP_v8_mismatch_spectrum_NTC.tsv"), sep = "\t")
fwrite(sc_w[, !c("tier_lab", "gate_lab", "caller_lab", "strand_lab")],
       file.path(REPORT_DIR, "DEDUP_v8_scores_umi_vs_coord.tsv"), sep = "\t")
fwrite(ov, file.path(REPORT_DIR, "DEDUP_v8_truth_overlap.tsv"), sep = "\t")

# ── Compose ─────────────────────────────────────────────────────────────────
leg_dd <- cowplot::get_plot_component(p.b + theme(legend.position = "bottom"), "guide-box-bottom", return_all = TRUE)
nol <- theme(legend.position = "none")
row1 <- cowplot::plot_grid(p.a + nol, p.b + nol, p.d, labels = c("a", "b", "d"), nrow = 1, label_size = 16,
                           align = "h", axis = "t", rel_widths = c(1, 1.15, 0.7))
row2 <- cowplot::plot_grid(p.c, p.e + nol, labels = c("c", "e"), nrow = 1, label_size = 16, rel_widths = c(1.35, 1))
p_final <- cowplot::plot_grid(row1, leg_dd, row2, ncol = 1, rel_heights = c(1, 0.06, 1.35))
STEM <- "TERATOMA_v8_dedup_compare"
srcs <- c(READS, DEPTH, SCORE, CALLS, DISC, V8_TRUTH_SITES)
save_figure(p_final, STEM, width = 18, height = 13.5, sources = srcs, script = SCRIPT)
save_figure(p.b, paste0(STEM, "_b_truth_coverage"), 7.5, 4.8, sources = DEPTH, data = dep, grid = "y", script = SCRIPT)
save_figure(p.c, paste0(STEM, "_c_f1"), 11, 7.5, sources = SCORE, data = sc_w, grid = "both", script = SCRIPT)
cat("\n=== v8 truth sites >= 10 reads, mean over teratomas ===\n")
print(dcast(dep_arm, arm ~ dedup, value.var = "n_cov"))
cat("\n=== F1 UMI - coordinate (t50), median per tier x gate ===\n")
print(sc_w[, .(median_delta = median(f1_umi_dedup - f1_coord_dedup), max_abs = max(abs(f1_umi_dedup - f1_coord_dedup))), by = .(tier, gate)])
message("Done -> figures/", DATASET, "/", STEM, "*.svg")
