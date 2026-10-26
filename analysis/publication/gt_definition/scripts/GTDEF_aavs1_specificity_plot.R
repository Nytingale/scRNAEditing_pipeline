#!/usr/bin/env Rscript
# =============================================================================
# GTDEF_aavs1_specificity_plot.R
#
# Is the v5 TRUE set enzyme-specific, or does it just pick up a generic
# CRISPR / cutting-associated drop in editing?
#
# AAVS1_ctrl is an isogenic cutting control (Cas9 + AAVS1 gRNA) that should
# NOT deplete ADAR substrates. For each tier's truth sites we plot
#   x = edit(AAVS1) / edit(NTC)
#   y = edit(KO)    / edit(NTC)
# on the same UMI merged-arm pileups used for Phase 1. Enzyme-specific TRUE
# sites sit at x ≈ 1 and y ≤ 0.5; a generic-perturbation artefact would sit
# near the diagonal y = x (both ratios depressed).
#
# Inputs
#   true_def_v5_full/{sample,dataset,celltype}/gt_truth_sites.tsv.gz
#   truth_combined_v5_umi/pileup_merged_{NTC_ctrl,AAVS1_ctrl,ADAR1_KO,
#                                       ADARB1_KO,ADARB2_KO}.tsv
#
# Outputs (publication/gt_definition/)
#   figures/pm_gt/GTDEF_hex_aavs1_vs_adar1_ratio_TRUE_by_tier.svg
#   figures/pm_gt/GTDEF_hex_aavs1_vs_ko_ratio_TRUE_by_arm_tier.svg
#   figures/pm_gt/GTDEF_hex_aavs1_vs_adar1_ratio_by_class_tier.svg
#   report/pm_gt/GTDEF_aavs1_specificity_summary.tsv
#
# Usage: Rscript GTDEF_aavs1_specificity_plot.R [DATASET]
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "GTDEF_common.R"))

argv <- commandArgs(trailingOnly = TRUE)
gtdef_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT <- caller_script()

PILE_DIR <- file.path(GT_DIR, "truth_combined_v5_umi")
V5_DIR   <- file.path(GT_DIR, "true_def_v5_full")
TIERS    <- c("sample", "dataset", "celltype")
KO_ARMS  <- c("ADAR1_KO", "ADARB1_KO", "ADARB2_KO")
min_cov  <- GT_MIN_COV
dep_cut  <- 1 - GT_DEPLETION   # 0.50 primary -> KO/NTC <= 0.5

hex_layer <- function(bins = 50, dpi = 900) {
  g <- geom_hex(bins = bins)
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::rasterize(g, dpi = dpi)
  } else {
    message("NOTE: ggrastr not installed; writing non-rasterized hex")
    g
  }
}

tier_title <- function(x) {
  factor(relabel(as.character(x), label_gt_tier),
         levels = c("Sample", "Arm", "Cell type"))
}

read_pileup <- function(arm) {
  f <- file.path(PILE_DIR, paste0("pileup_merged_", arm, ".tsv"))
  if (!file.exists(f)) stop("missing pileup: ", f)
  p <- fread(f, showProgress = FALSE)
  p[, site := paste0(chr, ":", pos)]
  p[, .(site, cov = as.numeric(coverage), edit = as.numeric(editing_level))]
}

message("reading merged UMI pileups from ", PILE_DIR)
pile <- list(
  NTC_ctrl   = read_pileup("NTC_ctrl"),
  AAVS1_ctrl = read_pileup("AAVS1_ctrl"),
  ADAR1_KO   = read_pileup("ADAR1_KO"),
  ADARB1_KO  = read_pileup("ADARB1_KO"),
  ADARB2_KO  = read_pileup("ADARB2_KO")
)
for (a in names(pile))
  message("  ", a, ": ", nrow(pile[[a]]), " covered sites")

# Per-tier truth classes + ratios -------------------------------------------
rows <- list()
for (tier in TIERS) {
  tf <- file.path(V5_DIR, tier, "gt_truth_sites.tsv.gz")
  if (!file.exists(tf)) stop("missing ", tf)
  tr <- fread(cmd = paste("zcat", shQuote(tf)), showProgress = FALSE)
  tr <- tr[, .(site, class, edit_NTC = edit_NTC_ctrl, cov_NTC = cov_NTC_ctrl,
               edit_ADAR1 = edit_ADAR1_KO, cov_ADAR1 = cov_ADAR1_KO,
               edit_ADARB1 = edit_ADARB1_KO, cov_ADARB1 = cov_ADARB1_KO)]
  aavs <- pile$AAVS1_ctrl[, .(site, cov_AAVS1 = cov, edit_AAVS1 = edit)]
  ad2  <- pile$ADARB2_KO[, .(site, cov_ADARB2 = cov, edit_ADARB2 = edit)]
  d <- merge(tr, aavs, by = "site", all.x = TRUE)
  d <- merge(d, ad2, by = "site", all.x = TRUE)
  for (j in c("cov_AAVS1", "edit_AAVS1", "cov_ADARB2", "edit_ADARB2"))
    set(d, which(is.na(d[[j]])), j, 0)
  d[, tier := tier]
  rows[[tier]] <- d
  message("tier ", tier, ": ", nrow(d), " sites, TRUE=", sum(d$class == "TRUE"))
}
sites <- rbindlist(rows, use.names = TRUE, fill = TRUE)

# Long ratios: one row per site × tier × ko_arm when both ratios are defined
mk_ratio <- function(edit_num, cov_num, edit_den, cov_den) {
  ok <- cov_num >= min_cov & cov_den >= min_cov & edit_den > 0
  fifelse(ok, edit_num / edit_den, NA_real_)
}
long <- rbindlist(lapply(KO_ARMS, function(arm) {
  if (arm == "ADAR1_KO") {
    en <- sites$edit_ADAR1; cn <- sites$cov_ADAR1
  } else if (arm == "ADARB1_KO") {
    en <- sites$edit_ADARB1; cn <- sites$cov_ADARB1
  } else {
    en <- sites$edit_ADARB2; cn <- sites$cov_ADARB2
  }
  data.table(
    site = sites$site, tier = sites$tier, class = sites$class, ko_arm = arm,
    ratio_aavs1 = mk_ratio(sites$edit_AAVS1, sites$cov_AAVS1,
                           sites$edit_NTC, sites$cov_NTC),
    ratio_ko = mk_ratio(en, cn, sites$edit_NTC, sites$cov_NTC)
  )
}))
long <- long[is.finite(ratio_aavs1) & is.finite(ratio_ko)]
long[, class_group := fcase(
  class == "TRUE", "TRUE",
  class %in% c("FP_known_not_depleted", "FP_novel", "FP_unsplit"),
    "Not depleted in KO",
  default = "Ambiguous"
)]
long[, class_group := factor(class_group,
  levels = c("TRUE", "Not depleted in KO", "Ambiguous"))]
long[, tier_lab := tier_title(tier)]
long[, ko_lab := factor(ko_arm, levels = KO_ARMS)]

# Summary --------------------------------------------------------------------
summ <- long[, {
  n <- .N
  .(n_sites = n,
    median_aavs1_ntc = round(median(ratio_aavs1), 3),
    median_ko_ntc    = round(median(ratio_ko), 3),
    pct_aavs1_ge_0.5 = round(100 * mean(ratio_aavs1 >= 0.5), 1),
    pct_ko_le_0.5    = round(100 * mean(ratio_ko <= dep_cut), 1),
    pct_enzyme_spec  = round(100 * mean(ratio_aavs1 >= 0.5 & ratio_ko <= dep_cut), 1),
    pct_generic_diag = round(100 * mean(ratio_aavs1 <= dep_cut & ratio_ko <= dep_cut), 1))
}, by = .(tier, ko_arm, class_group)]
setorder(summ, tier, ko_arm, class_group)
fwrite(summ, file.path(REPORT_DIR, "GTDEF_aavs1_specificity_summary.tsv"), sep = "\t")
message("wrote GTDEF_aavs1_specificity_summary.tsv")
print(summ[ko_arm == "ADAR1_KO" & class_group == "TRUE"])

src_pile <- file.path(PILE_DIR, paste0("pileup_merged_",
  c("NTC_ctrl", "AAVS1_ctrl", KO_ARMS), ".tsv"))
src_truth <- file.path(V5_DIR, TIERS, "gt_truth_sites.tsv.gz")
src_all <- c(src_pile, src_truth)

ratio_theme <- function() {
  theme_manuscript(facets = TRUE) +
    theme(panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          panel.spacing.x = unit(1.15, "lines"),
          legend.title.position = "top",
          aspect.ratio = 1)
}

# Annotate panels with n + enzyme-specific %
ann_lab <- function(dt) {
  dt[, .(tier_lab, ko_lab, class_group,
         lab = sprintf("n = %s\nenzyme-spec. %.0f%%",
                       formatC(.N, big.mark = ","),
                       100 * mean(ratio_aavs1 >= 0.5 & ratio_ko <= dep_cut))),
     by = .(tier_lab, ko_lab, class_group)]
}

# (1) TRUE sites, ADAR1 only, by tier ----------------------------------------
STEM <- "GTDEF_hex_aavs1_vs_adar1_ratio_TRUE_by_tier"
pts <- long[ko_arm == "ADAR1_KO" & class_group == "TRUE"]
ann <- ann_lab(pts)[, .(tier_lab, lab)]
p1 <- ggplot(pts, aes(ratio_aavs1, ratio_ko)) +
  hex_layer(bins = 60) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.3) +
  geom_hline(yintercept = dep_cut, linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40", linewidth = 0.3) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  facet_wrap(~ tier_lab, nrow = 1) +
  geom_text(data = ann, aes(x = 1.95, y = 0.08, label = lab),
            hjust = 1, vjust = 0, size = 2.8, lineheight = 0.95,
            inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
  labs(x = "AAVS1 / NTC editing ratio",
       y = "ADAR1 KO / NTC editing ratio",
       title = "TRUE sites: cutting control preserves editing; ADAR1 KO depletes it",
       subtitle = "Dashed lines mark the 50% depletion threshold; red = y = x (generic perturbation)") +
  ratio_theme()
save_figure(p1, STEM, width = 9.5, height = 3.8, sources = src_all,
            data = pts[, .(tier, site, ratio_aavs1, ratio_ko)],
            grid = "both", script = SCRIPT,
            rasterized = points_were_rasterized(nrow(pts)))

# (2) TRUE sites, all three KO arms × tier -----------------------------------
STEM <- "GTDEF_hex_aavs1_vs_ko_ratio_TRUE_by_arm_tier"
pts2 <- long[class_group == "TRUE"]
ann2 <- ann_lab(pts2)[, .(tier_lab, ko_lab, lab)]
p2 <- ggplot(pts2, aes(ratio_aavs1, ratio_ko)) +
  hex_layer(bins = 50) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.3) +
  geom_hline(yintercept = dep_cut, linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40", linewidth = 0.3) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  facet_grid(ko_lab ~ tier_lab) +
  geom_text(data = ann2, aes(x = 1.95, y = 0.08, label = lab),
            hjust = 1, vjust = 0, size = 2.4, lineheight = 0.95,
            inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
  labs(x = "AAVS1 / NTC editing ratio",
       y = "KO / NTC editing ratio",
       title = "TRUE sites vs cutting control: catalytic KOs deplete; ADARB2 does not",
       subtitle = "Dashed = 50% depletion; red = y = x") +
  ratio_theme()
save_figure(p2, STEM, width = 9.5, height = 8.2, sources = src_all,
            data = pts2[, .(tier, ko_arm, site, ratio_aavs1, ratio_ko)],
            grid = "both", script = SCRIPT,
            rasterized = points_were_rasterized(nrow(pts2)))

# (3) ADAR1: TRUE vs not-depleted FP, by tier --------------------------------
STEM <- "GTDEF_hex_aavs1_vs_adar1_ratio_by_class_tier"
pts3 <- long[ko_arm == "ADAR1_KO" & class_group %in% c("TRUE", "Not depleted in KO")]
ann3 <- ann_lab(pts3)[, .(tier_lab, class_group, lab)]
p3 <- ggplot(pts3, aes(ratio_aavs1, ratio_ko)) +
  hex_layer(bins = 50) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.3) +
  geom_hline(yintercept = dep_cut, linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40", linewidth = 0.3) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  facet_grid(class_group ~ tier_lab) +
  geom_text(data = ann3, aes(x = 1.95, y = 0.08, label = lab),
            hjust = 1, vjust = 0, size = 2.4, lineheight = 0.95,
            inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) +
  labs(x = "AAVS1 / NTC editing ratio",
       y = "ADAR1 KO / NTC editing ratio",
       title = "TRUE sites are ADAR1-depleted without AAVS1 depletion; FP sites are not",
       subtitle = "Dashed = 50% depletion; red = y = x (generic perturbation)") +
  ratio_theme()
save_figure(p3, STEM, width = 9.5, height = 6.2, sources = src_all,
            data = pts3[, .(tier, class_group, site, ratio_aavs1, ratio_ko)],
            grid = "both", script = SCRIPT,
            rasterized = points_were_rasterized(nrow(pts3)))

message("done")
