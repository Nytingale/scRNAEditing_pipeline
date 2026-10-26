#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_bam_processing_plot.R   (2026-10-13)
#
# GTDEF_bam_processing_plot.R on ground truth v8 (human/pm/gt_eval_v8/README.md):
#   mismatch types among sites called in NTC on each BAM type, before vs after restricting
#   to the 12-type ground truth, and the overlap of v8 truth sites across BAM types.
# Changes from the v6 script:
#   calls   v8 call state (nearby-mismatch pass; gt_eval_v8/dedup_compare/discovery from
#           gt_v8_dedup_inputs.R). varType is RNA frame for every caller (no REDItools3 flip).
#   12-type truth  expected-gate 12-type gate sites with d_max >= 0.5 that pass the v8 cleaning
#           (nearby-mismatch pass with PGP1 removed, not a PGP1 position; truth12_site_flags)
#   overlap v8 truth positives (expected, t = 50%) with a correct call, any reported configuration
#           of that BAM type (gt_eval_v8/scoring/mismatch/<tier>/gt_v7_calls_NTC.tsv.gz)
# The coordinate-vs-UMI hex panels are truth-independent (GTDEF_bam_processing_plot.R) and not redrawn.
# Unprocessed BAMs (bam_type `sample`) are kept in the report tables but not plotted.
# Outputs: GTDEF_v8_bar_mismatch_spectrum{,6}_by_bamtype_vs_TRUE, GTDEF_v8_upset_TRUE_sites_by_bamtype{,_<tier>}
#          (no Euler diagrams: user 2026-10-13)
#
# Usage: Rscript GTDEF_v8_bam_processing_plot.R [DATASET]
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggh4x)
  library(egg)
  library(forcats)
})
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "GTDEF_common.R"))

argv <- commandArgs(trailingOnly = TRUE)
gtdef_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT <- caller_script()

DISC_DIR   <- file.path(V8_DIR, "dedup_compare", "discovery")
TRUTH_GATE <- "expected"
TRUTH_D    <- 0.5
TRUTH12    <- file.path(V8_SUPPORT, "truth12_site_flags.tsv.gz")
# tier token -> discovery subdir
TIERS <- c(dataset = "arm", sample = "sample", celltype = "celltype")

# Strand-complement collapse of the 12 mismatch types -> 6 pairs.
# Order follows the first member of each pair in all_vartypes().
varType6_map <- c(
  AC = "AC", TG = "AC",
  AG = "AG", TC = "AG",
  AT = "AT", TA = "AT",
  CA = "CA", GT = "CA",
  CG = "CG", GC = "CG",
  CT = "CT", GA = "CT"
)
varType6_levels <- c("AC", "AG", "AT", "CA", "CG", "CT")

BT_ALL      <- c("sample", "sampleRmDup", "sampleUmiDup",
                 "recal", "recalq20", "recaldefault")
BT_REPORTED <- setdiff(BT_ALL, "sample")
BT_DEDUP    <- c("sampleRmDup", "sampleUmiDup")

tier_title <- function(x) {
  x <- as.character(x)
  # Pipeline token `dataset` → display "Arm" (arm-level merge across teratomas).
  x <- ifelse(x %in% c("dataset", "Dataset"), "Arm",
              ifelse(x %in% c("celltype", "Celltype"), "Cell type",
                     stringr::str_to_title(x)))
  factor(x, levels = c("Arm", "Sample", "Cell type"))
}

set_levels <- c("Sites called\nin NTC",
                sprintf("Ground-truth sites\n(\u2265 %d%% editing\nreduction in KO)", round(100 * TRUTH_D)))
set_title <- function(x) {
  factor(c(all_discovered = set_levels[1], pass_GT_criteria = set_levels[2])[as.character(x)],
         levels = set_levels)
}

bam_title <- function(x) {
  labs <- unname(label_bamType[as.character(x)])
  labs <- gsub(" (", "\n(", labs, fixed = T)
  factor(labs, levels = gsub(" (", "\n(", unname(label_bamType[BT_REPORTED]), fixed = T))
}

fill_spectrum <- function(spec, all_levels, key_cols) {
  # Ensure every group has every varType level (percent 0 if absent).
  ns <- unique(spec[, c(key_cols, "n_sites"), with = FALSE])
  grid <- unique(spec[, ..key_cols])[, .(varType = all_levels), by = key_cols]
  grid <- merge(grid, ns, by = key_cols, all.x = TRUE)
  m <- merge(grid, spec[, c(key_cols, "varType", "n"), with = FALSE],
             by = c(key_cols, "varType"), all.x = TRUE)
  m[is.na(n), n := 0L]
  m[is.na(n_sites), n_sites := 0L]
  m[, percent := fifelse(n_sites > 0L, as.numeric(n) / n_sites, 0)]
  m
}

# 12-type ground truth ----
if (!file.exists(TRUTH12)) stop("missing ", TRUTH12)
truth <- fread(cmd = paste("gzip -cd", shQuote(TRUTH12)), showProgress = FALSE)
truth <- unique(truth[gate == TRUTH_GATE & d_max >= TRUTH_D & nomm_pass_pgp1 & !pgp1,
                      .(site = paste0(chr, ":", pos), varType = type)])
message(sprintf("ground truth (%s, d >= %.2f): %d sites, %d A-to-G",
                TRUTH_GATE, TRUTH_D, nrow(truth), truth[varType == "AG", .N]))

# Per-bam discovered vs ground-truth spectra (12-type + collapsed 6-type) ----
# Before = unique sites called in NTC on that bam_type (majority varType)
# After  = those sites whose (site, varType) is in the ground truth
message("building per-bam discovered-vs-TRUE mismatch spectra")
all12 <- all_vartypes()

spec_rows <- list()
true_by_bam <- list()

for (tier_i in names(TIERS)) {
  disc_f <- file.path(DISC_DIR, TIERS[[tier_i]], "gt_eval_discovery.tsv.gz")
  if (!file.exists(disc_f)) stop("missing ", disc_f)
  d <- fread(cmd = paste("gzip -cd", shQuote(disc_f)),
             select = c("caller", "bam_type", "strand", "site", "varType"),
             showProgress = FALSE)
  keep_bt <- intersect(BT_ALL, unique(d$bam_type))
  true_by_bam[[tier_i]] <- list()

  for (bt in keep_bt) {
    x <- d[bam_type == bt & varType %in% all12]
    vt <- x[, .N, by = .(site, varType)]
    setorder(vt, site, -N, varType)
    vt <- vt[, .SD[1], by = site]
    n_disc <- nrow(vt)
    if (n_disc < 1L) next

    spec_d <- vt[, .(n = .N), by = varType]
    spec_d[, `:=`(tier = tier_i, bam_type = bt, set = "all_discovered",
                  n_sites = n_disc)]
    spec_rows[[length(spec_rows) + 1L]] <- spec_d

    vt_t <- vt[truth, on = .(site, varType), nomatch = NULL]
    n_true <- nrow(vt_t)
    true_by_bam[[tier_i]][[bt]] <- vt_t$site
    if (n_true < 1L) {
      spec_t <- data.table(varType = all12, n = 0L, tier = tier_i,
                           bam_type = bt, set = "pass_GT_criteria",
                           n_sites = 0L)
    } else {
      spec_t <- vt_t[, .(n = .N), by = varType]
      spec_t[, `:=`(tier = tier_i, bam_type = bt, set = "pass_GT_criteria",
                    n_sites = n_true)]
    }
    spec_rows[[length(spec_rows) + 1L]] <- spec_t
    message(sprintf("  %-8s %-14s disc=%7d AG=%.1f%% | TRUE∩=%5d AG=%.1f%%",
                    tier_i, bt, n_disc,
                    100 * spec_d[varType == "AG", n] / n_disc,
                    n_true,
                    if (n_true) 100 * spec_t[varType == "AG", sum(n)] / n_true else 0))
  }
  rm(d); gc(verbose = F)
}

both_bam <- rbindlist(spec_rows, fill = TRUE)
both_bam <- fill_spectrum(both_bam, all12, c("tier", "bam_type", "set"))
both_bam[, reported := bam_type %in% BT_REPORTED]
rp_bam_vs_true <- file.path(REPORT_DIR, "GTDEF_v8_mismatch_spectrum_by_bamtype_vs_TRUE.tsv")
fwrite(both_bam, rp_bam_vs_true, sep = "\t")
srcs <- c(rp_bam_vs_true, file.path(DISC_DIR, TIERS, "gt_eval_discovery.tsv.gz"), TRUTH12)

plot_bam_spectrum <- function(dt, stem, title, x_lab, fill_vals, width, height) {
  d <- copy(dt)
  d[, varType := factor(as.character(varType), levels = names(fill_vals))]
  d[, tier := tier_title(tier)]
  d[, set := set_title(set)]
  d[, bam_label := bam_title(bam_type)]
  d[, pct_label := paste0(round(percent*100, digits = 1), "%")]
  d[varType != "AG", pct_label := NA_character_]

  p <- ggplot(d, aes(varType, percent, fill = varType)) +
    geom_col(colour = "grey20", linewidth = 0.2, width = 0.85) +
    geom_text(aes(label = pct_label), vjust = -1, hjust = 0.1) +
    ggh4x::facet_nested(tier ~ bam_label + set) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_y_continuous(labels = function(x) round(100 * x, 1),
                       expand = expansion(mult = c(0, 0.25)),
                       breaks = seq(0, 1, 0.25),
                       limits = c(0, 1)) +
    labs(title = title, x = x_lab, y = "Percentage of sites (%)") +
    theme_manuscript(facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          panel.spacing.x = unit(0.6, "lines"),
          legend.title.position = "none")
  save_figure(p, stem, width = width, height = height,
              sources = srcs, data = d, grid = "y", script = SCRIPT)
  invisible(p)
}

STEM <- "GTDEF_v8_bar_mismatch_spectrum_by_bamtype_vs_TRUE"
plot_bam_spectrum(
  both_bam[bam_type %in% BT_REPORTED],
  STEM,
  title = "Mismatch types by BAM processing: all sites called in NTC vs those that pass the ground-truth criteria (all 12 mismatch types)",
  x_lab = "Mismatch type (RNA strand)",
  fill_vals = cmap_varType,
  width = 19, height = 7.2
)

# Collapsed 6-type (strand complements); red bar = A>G + T>C
collapse6 <- function(dt) {
  x <- copy(dt)
  x[, varType6 := varType6_map[as.character(varType)]]
  x <- x[!is.na(varType6),
         .(n = sum(as.numeric(n)), n_sites = max(n_sites)),
         by = .(tier, bam_type, set, varType6)]
  x[, percent := fifelse(n_sites > 0L, n / n_sites, 0)]
  setnames(x, "varType6", "varType")
  x[, reported := bam_type %in% BT_REPORTED]
  x
}
both6 <- collapse6(both_bam)
rp_bam6 <- file.path(REPORT_DIR, "GTDEF_v8_mismatch_spectrum6_by_bamtype_vs_TRUE.tsv")
fwrite(both6, rp_bam6, sep = "\t")

plot_bam_spectrum6 <- function(dt, stem, title, width, height) {
  d <- copy(dt)
  d[, varType := factor(as.character(varType), levels = varType6_levels)]
  d[, tier := tier_title(tier)]
  d[, set := set_title(set)]
  d[, bam_label := bam_title(bam_type)]
  d[, pct_label := paste0(round(percent*100, digits = 1), "%")]
  d[as.character(varType) != "AG", pct_label := NA_character_]

  p <- ggplot(d, aes(varType, percent, fill = varType)) +
    geom_col(colour = "grey20", linewidth = 0.2, width = 0.85) +
    geom_text(aes(label = pct_label), vjust = -1) +
    ggh4x::facet_nested(tier ~ bam_label + set) +
    scale_fill_manual(values = cmap_varType, guide = "none") +
    scale_y_continuous(labels = function(x) round(100 * x, 1),
                       expand = expansion(mult = c(0, 0.25)),
                       breaks = seq(0, 1, 0.25),
                       limits = c(0, 1)) +
    labs(title = title,
         x = "Mismatch type (strand complements collapsed)",
         y = "Percentage of sites (%)") +
    theme_manuscript(facets = TRUE) + grid_y +
    theme(axis.text.x = element_text(size = 8),
          panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1, "lines"),
          panel.spacing.x = unit(0.6, "lines"))
  save_figure(p, stem, width = width, height = height,
              sources = c(rp_bam6, srcs[-1]),
              data = d, grid = "y", script = SCRIPT)
}

STEM <- "GTDEF_v8_bar_mismatch_spectrum6_by_bamtype_vs_TRUE"
plot_bam_spectrum6(
  both6[bam_type %in% BT_REPORTED],
  STEM,
  title = "Mismatch types by BAM processing, strand complements collapsed: all sites called in NTC vs those that pass the ground-truth criteria",
  width = 19, height = 7.2
)

# Ground-truth site overlap across BAM processing methods ----
message("ground-truth site overlap across BAM methods (Euler + Upset)")

POS50 <- fread(V8_TRUTH_SITES)[gate == TRUTH_GATE & class_t50 == "positive", site]
CALLS8 <- setNames(file.path(v8_scoring_dir(), TIERS, "gt_v7_calls_NTC.tsv.gz"), names(TIERS))
REP8 <- setNames(file.path(v8_scoring_dir(), TIERS, "gt_v7_scores_reported.tsv"), names(TIERS))
memb <- rbindlist(lapply(names(TIERS), function(tier_i) {
  rep_cfg <- unique(fread(REP8[[tier_i]])$config)
  d <- fread(cmd = paste("gzip -cd", shQuote(CALLS8[[tier_i]])), showProgress = FALSE)[
    gate == TRUTH_GATE & correct == TRUE & site %chin% POS50 & config %in% rep_cfg]
  d[, bam_type := sub("^[^_]+__([^_]+)__.*$", "\\1", config)]
  unique(d[, .(tier = tier_i, bam_type, site)])
}))
true_by_bam <- lapply(split(memb, by = "tier", keep.by = FALSE), function(x) split(x$site, x$bam_type))
memb[, reported := bam_type %in% BT_REPORTED]
rp_memb <- file.path(REPORT_DIR, "GTDEF_v8_TRUE_sites_by_bamtype.tsv")
fwrite(memb, rp_memb, sep = "\t")
srcs_memb <- c(rp_memb, CALLS8, REP8, V8_TRUTH_SITES)

# Upset across the reported BAM methods (ggplot; UpSetR not in env)
upset_from_membership <- function(memb_tier, bt_levels, title, max_combos = 40L) {
  # Presence matrix for sites in any set
  wide <- dcast(unique(memb_tier[, .(site, bam_type)]),
                site ~ bam_type, fun.aggregate = length)
  for (bt in bt_levels) if (!bt %in% names(wide)) wide[, (bt) := 0L]
  mat <- as.matrix(wide[, ..bt_levels])
  mat <- (mat > 0L) * 1L
  rownames(mat) <- wide$site
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  comb_key <- apply(mat, 1L, paste0, collapse = "")
  comb_n <- as.data.table(table(comb_key))
  setnames(comb_n, c("comb", "n"))
  setorder(comb_n, -n)
  if (nrow(comb_n) > max_combos) comb_n <- comb_n[seq_len(max_combos)]
  comb_n[, idx := .I]
  decode <- function(key) as.integer(strsplit(key, "")[[1]])
  bar_dt <- comb_n[, .(idx, n, comb)]
  mat_long <- rbindlist(lapply(seq_len(nrow(comb_n)), function(i) {
    bits <- decode(comb_n$comb[i])
    data.table(idx = i, bam_type = bt_levels,
               present = factor(bits == 1L, levels = c(FALSE, TRUE)))
  }))
  mat_long[, bam_label := factor(unname(label_bamType[bam_type]),
                                 levels = unname(label_bamType[bt_levels]))]
  set_n <- data.table(bam_type = bt_levels,
                      n = as.integer(colSums(mat)[bt_levels]))
  set_n[, bam_label := factor(unname(label_bamType[bam_type]),
                              levels = unname(label_bamType[bt_levels]))]

  p_bar <- ggplot(bar_dt, aes(idx, n)) +
    geom_col(fill = "#e6e9f5", width = 0.7, color = "black", linewidth = 0.25) +
    scale_x_continuous(breaks = bar_dt$idx, expand = expansion(add = 0.3)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       n.breaks = 6,
                       labels = shorten_axis) +
    labs(x = NULL, y = "Intersection size", title = title) +
    theme_manuscript() + grid_y +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = 14, hjust = 0.5))

  p_mat <- ggplot() +
    geom_point(data = mat_long[present == FALSE],
               aes(idx, forcats::fct_rev(bam_label),
                   colour = present, size = present)) +
    geom_line(data = mat_long[present == TRUE],
              aes(idx, forcats::fct_rev(bam_label),
                  group = idx),
              colour = "grey20", linewidth = 0.4) +
    geom_point(data = mat_long[present == TRUE],
               aes(idx, forcats::fct_rev(bam_label),
                   colour = present, size = present)) +
    scale_colour_manual(values = c("FALSE" = "grey85", "TRUE" = "black"),
                        guide = "none") +
    scale_size_manual(values = c("FALSE" = 1.6, "TRUE" = 2.8), guide = "none") +
    scale_x_continuous(breaks = bar_dt$idx, expand = expansion(add = 0.6)) +
    labs(x = "Intersection", y = NULL) +
    theme_manuscript() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "black"),
          panel.grid.major.x = element_blank())

  p_set <- ggplot(set_n, aes(n, forcats::fct_rev(bam_label), fill = bam_label)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.25, show.legend = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = shorten_axis) +
    scale_fill_manual(values = cmap_proc) +
    labs(x = "Set size", y = NULL) +
    theme_manuscript() + grid_x +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(5.5, 22, 5.5, 5.5))

  egg::ggarrange(p_bar, ggplot() + theme_void(),
                 p_mat, p_set, byrow = TRUE, nrow = 2,
                 heights = c(0.6, 0.4), widths = c(0.75, 0.25))
}

upset_plots <- list()
upset_src <- list()
for (tier_i in names(TIERS)) {
  mt <- memb[tier == tier_i & bam_type %in% BT_REPORTED]
  if (!nrow(mt)) next
  bt_lev <- intersect(BT_REPORTED, unique(mt$bam_type))
  wide_n <- length(unique(mt$site))
  upset_plots[[tier_i]] <- upset_from_membership(
    mt, bt_lev,
    title = paste0("Ground-truth sites called in NTC by BAM processing method: ",
                   as.character(tier_title(tier_i)[1]),
                   " (n = ", format(wide_n, big.mark = ","), ")"))
  upset_src[[tier_i]] <- mt[, .(n_sites = uniqueN(site)), by = bam_type]
}
for (tier_i in names(upset_plots)) {
  stem <- paste0("GTDEF_v8_upset_TRUE_sites_by_bamtype_", tier_file_stem(tier_i))
  save_figure(upset_plots[[tier_i]], stem, width = 10, height = 6,
              sources = srcs_memb,
              data = upset_src[[tier_i]], grid = "y", script = SCRIPT)
}

# Combined upset panel (tiers as rows)
upset_all_plots <- lapply(names(TIERS), function(tier_i) {
  mt <- memb[tier == tier_i & bam_type %in% BT_REPORTED]
  if (!nrow(mt)) return(NULL)
  bt_lev <- intersect(BT_REPORTED, unique(mt$bam_type))
  upset_from_membership(
    mt, bt_lev,
    title = as.character(tier_title(tier_i)[1]))
})
upset_all_plots <- Filter(Negate(is.null), upset_all_plots)
if (length(upset_all_plots)) {
  p_up_panel <- cowplot::plot_grid(plotlist = upset_all_plots, ncol = 1,
                                   labels = "auto", label_size = 14)
  save_figure(p_up_panel, "GTDEF_v8_upset_TRUE_sites_by_bamtype",
              width = 9.5, height = 12.5,
              sources = srcs_memb,
              data = memb[bam_type %in% BT_REPORTED,
                          .(n_sites = uniqueN(site)), by = .(tier, bam_type)],
              grid = "y", script = SCRIPT)
}

message("done")
