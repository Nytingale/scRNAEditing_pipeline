#!/usr/bin/env Rscript
# gt_eval_classify.R  (dataset-agnostic GT evaluator, step 3/3)
#
# Site-matched, coverage-controlled evaluation of discovered A>G sites for an
# ADAR KO/KD-vs-control ground truth. Editing levels are taken from the
# INDEPENDENT pileups (gt_site_pileup.py), never from a caller's own counts.
#
# For every method combination (caller x bam_type x strand_mode) and every
# treatment-vs-control contrast:
#   * reference set = A>G sites discovered by that combo in either arm
#   * evaluable set = reference sites with total coverage >= min_cov in BOTH arms
#     (sufficiently covered in both conditions)
#   * per site, editing level = edited/(edited+unedited), strand-aware:
#       '+': edited=G, unedited=A     '-': edited=C, unedited=T
#   * classify:
#       true_like  : edited in control (>= min_ctrl_edit) AND editing DECREASED
#                    in KO/KD (delta >= min_abs_drop AND rel_drop >= min_drop_frac)
#                    -> ADAR-dependent editing == very likely a real site
#       fp_like    : substantial editing in KO/KD (>= min_ctrl_edit) that is NOT
#                    decreased (delta < min_abs_drop) -> likely false positive
#                    (incl. sites edited at high level in KO/KD)
#       ambiguous  : neither (both low, or small/borderline change)
#
# The headline number is: of sites a method discovered that are evaluable
# (covered in both), how many are true_like vs fp_like.
#
# Outputs (in --outdir):
#   gt_eval_persite.tsv.gz    per site x combo x contrast with levels + category
#   gt_eval_summary.tsv       counts/fractions per combo x contrast
#   plots/gt_eval_scatter_<treat>_vs_<ctrl>.svg
#   plots/gt_eval_bench_<treat>_vs_<ctrl>.svg
#   plots/gt_eval_fraction_true.svg
#
# Usage:
#   Rscript gt_eval_classify.R --discovery D/gt_eval_discovery.tsv.gz \
#       --pileup_dir D --control_arms NTC_ctrl --treatment_arms ADAR_KD,ADARB1_KD \
#       --outdir D --dataset_label replogle_K562 \
#       [--min_cov 10 --min_ctrl_edit 0.05 --min_abs_drop 0.02 --min_drop_frac 0.25] \
#       [--plotting_functions .../plotting/plotting_functions.R]

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

`%||%` <- function(x, y) if (is.null(x) || is.na(x) || !length(x)) y else x
parse_cli <- function() {
  a <- commandArgs(trailingOnly = TRUE); out <- list(); i <- 1L
  while (i <= length(a)) {
    if (grepl("^--", a[i])) { out[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L } else i <- i + 1L
  }
  out
}
args <- parse_cli()
discovery_f <- args[["discovery"]]  %||% stop("--discovery required")
pileup_dir  <- args[["pileup_dir"]] %||% stop("--pileup_dir required")
outdir      <- args[["outdir"]]     %||% stop("--outdir required")
ctrl_arms   <- strsplit(args[["control_arms"]]   %||% stop("--control_arms required"), ",")[[1]]
treat_arms  <- strsplit(args[["treatment_arms"]] %||% stop("--treatment_arms required"), ",")[[1]]
min_cov       <- as.numeric(args[["min_cov"]]       %||% 10)
min_ctrl_edit <- as.numeric(args[["min_ctrl_edit"]] %||% 0.05)
min_abs_drop  <- as.numeric(args[["min_abs_drop"]]  %||% 0.02)
min_drop_frac <- as.numeric(args[["min_drop_frac"]] %||% 0.25)
dataset_label <- args[["dataset_label"]] %||% "GT"
pf <- args[["plotting_functions"]] %||%
  file.path(Sys.getenv("SCRNA_DIR", "."), "plotting/plotting_functions.R")
plotdir <- file.path(outdir, "plots")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

have_plotfns <- file.exists(pf)
if (have_plotfns) source(pf) else message("NOTE: plotting_functions.R not found; using plain ggplot themes.")
script_id <- tryCatch({ ca <- commandArgs(FALSE); f <- sub("^--file=", "", ca[grep("^--file=", ca)]); if (length(f)) normalizePath(f[1], mustWork = FALSE) else NA_character_ }, error = function(e) NA_character_)

cat_cols <- c(true_like = "#009E73", fp_like = "#D55E00", ambiguous = "grey70")

read_gz <- function(f) tryCatch(fread(cmd = paste("gzip -cd", shQuote(f)), sep = "\t", header = TRUE),
                                error = function(e) fread(f, sep = "\t", header = TRUE))

# ── Discovery + consensus strand ────────────────────────────────────────────
disc <- read_gz(discovery_f)
stopifnot(all(c("combo", "arm", "site", "chr", "pos") %in% names(disc)))
if (!"strand" %in% names(disc)) disc[, strand := NA_character_]
sv <- disc[strand %in% c("+", "-"), .N, by = .(site, strand)]
consensus <- sv[order(site, -N)][, .(strand = strand[1]), by = site]
site_strand <- setNames(consensus$strand, consensus$site)

all_arms <- unique(c(ctrl_arms, treat_arms))

# ── Per-arm editing levels from independent pileups ─────────────────────────
# Two accepted pileup formats:
#  (A) parse_pileup_stranded.pl (preferred): has columns 'coverage','editing_level'
#      (read-strand-aware, direction from genomic ref base + library type).
#  (B) gt_site_pileup.py (fallback): has 'A','C','G','T','cov' base counts; the
#      editing direction is taken from the site's consensus strand.
arm_edit <- list()
for (arm in all_arms) {
  pf_arm <- file.path(pileup_dir, paste0("pileup_", arm, ".tsv"))
  if (!file.exists(pf_arm)) stop("Missing pileup: ", pf_arm)
  p <- fread(pf_arm, sep = "\t", header = TRUE)
  if (all(c("editing_level", "coverage") %in% names(p)) && all(c("chr", "pos") %in% names(p))) {
    # Perl parser output
    p[, site := paste0(chr, ":", pos)]
    p[, edit_level := suppressWarnings(as.numeric(editing_level))]
    setnames(p, "coverage", "cov")
    arm_edit[[arm]] <- p[, .(site, cov = as.numeric(cov), edit_level)]
  } else {
    # pysam base-count output; strand from site consensus
    stopifnot(all(c("site", "A", "C", "G", "T", "cov") %in% names(p)))
    p[, strand := site_strand[site]]
    p[is.na(strand), strand := "+"]
    p[, edited   := fifelse(strand == "+", G, C)]
    p[, unedited := fifelse(strand == "+", A, T)]
    p[, edit_level := fifelse((edited + unedited) > 0, edited / (edited + unedited), NA_real_)]
    arm_edit[[arm]] <- p[, .(site, cov = as.numeric(cov), edit_level)]
  }
}

# combo -> arm -> discovered site set
combos <- unique(disc$combo)
disc_key <- unique(disc[, .(combo, arm, site)])

persite_all <- list(); summ_all <- list()
for (this_combo in combos) {
  cc <- strsplit(this_combo, "__")[[1]]
  caller <- cc[1]; bam_type <- cc[2]; strand_mode <- cc[3]
  for (treat in treat_arms) {
    if (is.null(arm_edit[[treat]])) next
    for (ctrl in ctrl_arms) {
      if (is.null(arm_edit[[ctrl]])) next
      # reference = sites discovered by this combo in either arm of the contrast
      ref_sites <- unique(disc_key[combo == this_combo & arm %in% c(treat, ctrl), site])
      if (!length(ref_sites)) next
      ce <- arm_edit[[ctrl]][site %in% ref_sites]
      te <- arm_edit[[treat]][site %in% ref_sites]
      setnames(ce, c("cov", "edit_level"), c("ctrl_cov", "ctrl_edit"))
      setnames(te, c("cov", "edit_level"), c("treat_cov", "treat_edit"))
      m <- merge(ce, te, by = "site")
      m <- m[ctrl_cov >= min_cov & treat_cov >= min_cov &
               !is.na(ctrl_edit) & !is.na(treat_edit)]
      if (!nrow(m)) next
      m[, delta := ctrl_edit - treat_edit]
      m[, rel_drop := fifelse(ctrl_edit > 0, 1 - treat_edit / ctrl_edit, NA_real_)]
      m[, category := "ambiguous"]
      m[ctrl_edit >= min_ctrl_edit & delta >= min_abs_drop &
          !is.na(rel_drop) & rel_drop >= min_drop_frac, category := "true_like"]
      m[category != "true_like" & treat_edit >= min_ctrl_edit & delta < min_abs_drop,
        category := "fp_like"]
      m[, `:=`(combo = this_combo, caller = caller, bam_type = bam_type,
               strand_mode = strand_mode, control = ctrl, treatment = treat,
               strand = site_strand[site])]
      persite_all[[length(persite_all) + 1L]] <- m

      n_true <- sum(m$category == "true_like"); n_fp <- sum(m$category == "fp_like")
      n_amb  <- sum(m$category == "ambiguous")
      summ_all[[length(summ_all) + 1L]] <- data.table(
        combo = this_combo, caller, bam_type, strand_mode, control = ctrl, treatment = treat,
        n_reference = length(ref_sites), n_covered_both = nrow(m),
        n_true_like = n_true, n_fp_like = n_fp, n_ambiguous = n_amb,
        frac_true_of_called = ifelse((n_true + n_fp) > 0, n_true / (n_true + n_fp), NA_real_),
        median_delta = median(m$delta, na.rm = TRUE),
        median_edit_ctrl = median(m$ctrl_edit, na.rm = TRUE),
        median_edit_treat = median(m$treat_edit, na.rm = TRUE)
      )
    }
  }
}
if (!length(summ_all)) stop("No evaluable sites (covered >= ", min_cov, " in both arms). ",
                            "Check pileups / arm names.")
persite <- rbindlist(persite_all, fill = TRUE)
summ <- rbindlist(summ_all, fill = TRUE)
setorder(summ, treatment, control, -frac_true_of_called)
fwrite(persite, file.path(outdir, "gt_eval_persite.tsv.gz"), sep = "\t")
fwrite(summ, file.path(outdir, "gt_eval_summary.tsv"), sep = "\t")
cat("\n=== GT evaluation summary (", dataset_label, ") ===\n", sep = "")
print(summ[, .(combo, control, treatment, n_covered_both, n_true_like, n_fp_like,
               n_ambiguous, frac_true_of_called)])

# ── Plots ────────────────────────────────────────────────────────────────────
base_theme <- if (have_plotfns) theme_manuscript() else theme_bw(base_size = 11)
gy <- if (have_plotfns) grid_y else theme()
save_fig <- function(p, path, w, h, data, srcs) {
  if (have_plotfns) {
    save_panel(p, path, width = w, height = h, sources = srcs, data = data,
               grid = "y", script = script_id)
  } else {
    ggsave(sub("\\.svg$", ".svg", path), p, width = w, height = h, bg = "white")
  }
}
srcs <- c(discovery_f, list.files(pileup_dir, "^pileup_.*\\.tsv$", full.names = TRUE))

persite[, category := factor(category, levels = names(cat_cols))]
summ_long <- melt(summ, id.vars = c("combo", "control", "treatment"),
                  measure.vars = c("n_true_like", "n_fp_like", "n_ambiguous"),
                  variable.name = "category", value.name = "n")
summ_long[, category := sub("^n_", "", category)]
summ_long[, category := factor(category, levels = names(cat_cols))]

for (treat in unique(persite$treatment)) {
  for (ctrl in unique(persite[treatment == treat, control])) {
    sub <- persite[treatment == treat & control == ctrl]
    if (!nrow(sub)) next
    tag <- paste0(treat, "_vs_", ctrl)

    p_sc <- ggplot(sub, aes(ctrl_edit, treat_edit, colour = category)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(size = 0.5, alpha = 0.4) +
      scale_colour_manual(values = cat_cols, drop = FALSE, name = NULL) +
      facet_wrap(~ combo) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = paste0("editing level — ", ctrl),
           y = paste0("editing level — ", treat),
           title = paste0(dataset_label, ": site-matched editing  (", tag, ")"),
           subtitle = "each point = discovered A>G site covered >=cov in both; below diagonal = decreased in KO/KD") +
      base_theme +
      theme(legend.position = "bottom")
    save_fig(p_sc, file.path(plotdir, paste0("gt_eval_scatter_", tag, ".svg")),
             w = 11, h = 9, data = sub[, .(combo, site, ctrl_edit, treat_edit, category)], srcs = srcs)

    sl <- summ_long[treatment == treat & control == ctrl]
    p_b <- ggplot(sl, aes(reorder(combo, -n), n, fill = category)) +
      geom_col(colour = "grey20", linewidth = 0.2) +
      scale_fill_manual(values = cat_cols, drop = FALSE, name = NULL) +
      labs(x = NULL, y = "discovered sites (covered in both)",
           title = paste0(dataset_label, ": discovered-site classification  (", tag, ")")) +
      base_theme + gy +
      theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")
    save_fig(p_b, file.path(plotdir, paste0("gt_eval_bench_", tag, ".svg")),
             w = 10, h = 6, data = sl, srcs = srcs)
  }
}

# Overall fraction-true benchmark across combos x contrasts
p_f <- ggplot(summ, aes(reorder(combo, frac_true_of_called), frac_true_of_called,
                        fill = treatment)) +
  geom_col(position = position_dodge(width = 0.8), colour = "grey20", linewidth = 0.2) +
  facet_wrap(~ paste0(treatment, " vs ", control), scales = "free_x") +
  labs(x = NULL, y = "fraction ADAR-dependent  [true / (true+FP)]",
       title = paste0(dataset_label, ": editing-dependent fraction of discovered sites")) +
  coord_flip(ylim = c(0, 1)) +
  base_theme +
  theme(legend.position = "bottom")
save_fig(p_f, file.path(plotdir, "gt_eval_fraction_true.svg"),
         w = 11, h = 8, data = summ, srcs = srcs)

message("Wrote ", file.path(outdir, "gt_eval_summary.tsv"))
message("Wrote ", file.path(outdir, "gt_eval_persite.tsv.gz"))
message("Wrote plots under ", plotdir)
