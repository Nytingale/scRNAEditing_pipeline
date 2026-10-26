#!/usr/bin/env Rscript
# =============================================================================
# TERATOMA_dedup_audit_plot.R
#
# PI concern 2 — is dedup needed at TRUE editing sites?
#
# ONE figure (cowplot / manuscript theme; not patchwork): per CRISPR arm, mean
# ± SD across T1–T4 of samtools depth at the frozen v3 TRUE BED, for
# Unprocessed / De-duplicated (coordinate) / De-duplicated (UMI).
#
# Depth table (from SGE depth_truth array + collect):
#   report/pm_gt/DEDUP_persample_truth_depth.tsv
#
# Supplementary tables retained from audit_dedup_editing_relevance.R (global
# MarkDuplicates % + F1 by bam_type). Do NOT mix idxstats library size with
# site depth on the same axes.
#
# Usage: Rscript TERATOMA_dedup_audit_plot.R [DATASET]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "TERATOMA_common.R"))

argv <- commandArgs(trailingOnly = TRUE)
teratoma_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT <- caller_script()

DEPTH <- file.path(REPORT_DIR, "DEDUP_persample_truth_depth.tsv")
GLOB  <- file.path(REPORT_DIR, "DEDUP_global_rates.tsv")
F1DEL <- file.path(REPORT_DIR, "DEDUP_f1_delta_by_bamtype.tsv")
if (!file.exists(DEPTH))
  stop("missing ", DEPTH, " — wait for depth_truth collect (DEDUP_persample_truth_depth.tsv)")

# Manuscript labels for the three BAMs measured at TRUE sites (depth-only).
# Colours from cmap_proc (first 3 entries: Unprocessed, coord, UMI).
label_depth_bam <- c(
  unprocessed = "Unprocessed",
  coord_dedup = "De-duplicated (coordinate)",
  umi_dedup   = "De-duplicated (UMI)"
)
cmap_depth_bam <- cmap_proc[1:3]
arm_levels <- c("NTC_ctrl", "AAVS1_ctrl", "ADAR1_KO", "ADARB1_KO", "ADARB2_KO")

d <- fread(DEPTH)
need <- c("teratoma", "arm", "bam_type", "mean_depth", "sum_depth", "n_sites")
miss <- setdiff(need, names(d))
if (length(miss)) stop("DEDUP_persample_truth_depth.tsv missing columns: ",
                       paste(miss, collapse = ", "))
if (!all(d$bam_type %in% names(label_depth_bam)))
  stop("unexpected bam_type values: ", paste(unique(d$bam_type), collapse = ", "))

d[, layer := factor(relabel(bam_type, label_depth_bam),
                    levels = names(cmap_depth_bam))]
d[, arm := factor(arm, levels = intersect(arm_levels, unique(arm)))]

# Mean ± SD across the four teratomas (biological axis for this figure).
agg <- d[, .(
  n_teratoma = uniqueN(teratoma),
  n_sites    = as.integer(median(n_sites)),
  mean_depth = mean(mean_depth),
  sd_depth   = sd(mean_depth),
  mean_sum   = mean(sum_depth),
  sd_sum     = sd(sum_depth)
), by = .(arm, layer)]
if (any(agg$n_teratoma != 4L))
  warning("expected 4 teratomas per arm x layer; got: ",
          paste(unique(agg$n_teratoma), collapse = ","))

# ── Supp: per-arm depth summary (depth-consistent; replaces draft site table) ─
supp_depth <- agg[, .(
  arm,
  bam_processing = as.character(layer),
  n_true_sites = n_sites,
  n_teratoma,
  mean_of_mean_depth = round(mean_depth, 4),
  sd_of_mean_depth   = round(sd_depth, 4),
  mean_of_sum_depth  = round(mean_sum, 1),
  sd_of_sum_depth    = round(sd_sum, 1)
)]
fwrite(supp_depth, file.path(REPORT_DIR, "TERATOMA_supp_dedup_true_site_depth.tsv"),
       sep = "\t")
message("  table       : report/", DATASET, "/TERATOMA_supp_dedup_true_site_depth.tsv")

# Keep global rates + F1 tables if the earlier audit already wrote them.
if (file.exists(GLOB)) {
  g <- fread(GLOB)
  label_dataset <- c(teratoma_merged = "Teratoma", replogle = "Replogle")
  supp_reads <- g[, .(
    dataset = relabel(dataset, label_dataset),
    arm,
    n_input_reads = as.integer(unpaired_reads),
    n_duplicate_reads = as.integer(unpaired_dups),
    pct_duplication_MarkDuplicates = round(100 * percent_duplication, 2),
    n_after_coordinate = as.integer(n_rmDup),
    n_after_UMI = as.integer(n_umiDup),
    pct_removed_coordinate = round(100 * frac_removed_coord, 2),
    pct_removed_UMI = round(100 * frac_removed_umi, 2)
  )]
  fwrite(supp_reads, file.path(REPORT_DIR, "TERATOMA_supp_dedup_read_counts.tsv"),
         sep = "\t")
  message("  table       : report/", DATASET, "/TERATOMA_supp_dedup_read_counts.tsv")
} else {
  message("  WARNING: ", GLOB, " absent — skipping global-rate supp table")
}

if (file.exists(F1DEL)) {
  f1 <- fread(F1DEL)
  supp_f1 <- copy(f1)
  supp_f1[, `:=`(
    caller       = relabel(caller, label_caller),
    strand_mode  = relabel(strand_mode, label_strandMode),
    snp_subtract = relabel(snp_subtract, label_snpSubtract)
  )]
  setnames(supp_f1,
           c("sample", "sampleRmDup", "sampleUmiDup",
             "delta_sample_minus_RmDup", "delta_sample_minus_Umi"),
           c("F1_Unprocessed", "F1_coordinate_dedup", "F1_UMI_dedup",
             "delta_F1_Unprocessed_minus_coordinate",
             "delta_F1_Unprocessed_minus_UMI"),
           skip_absent = TRUE)
  fwrite(supp_f1, file.path(REPORT_DIR, "TERATOMA_supp_dedup_f1_by_bamtype.tsv"),
         sep = "\t")
  message("  table       : report/", DATASET, "/TERATOMA_supp_dedup_f1_by_bamtype.tsv")
} else {
  message("  WARNING: ", F1DEL, " absent — skipping F1 supp table")
}

# ── Single figure: mean ± SD depth at TRUE sites ────────────────────────────
STEM <- "TERATOMA_bar_dedup_true_site_depth"
ymax <- max(agg$mean_depth + agg$sd_depth, na.rm = TRUE)

p <- ggplot(agg, aes(arm, mean_depth, fill = layer)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "grey20", linewidth = 0.15) +
  geom_errorbar(aes(ymin = pmax(0, mean_depth - sd_depth),
                    ymax = mean_depth + sd_depth),
                position = position_dodge(width = 0.75),
                width = 0.22, linewidth = 0.4, colour = "grey20", na.rm = TRUE) +
  scale_fill_manual_cvd(values = cmap_depth_bam, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Mean depth at TRUE sites",
       title = paste0("Coverage at TRUE editing sites by dedup (mean \u00b1 SD, n = ",
                      unique(agg$n_sites)[1], " sites)")) +
  theme_manuscript() + grid_y +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

save_figure(p, STEM, width = 9.0, height = 4.8,
            sources = DEPTH, data = agg, grid = "y", script = SCRIPT)

message("Done. Single depth figure -> figures/", DATASET, "/", STEM, ".svg")
message("Draft multi-panel SVGs (read_counts / percent / panel_dedup_audit) are obsolete.")
