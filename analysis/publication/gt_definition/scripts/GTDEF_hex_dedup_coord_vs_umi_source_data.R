#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_hex_dedup_coord_vs_umi_source_data.R   (2026-10-25; the hex section of GTDEF_bam_processing_plot.R)
#
# Three de-duplication methods (2026-11-16, v9; seed-8787 umi_tools BAMs): coordinate-based (cell-aware
# MarkDuplicates), umi_tools and UMI-aware MarkDuplicates, compared in pairs. A-to-G editing level and KO/NTC
# depletion ratio of each candidate site, one point per site and tier (mean over the tier's groups: merged arm /
# teratomas / cell types; UMI-aware MarkDuplicates cell types: NTC only).
#   el_points / el_rho    editing level, x vs y method, per pair and arm; pooled Spearman's rho
#   dep_points / dep_rho  KO/NTC editing ratio, x vs y method, per pair and KO arm; Spearman's rho and agreement
#                         on 50% depletion
# Site universe: v6 candidate floor on each of the three merged NTC BAMs, union, minus dbSNP common.
#
# INPUT : human/pm/gt_eval_v9/dedup_compare/editlevel/{editlevel_pooled_points.tsv.gz,
#         editlevel_correlation_pooled.tsv, depletion_ratio_points.tsv.gz, depletion_concordance.tsv}
#         (human/pm/scripts/run_gt_editlevel_by_bamtype_v9.sh)
# OUTPUT: source_data/pm_gt/GTDEF_hex_dedup_coord_vs_umi.{el_points,dep_points}.tsv.gz, .{el_rho,dep_rho}.tsv
#         (drawn by GTDEF_hex_dedup_coord_vs_umi_plot.R: response figure 11);
#         report/pm_gt/GTDEF_editlevel_{correlation_pooled,correlation_by_group,depletion_concordance}.tsv
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_hex_dedup_coord_vs_umi"

EL_DIR <- file.path(REPO_DIR, "human", "pm", "gt_eval_v9", "dedup_compare", "editlevel")
el_files <- file.path(EL_DIR, c("editlevel_pooled_points.tsv.gz", "editlevel_correlation_pooled.tsv",
                                "depletion_ratio_points.tsv.gz", "depletion_concordance.tsv"))
stopifnot(all(file.exists(el_files)))
file.copy(el_files[2], file.path(REPORT_DIR, "GTDEF_editlevel_correlation_pooled.tsv"), overwrite = TRUE)
file.copy(el_files[4], file.path(REPORT_DIR, "GTDEF_editlevel_depletion_concordance.tsv"), overwrite = TRUE)
file.copy(file.path(EL_DIR, "editlevel_correlation_by_group.tsv"),
          file.path(REPORT_DIR, "GTDEF_editlevel_correlation_by_group.tsv"), overwrite = TRUE)

el_tier <- function(x) factor(c(arm = "Arm", sample = "Sample", celltype = "Cell type")[x],
                              levels = c("Arm", "Sample", "Cell type"))
arm_title <- function(x) factor(c(NTC_ctrl = "NTC", ADAR1_KO = "ADAR1 KO",
                                  ADARB1_KO = "ADARB1 KO")[as.character(x)],
                                levels = c("NTC", "ADAR1 KO", "ADARB1 KO"))
PAIR_LAB <- c(coord_dedup_vs_umi_dedup = "x: coordinate\ny: umi_tools",
              coord_dedup_vs_umimd_dedup = "x: coordinate\ny: UMI-aware\nMarkDuplicates",
              umi_dedup_vs_umimd_dedup = "x: umi_tools\ny: UMI-aware\nMarkDuplicates")
el_pair <- function(x) factor(PAIR_LAB[x], levels = PAIR_LAB)

pts  <- fread(el_files[1], showProgress = FALSE)
pool <- fread(el_files[2])
pts[, `:=`(pair = el_pair(pair), tier = el_tier(tier), arm = arm_title(arm))]
rho <- pool[, .(pair = el_pair(pair), tier = el_tier(tier), arm = arm_title(arm),
                lab = sprintf("\u03C1 = %.3f\nn = %s", spearman, format(n_sites, big.mark = ",")))]
grid <- CJ(pair = levels(rho$pair), tier = levels(rho$tier), arm = levels(rho$arm))
grid <- grid[!rho, on = .(pair, tier, arm)][, lab := "Not run\n(NTC only)"]
rho <- rbind(rho, grid[, .(pair = factor(pair, levels(rho$pair)), tier = factor(tier, levels(rho$tier)),
                           arm = factor(arm, levels(rho$arm)), lab)])
write_sd(pts, SRC_DIR, STEM, "el_points", gzip = TRUE)
write_sd(rho, SRC_DIR, STEM, "el_rho")

dep  <- fread(el_files[3], showProgress = FALSE)
conc <- fread(el_files[4])
dep[, `:=`(pair = el_pair(pair), tier = el_tier(tier), ko_arm = arm_title(ko_arm))]
rho_d <- conc[, .(pair = el_pair(pair), tier = el_tier(tier), ko_arm = arm_title(ko_arm),
                  lab = sprintf("\u03C1 = %.3f\nagree (50%%) = %.3f", spearman_ratio, agree_50))]
grid <- CJ(pair = levels(rho_d$pair), tier = levels(rho_d$tier), ko_arm = unique(as.character(rho_d$ko_arm)))
grid <- grid[!rho_d, on = .(pair, tier, ko_arm)][, lab := "Not run\n(NTC only)"]
rho_d <- rbind(rho_d, grid[, .(pair = factor(pair, levels(rho_d$pair)), tier = factor(tier, levels(rho_d$tier)),
                               ko_arm = factor(ko_arm, levels(rho_d$ko_arm)), lab)])
write_sd(dep, SRC_DIR, STEM, "dep_points", gzip = TRUE)
write_sd(rho_d, SRC_DIR, STEM, "dep_rho")
