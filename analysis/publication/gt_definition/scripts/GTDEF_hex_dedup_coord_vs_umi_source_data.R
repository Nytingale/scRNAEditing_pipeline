#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_hex_dedup_coord_vs_umi_source_data.R   (2026-10-25; the hex section of GTDEF_bam_processing_plot.R)
#
# Coordinate vs UMI de-duplication (v6): A-to-G editing level and KO/NTC depletion ratio of each candidate site,
# one point per site and tier (mean over the tier's groups: merged arm / teratomas / cell types).
#   el_points / el_rho    editing level, coordinate vs UMI BAM, per arm; pooled Spearman's rho
#   dep_points / dep_rho  KO/NTC editing ratio, coordinate vs UMI BAM, per KO arm; Spearman's rho and agreement
#                         on 50% depletion
# Site universe: v6 candidate floor on the UMI and coordinate merged NTC, union, minus dbSNP common.
#
# INPUT : human/pm/gt_eval_v6/dedup_compare/editlevel/{editlevel_pooled_points.tsv.gz,
#         editlevel_correlation_pooled.tsv, depletion_ratio_points.tsv.gz, depletion_concordance.tsv}
#         (human/pm/scripts/run_gt_editlevel_by_bamtype_v6.sh)
# OUTPUT: source_data/pm_gt/GTDEF_hex_dedup_coord_vs_umi.{el_points,dep_points}.tsv.gz, .{el_rho,dep_rho}.tsv
#         (drawn by GTDEF_hex_dedup_coord_vs_umi_plot.R: response figure 11);
#         report/pm_gt/GTDEF_editlevel_{correlation_pooled,correlation_by_group,depletion_concordance}.tsv
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_hex_dedup_coord_vs_umi"

EL_DIR <- file.path(GT_V6_DIR, "dedup_compare", "editlevel")
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

pts  <- fread(el_files[1], showProgress = FALSE)
pool <- fread(el_files[2])
pts[, `:=`(tier = el_tier(tier), arm = arm_title(arm))]
rho <- pool[, .(tier = el_tier(tier), arm = arm_title(arm),
                lab = sprintf("\u03C1 = %.3f\nn = %s", spearman, format(n_sites, big.mark = ",")))]
write_sd(pts, SRC_DIR, STEM, "el_points", gzip = TRUE)
write_sd(rho, SRC_DIR, STEM, "el_rho")

dep  <- fread(el_files[3], showProgress = FALSE)
conc <- fread(el_files[4])
dep[, `:=`(tier = el_tier(tier), ko_arm = arm_title(ko_arm))]
rho_d <- conc[, .(tier = el_tier(tier), ko_arm = arm_title(ko_arm),
                  lab = sprintf("\u03C1 = %.3f\nagree (50%%) = %.3f", spearman_ratio, agree_50))]
write_sd(dep, SRC_DIR, STEM, "dep_points", gzip = TRUE)
write_sd(rho_d, SRC_DIR, STEM, "dep_rho")
