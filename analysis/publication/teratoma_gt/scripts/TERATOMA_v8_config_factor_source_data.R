#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_config_factor_source_data.R   (2026-10-15; split from the _plot.R script 2026-10-25;
#   v9 truths, 76 configurations, caller column and held-out aorta 2026-11-14)
# Source data for response figure 13:
# Effect of each configuration choice on F1 (review: R1.3 GATK default / strand modes / de-duplication, R2.5
# stranding, R3 per-sample SNP blacklist, BQSR): paired F1 difference between two configurations that differ in one
# choice only, per tier and truth set: teratoma v9 KO truths (6) and REDIportal truths (3), held-out aorta v9 truths (6);
# main call state (mismatch).
#   Caller: each caller - BCFtools, at every BAM / stranding / blacklist setting the pair shares (e.g. GATK and BCFtools
#     with the sample blacklist on)
#   Strand: hybrid - annotation (BCFtools, GATK, Cellsnp-lite); caller's own - annotation (JACUSA2, REDItools3, scAllele)
#   De-duplication: umi_tools - coordinate; UMI-aware MarkDuplicates - coordinate; UMI-aware MarkDuplicates - umi_tools
#   SNP removal: Database + sample blacklist - Database (BCFtools, GATK)
#   GATK on the recalibrated BAM: our Q20 run (mapping and base quality 20) - HaplotypeCaller defaults (mapping
#     quality 20, base quality 10), i.e. base-quality threshold 20 vs 10 at equal mapping quality. Both keep the default
#     read filters (only the non-recalibrated BAMs run GATK with 8 read filters disabled).
#   BQSR: recalibrated - coordinate-based de-duplicated; recalibrated (Q20) - recalibrated (Q30)
#   One column per contrast (nested under its choice), truth sets on the y-axis (row strip: truth family and dataset);
#   points pool the tiers.
# Inputs : human/pm/gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv (teratoma KO truths, f1_<tier>)
#          human/pm/gt_eval_v9/rediportal_v9/scoring_mismatch76/gt_tier_scoring_<set>.tsv
#          mouse/aorta_gse255661/gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv (aorta, arm and sample tiers)
# Output : source_data/pm_gt/TERATOMA_v8_config_factor_delta_f1.{pairs,summary}.tsv
#          (drawn by TERATOMA_v8_config_factor_plot.R)
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
STEM <- "TERATOMA_v8_config_factor_delta_f1"

TIER_LAB <- c(arm = "Arm", sample = "Sample", celltype = "Cell type")
KO_IDS <- c("expected_t25", "expected_t50", "expected_t75", "depth_t25", "depth_t50", "depth_t75")
REDI_IDS <- c("redi_cov10_alt1", "redi_alt3", "redi_alt3_e05")
KO_LAB <- setNames(rep(sprintf("\u2265 %d%% reduction", c(25, 50, 75)), 2), KO_IDS)
TRUTH_LAB <- c(KO_LAB, setNames(c("\u2265 1 variant read", "\u2265 3 variant reads",
                                  "\u2265 3 variant reads,\n\u2265 5% editing level"), REDI_IDS))
FAMILY_LAB <- c(ko_expected = "Power-gated\n(Human teratoma)", ko_depth = "Depth-gated\n(Human teratoma)",
                redi = "REDIportal-based\n(Human teratoma)",
                aorta_expected = "Power-gated\n(Mouse aorta)", aorta_depth = "Depth-gated\n(Mouse aorta)")
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
DEDUP <- c("sampleRmDup", "sampleUmiDup", "sampleUmiMd")
N_CONFIG <- 76L
KO_SRC <- file.path(REPO_DIR, "human/pm/gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv")
AO_SRC <- file.path(REPO_DIR, "mouse/aorta_gse255661/gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv")
REDI_SRC <- setNames(file.path(REPO_DIR, "human/pm/gt_eval_v9/rediportal_v9/scoring_mismatch76",
                               paste0("gt_tier_scoring_", sub("^redi_", "", REDI_IDS), ".tsv")), REDI_IDS)
stopifnot(file.exists(KO_SRC), file.exists(AO_SRC), all(file.exists(REDI_SRC)))

rank_long <- function(f, fam) {
  x <- fread(f)
  tiers <- intersect(names(TIER_LAB), sub("^f1_", "", grep("^f1_", names(x), value = TRUE)))
  x <- melt(x[, c("gate", "t", COMBO, paste0("f1_", tiers)), with = FALSE], measure.vars = paste0("f1_", tiers),
            variable.name = "tier", value.name = "f1")
  x[, `:=`(tier = sub("^f1_", "", tier), truth = paste(gate, t, sep = "_"), fam = fam)][, c("gate", "t") := NULL]
}
redi <- rbindlist(lapply(names(REDI_SRC), function(n) fread(REDI_SRC[[n]])[
  , .(tier = fifelse(tier == "dataset", "arm", tier), caller, bam_type, strand_mode, snp_subtract, f1, truth = n, fam = "redi")]))
f1 <- rbind(rank_long(KO_SRC, "ko"), redi, rank_long(AO_SRC, "aorta"), use.names = TRUE)
f1[, ds_truth := paste(fam, truth)]
stopifnot(f1[, uniqueN(do.call(paste, .SD)), .SDcols = COMBO, by = fam]$V1 == N_CONFIG,
          f1[, .N, by = .(ds_truth, tier)]$N == N_CONFIG)

# one row per contrast: the configuration pair differs in exactly one column
contrast <- function(label, col, a, b, callers = NULL) {
  x <- if (is.null(callers)) f1 else f1[caller %in% callers]
  keys <- setdiff(c(COMBO, "tier", "truth", "fam", "ds_truth"), col)
  w <- merge(x[get(col) == a, c(keys, "f1"), with = FALSE], x[get(col) == b, c(keys, "f1"), with = FALSE],
             by = keys, suffixes = c("_a", "_b"))
  w[, `:=`(contrast = label, delta = f1_a - f1_b)]
}
# caller: every other caller against BCFtools wherever both have the same BAM, stranding and blacklist setting;
# the point keeps the other caller's colour
cal <- merge(f1[caller != "bcftools"], f1[caller == "bcftools", .(bam_type, strand_mode, snp_subtract, tier, ds_truth, f1_b = f1)],
             by = c("bam_type", "strand_mode", "snp_subtract", "tier", "ds_truth"))
cal[, `:=`(f1_a = f1, contrast = "Caller: Each caller - BCFtools", delta = f1 - f1_b)][, f1 := NULL]
print(cal[, .N, by = .(fam, caller, bam_type, strand_mode, snp_subtract)][, .(n_settings = .N), by = .(fam, caller)])
d <- rbind(
  cal,
  contrast("Strand: Hybrid - Annotation", "strand_mode", "hybrid", "annotation"),
  contrast("Strand: Caller's own - Annotation", "strand_mode", "caller", "annotation"),
  contrast("De-duplication: umi_tools - Coordinate-based", "bam_type", "sampleUmiDup", "sampleRmDup"),
  contrast("De-duplication: UMI-aware MarkDuplicates - Coordinate-based", "bam_type", "sampleUmiMd", "sampleRmDup"),
  contrast("De-duplication: UMI-aware MarkDuplicates - umi_tools", "bam_type", "sampleUmiMd", "sampleUmiDup"),
  contrast("SNP removal: Database + sample blacklist - Database", "snp_subtract", "on", "off"),
  contrast("GATK base quality: Recalibrated (Q20) - Default (Q10)", "bam_type", "recalq20", "recaldefault", callers = "gatk"),
  contrast("Recalibration: Recalibrated - Coordinate-based", "bam_type", "recal", "sampleRmDup"),
  contrast("Recalibration: Recalibrated (Q20) - Recalibrated", "bam_type", "recalq20", "recal"), use.names = TRUE, fill = TRUE)
d <- d[!is.na(delta)]
# column label: "<first> -\n<second>", each side wrapped at 15 characters (SNP removal: unwrapped)
detail_lab <- function(x) vapply(x, function(s) {
  p <- strsplit(sub("^[^:]*: ", "", s), " - ", fixed = TRUE)[[1]]
  w <- if (startsWith(s, "SNP removal")) 100 else 15
  paste0(paste(strwrap(p[1], w), collapse = "\n"), " \u2212\n", paste(strwrap(p[2], w), collapse = "\n"))
}, "", USE.NAMES = FALSE)
d[, `:=`(Tier = factor(TIER_LAB[tier], TIER_LAB), method = caller_factor(caller),
         family = fam,
         Family = factor(FAMILY_LAB[fifelse(fam == "redi", "redi", paste0(fam, "_", sub("_.*$", "", truth)))], FAMILY_LAB),
         Truth = factor(TRUTH_LAB[truth], rev(unique(TRUTH_LAB))),
         Choice = factor(sub(":.*$", "", contrast), unique(sub(":.*$", "", contrast))),
         Detail = factor(detail_lab(contrast), unique(detail_lab(contrast))))]
d[, c("fam", "ds_truth") := NULL]
summ <- d[, .(n = .N, median_delta = round(median(delta), 3), q25 = round(quantile(delta, 0.25), 3),
              q75 = round(quantile(delta, 0.75), 3), frac_positive = round(mean(delta > 0), 2)),
          by = .(contrast, tier, family)][order(contrast, family, tier)]
print(summ, nrows = 300)
print(d[, .(median_delta = round(median(delta), 3), frac_positive = round(mean(delta > 0), 2)), by = .(contrast, family, caller)][
  order(contrast, family, caller)], nrows = 300)

write_sd(d, SRC_DIR, STEM, "pairs")
write_sd(summ, SRC_DIR, STEM, "summary")
