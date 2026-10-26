#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_config_factor_source_data.R   (2026-10-15; split from the _plot.R script 2026-10-25)
# Source data for response figure 10:
# Effect of each configuration choice on F1 (review: R1.3 GATK default / strand modes, R2.5 stranding,
# R3 per-sample SNP blacklist, BQSR): paired F1 difference between two configurations that differ in one
# choice only, per tier and v8 truth set (6 KO + 3 REDIportal), v8 main call state (mismatch).
#   Strand: hybrid - annotation (BCFtools, GATK); caller's own - annotation (JACUSA2, REDItools3, scAllele)
#   De-duplication: UMI - coordinate
#   SNP removal: dbSNP + sample blacklist - dbSNP only (BCFtools, GATK)
#   GATK on the recalibrated BAM: HaplotypeCaller defaults (mapping quality 20, base quality 10) - our Q20 run (mapping
#     and base quality 20), i.e. base-quality threshold 10 vs 20 at equal mapping quality. Both keep the default read
#     filters (only the non-recalibrated BAMs run GATK with 8 read filters disabled).
#   BQSR: recalibrated - coordinate de-duplicated; recalibrated (Q20) - recalibrated (Q30)
#   One column per contrast (nested under its choice), truth sets on the y-axis; points pool the three tiers.
# Inputs : human/pm/gt_eval_v8/scoring/mismatch/config_rank.tsv (KO truths, f1_<tier>)
#          human/pm/gt_eval_v8/rediportal/scoring_mismatch/gt_tier_scoring_<set>.tsv
# Output : source_data/pm_gt/TERATOMA_v8_config_factor_delta_f1.{pairs,summary}.tsv
#          (drawn by TERATOMA_v8_config_factor_plot.R)
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
STEM <- "TERATOMA_v8_config_factor_delta_f1"

TIER_LAB <- c(arm = "Arm", sample = "Sample", celltype = "Cell type")
TRUTH_STRIP <- setNames(c(sprintf("KO, \u2265 %d%%\npower-gated", c(25, 50, 75)),
                          sprintf("KO, \u2265 %d%%\ndepth-gated", c(25, 50, 75)),
                          "REDIportal\n\u2265 1 variant\nread", "REDIportal\n\u2265 3 variant\nreads",
                          "REDIportal\n\u2265 3 variant reads,\n\u2265 5% editing level"), V8_TRUTH_IDS)
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
KO_SRC <- v8_config_rank(V8_STATE)
REDI_SRC <- v8_redi_files(V8_STATE)
stopifnot(file.exists(KO_SRC), all(file.exists(REDI_SRC)))

ko <- melt(fread(KO_SRC)[, c("gate", "t", COMBO, paste0("f1_", names(TIER_LAB))), with = FALSE],
           measure.vars = paste0("f1_", names(TIER_LAB)), variable.name = "tier", value.name = "f1")
ko[, `:=`(tier = sub("^f1_", "", tier), truth = paste(gate, t, sep = "_"))][, c("gate", "t") := NULL]
redi <- rbindlist(lapply(names(REDI_SRC), function(n) fread(REDI_SRC[[n]])[
  , .(tier = fifelse(tier == "dataset", "arm", tier), caller, bam_type, strand_mode, snp_subtract, f1, truth = n)]))
f1 <- rbind(ko, redi, use.names = TRUE)
stopifnot(setequal(f1$truth, V8_TRUTH_IDS), f1[, uniqueN(do.call(paste, .SD)), .SDcols = COMBO] == 58L)

# one row per contrast: the configuration pair differs in exactly one column
contrast <- function(label, col, a, b, callers = NULL) {
  x <- if (is.null(callers)) f1 else f1[caller %in% callers]
  keys <- setdiff(c(COMBO, "tier", "truth"), col)
  w <- merge(x[get(col) == a, c(keys, "f1"), with = FALSE], x[get(col) == b, c(keys, "f1"), with = FALSE],
             by = keys, suffixes = c("_a", "_b"))
  w[, `:=`(contrast = label, delta = f1_a - f1_b)]
}
d <- rbind(
  contrast("Strand: hybrid - annotation", "strand_mode", "hybrid", "annotation"),
  contrast("Strand: caller's own - annotation", "strand_mode", "caller", "annotation"),
  contrast("De-duplication: UMI - coordinate", "bam_type", "sampleUmiDup", "sampleRmDup"),
  contrast("SNP removal: + sample blacklist - dbSNP only", "snp_subtract", "on", "off"),
  contrast("GATK base quality: Default (BQ 10) - recalibrated (Q20)", "bam_type", "recaldefault", "recalq20", callers = "gatk"),
  contrast("Recalibration: recalibrated - coordinate de-dup.", "bam_type", "recal", "sampleRmDup"),
  contrast("Recalibration: recalibrated (Q20) - recalibrated", "bam_type", "recalq20", "recal"), fill = TRUE)
d <- d[!is.na(delta)]
d[, `:=`(Tier = factor(TIER_LAB[tier], TIER_LAB), method = caller_factor(caller),
         family = fifelse(grepl("^redi_", truth), "REDIportal truths", "KO truths"),
         Truth = factor(gsub("\n", " ", TRUTH_STRIP[truth]), rev(gsub("\n", " ", TRUTH_STRIP))),
         Choice = factor(sub(":.*$", "", contrast), unique(sub(":.*$", "", contrast))),
         Detail = factor(sub(" - ", " \u2212\n", sub("^[^:]*: ", "", contrast)),
                         unique(sub(" - ", " \u2212\n", sub("^[^:]*: ", "", contrast)))))]
d[, Family := factor(fifelse(family == "KO truths", "Editing reduction\nin KO", "REDIportal"),
                     c("Editing reduction\nin KO", "REDIportal"))]
summ <- d[, .(n = .N, median_delta = round(median(delta), 3), q25 = round(quantile(delta, 0.25), 3),
              q75 = round(quantile(delta, 0.75), 3), frac_positive = round(mean(delta > 0), 2)),
          by = .(contrast, tier, family)][order(contrast, family, tier)]
print(summ, nrows = 200)
print(d[, .(median_delta = round(median(delta), 3), frac_positive = round(mean(delta > 0), 2)), by = .(contrast, caller)][
  order(contrast, caller)], nrows = 200)

write_sd(d, SRC_DIR, STEM, "pairs")
write_sd(summ, SRC_DIR, STEM, "summary")
