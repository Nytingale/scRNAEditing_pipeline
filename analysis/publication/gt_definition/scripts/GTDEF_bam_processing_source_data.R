#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_bam_processing_source_data.R   (split from GTDEF_bam_processing_plot.R 2026-10-25)
#
# Ground-truth DEFINITION on BAM processing (v6): mismatch types among sites called in NTC on each BAM type, before
# vs after restricting to the 12-type ground truth, and the overlap of ground-truth sites across BAM types.
#   spectrum12  percentage of sites per mismatch type, per tier x BAM type x set (called in NTC / ground truth)
#   spectrum6   the same with strand complements collapsed
#   euler       ground-truth sites in coordinate-only, UMI-only and both de-duplicated BAMs, per tier
#   upset_combos  ground-truth sites per combination of the reported BAM types (presence bits in BAM type order)
#   upset_sets    ground-truth sites per BAM type
#
# INPUT (human/pm/gt_eval_v6/):
#   dedup_compare/discovery/{arm,sample,celltype}/gt_eval_discovery.tsv.gz
#       all-varType pre100 discovery, NTC units only (run_dedup_compare_v6.sh)
#   spectrum12/<gate>/truth12_sites.tsv.gz
#       12-type truth (RNA strand); ground truth = d_max >= TRUTH_D
# OUTPUT: source_data/pm_gt/GTDEF_bam_processing.<part>.tsv (drawn by GTDEF_bam_processing_plot.R) and the report
#         tables report/pm_gt/GTDEF_{mismatch_spectrum,mismatch_spectrum6}_by_bamtype_vs_TRUE.tsv,
#         GTDEF_TRUE_sites_by_bamtype.tsv, GTDEF_euler_TRUE_sites_by_bamtype_counts.tsv
# ============================================================================
suppressPackageStartupMessages(library(data.table))
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_bam_processing"

DISC_DIR   <- file.path(GT_V6_DIR, "dedup_compare", "discovery")
TRUTH_GATE <- "eval"   # expected-edited-reads gate (m = 3, frozen); "eval_depth" = depth gate
TRUTH_D    <- 0.5
TRUTH12    <- file.path(GT_V6_DIR, "spectrum12", TRUTH_GATE, "truth12_sites.tsv.gz")
# tier token -> discovery subdir
TIERS <- c(dataset = "arm", sample = "sample", celltype = "celltype")

# Strand-complement collapse of the 12 mismatch types -> 6 pairs.
varType6_map <- c(
  AC = "AC", TG = "AC",
  AG = "AG", TC = "AG",
  AT = "AT", TA = "AT",
  CA = "CA", GT = "CA",
  CG = "CG", GC = "CG",
  CT = "CT", GA = "CT"
)

# Unprocessed BAMs (bam_type `sample`) are kept in the report tables but not plotted.
BT_ALL      <- c("sample", "sampleRmDup", "sampleUmiDup",
                 "recal", "recalq20", "recaldefault")
BT_REPORTED <- setdiff(BT_ALL, "sample")
BT_DEDUP    <- c("sampleRmDup", "sampleUmiDup")

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
truth <- unique(truth[d_max >= TRUTH_D, .(site = paste0(chr, ":", pos), varType = type)])
message(sprintf("ground truth (%s, d >= %.2f): %d sites, %d A-to-G",
                TRUTH_GATE, TRUTH_D, nrow(truth), truth[varType == "AG", .N]))

# Per-bam discovered vs ground-truth spectra (12-type + collapsed 6-type) ----
# Before = unique sites called in NTC on that bam_type (majority varType)
# After  = those sites whose (site, varType) is in the ground truth
all12 <- all_vartypes()
spec_rows <- list()
true_by_bam <- list()
for (tier_i in names(TIERS)) {
  disc_f <- file.path(DISC_DIR, TIERS[[tier_i]], "gt_eval_discovery.tsv.gz")
  if (!file.exists(disc_f)) stop("missing ", disc_f)
  d <- fread(cmd = paste("gzip -cd", shQuote(disc_f)),
             select = c("caller", "bam_type", "strand", "site", "varType"),
             showProgress = FALSE)
  # REDItools3 reports the genomic-strand type on - strand rows; other callers the RNA strand
  d[caller == "reditools3" & strand == "-", varType := chartr("ACGT", "TGCA", varType)]
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
    spec_d[, `:=`(tier = tier_i, bam_type = bt, set = "all_discovered", n_sites = n_disc)]
    spec_rows[[length(spec_rows) + 1L]] <- spec_d
    vt_t <- vt[truth, on = .(site, varType), nomatch = NULL]
    n_true <- nrow(vt_t)
    true_by_bam[[tier_i]][[bt]] <- vt_t$site
    if (n_true < 1L) {
      spec_t <- data.table(varType = all12, n = 0L, tier = tier_i, bam_type = bt, set = "pass_GT_criteria",
                           n_sites = 0L)
    } else {
      spec_t <- vt_t[, .(n = .N), by = varType]
      spec_t[, `:=`(tier = tier_i, bam_type = bt, set = "pass_GT_criteria", n_sites = n_true)]
    }
    spec_rows[[length(spec_rows) + 1L]] <- spec_t
    message(sprintf("  %-8s %-14s disc=%7d AG=%.1f%% | TRUE∩=%5d AG=%.1f%%",
                    tier_i, bt, n_disc, 100 * spec_d[varType == "AG", n] / n_disc, n_true,
                    if (n_true) 100 * spec_t[varType == "AG", sum(n)] / n_true else 0))
  }
  rm(d); gc(verbose = FALSE)
}

both_bam <- rbindlist(spec_rows, fill = TRUE)
both_bam <- fill_spectrum(both_bam, all12, c("tier", "bam_type", "set"))
both_bam[, reported := bam_type %in% BT_REPORTED]
fwrite(both_bam, file.path(REPORT_DIR, "GTDEF_mismatch_spectrum_by_bamtype_vs_TRUE.tsv"), sep = "\t")
write_sd(both_bam[bam_type %in% BT_REPORTED], SRC_DIR, STEM, "spectrum12")

both6 <- copy(both_bam)
both6[, varType6 := varType6_map[as.character(varType)]]
both6 <- both6[!is.na(varType6), .(n = sum(as.numeric(n)), n_sites = max(n_sites)),
               by = .(tier, bam_type, set, varType6)]
both6[, percent := fifelse(n_sites > 0L, n / n_sites, 0)]
setnames(both6, "varType6", "varType")
both6[, reported := bam_type %in% BT_REPORTED]
fwrite(both6, file.path(REPORT_DIR, "GTDEF_mismatch_spectrum6_by_bamtype_vs_TRUE.tsv"), sep = "\t")
write_sd(both6[bam_type %in% BT_REPORTED], SRC_DIR, STEM, "spectrum6")

# Ground-truth site overlap across BAM processing methods ----
memb <- rbindlist(lapply(names(true_by_bam), function(tier_i)
  rbindlist(lapply(names(true_by_bam[[tier_i]]), function(bt) {
    ss <- true_by_bam[[tier_i]][[bt]]
    if (length(ss)) data.table(tier = tier_i, bam_type = bt, site = ss)
  }))))
memb[, reported := bam_type %in% BT_REPORTED]
fwrite(memb, file.path(REPORT_DIR, "GTDEF_TRUE_sites_by_bamtype.tsv"), sep = "\t")

euler <- rbindlist(lapply(names(TIERS), function(tier_i) {
  sets <- true_by_bam[[tier_i]]
  if (is.null(sets) || !all(BT_DEDUP %in% names(sets))) return(NULL)
  R <- sets[["sampleRmDup"]]
  U <- sets[["sampleUmiDup"]]
  data.table(tier = tier_i, region = c("Coord", "UMI", "Coord&UMI"),
             n_sites = c(length(setdiff(R, U)), length(setdiff(U, R)), length(intersect(R, U))))
}))
fwrite(euler, file.path(REPORT_DIR, "GTDEF_euler_TRUE_sites_by_bamtype_counts.tsv"), sep = "\t")
write_sd(euler, SRC_DIR, STEM, "euler")

# UpSet: sites per combination of the reported BAM types present in the tier (bits in BT_REPORTED order), counted as
# table() does (combinations in key order), and set sizes.
upset_combos <- list()
upset_sets <- list()
for (tier_i in names(TIERS)) {
  mt <- memb[tier == tier_i & bam_type %in% BT_REPORTED]
  if (!nrow(mt)) next
  bt_lev <- intersect(BT_REPORTED, unique(mt$bam_type))
  wide <- dcast(unique(mt[, .(site, bam_type)]), site ~ bam_type, fun.aggregate = length)
  for (bt in bt_lev) if (!bt %in% names(wide)) wide[, (bt) := 0L]
  mat <- (as.matrix(wide[, ..bt_lev]) > 0L) * 1L
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  comb_n <- as.data.table(table(apply(mat, 1L, paste0, collapse = "")))
  setnames(comb_n, c("comb", "n"))
  upset_combos[[tier_i]] <- comb_n[, .(tier = tier_i, bam_types = paste(bt_lev, collapse = ","), comb, n)]
  upset_sets[[tier_i]] <- data.table(tier = tier_i, bam_type = bt_lev, n = as.integer(colSums(mat)[bt_lev]))
}
write_sd(rbindlist(upset_combos), SRC_DIR, STEM, "upset_combos")
write_sd(rbindlist(upset_sets), SRC_DIR, STEM, "upset_sets")
message("done")
