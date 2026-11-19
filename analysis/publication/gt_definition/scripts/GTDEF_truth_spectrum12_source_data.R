#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_truth_spectrum12_source_data.R
#
# Mismatch spectrum (12 types, read orientation) of the KO-depleted gate sites behind each v9 ground truth,
# before and after the truth cleaning, for both gates and t = 25 / 50 / 75% editing reduction in KO.
#   teratoma  v9 build (gt_v9_build.sh truth: gt_truth_v7_variants.R, GT_V7_VARIANT=AGTC)
#             before = none, after = nomm_pgp1 (removal of sites near other mismatch types + PGP1 SNV positions)
#   aorta     aorta_truth_eval.R, SNP mode excluded (MGP SNP positions removed before anything else)
#             before = none, after = nomm (removal of sites near other mismatch types)
# Truths: the headline (umi_tools, --random-seed 8787) and the read-set variants of the response
# (umi_tools seeds 1 and 2, UMI-aware MarkDuplicates, coordinate MarkDuplicates; human/pm/scripts/umi_seed_reps.sh).
#
# INPUT : human/pm/gt_eval_v9/build/variants/spectrum12_cleaning.tsv
#         human/pm/umi_seed_reps/<seed1|seed2|umimd|coord>/v9/build/variants/spectrum12_cleaning.tsv
#         mouse/aorta_gse255661/truth/eval/excluded/spectrum12_by_t.tsv
#         mouse/aorta_gse255661/umi_seed_reps/<seed1|seed2|umimd|coord>/eval/excluded/spectrum12_by_t.tsv
# OUTPUT: source_data/pm_gt/GTDEF_truth_spectrum12.spectrum.tsv (drawn by GTDEF_truth_spectrum12_plot.R)
#         dataset truth truth_order gate t stage type N pct w1 source
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_truth_spectrum12"

TRUTHS <- data.table(
  truth = c("umi_tools seed 8787", "umi_tools seed 1", "umi_tools seed 2", "UMI-aware MarkDuplicates",
            "Coordinate MarkDuplicates"),
  rep = c(NA, "seed1", "seed2", "umimd", "coord"))
TER <- file.path(REPO_DIR, "human", "pm")
AOR <- file.path(REPO_DIR, "mouse", "aorta_gse255661")
ter_file <- function(r) if (is.na(r)) file.path(TER, "gt_eval_v9", "build", "variants", "spectrum12_cleaning.tsv") else
  file.path(TER, "umi_seed_reps", r, "v9", "build", "variants", "spectrum12_cleaning.tsv")
aor_file <- function(r) if (is.na(r)) file.path(AOR, "truth", "eval", "excluded", "spectrum12_by_t.tsv") else
  file.path(AOR, "umi_seed_reps", r, "eval", "excluded", "spectrum12_by_t.tsv")
STAGES <- list(Teratoma = c(before = "none", after = "nomm_pgp1"), Aorta = c(before = "none", after = "nomm"))

read_one <- function(ds, truth, f) {
  if (!file.exists(f)) stop("missing ", f)
  x <- fread(f)
  if (!"w1" %in% names(x)) x[, w1 := NA_integer_]
  st <- STAGES[[ds]]
  x <- x[cleaning %chin% st, .(gate, t, stage = names(st)[match(cleaning, st)], type, N, w1)]
  stopifnot(nrow(x) > 0, all(x$type %chin% all_vartypes()))
  # every (gate, t, stage) gets all 12 types, absent ones at 0
  g <- unique(x[, .(gate, t, stage)])
  x <- merge(g[, .(type = all_vartypes()), by = .(gate, t, stage)], x, by = c("gate", "t", "stage", "type"), all.x = TRUE)
  x[is.na(N), N := 0L]
  x[, w1 := w1[!is.na(w1)][1], by = .(gate, t, stage)]
  if (ds == "Aorta") {
    # aorta spectrum tables carry no window; it is in the v9 truth counts of the same build
    w <- unique(fread(sub("/eval/excluded/spectrum12_by_t\\.tsv$", "/eval_v9/excluded/truth_counts.tsv", f))[, .(gate, w1a = w1)])
    x <- merge(x, w, by = "gate", all.x = TRUE)
    x[stage == "after", w1 := w1a][, w1a := NULL]
  }
  x[, pct := 100 * N / sum(N), by = .(gate, t, stage)]
  x[, `:=`(dataset = ds, truth = truth, truth_order = match(truth, TRUTHS$truth),
           source = sub(paste0("^", REPO_DIR, "/"), "", f))]
  x
}
sp <- rbindlist(c(
  lapply(seq_len(nrow(TRUTHS)), function(i) read_one("Teratoma", TRUTHS$truth[i], ter_file(TRUTHS$rep[i]))),
  lapply(seq_len(nrow(TRUTHS)), function(i) read_one("Aorta", TRUTHS$truth[i], aor_file(TRUTHS$rep[i])))))
setcolorder(sp, c("dataset", "truth", "truth_order", "gate", "t", "stage", "type", "N", "pct", "w1", "source"))
setorder(sp, dataset, truth_order, gate, t, stage, type)
write_sd(sp, SRC_DIR, STEM, "spectrum")
print(sp[stage == "after" & t == "t50", .(sites = sum(N), AG_pct = round(pct[type == "AG"], 1),
                                          TC_pct = round(pct[type == "TC"], 1)), by = .(dataset, truth, gate)])
