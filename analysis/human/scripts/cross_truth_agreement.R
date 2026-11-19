#!/usr/bin/env Rscript
# cross_truth_agreement.R DS
# Agreement truth: a deduplication-neutral truth from the four truths of the cross-truth test (v9 = unseeded
# umi_tools, seed1, seed2, umimd). Per (gate, site) and threshold, the class is kept when all four truths give the same
# class; otherwise, or when the site is missing from any truth, it is "unresolved" (ignored by gt_score_v7.R).
# Unit coverage (recall denominator, editing level for the expected gate): per (unit, site), the minimum coverage and
# minimum variant coverage over the four unit-coverage tables, so a positive enters the denominator only if it passes
# under every deduplication.
# Output: human/pm/gt_eval_v9/cross_truth/truth_agree/<DS>/{truth_sites,unit_cov}.tsv.gz, agreement_counts.tsv
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("SCRNA_DIR", ".")
setwd(B)
DS <- commandArgs(TRUE)[1]
stopifnot(DS %in% c("teratoma", "aorta"))
TRUTHS <- c("v9", "seed1", "seed2", "umimd")
if (DS == "teratoma") {
  dir <- c(v9 = "human/pm/gt_eval_v9", setNames(file.path("human/pm/umi_seed_reps", TRUTHS[-1], "v9"), TRUTHS[-1]))
  tf <- file.path(dir, "truth/gt_v9_truth_sites.tsv.gz")
  uf <- file.path(dir, "unit_cov/ko_truth_unit_cov.tsv.gz")
} else {
  A <- "mouse/aorta_gse255661"
  dir <- c(v9 = file.path(A, "truth/eval_v9/excluded"),
           setNames(file.path(A, "umi_seed_reps", TRUTHS[-1], "eval_v9/excluded"), TRUTHS[-1]))
  tf <- file.path(dir, "gt_aorta_truth_sites.tsv.gz")
  uf <- file.path(dir, "unit_cov.tsv.gz")
}
names(tf) <- names(uf) <- TRUTHS
CL <- paste0("class_", c("t25", "t50", "t75"))
O <- file.path("human/pm/gt_eval_v9/cross_truth/truth_agree", DS)
dir.create(O, recursive = TRUE, showWarnings = FALSE)

tr <- lapply(TRUTHS, function(x) fread(tf[[x]]))
names(tr) <- TRUTHS
key <- c("gate", "site")
base <- unique(rbindlist(lapply(tr, function(x) x[, ..key])))
ag <- copy(tr$v9)[base, on = key]
for (x in TRUTHS[-1]) {
  y <- tr[[x]][, c(key, "chr", "pos", "genomic_ref", "a_to_i_strand", CL), with = FALSE]
  ag[y, on = key, `:=`(chr = fcoalesce(chr, i.chr), pos = fcoalesce(pos, i.pos),
                       genomic_ref = fcoalesce(genomic_ref, i.genomic_ref),
                       a_to_i_strand = fcoalesce(a_to_i_strand, i.a_to_i_strand))]
}
cnt <- list()
for (cl in CL) {
  m <- Reduce(function(a, b) merge(a, b, by = key, all = TRUE),
              lapply(TRUTHS, function(x) setnames(tr[[x]][, c(key, cl), with = FALSE], cl, x)))
  m[, agree := fifelse(fcoalesce(v9 == seed1 & v9 == seed2 & v9 == umimd, FALSE), v9, "unresolved")]
  stopifnot(!anyNA(m$agree))
  ag[m, on = key, (cl) := i.agree]
  cnt[[cl]] <- rbind(
    m[, .(truth = "agree", n = sum(agree == "positive")), by = gate],
    rbindlist(lapply(TRUTHS, function(x) m[, .(truth = x, n = sum(get(x) == "positive", na.rm = TRUE)), by = gate])))[
    , class := cl]
}
cnt <- dcast(rbindlist(cnt), gate + class ~ truth, value.var = "n")[, c("gate", "class", TRUTHS, "agree"), with = FALSE]
cnt[, lost_vs_v9 := v9 - agree]
fwrite(cnt, file.path(O, "agreement_counts.tsv"), sep = "\t")
print(cnt)
cat("negatives (informative_negative) per gate x threshold:\n")
print(melt(ag[, lapply(.SD, function(v) sum(v == "informative_negative")), by = gate, .SDcols = CL], id.vars = "gate"))
fwrite(ag, file.path(O, "truth_sites.tsv.gz"), sep = "\t")

uc <- rbindlist(lapply(TRUTHS, function(x) fread(uf[[x]])[, truth := x]))
u <- uc[, .(coverage = min(coverage), variant_cov = min(variant_cov), n = .N), by = .(unit, site)][n == length(TRUTHS)]
psites <- unique(ag[get(CL[1]) == "positive" | get(CL[2]) == "positive" | get(CL[3]) == "positive", site])
cat("agreement positives:", length(psites), "; with unit coverage in all four tables:",
    sum(psites %chin% u$site), "\n")
stopifnot(all(psites %chin% u$site))
fwrite(u[, .(unit, site, coverage, variant_cov)], file.path(O, "unit_cov.tsv.gz"), sep = "\t")
