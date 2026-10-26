#!/usr/bin/env Rscript
# =============================================================================
# gt_multiallelic_counts_v6.R   (2026-10-01, nine-truth suite)
#
# Records how many sites of each of the nine truth sets are multiallelic in the merged NTC
# pileup. Record only: no truth set or score is changed.
#
# Multiallelic: a third allele (not the reference, not the edited base) on the transcribed
# strand has >= 3 reads and >= 5% of the transcribed-strand depth. With the stranded parser
# (forward library) the transcribed strand is set by the reference base:
#   ref A: plus-strand reads; edited base G; third alleles C, T (uppercase counts)
#   ref T: minus-strand reads; edited base c; third alleles a, g (lowercase counts)
#   depth = coverage (ref + edited reads) + third-allele reads
#
# Sets
#   rediportal_cov10_alt1 / _alt3 / _alt3_e05   positives only (REDIportal has no negatives)
#   gt_expected_t25 / t50 / t75   ground truth, KO arm covered if KO cov x NTC edit >= 3
#   gt_depth_t25 / t50 / t75      ground truth, KO arm covered if KO cov >= 10
#
# Inputs : gt_eval_v6/rediportal/scoring/positive_sites_<variant>.tsv.gz
#          gt_eval_v6/rediportal/pileup_full/merged_NTC_ctrl/chr*.tsv.gz
#          gt_eval_v6/scoring/sweep_infneg{,_depth}/truth_sweep_classes.tsv.gz
#          gt_eval_v6/pileups/bq30/pileup_merged_NTC_ctrl.tsv, gt_eval_v6/truth/gt_truth_sites.tsv.gz
# Outputs: gt_eval_v6/multiallelic/
#          sites_rediportal_<variant>.tsv.gz, sites_gt_{expected,depth}.tsv.gz  per-site flag tables
#          multiallelic_counts.tsv   set x class: n_sites, n_multiallelic, share
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("SCRNA_DIR", ".")
V6 <- file.path(B, "human/pm/gt_eval_v6")
OUT <- file.path(V6, "multiallelic"); dir.create(OUT, showWarnings = FALSE)
MIN_THIRD_READS <- 3; MIN_THIRD_FRAC <- 0.05
CHRS <- c(paste0("chr", 1:22), "chrX", "chrY")
REDI <- c("cov10_alt1", "alt3", "alt3_e05")
GATES <- c(expected = "sweep_infneg", depth = "sweep_infneg_depth")
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

flag <- function(p) {
  isA <- toupper(p$ref) == "A"
  p[, `:=`(third_1 = fifelse(isA, C, a), third_2 = fifelse(isA, T, g))]
  p[, `:=`(third_max = pmax(third_1, third_2), third_total = third_1 + third_2)]
  p[, depth := coverage + third_total]
  p[, multiallelic := third_max >= MIN_THIRD_READS & third_max >= MIN_THIRD_FRAC * depth]
  p[, .(site = paste0(chr, ":", pos), coverage, variant_cov, depth, third_max, third_total, multiallelic)]
}
count <- function(d, set, cls_col) {
  d[, .(n_sites = .N, n_multiallelic = sum(multiallelic)), by = c(class = cls_col)][
    , `:=`(set = set, share = n_multiallelic / n_sites)][, .(set, class, n_sites, n_multiallelic, share)]
}
res <- list()

# ---- REDIportal sets --------------------------------------------------------------------
pf <- file.path(V6, "rediportal", "pileup_full", "merged_NTC_ctrl", paste0(CHRS, ".tsv.gz"))
if (!all(file.exists(pf))) stop("missing REDIportal full pileup: ", sum(!file.exists(pf)), " chromosomes")
pr <- rbindlist(lapply(pf, function(f) fread(f, showProgress = FALSE)[coverage >= 10 & variant_cov >= 1]))
pr <- flag(pr)
for (v in REDI) {
  s <- fread(file.path(V6, "rediportal", "scoring", paste0("positive_sites_", v, ".tsv.gz")))
  d <- pr[s, on = "site"]
  if (anyNA(d$coverage)) stop(v, ": ", sum(is.na(d$coverage)), " positives missing from the full pileup")
  d[, class := "positive"]
  fwrite(d, file.path(OUT, paste0("sites_rediportal_", v, ".tsv.gz")), sep = "\t")
  res[[v]] <- count(d, paste0("rediportal_", v), "class")
  msg("rediportal ", v, ": ", nrow(d), " positives, ", sum(d$multiallelic), " multiallelic")
}
rm(pr)

# ---- ground-truth sets ------------------------------------------------------------------
pv <- fread(file.path(V6, "pileups", "bq30", "pileup_merged_NTC_ctrl.tsv"), showProgress = FALSE)
pv <- flag(pv)
tr <- fread(file.path(V6, "truth", "gt_truth_sites.tsv.gz"), select = c("site", "cov_NTC_ctrl", "alt_NTC_ctrl"))
for (g in names(GATES)) {
  cl <- fread(file.path(V6, "scoring", GATES[[g]], "truth_sweep_classes.tsv.gz"))
  d <- merge(cl, pv, by = "site", all.x = TRUE)
  chk <- merge(d[, .(site, coverage, variant_cov)], tr, by = "site")
  if (anyNA(d$coverage) || nrow(chk) != nrow(d) || any(chk$coverage != chk$cov_NTC_ctrl | chk$variant_cov != chk$alt_NTC_ctrl))
    stop(g, ": pileup does not match the truth table at eligible sites")
  fwrite(d, file.path(OUT, paste0("sites_gt_", g, ".tsv.gz")), sep = "\t")
  for (tc in grep("^class_t", names(d), value = TRUE)) {
    set <- paste0("gt_", g, "_", sub("^class_", "", tc))
    x <- count(d, set, tc)
    res[[set]] <- rbind(x, data.table(set = set, class = "eligible", n_sites = nrow(d),
                                      n_multiallelic = sum(d$multiallelic), share = mean(d$multiallelic)))
  }
  msg("ground truth ", g, ": ", nrow(d), " eligible, ", sum(d$multiallelic), " multiallelic")
}

out <- rbindlist(res)
out[, class := factor(class, c("positive", "informative_negative", "unresolved", "eligible"))]
setorder(out, set, class)
fwrite(out, file.path(OUT, "multiallelic_counts.tsv"), sep = "\t")
print(out)
msg("done -> ", OUT)
