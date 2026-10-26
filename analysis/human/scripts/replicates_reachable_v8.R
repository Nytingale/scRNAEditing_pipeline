#!/usr/bin/env Rscript
# replicates_reachable_v8.R   (2026-10-15; Reviewer 3 replicates, fairness check)
# The truths and the pooled recall denominator use the depth of the merged NTC arm, which a replicate-aware caller
# never sees. Here recall is restricted to "reachable" positives: >= MIN_COV reads in each of T1-T4 NTC, or in at least
# REACH_MIN_TER of them (env; default 4; outputs then in replicates/reachable_min<n>/).
#   mode prep : writes the reachable KO coverage table (unit NTC_ctrl; for gt_score_v7.R GT_UNIT_COV)
#   mode score: REDIportal recall at reachable positives (AG calls, as gt_score_rediportal_v8.R), precision unchanged
#               (taken from the existing unit_scores); KO scores from the rerun gt_score_v7.R (reachable/ko_<s>);
#               summary of pooled vs replicate strategies.
# Inputs : gt_eval_v8/unit_cov/ko_truth_unit_cov.tsv.gz; gt_eval_v6/rediportal/pileup/T{1..4}_NTC_ctrl/chr*.tsv.gz;
#          gt_eval_v8/rediportal/scoring_mismatch/positive_sites_<v>.tsv.gz; call sets *_step5_nomm.tsv.gz
# Output : gt_eval_v8/replicates/reachable/{ko_truth_unit_cov_reachable.tsv.gz, reachable_counts.tsv,
#          reachable_scores.tsv, reachable_summary.tsv}
suppressPackageStartupMessages(library(data.table))
mode <- commandArgs(trailingOnly = TRUE)[1]
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
MIN_COV <- 10L
TS <- paste0("T", 1:4, "_NTC_ctrl")
MIN_TER <- as.integer(Sys.getenv("REACH_MIN_TER", length(TS)))
V8 <- file.path(PM, "gt_eval_v8"); RR <- file.path(V8, "replicates")
O <- file.path(RR, if (MIN_TER == length(TS)) "reachable" else paste0("reachable_min", MIN_TER))
dir.create(O, showWarnings = FALSE)
STRATS <- c("pooled", "k2", "k3", "k4", "jacusa2rep")
res_dir <- function(s) if (s == "pooled") file.path(PM, "results_gt_v8") else file.path(RR, paste0("results_", s, "_arm"))

if (mode == "prep") {
  uc <- fread(cmd = paste("zcat", file.path(V8, "unit_cov", "ko_truth_unit_cov.tsv.gz")))
  r <- uc[unit %in% TS & coverage >= MIN_COV, .(n = .N, coverage = min(coverage)), by = site][n >= MIN_TER]
  fwrite(r[, .(unit = "NTC_ctrl", site, coverage)], file.path(O, "ko_truth_unit_cov_reachable.tsv.gz"), sep = "\t")
  cat("KO truth sites reachable in >=", MIN_TER, "teratomas:", nrow(r), "of", uniqueN(uc$site), "\n")
  quit(save = "no")
}
stopifnot(mode == "score")

CHRS <- paste0("chr", c(1:22, "X", "Y"))
key <- function(chr, pos) match(chr, CHRS) * 1e9 + as.numeric(pos)
site_key <- function(s) { p <- tstrsplit(s, ":", fixed = TRUE); key(p[[1]], p[[2]]) }
RDP <- file.path(PM, "gt_eval_v6", "rediportal", "pileup")
VARS <- c("cov10_alt1", "alt3", "alt3_e05")
POS <- lapply(setNames(VARS, VARS), function(v)
  site_key(fread(cmd = paste("zcat", file.path(V8, "rediportal", "scoring_mismatch", paste0("positive_sites_", v, ".tsv.gz"))))$site))
pos_any <- unique(unlist(POS))
ter_keys <- unlist(lapply(TS, function(u) unique(unlist(lapply(CHRS, function(c) {
  f <- file.path(RDP, u, paste0(c, ".tsv.gz")); if (!file.exists(f)) return(numeric())
  p <- fread(cmd = paste("zcat", f), select = c("chr", "pos", "coverage"), showProgress = FALSE)[coverage >= MIN_COV]
  k <- key(p$chr, p$pos); k[k %in% pos_any]
})))))
reach <- data.table(k = ter_keys)[, .N, by = k][N >= MIN_TER, k]
POS_R <- lapply(POS, function(p) p[p %in% reach])
ko_r <- fread(cmd = paste("zcat", file.path(O, "ko_truth_unit_cov_reachable.tsv.gz")))
truth <- fread(cmd = paste("zcat", file.path(V8, "truth", "gt_v7_truth_sites.tsv.gz")))
kp <- melt(truth, id.vars = c("gate", "site"), measure.vars = patterns("^class_t"), variable.name = "t")[value == "positive"]
kp[, `:=`(t = sub("class_", "", t), reachable = site %in% ko_r$site)]
counts <- rbind(data.table(truth = paste("REDIportal", VARS), n_pos = lengths(POS), n_reachable = lengths(POS_R)),
                kp[, .(n_pos = .N, n_reachable = sum(reachable)), by = .(truth = paste("KO", gate, t))])
counts[, pct_reachable := round(100 * n_reachable / n_pos, 1)]
print(counts); fwrite(counts, file.path(O, "reachable_counts.tsv"), sep = "\t")

redi_prec <- function(s) {
  f <- if (s == "pooled") file.path(V8, "rediportal", "scoring_mismatch", "unit_scores.tsv") else
    file.path(RR, "scoring", paste0("rediportal_", s), "unit_scores.tsv")
  fread(f)[tier == "arm", .(config, variant, precision)]
}
redi <- rbindlist(lapply(STRATS, function(s) {
  base <- file.path(res_dir(s), "NTC_ctrl", "04_filtered")
  f <- list.files(base, "_step5_nomm\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
  p <- tstrsplit(sub(paste0("^", base, "/"), "", f), "/", fixed = TRUE)
  cfg <- paste(p[[1]], p[[2]], p[[3]], sub("^snp_", "", p[[4]]), sep = "__")
  pr <- redi_prec(s)
  rbindlist(lapply(which(cfg %in% pr$config), function(i) {
    ck <- if (file.size(f[i]) == 0) numeric() else {
      x <- fread(cmd = paste("zcat", shQuote(f[i])), select = c("site", "varType"), showProgress = FALSE)[varType == "AG"]
      unique(site_key(x$site))
    }
    data.table(strategy = s, config = cfg[i], variant = VARS,
               recall = sapply(VARS, function(v) mean(POS_R[[v]] %in% ck)))
  }))[pr, on = .(config, variant), nomatch = 0L]
}))
redi[, `:=`(truth = paste("REDIportal", variant), variant = NULL)]

ko <- rbindlist(lapply(STRATS, function(s)
  fread(file.path(O, paste0("ko_", s), "gt_v7_scores_reported.tsv"))[, strategy := s]), fill = TRUE)
ko <- ko[, .(strategy, config, truth = paste("KO", gate, t), precision, recall)]

d <- rbind(redi[, .(strategy, config, truth, precision, recall)], ko)
d[, f1 := fifelse(is.na(precision) | precision + recall == 0, 0, 2 * precision * recall / (precision + recall))]
p <- tstrsplit(d$config, "__", fixed = TRUE); d[, `:=`(caller = p[[1]], bam_type = p[[2]])]
fwrite(d, file.path(O, "reachable_scores.tsv"), sep = "\t")

DEDUP <- c("sampleRmDup", "sampleUmiDup")
d[, fam := fifelse(grepl("^REDI", truth), "REDIportal", "KO")]
j <- dcast(d[caller == "jacusa2" & strategy %in% c("pooled", "jacusa2rep")], fam + truth + config + bam_type ~ strategy,
           value.var = c("precision", "recall", "f1"))
cat("\nJACUSA2, reachable positives: replicate-aware vs pooled\n")
print(j[, .(n = .N, rep_below_pooled_f1 = sum(f1_jacusa2rep < f1_pooled),
            f1_pooled = paste(round(range(f1_pooled), 2), collapse = "-"), f1_rep = paste(round(range(f1_jacusa2rep), 2), collapse = "-"),
            rec_ratio = paste(round(range(recall_pooled / recall_jacusa2rep), 1), collapse = "-")), by = .(fam, dedup = bam_type %in% DEDUP)])
m <- d[bam_type %in% DEDUP, .(precision = mean(precision, na.rm = TRUE), recall = mean(recall), f1 = mean(f1)),
       by = .(fam, strategy, config)]
sm <- m[, .(n_config = .N, median_precision = median(precision, na.rm = TRUE), median_recall = median(recall),
            median_f1 = median(f1)), by = .(fam, strategy)]
print(sm); fwrite(sm, file.path(O, "reachable_summary.tsv"), sep = "\t")
fwrite(j, file.path(O, "reachable_jacusa2.tsv"), sep = "\t")
