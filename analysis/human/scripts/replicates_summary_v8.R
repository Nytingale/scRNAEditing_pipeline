#!/usr/bin/env Rscript
# =============================================================================
# replicates_summary_v8.R   (2026-10-15; Reviewer 3 replicates; from _v6)
# Arm-level NTC call sets (main state, _step5_nomm) per replicate strategy (pooled, k2, k3, k4, jacusa2rep):
# precision / recall / F1 against the three REDIportal truths and the six KO truths, and A-to-G % of the call set.
# Reported configurations only (as in the scorers); pooled = the existing v8 arm-tier scoring.
# Output: gt_eval_v8/replicates/summary/{replicate_scores, replicate_summary, jacusa2_compare}.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
V8 <- file.path(PM, "gt_eval_v8")
R <- file.path(V8, "replicates"); O <- file.path(R, "summary"); dir.create(O, showWarnings = FALSE)
STRAT <- c(pooled = "Pooled", k2 = ">= 2 teratomas", k3 = ">= 3 teratomas", k4 = "All 4 teratomas",
           jacusa2rep = "JACUSA2 replicate-aware")
RE <- "_step5_nomm\\.tsv\\.gz$"
redi_dir <- function(s) if (s == "pooled") file.path(V8, "rediportal", "scoring_mismatch") else file.path(R, "scoring", paste0("rediportal_", s))
ko_dir <- function(s) if (s == "pooled") file.path(V8, "scoring", "mismatch", "arm") else file.path(R, "scoring", paste0("ko_", s))
res_dir <- function(s) if (s == "pooled") file.path(PM, "results_gt_v8") else file.path(R, paste0("results_", s, "_arm"))

redi <- rbindlist(lapply(names(STRAT), function(s)
  fread(file.path(redi_dir(s), "unit_scores.tsv"))[tier == "arm"][, strategy := s]))
redi <- redi[, .(strategy, truth = paste0("REDIportal ", variant), config, caller, bam_type, strand_mode,
                 n_calls, precision, recall, f1)]
ko <- rbindlist(lapply(names(STRAT), function(s)
  fread(file.path(ko_dir(s), "gt_v7_scores_reported.tsv"))[, strategy := s]))
ko <- ko[, .(strategy, truth = paste0("KO ", gate, " ", t), config, caller, bam_type, strand_mode,
             n_calls = n_tp + n_fp, precision, recall, f1)]
sc <- rbind(redi, ko)

ag <- rbindlist(lapply(names(STRAT), function(s) {
  root <- file.path(res_dir(s), "NTC_ctrl", "04_filtered")
  f <- list.files(root, RE, recursive = TRUE, full.names = TRUE)
  rbindlist(lapply(f, function(p) {
    x <- if (file.size(p) == 0) data.table(varType = character()) else fread(p, select = "varType")
    q <- tstrsplit(sub(paste0("^", root, "/"), "", p), "/", fixed = TRUE)
    data.table(strategy = s, config = paste(q[[1]], q[[2]], q[[3]], sub("^snp_", "", q[[4]]), sep = "__"),
               n_sites = nrow(x), pct_ag = if (nrow(x)) 100 * mean(x$varType == "AG") else NA_real_)
  }))
}))
d <- merge(sc, ag, by = c("strategy", "config"), all.x = TRUE)
d[, Strategy := factor(STRAT[strategy], STRAT)]
fwrite(d, file.path(O, "replicate_scores.tsv"), sep = "\t")

sm <- d[, .(n_config = .N, median_precision = median(precision, na.rm = TRUE), median_recall = median(recall, na.rm = TRUE),
            median_f1 = median(f1, na.rm = TRUE), median_pct_ag = median(pct_ag, na.rm = TRUE),
            median_sites = median(n_sites)), by = .(truth, Strategy)][order(truth, Strategy)]
fwrite(sm, file.path(O, "replicate_summary.tsv"), sep = "\t")
cat("all reported configurations (arm level):\n"); print(sm, digits = 3)

# JACUSA2: same BAM type x strand mode, SNP blacklist off (the replicate-aware run's setting)
j <- d[caller == "jacusa2" & sub(".*__", "", config) == "off"]
jc <- dcast(j, truth + bam_type + strand_mode ~ strategy, value.var = c("precision", "recall", "f1", "pct_ag", "n_sites"))
fwrite(jc, file.path(O, "jacusa2_compare.tsv"), sep = "\t")
cat("\nJACUSA2 (SNP blacklist off) by strategy, range over BAM types x strand modes:\n")
print(j[, .(n = .N, precision = sprintf("%.2f-%.2f", min(precision), max(precision)),
            recall = sprintf("%.3f-%.3f", min(recall), max(recall)), f1 = sprintf("%.2f-%.2f", min(f1), max(f1)),
            pct_ag = sprintf("%.0f-%.0f", min(pct_ag), max(pct_ag)), sites = median(n_sites)),
        by = .(truth, Strategy)][order(truth, Strategy)])
best <- d[, .SD[which.max(f1)], by = .(truth, Strategy)][, .(truth, Strategy, config, f1, precision, recall, pct_ag)]
cat("\nbest configuration per strategy:\n"); print(best[order(truth, Strategy)], digits = 3)
