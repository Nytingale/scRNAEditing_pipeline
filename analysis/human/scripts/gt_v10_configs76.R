#!/usr/bin/env Rscript
# gt_v10_configs76.R
# Seeded (v10) replacement of tmp/meeting_v9_2026-11-07/configs60_v9.tsv for the 76 ranked configurations:
#   ag_*, pct_*  A-to-G sites and A-to-G % of the unique (site, varType) of each configuration's step5_nomm call set in
#                the arm unit (teratoma NTC_ctrl, aorta CTRL); call-set roots as human/pm/scripts/cross_truth_score.sh
#                (published, Cellsnp-lite hybrid, UMI-aware MarkDuplicates)
#   f1_*, precision_*, recall_*  mean over the six KO truths of the cross-tier means (teratoma arm / sample / cell type,
#                aorta arm / sample) from gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv of each dataset
# Output: human/pm/gt_eval_v9/configs76/configs76_seeded.tsv
# Usage : qsub human/pm/scripts/gt_v10_configs76.sh
suppressPackageStartupMessages(library(data.table))
setwd(Sys.getenv("SCRNA_DIR", "."))
A <- "mouse/aorta_gse255661"
DS <- list(
  te = list(unit = "NTC_ctrl", tiers = c("arm", "sample", "celltype"),
            rank = "human/pm/gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv",
            roots = c("human/pm/results_gt_v8", "human/pm/gt_eval_v8/cellsnp_hybrid/results/centre/results_gt_v8",
                      "human/pm/results_gt_v8_umimd")),
  ao = list(unit = "CTRL", tiers = c("arm", "sample"),
            rank = file.path(A, "gt_eval_v9/scoring_pwcov76/mismatch/config_rank.tsv"),
            roots = file.path(A, c("results_v8", "gt_eval_v8/cellsnp_hybrid/results/results_v8", "results_v8_umimd"))))

out <- NULL
for (ds in names(DS)) {
  s <- DS[[ds]]
  rk <- fread(s$rank)
  stopifnot(rk[, .N, by = .(gate, t)][, all(N == 76L)], uniqueN(rk[, .(gate, t)]) == 6L)
  cfg <- unique(rk$config)
  met <- rk[, .(f1 = mean(mean_f1),
                precision = mean(rowMeans(.SD[, paste0("precision_", s$tiers), with = FALSE])),
                recall = mean(rowMeans(.SD[, paste0("recall_", s$tiers), with = FALSE]))), by = .(config, caller)]
  calls <- rbindlist(lapply(s$roots, function(r) {
    d <- file.path(r, s$unit, "04_filtered")
    f <- list.files(d, "_step5_nomm\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
    p <- tstrsplit(sub(paste0("^", d, "/"), "", dirname(f)), "/", fixed = TRUE)
    data.table(file = f, config = paste(p[[1]], p[[2]], p[[3]], sub("^snp_", "", p[[4]]), sep = "__"))
  }))[config %in% cfg]
  stopifnot(nrow(calls) == 76L, setequal(calls$config, cfg))
  calls[, c("ag", "pct") := {
    v <- if (file.size(file) == 0) data.table(varType = character()) else unique(fread(file, select = c("site", "varType")))
    list(sum(v$varType == "AG"), 100 * mean(v$varType == "AG"))
  }, by = file]
  x <- merge(met, calls[, .(config, ag, pct)], by = "config")
  setnames(x, c("ag", "pct", "f1", "precision", "recall"), paste0(c("ag", "pct", "f1", "precision", "recall"), "_", ds))
  out <- if (is.null(out)) x else merge(out, x, by = c("config", "caller"))
}
stopifnot(nrow(out) == 76L, !anyNA(out))
fwrite(out, "human/pm/gt_eval_v9/configs76/configs76_seeded.tsv", sep = "\t")
# comparison with the unseeded 60-configuration table, where it still exists (not part of the repository)
OLD <- "tmp/meeting_v9_2026-11-07/configs60_v9.tsv"
if (file.exists(OLD)) {
  cmp <- merge(fread(OLD), out, by = "config", suffixes = c("_v9", "_v10"))
  cat("configs shared with configs60_v9:", nrow(cmp), "\n")
  print(cmp[, .(max_abs_d_ag_te = max(abs(ag_te_v9 - ag_te_v10)), max_abs_d_f1_te = max(abs(f1_te_v9 - f1_te_v10)),
                max_abs_d_ag_ao = max(abs(ag_ao_v9 - ag_ao_v10)), max_abs_d_f1_ao = max(abs(f1_ao_v9 - f1_ao_v10)))])
}
