#!/usr/bin/env Rscript
# gt_v8_state_caller_summary.R   (2026-10-15)
# Which caller / configuration is on top, per call state, under all nine v8 truths (6 KO + 3 REDIportal)?
# Ranking as TERATOMA_v8_rank_robustness_plot.R: per truth, F1 rank per tier (min ties), rank sum over arm /
# sample / cell type, overall rank = rank of (rank sum, -mean F1). KO ranks are taken from scoring/<state>/config_rank.tsv.
# Output gt_eval_v8/state_summary/:
#   config_ranks.tsv   state x truth x configuration: overall_rank, mean_f1
#   config_robust.tsv  state x family (ko / redi / all) x configuration: best, worst, mean rank
#   top_by_truth.tsv   state x truth: rank-1 configuration and the best rank of each caller
#   family_rho.tsv     state: Spearman's rho of KO-family vs REDIportal-family mean ranks
# SCN=scoring_pwcov (headline): scoring_pwcov/, rediportal_pwcov/ -> state_summary_pwcov/
suppressMessages(library(data.table))
V8 <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8")
SCN <- Sys.getenv("SCN", "scoring")
SFX <- sub("^scoring", "", SCN)
OUT <- file.path(V8, paste0("state_summary", SFX)); dir.create(OUT, showWarnings = FALSE)
STATES <- c("none", "mismatch", "clustered", "both")
REDI <- c("cov10_alt1", "alt3", "alt3_e05")
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
TIER_MAP <- c(dataset = "arm", sample = "sample", celltype = "celltype")

rank_state <- function(s) {
  ko <- fread(file.path(V8, SCN, s, "config_rank.tsv"))
  ko <- ko[, c(COMBO, "gate", "t", "overall_rank", "mean_f1"), with = FALSE][, truth := paste0(gate, "_", t)][, c("gate", "t") := NULL]
  rd <- rbindlist(lapply(REDI, function(v) {
    x <- fread(file.path(V8, paste0("rediportal", SFX), paste0("scoring_", s), paste0("gt_tier_scoring_", v, ".tsv")))
    stopifnot(all(x$tier %in% names(TIER_MAP)))
    x[, f1_rank := frank(-f1, ties.method = "min"), by = tier]
    r <- x[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = COMBO]
    r[, overall_rank := frank(list(-mean_f1, rank_sum), ties.method = "min")]
    r[, c(COMBO, "overall_rank", "mean_f1"), with = FALSE][, truth := paste0("redi_", v)]
  }))
  out <- rbind(ko, rd)[, state := s]
  out[, family := ifelse(grepl("^redi_", truth), "redi", "ko")]
  out[, config := paste(caller, bam_type, strand_mode, snp_subtract, sep = "__")]
  stopifnot(out[, uniqueN(config), by = truth][, all(V1 == V1[1])])
  out
}
rk <- rbindlist(lapply(STATES[file.exists(file.path(V8, SCN, STATES, "config_rank.tsv")) &
                              dir.exists(file.path(V8, paste0("rediportal", SFX), paste0("scoring_", STATES)))], rank_state))
fwrite(rk, file.path(OUT, "config_ranks.tsv"), sep = "\t")

rob <- rbind(rk[, .(best = min(overall_rank), worst = max(overall_rank), mean_rank = mean(overall_rank)), by = .(state, family, config, caller)],
             rk[, .(best = min(overall_rank), worst = max(overall_rank), mean_rank = mean(overall_rank)), by = .(state, config, caller)][, family := "all"])
setorder(rob, state, family, mean_rank)
fwrite(rob, file.path(OUT, "config_robust.tsv"), sep = "\t")

top <- rk[, c(list(rank1 = config[which.min(overall_rank)], rank1_caller = caller[which.min(overall_rank)]),
              as.list(tapply(overall_rank, caller, min))), by = .(state, truth)]
fwrite(top, file.path(OUT, "top_by_truth.tsv"), sep = "\t")

rho <- rob[family != "all", .(mean_rank = mean_rank), by = .(state, family, config)]
rho <- dcast(rho, state + config ~ family, value.var = "mean_rank")[, .(rho_ko_vs_redi = cor(ko, redi, method = "spearman")), by = state]
fwrite(rho, file.path(OUT, "family_rho.tsv"), sep = "\t")

for (s in unique(rk$state)) {
  cat("\n=====", s, "=====\n")
  print(top[state == s, !"state"])
  print(rob[state == s & family == "all"][1:5, .(config, best, worst, mean_rank = round(mean_rank, 1))])
}
print(rho)
