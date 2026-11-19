# gt_v9_redi_ko_rho76.R   (2026-11-16; was tmp/s5/redi_rho76/redi_ko_rho76.R)
# Spearman rho between each REDIportal truth set and each KO truth set, 76 configurations (seed-8787 scoring)
# Output: human/pm/gt_eval_v9/redi_rho76/redi_ko_rho76.tsv
# Usage : Rscript human/scripts/gt_v9_redi_ko_rho76.R   (from the repo root)
suppressPackageStartupMessages(library(data.table))
PM = "human/pm/gt_eval_v9"
COMBO = c("caller", "bam_type", "strand_mode", "snp_subtract")
ko = fread(file.path(PM, "scoring_pwcov76/mismatch/config_rank.tsv"))[, c(COMBO, "gate", "t", "overall_rank"), with = FALSE]
ko[, truth := paste(gate, t, sep = "_")]
redi = rbindlist(lapply(c("cov10_alt1", "alt3", "alt3_e05"), function(v) {
  s = fread(file.path(PM, "rediportal_v9/scoring_mismatch76", paste0("gt_tier_scoring_", v, ".tsv")))
  s[, f1_rank := frank(-f1, ties.method = "min"), by = tier]
  r = s[, .(rank_sum = sum(f1_rank), mean_f1 = mean(f1, na.rm = TRUE)), by = COMBO]
  r[, overall_rank := frank(list(-mean_f1, rank_sum), ties.method = "min")][, truth := v]
}))
stopifnot(redi[, .N, by = truth]$N == 76, ko[, .N, by = truth]$N == 76)
res = CJ(redi_set = unique(redi$truth), ko_set = unique(ko$truth))
res[, rho := mapply(function(a, b) {
  m = merge(redi[truth == a, c(COMBO, "overall_rank"), with = FALSE], ko[truth == b, c(COMBO, "overall_rank"), with = FALSE], by = COMBO)
  stopifnot(nrow(m) == 76); cor(m$overall_rank.x, m$overall_rank.y, method = "spearman") }, redi_set, ko_set)]
print(res[, .(redi_set, ko_set, rho = round(rho, 3))])
cat("t50 range:", res[grepl("t50", ko_set), sprintf("%.3f-%.3f", min(rho), max(rho))], "; all t:", res[, sprintf("%.3f-%.3f", min(rho), max(rho))], "\n")
print(redi[overall_rank <= 3][order(truth, overall_rank), .(truth, overall_rank, caller, bam_type, strand_mode, snp_subtract, mean_f1 = round(mean_f1, 3))])
fwrite(res, file.path(PM, "redi_rho76", "redi_ko_rho76.tsv"), sep = "\t")
