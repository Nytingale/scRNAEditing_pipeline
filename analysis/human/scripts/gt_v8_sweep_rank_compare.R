#!/usr/bin/env Rscript
# Rank stability across the v8 truth-parameter sweep (human/pm/scripts/run_gt_v8_sweep_rank.sh).
# For each neighbour truth (one parameter one step from the centre) x gate x depletion
# threshold: Spearman rho of overall_rank vs the centre truth, rank-1 config, and the
# rank of the centre's rank-1 config.
# Out: <W>/sweep_rank_compare.tsv, <W>/sweep_rank_top.tsv
#      (figure: publication/teratoma_gt/scripts/TERATOMA_v8_sweep_rank_plot.R)
suppressPackageStartupMessages(library(data.table))

W   = Sys.getenv("W", file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8/sweep_rank"))

# SCN=scoring_pwcov: headline recall denominator (gt_score_v7.R GT_UNIT_RULE=gate_cov); outputs suffixed _pwcov
SCN = Sys.getenv("SCN", "scoring")
SFX = sub("^scoring", "", SCN)
truths = basename(list.dirs(file.path(W, SCN), recursive = FALSE))
rk = rbindlist(lapply(truths, function(tr) {
  f = file.path(W, SCN, tr, "mismatch", "config_rank.tsv")
  if (!file.exists(f)) return(NULL)
  d = fread(f, select = c("gate", "t", "config", "caller", "overall_rank", "mean_f1"))
  d[, truth := tr]
}))
stopifnot("centre" %in% rk$truth)

cen = rk[truth == "centre", .(gate, t, config, rank_centre = overall_rank)]
cmp = merge(rk[truth != "centre"], cen, by = c("gate", "t", "config"))
top_cen = cen[rank_centre == 1, .(gate, t, top_centre = config)]

res = cmp[, .(rho = cor(overall_rank, rank_centre, method = "spearman"),
              n_config = .N,
              top = config[which.min(overall_rank)]), by = .(truth, gate, t)]
res = merge(res, top_cen, by = c("gate", "t"))
res[, centre_top_rank := cmp[res, on = .(truth, gate, t, config = top_centre), x.overall_rank]]
res[, parameter := sub("_[^_]+$", "", truth)]
setorder(res, parameter, truth, gate, t)
fwrite(res, file.path(W, paste0("sweep_rank_compare", SFX, ".tsv")), sep = "\t")

top = rk[overall_rank <= 3, .(truth, gate, t, overall_rank, config, mean_f1)]
setorder(top, truth, gate, t, overall_rank)
fwrite(top, file.path(W, paste0("sweep_rank_top", SFX, ".tsv")), sep = "\t")

print(res[, .(truth, gate, t, rho = round(rho, 3), centre_top_rank)], nrows = 200)
print(res[, .N, by = top])