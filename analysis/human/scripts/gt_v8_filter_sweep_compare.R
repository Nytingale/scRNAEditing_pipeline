#!/usr/bin/env Rscript
# Rank stability across the downstream-filter sweep (human/pm/scripts/filter_sweep_v8.sh), as
# gt_v8_sweep_rank_compare.R for the truth sweep. For each setting (one downstream filter moved or turned off)
# x gate x depletion threshold: Spearman rho of overall_rank vs the centre (main filters), rank-1 config, and the
# rank of the centre's rank-1 config. The sweep centre must reproduce the main ranking (scoring/mismatch).
# Out: <W>/filter_sweep_compare.tsv, <W>/filter_sweep_top.tsv, <W>/filter_sweep_sites.tsv
#      (figure: publication/teratoma_gt/scripts/TERATOMA_v8_filter_sweep_rank_plot.R)
suppressPackageStartupMessages(library(data.table))

V8 = file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8")
W  = Sys.getenv("W", file.path(V8, "filter_sweep"))
# SCN=scoring_pw: sweep scored with the gate-matched recall denominator, compared with V8/scoring_pw; outputs suffixed _pw
SCN = Sys.getenv("SCN", "scoring")
SFX = sub("^scoring", "", SCN)

settings = basename(list.dirs(file.path(W, SCN), recursive = FALSE))
rk = rbindlist(lapply(settings, function(s) {
  f = file.path(W, SCN, s, "mismatch", "config_rank.tsv")
  if (!file.exists(f)) return(NULL)
  fread(f, select = c("gate", "t", "config", "caller", "overall_rank", "mean_f1"))[, setting := s]
}))
stopifnot("centre" %in% rk$setting)

main = fread(file.path(V8, SCN, "mismatch", "config_rank.tsv"), select = c("gate", "t", "config", "overall_rank", "mean_f1"))
chk = merge(rk[setting == "centre"], main, by = c("gate", "t", "config"), suffixes = c("", "_main"))
stopifnot(nrow(chk) == nrow(main), nrow(chk) == nrow(rk[setting == "centre"]))
cat("centre vs main ranking: identical ranks", chk[, all(overall_rank == overall_rank_main)],
    "; max |mean F1 diff|", chk[, signif(max(abs(mean_f1 - mean_f1_main)), 3)], "\n")
stopifnot(rk[, uniqueN(config), by = .(setting, gate, t)][, all(V1 == V1[1])])

cen = rk[setting == "centre", .(gate, t, config, rank_centre = overall_rank)]
cmp = merge(rk[setting != "centre"], cen, by = c("gate", "t", "config"))
top_cen = cen[rank_centre == 1, .(top_centre = config[1]), by = .(gate, t)]
res = cmp[, .(rho = cor(overall_rank, rank_centre, method = "spearman"), n_config = .N,
              top = config[which.min(overall_rank)]), by = .(setting, gate, t)]
res = merge(res, top_cen, by = c("gate", "t"))
res[, centre_top_rank := cmp[res, on = .(setting, gate, t, config = top_centre), x.overall_rank]]
setorder(res, setting, gate, t)
fwrite(res, file.path(W, paste0("filter_sweep_compare", SFX, ".tsv")), sep = "\t")

top = rk[overall_rank <= 3, .(setting, gate, t, overall_rank, config, mean_f1)]
setorder(top, setting, gate, t, overall_rank)
fwrite(top, file.path(W, paste0("filter_sweep_top", SFX, ".tsv")), sep = "\t")

print(res[t == "t50", .(setting, gate, rho = round(rho, 3), centre_top_rank, top)], nrows = 200)
