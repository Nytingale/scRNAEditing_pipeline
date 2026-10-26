#!/usr/bin/env Rscript
# ============================================================================
# gt_v8_sites_inflation.R   (2026-10-21)
# Can A-to-G site count be inflated by loosening filters without F1 gain? Uses the old filter sweep
# (gt_eval_v8/filter_sweep, old scoring, mismatch state). A-to-G sites / percentage from the arm tier
# (NTC_ctrl) calls of each setting; F1 / precision / recall = config_rank.tsv mean over the 6 truths.
# Settings restricted to some callers are used only for those callers.
# Output: gt_eval_v8/filter_sweep/sites_inflation.tsv (+ summary printed)
# ============================================================================
suppressPackageStartupMessages(library(data.table))
W <- "human/pm/gt_eval_v8/filter_sweep"
st <- fread(file.path(W, "settings.tsv"))
SET <- c("centre", "alt_off", "edit_off", "dbsnp_off", "cov_off", "multimm_off", "jscore_off", "window_off")
callers_of <- function(s) {
  if (s == "window_off") return(NULL)
  x <- st[setting == s, callers]; if (x == "all") NULL else strsplit(x, ",")[[1]]
}
cfg_of <- function(f) {
  p <- strsplit(f, "/")[[1]]; i <- which(p == "04_filtered")
  paste(p[i + 1], p[i + 2], p[i + 3], sub("snp_", "", p[i + 4]), sep = "__")
}
cnt <- rbindlist(lapply(SET, function(s) {
  src <- if (s == "window_off") "centre" else s
  re <- if (s == "window_off") "_step4_cov_edit_filtered\\.tsv\\.gz$" else "_step5_nomm\\.tsv\\.gz$"
  fs <- list.files(file.path(W, "results", src, "results_gt_v8", "NTC_ctrl"), re, recursive = TRUE, full.names = TRUE)
  rbindlist(lapply(fs, function(f) {
    v <- fread(cmd = paste("zcat", shQuote(f)), select = "varType")$varType
    data.table(setting = s, config = cfg_of(f), n = length(v), ag = sum(v == "AG"))
  }))
}))
cnt[, caller := sub("__.*", "", config)]
cnt <- cnt[mapply(function(s, c) { cl <- callers_of(s); is.null(cl) || c %in% cl }, setting, caller)]
sc <- rbindlist(lapply(SET, function(s) {
  r <- fread(file.path(W, "scoring", s, "mismatch", "config_rank.tsv"))
  r[, .(setting = s, mean_f1 = mean(mean_f1),
        precision = mean((precision_arm + precision_sample + precision_celltype) / 3),
        recall = mean((recall_arm + recall_sample + recall_celltype) / 3)), by = config]
}))
d <- merge(cnt[, .(setting, config, caller, n, ag, pct = 100 * ag / n)], sc, by = c("setting", "config"))
d <- d[config %in% d[setting == "centre", config]]
ce <- d[setting == "centre", .(config, ag0 = ag, pct0 = pct, f10 = mean_f1, p0 = precision, r0 = recall)]
d <- merge(d, ce, by = "config")
d[, `:=`(ag_fold = ag / ag0, d_pct = pct - pct0, d_f1 = mean_f1 - f10, d_prec = precision - p0, d_rec = recall - r0)]
fwrite(d, file.path(W, "sites_inflation.tsv"), sep = "\t")

cat("\n== per setting vs centre (median over configs) ==\n")
print(d[setting != "centre", .(n_cfg = .N, ag_fold = round(median(ag_fold), 2), d_pct = round(median(d_pct), 1),
                               d_f1 = round(median(d_f1), 3), d_prec = round(median(d_prec), 3),
                               d_rec = round(median(d_rec), 3),
                               more_ag_lower_f1 = sum(ag_fold > 1.05 & d_f1 < -0.005)), by = setting])

cat("\n== pooled over config x setting: which metric tracks F1 ==\n")
d[, `:=`(r_ag = frank(-ag), r_pct = frank(-pct))][, comp := (r_ag + r_pct) / 2]
rho <- function(x) round(cor(x, d$mean_f1, method = "spearman"), 3)
print(c(n = nrow(d), sites = rho(d$ag), pct = rho(d$pct), composite = rho(-d$comp)))
cat("\n== top 5 picked by each metric (pooled) ==\n")
for (m in c("ag", "pct", "comp", "mean_f1")) {
  o <- switch(m, ag = order(-d$ag), pct = order(-d$pct), comp = order(d$comp), mean_f1 = order(-d$mean_f1))
  cat(m, "\n"); print(d[o][1:5, .(setting, config, ag, pct = round(pct, 1), mean_f1 = round(mean_f1, 3),
                                  precision = round(precision, 3), recall = round(recall, 3))])
}
