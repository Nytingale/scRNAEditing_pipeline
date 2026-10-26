#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_fisher_compare.R   (2026-10-18)
#
# Config and caller ranks on the Fisher-test truth sets (run_gt_v8_fisher.sh: gates depth1, fisher_bonf,
# fisher_bh) against the v8 truth sets (expected = Set A, m = 3; depth = KO cov >= 10), mismatch state,
# reported configurations, per t.
# OUTPUT: gt_eval_v8/fisher/{fisher_truth_counts,fisher_rank_compare,fisher_caller_rank}.tsv
# =============================================================================
suppressPackageStartupMessages(library(data.table))
V8 = file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8")
FT = file.path(V8, "fisher")
# SCN=scoring_pwcov: headline recall denominator; outputs suffixed _pwcov
SCN = Sys.getenv("SCN", "scoring")
SFX = sub("^scoring", "", SCN)
REFS = c("expected", "depth")

cnt = rbind(fread(file.path(V8, "truth/gt_v8_truth_counts.tsv"))[, src := "v8"],
            fread(file.path(FT, "truth_fisher/fisher_truth_counts.tsv"))[, src := "fisher"], fill = TRUE)
cnt = cnt[, .(gate, t, eligible, positive, informative_negative)]
fwrite(cnt, file.path(FT, "fisher_truth_counts.tsv"), sep = "\t")
print(dcast(cnt, gate ~ t, value.var = "positive"))

cr = rbind(fread(file.path(V8, file.path(SCN, "mismatch/config_rank.tsv"))),
           fread(file.path(FT, file.path(SCN, "mismatch/config_rank.tsv"))), fill = TRUE)
stopifnot(cr[, uniqueN(config), by = .(gate, t)][, uniqueN(V1) == 1])
new_gates = setdiff(unique(cr$gate), REFS)
cmp = rbindlist(lapply(new_gates, function(g) rbindlist(lapply(REFS, function(r) rbindlist(lapply(unique(cr$t), function(tn) {
  a = cr[gate == g & t == tn]; b = cr[gate == r & t == tn]
  m = merge(a[, .(config, rank_new = overall_rank, f1_new = mean_f1)], b[, .(config, rank_ref = overall_rank, f1_ref = mean_f1)],
            by = "config")
  data.table(gate = g, reference = r, t = tn, n_config = nrow(m),
             spearman = cor(m$rank_new, m$rank_ref, method = "spearman"),
             top_new = a[overall_rank == 1, config][1], top_ref = b[overall_rank == 1, config][1],
             top5_overlap = length(intersect(a[overall_rank <= 5, config], b[overall_rank <= 5, config])),
             top10_overlap = length(intersect(a[overall_rank <= 10, config], b[overall_rank <= 10, config])),
             median_abs_shift = median(abs(m$rank_new - m$rank_ref)), max_abs_shift = max(abs(m$rank_new - m$rank_ref)),
             mean_f1_new = mean(m$f1_new), mean_f1_ref = mean(m$f1_ref))
}))))))
fwrite(cmp, file.path(FT, paste0("fisher_rank_compare", SFX, ".tsv")), sep = "\t")
print(cmp, digits = 3)

cb = rbind(fread(file.path(V8, file.path(SCN, "mismatch/caller_rank.tsv"))),
           fread(file.path(FT, file.path(SCN, "mismatch/caller_rank.tsv"))), fill = TRUE)
cbw = dcast(cb, caller ~ gate + t, value.var = "caller_rank")
fwrite(cbw, file.path(FT, paste0("fisher_caller_rank", SFX, ".tsv")), sep = "\t")
print(cbw)
