#!/usr/bin/env Rscript
# =============================================================================
# aorta_rank_compare.R   (2026-10-26)
# Configuration ranks on the aorta truth (gt_eval_v8/scoring_pwcov/<state>/config_rank.tsv, tiers arm + sample)
# vs the teratoma headline (human/pm/gt_eval_v8/scoring_pwcov/<state>/config_rank.tsv, arm + sample + cell type)
# and the Replogle arm (human/replogle/gt_eval_v8/scoring/<state>/v7_nomm_k562/config_rank.tsv), on shared
# reported configurations, per gate x t:
#   Spearman's rho of overall mean F1 (aorta 2 tiers vs teratoma 3 tiers), and tier by tier (arm vs arm,
#   sample vs sample; Replogle arm vs aorta arm); each dataset's top configuration and its rank in the other.
# Output: gt_eval_v8/compare/{rank_rho.tsv, top_configs.tsv, caller_rank.tsv}
# Env: AORTA_SCORE_DIR (default scoring_pwcov), AORTA_CMP_DIR (default compare), REPLOGLE_TRUTH (default v7_nomm_k562)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
A = file.path(B, "mouse/aorta_gse255661/gt_eval_v8")
SCD = Sys.getenv("AORTA_SCORE_DIR", "scoring_pwcov")
RPT = Sys.getenv("REPLOGLE_TRUTH", "v7_nomm_k562")
OUT = file.path(A, Sys.getenv("AORTA_CMP_DIR", "compare")); dir.create(OUT, showWarnings = FALSE)
KEY = c("gate", "t")
rho = list(); top = list(); crk = list()
for (st in c("mismatch", "none")) {
  ao = fread(file.path(A, SCD, st, "config_rank.tsv"))
  te = fread(file.path(B, "human/pm/gt_eval_v8/scoring_pwcov", st, "config_rank.tsv"))
  rf = file.path(B, "human/replogle/gt_eval_v8/scoring", st, RPT, "config_rank.tsv")
  rp = if (file.exists(rf)) fread(rf)[reported == TRUE] else NULL
  m = merge(ao, te, by = c(KEY, "config"), suffixes = c("_aorta", "_teratoma"))
  cor_s = function(x, y) { ok = is.finite(x) & is.finite(y); if (sum(ok) < 3) NA_real_ else cor(x[ok], y[ok], method = "spearman") }
  r = m[, .(n_configs = .N,
            rho_overall = cor_s(mean_f1_aorta, mean_f1_teratoma),
            rho_arm = cor_s(f1_arm_aorta, f1_arm_teratoma),
            rho_sample = cor_s(f1_sample_aorta, f1_sample_teratoma),
            rho_aorta_arm_vs_teratoma_overall = cor_s(f1_arm_aorta, mean_f1_teratoma)), by = KEY]
  if (!is.null(rp)) {
    mr = merge(ao, rp[, .(gate, t, config, f1_replogle = f1)], by = c(KEY, "config"))
    r = merge(r, mr[, .(n_configs_replogle = .N, rho_arm_vs_replogle_arm = cor_s(f1_arm, f1_replogle)), by = KEY], by = KEY, all.x = TRUE)
  }
  rho[[st]] = r[, state := st]
  top[[st]] = rbind(
    m[overall_rank_aorta == 1, .(state = st, gate, t, which = "aorta top", config, rank_aorta = overall_rank_aorta,
                                rank_teratoma = overall_rank_teratoma, mean_f1_aorta, mean_f1_teratoma)],
    m[overall_rank_teratoma == 1, .(state = st, gate, t, which = "teratoma top", config, rank_aorta = overall_rank_aorta,
                                   rank_teratoma = overall_rank_teratoma, mean_f1_aorta, mean_f1_teratoma)])
  ca = fread(file.path(A, SCD, st, "caller_rank.tsv"))
  ct = fread(file.path(B, "human/pm/gt_eval_v8/scoring_pwcov", st, "caller_rank.tsv"))
  crk[[st]] = merge(ca[, .(gate, t, caller, caller_rank_aorta = caller_rank, best_aorta = config)],
                    ct[, .(gate, t, caller, caller_rank_teratoma = caller_rank, best_teratoma = config)],
                    by = c(KEY, "caller"))[, state := st]
}
rho = rbindlist(rho, fill = TRUE); top = rbindlist(top); crk = rbindlist(crk)
fwrite(rho, file.path(OUT, "rank_rho.tsv"), sep = "\t")
fwrite(top, file.path(OUT, "top_configs.tsv"), sep = "\t")
fwrite(crk, file.path(OUT, "caller_rank.tsv"), sep = "\t")
print(rho); print(top[state == "mismatch" & t == "t50"]); print(crk[state == "mismatch" & t == "t50"][order(gate, caller_rank_aorta)])
