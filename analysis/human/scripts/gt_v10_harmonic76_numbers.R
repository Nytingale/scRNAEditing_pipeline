# gt_v10_harmonic76_numbers.R   (2026-11-15; was tmp/s5/harm76/harm76_numbers.R)
# Harmonic score (A-to-G sites and A-to-G %, each scaled to its maximum) vs mean F1 over the 76 configurations, both
# datasets: Spearman rho, top configuration by each and the F1 gap (response: harmonic-score sentence).
# Input : human/pm/gt_eval_v9/configs76/configs76_seeded.tsv (gt_v10_configs76.R)
# Usage : Rscript human/scripts/gt_v10_harmonic76_numbers.R
suppressPackageStartupMessages(library(data.table))
setwd(Sys.getenv("SCRNA_DIR", "."))
x0 <- fread("human/pm/gt_eval_v9/configs76/configs76_seeded.tsv"); stopifnot(nrow(x0) == 76)
x <- rbind(x0[, .(ds = "teratoma", config, sites = ag_te, pct = pct_te, f1 = f1_te)],
           x0[, .(ds = "aorta", config, sites = ag_ao, pct = pct_ao, f1 = f1_ao)])
x[, `:=`(a = sites / max(sites), b = pct / max(pct)), by = ds][, harm := 2 * a * b / (a + b)]
x[, `:=`(r_f1 = frank(-f1, ties.method = "min"), r_harm = frank(-harm, ties.method = "min")), by = ds]
print(x[, .(rho_f1_harm = round(cor(f1, harm, method = "spearman"), 3),
            rho_f1_sites = round(cor(f1, sites, method = "spearman"), 3), rho_f1_pct = round(cor(f1, pct, method = "spearman"), 3),
            harm1 = config[which.max(harm)], harm1_f1 = round(f1[which.max(harm)], 4), harm1_f1rank = r_f1[which.max(harm)],
            f1_1 = config[which.max(f1)], f1_best = round(max(f1), 4), gap = round(max(f1) - f1[which.max(harm)], 4)), by = ds])
