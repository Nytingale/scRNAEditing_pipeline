#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_dedup_editlevel_pearson_source_data.R   (2026-11-19)
#
# Pearson's r (and Spearman's rho) of candidate-site A-to-G editing levels between two de-duplication methods,
# per method pair x tier x arm, from the same site points as the hex panels (GTDEF_hex_dedup_coord_vs_umi).
# p values: two-sided, asymptotic t approximation (as cor.test(exact = FALSE)).
#
# INPUT : source_data/pm_gt/GTDEF_hex_dedup_coord_vs_umi.el_points.tsv.gz (GTDEF_hex_dedup_coord_vs_umi_source_data.R)
# OUTPUT: source_data/pm_gt/GTDEF_dedup_editlevel_pearson.cor.tsv
#         (drawn by GTDEF_dedup_editlevel_pearson_plot.R: response figure 15)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_dedup_editlevel_pearson"

pts <- fread(file.path(SRC_DIR, "GTDEF_hex_dedup_coord_vs_umi.el_points.tsv.gz"))
cor_p <- function(r, n) {
  t <- r * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.xmin))
  2 * pt(-abs(t), df = n - 2)
}
cr <- pts[, .(n = .N, r = cor(edit_x, edit_y), rho = cor(edit_x, edit_y, method = "spearman")),
          by = .(pair, tier, arm)]
cr[, `:=`(p_r = cor_p(r, n), p_rho = cor_p(rho, n))]
full <- CJ(pair = unique(cr$pair), tier = unique(cr$tier), arm = unique(cr$arm))
cr <- merge(full, cr, by = c("pair", "tier", "arm"), all.x = TRUE)
print(cr[, .(pair = gsub("\n", " ", pair), tier, arm, n, r = round(r, 4), rho = round(rho, 4))])
write_sd(cr, SRC_DIR, STEM, "cor")
