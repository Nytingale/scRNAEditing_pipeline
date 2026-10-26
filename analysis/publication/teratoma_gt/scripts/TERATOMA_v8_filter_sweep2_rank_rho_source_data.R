#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_filter_sweep2_rank_rho_source_data.R   (2026-10-21; split from TERATOMA_v8_filter_sweep_rank_plot.R
# DESIGN=2 2026-10-25)
# Source data for response figure 6: one downstream filter at a time moved from the main value or turned off
# (coverage 5 / 10 / 15, variant reads 2 / 3 / 4, editing level 1 / 2 / 3%, removal of sites near other mismatch
# types on / off, all four off together; multi-allelic, dbSNP and JACUSA2-score removal stay on), both coverage
# gates, every configuration replayed with its own caller, BAM type, strand mode and SNP setting; scored with the
# power- and coverage-gated per-unit recall denominator (scoring_pwcov).
#   top     Spearman rho of the overall configuration ranking vs the ranking under the main filters
#   bottom  rank, under the moved filter, of the configuration ranked first under the main filters
# The first sweep (DESIGN=1) and its preliminary layout are in
# ../_superseded/pre_split_2026-10-25/scripts/TERATOMA_v8_filter_sweep_rank_plot.R.
# Input : human/pm/gt_eval_v8/filter_sweep2/filter_sweep_compare_pwcov.tsv
#         (human/pm/scripts/filter_sweep2_v8.sh -> human/scripts/gt_v8_filter_sweep_compare.R)
# Output: source_data/pm_gt/TERATOMA_v8_filter_sweep2_rank_rho_pwcov.{t50,all_t}.tsv
#         (t50 is drawn by TERATOMA_v8_filter_sweep2_rank_rho_plot.R; all_t adds t = 25% and 75%)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_filter_sweep2_rank_rho_pwcov"

SRC <- file.path(V8_DIR, "filter_sweep2", "filter_sweep_compare_pwcov.tsv")
stopifnot(file.exists(SRC))
# setting -> filter axis, x position, x label; centre = the main pipeline value of that filter
AX <- rbindlist(list(
  list("cov5",   "cov",  1, "5"),   list("centre", "cov",  2, "10"),  list("cov15", "cov",  3, "15"),
  list("alt2",   "alt",  1, "2"),   list("centre", "alt",  2, "3"),   list("alt4",  "alt",  3, "4"),
  list("centre", "edit", 1, "1%"),  list("edit2",  "edit", 2, "2%"),  list("edit3", "edit", 3, "3%"),
  list("centre", "window", 1, "on"), list("window_off", "window", 2, "off"),
  list("centre", "all", 1, "on"),    list("all_off", "all", 2, "off")))
setnames(AX, c("setting", "axis", "ord", "lab"))
AXIS_LAB <- c(cov = "Coverage\n(reads)", alt = "Variant\nreads", edit = "Editing\nlevel",
              window = "Remove sites with\nother mismatch\ntypes nearby", all = "All four\nfilters")
GATE_LAB <- sub(" \\(.*\\)$", "", V8_GATE_LAB)

all_t <- fread(SRC)[setting %in% AX$setting]
stopifnot(setequal(setdiff(AX$setting, "centre"), all_t$setting))
print(all_t[, .(min_rho = round(min(rho), 3), max_rank_of_top = max(centre_top_rank)), by = t])
cen <- unique(all_t[, .(gate, t)])[, `:=`(setting = "centre", rho = 1, centre_top_rank = 1L)]
d <- merge(AX, rbind(all_t[, .(gate, t, setting, rho, centre_top_rank)], cen, use.names = TRUE),
           by = "setting", allow.cartesian = TRUE)
d[, centre := setting == "centre"]
src_data <- copy(d)
d <- d[t == "t50"]
setorder(d, axis, ord)
d[, x := factor(paste(axis, lab), unique(paste(axis, lab)))]
d[, `:=`(Axis = factor(AXIS_LAB[axis], AXIS_LAB), Gate = factor(GATE_LAB[gate], GATE_LAB[c("expected", "depth")]))]
# rho label above the higher of the two gates at each x and below the lower (power-gated above on ties)
d[, vj := if (.N == 2 && !anyNA(rho) && rho[1] != rho[2]) fifelse(rho == max(rho), -0.8, 1.8) else fifelse(gate == "depth", 1.8, -0.8),
  by = .(x)]

write_sd(d, SRC_DIR, STEM, "t50")
write_sd(src_data[, .(gate, t, axis, value = gsub("\n", " ", lab), setting, centre, rho, centre_top_rank)],
         SRC_DIR, STEM, "all_t")
