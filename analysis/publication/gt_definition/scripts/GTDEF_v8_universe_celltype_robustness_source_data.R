#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_universe_celltype_robustness_source_data.R   (split from GTDEF_gtdefinition_plot.R 2026-10-25)
# Source data for response figure 2: pileup quality, site universe, truth sites per cell type and the ranking
# under one-step threshold moves.
#   a_spectrum     NTC: share of all mismatches by type at BQ/MQ 20 and 30 (every aligned base)
#   b_ntc_vs_ko    ADAR1 KO / NTC rate per mismatch type, BQ/MQ 20 and 30
#   c_af_hist      NTC A-to-G allele fraction of candidates in / not in dbSNP common (c_groups: group labels)
#   d_ko_response  ADAR1 KO / NTC A-to-G fraction pooled per NTC allele-fraction bin, Wilson 95% CI
#   e_celltype     truth sites in each cell type's recall denominator, per truth set
#   f_rank_rho     Spearman's rho with the chosen-truth ranking and rank of its top configuration, one threshold
#                  moved one step (t = 50%), both gates
# INPUT : report/pm_gt/GTDEF_step1_{spectrum,ntc_vs_ko}.tsv, GTDEF_step2_{summary,af_hist,ko_response}.tsv
#         (GTDEF_prep_source_data.R, from the merged-BAM pileups);
#         publication/teratoma_gt/source_data/pm_gt/TERATOMA_v8_dot_celltype_truth_sites_min<N>.source_data.tsv
#         (TERATOMA_v8_celltype_truth_sites_dot_plot.R), TERATOMA_v8_sweep_rank_rho.source_data.tsv
#         (TERATOMA_v8_sweep_rank_plot.R)
# OUTPUT: source_data/pm_gt/GTDEF_v8_universe_celltype_robustness.<part>.tsv
#         (drawn by GTDEF_v8_universe_celltype_robustness_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(tibble) })
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_v8_universe_celltype_robustness"
rp <- function(stem) file.path(REPORT_DIR, paste0(stem, ".tsv"))
rd <- function(stem) as_tibble(fread(rp(stem), sep = "\t", showProgress = FALSE))
count_label <- function(n) {
  out = formatC(n, format = "d", big.mark = ",")
  big = !is.na(n) & n >= 1000
  out[big] = vapply(n[big], shorten_num, character(1))
  out
}

# a-b: pileup quality, every aligned base of merged NTC and ADAR1 KO
q_lab <- function(q) paste0("BQ/MQ ", q)
Q_LEVELS <- q_lab(c(20, 30))
stopifnot(setequal(Q_LEVELS, names(cmap_bq_mq)))
s_lab <- c(merged_NTC_ctrl = "NTC", merged_ADAR1_KO = "ADAR1 KO")
ty <- rd("GTDEF_step1_spectrum")
ty_levels <- sort(unique(ty$type))
ty <- ty %>% mutate(q = factor(q_lab(bq_mq), levels = Q_LEVELS), type = factor(type, levels = ty_levels),
                    key = gsub(">", "", type), sample = factor(s_lab[sample], levels = s_lab))
write_sd(ty %>% filter(sample == "NTC"), SRC_DIR, STEM, "a_spectrum")
nk <- rd("GTDEF_step1_ntc_vs_ko") %>%
  mutate(q = factor(q_lab(bq_mq), levels = Q_LEVELS), type = factor(type, levels = ty_levels))
write_sd(nk, SRC_DIR, STEM, "b_ntc_vs_ko")

# c-d: site universe, every candidate with NTC >= 10 reads; sites in dbSNP = dbSNP156 common
sm <- rd("GTDEF_step2_summary")
grp_lab <- sm %>% transmute(group, lab = paste0(group, " (", count_label(n_sites), ")"))
relab_g <- function(g) setNames(grp_lab$lab, grp_lab$group)[g]
write_sd(grp_lab, SRC_DIR, STEM, "c_groups")
write_sd(rd("GTDEF_step2_af_hist") %>% mutate(g = factor(relab_g(group), levels = grp_lab$lab)), SRC_DIR, STEM, "c_af_hist")
bin_lab <- function(b) { x = as.numeric(regmatches(b, gregexpr("[0-9.]+", b))[[1]]) * 100
  paste0(round(x[1]), "-", pmin(100, round(x[2])), "%") }
kr <- rd("GTDEF_step2_ko_response") %>%
  mutate(g = factor(relab_g(group), levels = grp_lab$lab), bl = map_chr(af_bin, bin_lab),
         bl = factor(bl, levels = unique(bl[order(ntc_edit)])))
write_sd(kr, SRC_DIR, STEM, "d_ko_response")

# e: truth sites per cell type, one value per truth set; cell types ordered by their minimum over the nine truths
TSD <- file.path(REPO_DIR, "publication", "teratoma_gt", "source_data", "pm_gt")
DOT_SD <- file.path(TSD, sprintf("TERATOMA_v8_dot_celltype_truth_sites_min%d.source_data.tsv", V8_CT_MIN_SITES))
dn <- fread(DOT_SD)
stopifnot(setequal(dn$truth, V8_TRUTH_IDS), uniqueN(dn$unit) == 23,
          sum(unique(dn[, .(unit, kept)])$kept) == V8_N_CELLTYPES)
lev_ct <- unique(dn[order(-min_all, unit), unit])
SHORT_LAB <- setNames(c(sprintf("\u2265 %d%% reduction in KO, power-gated", c(25, 50, 75)),
                        sprintf("\u2265 %d%% reduction in KO, depth-gated", c(25, 50, 75)),
                        "REDIportal, \u2265 1 variant read", "REDIportal, \u2265 3 variant reads",
                        "REDIportal, \u2265 3 variant reads, \u2265 5% editing level"), V8_TRUTH_IDS)
dn[, `:=`(Truth = factor(SHORT_LAB[truth], levels = SHORT_LAB[V8_TRUTH_IDS]),
          unit = factor(unit, levels = lev_ct), y = log10(n_sites + 1))]
ct_min <- unique(dn[, .(unit, kept)])
stopifnot(all(ct_min[order(unit)]$kept == (seq_len(23) <= sum(ct_min$kept))))
write_sd(dn, SRC_DIR, STEM, "e_celltype")

# f-g: ranking under one-step threshold moves, t = 50%
sw <- fread(file.path(TSD, "TERATOMA_v8_sweep_rank_rho.source_data.tsv"))[t == "t50"]
SW_AX <- c(min_cov = "NTC coverage", min_edit_reads = "NTC edited reads", min_ctrl_edit = "NTC editing level",
           min_cov_teratoma = "Per-sample NTC coverage", rule = "Reproducibility across samples",
           ko_min_expected = "Expected edited KO reads (m)", ko_min_cov = "KO coverage")
sw[, value := as.character(value)]
sw_num <- suppressWarnings(as.numeric(sw$value))
sw[, ord := fifelse(parameter == "rule", c(strict2 = 1, evaluable = 2, strict3 = 3)[value], sw_num)]
sw[, lab := fcase(parameter == "min_ctrl_edit", paste0(round(100 * sw_num), "%"),
                  parameter == "rule" & value == "strict2", "\u22652\nof 4",
                  parameter == "rule" & value == "strict3", "\u22653\nof 4",
                  parameter == "rule", "all\nbut 1", default = value)]
setorder(sw, parameter, ord)
sw[, x := factor(paste(parameter, lab), unique(paste(parameter, lab)))]
sw[, `:=`(Axis = factor(SW_AX[parameter], SW_AX), g = factor(V8_GATE_LAB[gate], levels = V8_GATE_LAB))]
write_sd(sw, SRC_DIR, STEM, "f_rank_rho")
