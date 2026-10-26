#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_metric_case_source_data.R   (2026-10-21; split from the _plot.R script 2026-10-25)
# Source data for response figure 9 (why the A-to-G composite, R1 major 1): the two truth-free measures disagree,
# and the number of sites tracks F1 only because F1 is mostly recall under the standard filters.
# 58 configurations, standard downstream filters (published call sets, headline scoring). A-to-G sites /
# percentage = merged NTC arm calls; mean F1 / precision / recall = mean over the 6 KO-based truths of the
# cross-tier means (human/scripts/gt_v8_metric_case_data.R). Composite = mean of the two ranks.
#   a  rank by number of A-to-G sites vs rank by A-to-G percentage (a_ranks)
#   b  mean recall and mean precision vs mean F1 (a_ranks)
#   c  Spearman's rho of each measure with mean F1 and with the KO response of the whole call set (c_rho)
# all_configs adds the configurations without the >= 3 variant reads minimum (filter sweep 2, setting alt1);
# the layout with those rows (WITH_RELAXED=1) is in ../_superseded/pre_split_2026-10-25/scripts/.
# Input : human/pm/gt_eval_v8/metric_case/metric_case.tsv,
#         human/pm/gt_eval_v8/callset_ko/callset_ko_response.tsv (state mismatch; same call sets)
# Output: source_data/pm_gt/TERATOMA_v8_metric_case.{a_ranks,c_rho,all_configs}.tsv
#         (drawn by TERATOMA_v8_metric_case_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_metric_case"
SHIFT <- 20
SRC <- file.path(V8_DIR, "metric_case", "metric_case.tsv")
SRC_KO <- file.path(V8_DIR, "callset_ko", "callset_ko_response.tsv")

SET <- c(centre = "Main analysis", alt1 = "Relaxed")

d <- fread(SRC)[setting %in% names(SET)]
stopifnot(d[, uniqueN(config), by = setting][, all(V1 == 58L)])
d[, c("bam_type", "strand_mode", "snp_subtract") := tstrsplit(config, "__")[2:4]]
d[, `:=`(r_ag = frank(-ag), r_pct = frank(-pct)), by = setting]
d[, comp := frank((r_ag + r_pct) / 2), by = setting]
d[, `:=`(col = factor(SET[setting], SET), method = caller_factor(caller), proc = relabel(bam_type, label_bamType))]
PROC_LEVELS <- intersect(label_bamType[c("sampleRmDup", "sampleUmiDup", "recal", "recalq20", "recaldefault")], d$proc)
d[, proc := factor(proc, PROC_LEVELS)]
rho <- function(x, y) cor(x, y, method = "spearman")
print(d[, .(sites = rho(ag, mean_f1), pct = rho(pct, mean_f1), comp = rho(-comp, mean_f1),
            best = max(mean_f1)), by = setting])

s0 <- d[setting == "centre"]
s0[, shift := r_pct - r_ag]
cat("standard: rank sites vs rank pct rho", round(rho(s0$r_ag, s0$r_pct), 3), "; |shift| >= ", SHIFT, ":",
    s0[abs(shift) >= SHIFT, .N], "\n")
print(s0[abs(shift) >= SHIFT][order(-shift), .(config, r_ag, r_pct, shift)])
write_sd(s0, SRC_DIR, STEM, "a_ranks")

ko <- fread(SRC_KO)[state == "mismatch", .(config, n_ko = n_calls, ko = resp_share_all)]
s0 <- merge(s0, ko, by = "config")
stopifnot(nrow(s0) == 58L, s0[, all(n == n_ko)])
MEAS <- c(pct = "A-to-G\npercentage", ag = "Number of\nA-to-G sites", comp = "A-to-G\ncomposite rank")
TGT <- c(f1 = "With mean F1", ko = "With KO response")
rc <- rbindlist(lapply(c("f1", "ko"), function(tg) {
  y <- if (tg == "f1") s0$mean_f1 else s0$ko
  data.table(target = TGT[[tg]], measure = MEAS[c("ag", "pct", "comp")],
             rho = c(rho(s0$ag, y), rho(s0$pct, y), rho(-s0$comp, y)))
}))
rc[, `:=`(target = factor(target, TGT), measure = factor(measure, MEAS))]
print(rc)
write_sd(rc, SRC_DIR, STEM, "c_rho")
write_sd(d[, .(setting, config, caller, bam_type, strand_mode, snp_subtract, ag, pct, r_ag, r_pct, comp,
               mean_f1, precision, recall)], SRC_DIR, STEM, "all_configs")
