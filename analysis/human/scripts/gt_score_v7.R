#!/usr/bin/env Rscript
# =============================================================================
# gt_score_v7.R   (2026-10-08)
#
# Score every method configuration against the v7 A-to-I position truth
# (gt_eval_v7/truth/gt_v7_truth_sites.tsv.gz; gt_truth_v7.R). NTC calls only.
#
#   call            one row of <unit>/04_filtered/<caller>/<bam_type>/<strand_mode>/snp_<x>/*_pre100.tsv.gz,
#                   any varType; matched to the truth on chr:pos. varType is in the RNA frame of the
#                   assigned strand, so the genomic substitution is varType on + and its complement on -.
#   correct         genomic A>G on + or genomic T>C on - (= RNA-frame AG on either strand)
#   wrong           any other call, including + strand T>C and - strand A>G
#   TP              positives(t) with a correct call
#   FP              positives(t) with a wrong call (a miss as well, unless the site also has a correct call)
#                   + informative negatives with any call
#   recall          TP / positives(t) (global denominator); with GT_UNIT_COV, TP at positives covered by
#                   >= GT_UNIT_MIN_COV (10) stranded reads in the unit / those covered positives (n_pos_unit);
#                   GT_UNIT_RULE=gate: covered = unit coverage x merged NTC editing level >= GT_UNIT_MIN_EXP (3)
#                   for the expected (power) gate, >= GT_UNIT_MIN_COV reads for every other gate; GT_UNIT_RULE=gate_cov:
#                   as gate, and expected-gate positives also need >= GT_UNIT_MIN_COV reads (the downstream min_cov)
#   precision       TP / (TP + FP); NA if 0
#   F1              2PR / (P + R) = 2TP / (2TP + FP + FN); 0 when the unit has positives but no TP (including no
#                   scored calls, P = NA); NA only when the unit has no positives (recall NA)
#   unit            one NTC pseudobulk in which calls are made: arm = all NTC cells (T1-T4 merged), sample = the NTC
#                   cells of one teratoma, celltype = the NTC cells of one cell type (pooled over teratomas)
#   tier value      mean over units (precision ignores NA units; F1 and recall ignore units without positives)
# Configurations, units and the reported subset are as gt_score_v6.R.
#
# Usage: [GT_V7_TRUTH=<truth_sites.tsv.gz>] Rscript gt_score_v7.R <tier> [results_dir] [outdir] [call_file_regex]
#        GT_V7_TRUTH defaults to gt_eval_v7/truth/gt_v7_truth_sites.tsv.gz (truth variants: gt_truth_v7_variants.R)
#        GT_UNIT_COV=<unit site coverage table> (gt_v8_unit_cov.sh) switches recall to the per-unit denominator
#        GT_UNIT_RULE=cov (default) | gate | gate_cov: per-unit denominator rule (above); gate* take the editing level
#        from the GT_UNIT_EDIT_UNIT (NTC_ctrl) rows of GT_UNIT_COV, the merged-NTC pileup of the truth build
#        GT_UNIT_KEEP=<file, one cell-type unit per line> restricts the celltype tier to those units
# Outputs (<outdir>, default gt_eval_v7/scoring/<tier>/):
#   gt_v7_unit_scores.tsv      unit x configuration x gate x t
#   gt_v7_scores.tsv           configuration x gate x t (tier summary; all BAM types)
#   gt_v7_scores_reported.tsv  the same, reported configurations only
#   gt_v7_calls_NTC.tsv.gz     calls at truth positions: unit, config, gate, site, varType, strand, correct
# =============================================================================
suppressPackageStartupMessages(library(data.table))
args = commandArgs(trailingOnly = TRUE)
TIER = args[1]
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
DEF = list(arm = "results_gt_v6", sample = "results_gt_persample_v6", celltype = "results_gt_celltype_v6")
UNIT_RE = list(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")
stopifnot(TIER %in% names(DEF))
RES = if (length(args) >= 2) args[2] else file.path(PM, DEF[[TIER]])
OUT = if (length(args) >= 3) args[3] else file.path(PM, "gt_eval_v7", "scoring", TIER)
CALL_RE = if (length(args) >= 4) args[4] else "_pre100\\.tsv\\.gz$"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DS = c("t25", "t50", "t75")
COMBO = c("caller", "bam_type", "strand_mode", "snp_subtract")
UNREPORTED_BAM = "sample"
UNREPORTED_CALLER_BAM = list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))
unreported = function(caller, bam_type) bam_type %in% UNREPORTED_BAM |
  mapply(function(c, b) b %in% UNREPORTED_CALLER_BAM[[c]], caller, bam_type) &
  caller %in% names(UNREPORTED_CALLER_BAM)
COMP = c(A = "T", C = "G", G = "C", T = "A")
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

TRUTH = Sys.getenv("GT_V7_TRUTH", file.path(PM, "gt_eval_v7", "truth", "gt_v7_truth_sites.tsv.gz"))
msg("truth table: ", TRUTH)
truth = fread(TRUTH, select = c("gate", "site", "genomic_ref", paste0("class_", DS)))
GATES = unique(truth$gate)
n_pos = truth[, lapply(.SD, function(v) sum(v == "positive")), by = gate, .SDcols = paste0("class_", DS)]
msg("truth: ", paste(n_pos[, sprintf("%s t25/t50/t75 = %d/%d/%d", gate, class_t25, class_t50, class_t75)], collapse = "; "))
tsites = unique(truth$site)
UNIT_COV = Sys.getenv("GT_UNIT_COV", "")
MIN_COV = as.integer(Sys.getenv("GT_UNIT_MIN_COV", "10"))
UNIT_RULE = Sys.getenv("GT_UNIT_RULE", "cov")
MIN_EXP = as.numeric(Sys.getenv("GT_UNIT_MIN_EXP", "3"))
EDIT_UNIT = Sys.getenv("GT_UNIT_EDIT_UNIT", "NTC_ctrl")
stopifnot(UNIT_RULE %in% c("cov", "gate", "gate_cov"), UNIT_RULE == "cov" || nzchar(UNIT_COV))

units = grep(UNIT_RE[[TIER]], list.dirs(RES, recursive = FALSE, full.names = FALSE), value = TRUE)
if (nzchar(Sys.getenv("GT_UNIT_KEEP")) && TIER == "celltype") {
  keep = readLines(Sys.getenv("GT_UNIT_KEEP"))
  stopifnot(all(keep %in% units))
  units = intersect(units, keep)
}
if (!length(units)) stop("no NTC units matching ", UNIT_RE[[TIER]], " in ", RES)
files = rbindlist(lapply(units, function(u) {
  f = list.files(file.path(RES, u, "04_filtered"), CALL_RE, recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NULL)
  p = tstrsplit(sub(paste0("^", file.path(RES, u, "04_filtered"), "/"), "", f), "/", fixed = TRUE)
  data.table(unit = u, caller = p[[1]], bam_type = p[[2]], strand_mode = p[[3]],
             snp_subtract = sub("^snp_", "", p[[4]]), file = f)
}))
files[, config := do.call(paste, c(.SD, sep = "__")), .SDcols = COMBO]
msg("tier ", TIER, ": ", length(units), " unit(s), ", nrow(files), " call files, ", uniqueN(files$config), " configurations")

calls = rbindlist(lapply(seq_len(nrow(files)), function(i) {
  if (file.size(files$file[i]) == 0) return(NULL)
  x = fread(files$file[i], select = c("site", "varType", "strand"), showProgress = FALSE)
  x = unique(x[site %chin% tsites])
  if (!nrow(x)) return(NULL)
  x[, `:=`(unit = files$unit[i], config = files$config[i])]
}))
stopifnot(all(calls$strand %in% c("+", "-")), all(nchar(calls$varType) == 2L))
calls[, genomic_sub := fifelse(strand == "+", varType,
                               paste0(COMP[substr(varType, 1, 1)], COMP[substr(varType, 2, 2)]))]
calls[, correct := (genomic_sub == "AG" & strand == "+") | (genomic_sub == "TC" & strand == "-")]
stopifnot(identical(calls$correct, calls$varType == "AG"))
msg("calls at truth positions: ", nrow(calls))

cg = merge(calls, truth, by = "site", allow.cartesian = TRUE)
# a correct call is only possible on the truth position's genomic base
stopifnot(cg[correct == TRUE, all(substr(genomic_sub, 1, 1) == genomic_ref)])
fwrite(cg[, .(unit, config, gate, site, genomic_ref, varType, strand, genomic_sub, correct)],
       file.path(OUT, "gt_v7_calls_NTC.tsv.gz"), sep = "\t")
site_call = cg[, .(any_correct = any(correct), any_wrong = any(!correct),
                   wrong_TC = any(!correct & varType == "TC")), by = c("unit", "config", "gate", "site", paste0("class_", DS))]

if (nzchar(UNIT_COV)) {
  uc = fread(UNIT_COV, select = c("unit", "site", "coverage", "variant_cov"))
  if (!all(units %chin% uc$unit)) stop("units missing from ", UNIT_COV, ": ", paste(setdiff(units, uc$unit), collapse = ", "))
  if (UNIT_RULE != "cov") {
    ed = uc[unit == EDIT_UNIT & coverage > 0, .(site, edit = variant_cov / coverage)]
    if (!all(tsites %chin% ed$site)) stop(sum(!tsites %chin% ed$site), " truth sites without a ", EDIT_UNIT, " editing level")
    uc = merge(uc[unit %chin% units], ed, by = "site")
    ucov = rbindlist(lapply(GATES, function(g)
      uc[if (g != "expected") coverage >= MIN_COV
         else coverage * edit >= MIN_EXP & (UNIT_RULE == "gate" | coverage >= MIN_COV), .(unit, site, gate = g)]))
    msg("per-unit recall denominator by gate: expected, unit coverage x ", EDIT_UNIT, " editing level >= ", MIN_EXP,
        if (UNIT_RULE == "gate_cov") paste0(" and >= ", MIN_COV, " stranded reads"),
        "; depth, >= ", MIN_COV, " stranded reads in the unit (", UNIT_COV, ")")
  } else {
    ucov = unique(uc[unit %chin% units & coverage >= MIN_COV, .(unit, site)])[, .(gate = GATES), by = .(unit, site)]
    msg("per-unit recall denominator: positives with >= ", MIN_COV, " stranded reads in the unit (", UNIT_COV, ")")
  }
} else {
  ucov = CJ(unit = units, site = tsites, gate = GATES)
}
site_call[, covered := FALSE][ucov, covered := TRUE, on = c("unit", "site", "gate")]
tc = merge(truth, ucov, by = c("site", "gate"), allow.cartesian = TRUE)
npu = rbindlist(lapply(DS, function(tn)
  tc[, .(n_pos_unit = sum(get(paste0("class_", tn)) == "positive")), by = .(unit, gate)][, t := tn]))

grid = CJ(unit = units, config = unique(files$config), gate = GATES, t = DS)
grid = merge(grid, unique(files[, c("unit", "config", COMBO), with = FALSE]), by = c("unit", "config"))
us = rbindlist(lapply(DS, function(tn) {
  cl = paste0("class_", tn)
  s = site_call[, .(n_tp = sum(get(cl) == "positive" & any_correct),
                    n_tp_cov = sum(get(cl) == "positive" & any_correct & covered),
                    n_fp_wrong = sum(get(cl) == "positive" & any_wrong),
                    n_fp_wrong_tc = sum(get(cl) == "positive" & wrong_TC),
                    n_wrong_only = sum(get(cl) == "positive" & any_wrong & !any_correct),
                    n_fp_neg = sum(get(cl) == "informative_negative"),
                    n_unresolved = sum(get(cl) == "unresolved")), by = .(unit, config, gate)]
  s[, t := tn]
}))
us = merge(grid, us, by = c("unit", "config", "gate", "t"), all.x = TRUE)
cnt_cols = c("n_tp", "n_tp_cov", "n_fp_wrong", "n_fp_wrong_tc", "n_wrong_only", "n_fp_neg", "n_unresolved")
for (k in cnt_cols) set(us, which(is.na(us[[k]])), k, 0L)
npl = melt(n_pos, id.vars = "gate", variable.name = "t", value.name = "n_pos")[, t := sub("^class_", "", t)]
us = merge(us, npl, by = c("gate", "t"))
us = merge(us, npu, by = c("unit", "gate", "t"), all.x = TRUE)
us[is.na(n_pos_unit), n_pos_unit := 0L]
us[, n_fp := n_fp_wrong + n_fp_neg]
us[, `:=`(precision = fifelse(n_tp + n_fp > 0, n_tp / (n_tp + n_fp), NA_real_),
          recall = fifelse(n_pos_unit > 0, n_tp_cov / n_pos_unit, NA_real_))]
us[, f1 := fifelse(is.na(recall), NA_real_,
                   fifelse(is.na(precision) | precision + recall == 0, 0, 2 * precision * recall / (precision + recall)))]
setorder(us, gate, t, config, unit)
fwrite(us, file.path(OUT, "gt_v7_unit_scores.tsv"), sep = "\t")

se = function(x) { x = x[!is.na(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }
sdv = function(x) { x = x[!is.na(x)]; if (length(x) > 1) sd(x) else NA_real_ }
spread = if (TIER == "celltype") sdv else se
mn = function(x) { x = x[!is.na(x)]; if (length(x)) mean(x) else NA_real_ }
summ = us[, c(list(n_units = .N), lapply(.SD, mean),
              list(precision = mn(precision), recall = mn(recall), f1 = mn(f1),
                   precision_err = spread(precision), recall_err = spread(recall), f1_err = spread(f1))),
          by = c("gate", "t", "n_pos", "config", COMBO), .SDcols = c(cnt_cols, "n_fp", "n_pos_unit")]
summ[, err_type := fifelse(TIER == "celltype", "SD across cell types",
                    fifelse(TIER == "sample", "SE across teratomas", "none (1 unit)"))]
summ[, reported := !unreported(caller, bam_type)]
setorder(summ, gate, t, -f1)
fwrite(summ, file.path(OUT, "gt_v7_scores.tsv"), sep = "\t")
fwrite(summ[reported == TRUE], file.path(OUT, "gt_v7_scores_reported.tsv"), sep = "\t")
msg("wrote ", summ[gate == GATES[1] & t == "t50", .N], " configurations x ", length(GATES), " gates x ", length(DS),
    " thresholds (", summ[reported == TRUE & gate == GATES[1] & t == "t50", .N], " reported) -> ", OUT)
print(head(summ[reported == TRUE & gate == "expected" & t == "t50",
                .(config, n_tp = round(n_tp, 1), n_fp_wrong = round(n_fp_wrong, 1), n_fp_neg = round(n_fp_neg, 1),
                  precision = round(precision, 3), recall = round(recall, 3), f1 = round(f1, 3))], 12))
