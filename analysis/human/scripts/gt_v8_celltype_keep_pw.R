#!/usr/bin/env Rscript
# gt_v8_celltype_keep_pw.R   (2026-10-21)
# Cell-type tier under the gate-matched recall denominators of gt_score_v7.R:
#   gate      power-gated positives count in a cell type when its coverage x merged NTC editing level >= 3
#   gate_cov  as gate, and >= 10 reads in the cell type (the downstream min_cov)
# Depth-gated and REDIportal counts keep >= 10 reads (copied from celltype_covered_positives.tsv). Cell types with
# >= MIN sites under all nine truths form the tier.
# Usage : Rscript gt_v8_celltype_keep_pw.R [gate | gate_cov] [MIN (30)]
# Input : gt_eval_v8/truth/gt_v7_truth_sites.tsv.gz, gt_eval_v8/unit_cov/{ko_truth_unit_cov.tsv.gz,
#         celltype_covered_positives.tsv}
# Output: gt_eval_v8/unit_cov/celltype_covered_positives_<pw | pwcov>.tsv, celltype_keep_<pw | pwcov>.txt
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
RULE <- if (length(args) >= 1) args[1] else "gate"
MIN <- if (length(args) >= 2) as.integer(args[2]) else 30L
stopifnot(RULE %in% c("gate", "gate_cov"))
TAG <- c(gate = "pw", gate_cov = "pwcov")[[RULE]]
D <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/gt_eval_v8")
MIN_EXP <- 3
MIN_COV <- 10
tr <- fread(file.path(D, "truth", "gt_v7_truth_sites.tsv.gz"))
uc <- fread(file.path(D, "unit_cov", "ko_truth_unit_cov.tsv.gz"))
old <- fread(file.path(D, "unit_cov", "celltype_covered_positives.tsv"))
ed <- uc[unit == "NTC_ctrl" & coverage > 0, .(site, edit = variant_cov / coverage)]
stopifnot(all(unique(tr$site) %chin% ed$site))
ct <- merge(uc[startsWith(unit, "NTC_ctrl__")], ed, by = "site")[, unit := sub("^NTC_ctrl__", "", unit)]
stopifnot(setequal(ct$unit, old$unit))

cnt <- function(g, rule) rbindlist(lapply(c("t25", "t50", "t75"), function(tn) {
  pos <- tr[gate == g & get(paste0("class_", tn)) == "positive", site]
  ct[, .(n = sum(site %chin% pos & rule(coverage, edit))), by = unit][, truth := paste0(g, "_", tn)]
}))
# the >= 10-read counts must reproduce the existing table before its power-gated columns are replaced
chk <- dcast(rbind(cnt("expected", function(c, e) c >= MIN_COV), cnt("depth", function(c, e) c >= MIN_COV)),
             unit ~ truth, value.var = "n")
stopifnot(all(merge(old, chk, by = "unit", suffixes = c("", ".y"))[, sapply(names(chk)[-1], function(k) all(get(k) == get(paste0(k, ".y"))))]))

pw_rule <- if (RULE == "gate") function(c, e) c * e >= MIN_EXP else function(c, e) c * e >= MIN_EXP & c >= MIN_COV
pw <- dcast(cnt("expected", pw_rule), unit ~ truth, value.var = "n")
new <- copy(old)
for (k in names(pw)[-1]) new[pw, (k) := get(paste0("i.", k)), on = "unit"]
TR <- c("expected_t25", "expected_t50", "expected_t75", "depth_t25", "depth_t50", "depth_t75",
        "redi_cov10_alt1", "redi_alt3", "redi_alt3_e05")
new[, min_all := do.call(pmin, .SD), .SDcols = TR]
setorder(new, -min_all)
fwrite(new, file.path(D, "unit_cov", sprintf("celltype_covered_positives_%s.tsv", TAG)), sep = "\t")
keep <- paste0("NTC_ctrl__", new[min_all >= MIN, unit])
writeLines(keep, file.path(D, "unit_cov", sprintf("celltype_keep_%s.txt", TAG)))
print(new[, c("unit", TR, "min_all"), with = FALSE])
s <- sort(new$min_all)
cat(RULE, ": ", length(keep), " cell types with >= ", MIN, " sites under all nine truths; largest excluded ",
    max(s[s < MIN]), ", smallest kept ", min(s[s >= MIN]), "\n", sep = "")
