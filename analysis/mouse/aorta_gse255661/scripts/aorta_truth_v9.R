#!/usr/bin/env Rscript
# =============================================================================
# aorta_truth_v9.R   (user 2026-11-07)
#
# Aorta ground truth under the teratoma v9 rules (human/scripts/gt_truth_v9.R), from the 12-type gate tables of
# aorta_truth_eval.R (truth/eval/{excluded,included}/truth12_sites.tsv.gz; no new pileups):
#   orientation  per position (both gates; unit coverage is keyed by site) the read A>G or read T>C gate row (multiallelic removed) with more merged
#                CTRL reads (ties: A>G); genomic_ref A if (A>G on +) or (T>C on -), else T
#   positive(t)  chosen row d_max >= t, not an MGP SNP position, passes removal of sites near other mismatch
#                types (excluded-mode table; 12 types, w1 chosen by the A>G share as before)
#   negative     chosen row informative negative (not depleted; Clopper-Pearson rule), before removal of SNP
#                positions and of sites near other mismatch types; REDIportal negatives stay unscored
#                (unscored_negative_rediportal: with only an Adar1 arm, ADAR2 sites are not depleted)
#   No germline negatives: the mice have no genotype (MGP SNPs are strain positions, not these animals).
#   Scoring: gt_score_v7.R GT_FP_RULE=ag (calls are AG in the RNA frame; FP = AG call at a negative).
# Output: truth/eval_v9/excluded/{gt_aorta_truth_sites.tsv.gz, unit_cov.tsv.gz, truth_counts.tsv}
# =============================================================================
suppressPackageStartupMessages(library(data.table))
A = file.path(Sys.getenv("SCRNA_DIR", "."), "mouse/aorta_gse255661")
EV = Sys.getenv("AORTA_EVAL_DIR", file.path(A, "truth/eval"))
OUT = Sys.getenv("AORTA_V9_DIR", file.path(A, "truth/eval_v9/excluded")); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DS = c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
K = c("gate", "chr", "pos", "strand", "type")

inc = fread(file.path(EV, "included/truth12_sites.tsv.gz"))[type %in% c("AG", "TC") & multiallelic == FALSE]
exc = fread(file.path(EV, "excluded/truth12_sites.tsv.gz"), select = c(K, "nomm_pass", "w1_chosen"))[type %in% c("AG", "TC")]
x = exc[inc, on = K]
stopifnot(!anyNA(x[snp == FALSE, nomm_pass]), all(is.na(x[snp == TRUE, nomm_pass])))
x[is.na(nomm_pass), nomm_pass := FALSE]
ch = unique(x[, .(chr, pos, strand, type, ctrl_cov)])
setorder(ch, chr, pos, -ctrl_cov, type)
ch = ch[, .SD[1], by = .(chr, pos)]
x = x[ch[, .(chr, pos, strand, type)], on = .(chr, pos, strand, type), nomatch = NULL]
stopifnot(!anyDuplicated(x[, .(gate, chr, pos)]))
x[, `:=`(site = paste0(chr, ":", pos), genomic_ref = fifelse((type == "AG") == (strand == "+"), "A", "T"),
         read_type = type, read_strand = strand)]
x[, a_to_i_strand := fifelse(genomic_ref == "A", "+", "-")]
clean = x$snp == FALSE & x$nomm_pass
for (tn in names(DS)) set(x, j = paste0("class_", tn),
  value = fifelse(clean & x$d_max >= DS[[tn]], "positive",
          fifelse(x$informative_negative & x$redi, "unscored_negative_rediportal",
          fifelse(x$informative_negative, "informative_negative",
          fifelse(clean, "unresolved", NA_character_)))))
stopifnot(!any(x$informative_negative & x$d_max >= DS[["t25"]]))
tr = x[!is.na(class_t25)]

COLS = c("gate", "site", "chr", "pos", "genomic_ref", "a_to_i_strand", "read_type", "read_strand", paste0("class_", names(DS)))
setorder(tr, gate, chr, pos)
fwrite(tr[, ..COLS], file.path(OUT, "gt_aorta_truth_sites.tsv.gz"), sep = "\t")

CM = 1:3
us = unique(tr[, c("site", "gate", "ctrl_cov", "ctrl_alt", paste0("c", CM, "_cov"), paste0("c", CM, "_alt")), with = FALSE])
us = unique(us[, !"gate"])
stopifnot(!anyDuplicated(us$site))
design = fread(file.path(A, "truth/design.tsv")); CMN = design[arm == "CTRL", sample]
uc = rbindlist(c(list(us[, .(unit = "CTRL", site, coverage = ctrl_cov, variant_cov = ctrl_alt)]),
  lapply(CM, function(j) us[, .(unit = CMN[j], site, coverage = get(paste0("c", j, "_cov")), variant_cov = get(paste0("c", j, "_alt")))])))
fwrite(uc, file.path(OUT, "unit_cov.tsv.gz"), sep = "\t")

cnt = rbindlist(lapply(names(DS), function(tn) { cl = paste0("class_", tn)
  tr[, .(t = tn, positive = sum(get(cl) == "positive"), positive_AG = sum(get(cl) == "positive" & read_type == "AG"),
         positive_TC = sum(get(cl) == "positive" & read_type == "TC"),
         informative_negative = sum(get(cl) == "informative_negative"),
         negative_mgp_snp = sum(get(cl) == "informative_negative" & snp),
         unscored_negative_rediportal = sum(get(cl) == "unscored_negative_rediportal"),
         unresolved = sum(get(cl) == "unresolved"), w1 = w1_chosen[!is.na(w1_chosen)][1]), by = gate] }))
fwrite(cnt, file.path(OUT, "truth_counts.tsv"), sep = "\t")
print(cnt)

v8 = fread(file.path(EV, "excluded/gt_aorta_truth_sites.tsv.gz"))
cat("\nv8 vs v9 (t50 positives, by gate):\n")
print(rbindlist(lapply(c("expected", "depth"), function(g) {
  a = v8[gate == g & class_t50 == "positive", site]; b = tr[gate == g & class_t50 == "positive", site]
  bag = tr[gate == g & class_t50 == "positive" & read_type == "AG", site]
  data.table(gate = g, v8 = length(a), v9 = length(b), v9_AG = length(bag), v8_in_v9AG = sum(a %in% bag),
             v8_lost = sum(!a %in% b), v8_now_TC = sum(a %in% setdiff(b, bag)))
})))
