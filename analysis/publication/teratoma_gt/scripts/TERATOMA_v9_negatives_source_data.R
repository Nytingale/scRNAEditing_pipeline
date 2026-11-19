#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v9_negatives_source_data.R   (2026-11-10; v9 copy of TERATOMA_v8_negatives_source_data.R)
# Source data for response figure 3 under the v9 truth: pooled NTC editing level of the informative negatives vs
# the positives (t = 50%), per coverage gate, and the share of each class in REDIportal. Positives are read-orientation
# A-to-G or T-to-C rows (gt_truth_v9.R); an informative negative (neg_source infneg, incl. infneg+pgp1_ag) is keyed
# by position, and its editing level is that of the A-to-G / T-to-C row on its genomic base with more NTC reads
# (as the orientation choice of the truth build). Germline negatives (pgp1_ag / pgp1_hom only) are not drawn.
# Input : human/pm/gt_eval_v9/truth/gt_v9_truth_sites.tsv.gz,
#         human/pm/gt_eval_v6/spectrum12/{eval,eval_depth}/truth12_sites{,_rediportal}.tsv.gz
# Output: source_data/pm_gt/TERATOMA_v9_negatives_ntc_editing.{sites,labels}.tsv (TERATOMA_v9_negatives_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v9_negatives_ntc_editing"
PM <- file.path(REPO_DIR, "human", "pm")
SPEC <- c(expected = "eval", depth = "eval_depth")

v9 <- fread(file.path(PM, "gt_eval_v9", "truth", "gt_v9_truth_sites.tsv.gz"))
d <- rbindlist(lapply(names(SPEC), function(g) {
  dir <- file.path(PM, "gt_eval_v6", "spectrum12", SPEC[[g]])
  s <- fread(file.path(dir, "truth12_sites.tsv.gz"),
             select = c("type", "chr", "pos", "strand", "ntc_cov", "ntc_edit"))[type %in% c("AG", "TC")]
  s[, genomic_ref := fifelse((type == "AG") == (strand == "+"), "A", "T")]
  setorder(s, chr, pos, -ntc_cov)
  s1 <- s[, .SD[1], by = .(chr, pos, genomic_ref)]
  redi <- unique(fread(file.path(dir, "truth12_sites_rediportal.tsv.gz")))
  v <- v9[gate == g & (class_t50 == "positive" | grepl("infneg", neg_source))]
  pos <- s[v[class_t50 == "positive"], on = .(type = read_type, chr, pos, strand = read_strand), nomatch = NA,
           .(class = "positive_t50", site, chr, pos, edit = ntc_edit)]
  neg <- s1[v[class_t50 != "positive"], on = .(chr, pos, genomic_ref), nomatch = NA,
            .(class = "informative_negative", site, chr, pos, edit = ntc_edit)]
  x <- rbind(pos, neg)
  stopifnot(!anyNA(x$edit), nrow(x) == nrow(v))
  x <- redi[x, on = .(chr, pos)]
  x[, `:=`(gate = g, in_redi = in_rediportal %in% TRUE)]
}))
CLASS_LAB <- c(positive_t50 = "Positives (KO depletion \u2265 50%)", informative_negative = "Informative negatives")
d[, `:=`(Class = factor(CLASS_LAB[class], CLASS_LAB), Gate = factor(V8_GATE_LAB[gate], V8_GATE_LAB))]
stopifnot(!anyNA(d$Gate))
lab <- d[, .(n = .N, med = median(edit), redi = mean(in_redi)), by = .(Gate, Class)]
lab[, txt := sprintf("n = %s, median %.2f, %.0f%% are REDIportal sites", format(n, big.mark = ",", trim = TRUE), med, 100 * redi)]
print(lab)

write_sd(d[, .(Gate, Class, site, edit, in_redi)], SRC_DIR, STEM, "sites")
write_sd(lab, SRC_DIR, STEM, "labels")
