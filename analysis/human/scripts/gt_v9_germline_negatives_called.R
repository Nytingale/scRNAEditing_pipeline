#!/usr/bin/env Rscript
# =============================================================================
# gt_v9_germline_negatives_called.R   (user 2026-11-15: are the 1.31 million germline negatives covered?)
# The germline negatives of the v9 KO truth are all PGP1 PASS SNVs that rule out editing (gt_truth_v9.R), with
# no coverage requirement. Counts the negatives (t = 50%) that receive an A-to-G call in the NTC arm from at
# least one of the 76 ranked configurations, and the calls per configuration, by negative kind.
# Input : human/pm/gt_eval_v9/truth/gt_v9_truth_sites.tsv.gz
#         human/pm/gt_eval_v9/scoring_pwcov76/none/{arm/gt_v7_calls_NTC.tsv.gz,config_rank.tsv}
# Output: human/pm/gt_eval_v9/truth/germline_negatives_called.tsv (+ printed summary)
# Usage : srun -p slice --mem=32G Rscript human/scripts/gt_v9_germline_negatives_called.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
V9 <- "human/pm/gt_eval_v9"
tr <- fread(file.path(V9, "truth/gt_v9_truth_sites.tsv.gz"), select = c("gate", "site", "class_t50", "neg_source"))
neg <- tr[class_t50 == "informative_negative"]
neg[, kind := fifelse(neg_source == "infneg", "informative", fifelse(grepl("pgp1", neg_source) & !grepl("infneg", neg_source), "germline", "both"))]
print(neg[, .N, by = .(gate, kind)])
cl <- fread(file.path(V9, "scoring_pwcov76/none/arm/gt_v7_calls_NTC.tsv.gz"))
cl <- cl[config %chin% unique(fread(file.path(V9, "scoring_pwcov76/none/config_rank.tsv"))$config)]
cat("configs:", uniqueN(cl$config), " units:", paste(unique(cl$unit), collapse = ","), "\n")
x <- merge(cl[varType == "AG"], neg[, .(gate, site, kind, neg_source)], by = c("gate", "site"))
cat("\nDistinct negative sites with an A-to-G call from >= 1 of the configurations:\n")
print(x[, .(sites_called = uniqueN(site)), by = .(gate, kind)])
cat("\nPer configuration, calls at negatives (median, range):\n")
pc <- x[, .(n = uniqueN(site)), by = .(gate, kind, config)]
pc <- merge(CJ(gate = unique(cl$gate), kind = unique(pc$kind), config = unique(cl$config)), pc, all.x = TRUE)[is.na(n), n := 0L]
print(pc[, .(median = median(n), min = min(n), max = max(n)), by = .(gate, kind)])
cat("\nGermline negatives called, by PGP1 class:\n")
print(x[kind != "informative", .(sites_called = uniqueN(site)), by = .(gate, neg_source)])

# REDIportal truth sets (gt_score_rediportal_v9.R): germline negatives = neg_source exactly pgp1_ag / pgp1_hom,
# a call at a site that is a positive of the set is not a germline false positive
RP <- file.path(V9, "rediportal_v9/scoring_mismatch76")
g_only <- unique(x[neg_source %chin% c("pgp1_ag", "pgp1_hom"), site])
cat("\nREDIportal sets, germline negatives with an A-to-G call from >= 1 of the 76 configurations:\n")
for (v in c("cov10_alt1", "alt3", "alt3_e05")) {
  ps <- fread(file.path(RP, paste0("positive_sites_", v, ".tsv.gz")))$site
  cat(v, ": ", sum(!g_only %chin% ps), " (", sum(g_only %chin% ps), " are positives of the set)\n", sep = "")
}
cat("REDIportal germline negatives in total:", uniqueN(tr[neg_source %chin% c("pgp1_ag", "pgp1_hom"), site]), "\n")

# for the ground truth sets table (publication/gt_definition/scripts/GTDEF_truth_sets_table.R): germline negatives
# called by >= 1 configuration, not counting those that are also informative negatives
out <- rbind(x[neg_source %chin% c("pgp1_ag", "pgp1_hom"), .(truth = "KO", n_germline_called = uniqueN(site)), by = gate],
             data.table(gate = "any", truth = "REDIportal", n_germline_called = length(g_only)))
fwrite(out, file.path(V9, "truth/germline_negatives_called.tsv"), sep = "\t")
print(out)
