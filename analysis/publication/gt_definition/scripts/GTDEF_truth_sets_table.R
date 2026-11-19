#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_truth_sets_table.R
#
# One row per dataset x truth set x editing reduction threshold t (REDIportal sets: t = NA). Same numbers as the
# "Ground truth sets" table of the response (review/EDIT_SCRIPTS/edit_user_truth_table_2026-11-15.py).
#   positives  total = A-to-G + T-to-C (read orientation)
#   negatives  total = germline + other
#     germline  PGP1 germline SNVs with an A-to-G call from >= 1 of the 76 configurations in the NTC arm, not counting
#               those that are also informative negatives (teratoma only; the mice have no genotype)
#     other     informative negatives (KO-based sets) or NTC candidate sites not in REDIportal (REDIportal sets)
#
# INPUT : human/pm/gt_eval_v9/truth/gt_v9_truth_counts.tsv
#         human/pm/gt_eval_v9/truth/germline_negatives_called.tsv (human/scripts/gt_v9_germline_negatives_called.R)
#         human/pm/gt_eval_v9/rediportal_v9/scoring_mismatch76/positives.tsv
#         mouse/aorta_gse255661/truth/eval_v9/excluded/truth_counts.tsv
# OUTPUT: publication/gt_definition/supplementary_tables/ground_truth_sets.tsv
# ============================================================================
suppressPackageStartupMessages(library(data.table))
REPO <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "..", "..", ".."))
OUT <- file.path(REPO, "publication/gt_definition/supplementary_tables/ground_truth_sets.tsv")

GATE <- data.table(gate = c("expected", "depth"), truth_set = c("Power-gated", "Depth-gated"),
                   definition = c("KO coverage x NTC editing level >= 3",
                                  ">= 10 reads in the NTC and in each KO arm used"))
REDI <- data.table(variant = c("cov10_alt1", "alt3", "alt3_e05"),
                   truth_set = c("REDIportal coverage only", "REDIportal coverage and variant reads",
                                 "REDIportal coverage, variant reads and editing"),
                   definition = c("NTC coverage >= 10, variant reads >= 1, editing >= 1%",
                                  "NTC coverage >= 10, variant reads >= 3, editing >= 1%",
                                  "NTC coverage >= 10, variant reads >= 3, editing >= 5%"))
t_of <- function(t) as.integer(sub("^t", "", t))

germ <- fread(file.path(REPO, "human/pm/gt_eval_v9/truth/germline_negatives_called.tsv"))
ter <- merge(fread(file.path(REPO, "human/pm/gt_eval_v9/truth/gt_v9_truth_counts.tsv")),
             germ[truth == "KO", .(gate, n_germline_called)], by = "gate")
ter <- merge(ter, GATE, by = "gate")[, .(
  dataset = "Human teratoma", truth_set, definition, t_pct = t_of(t),
  positives_total = positive, positives_AG = positive_read_AG, positives_TC = positive_read_TC,
  negatives_germline = n_germline_called, negatives_other = neg_infneg, negatives_other_type = "informative negatives")]

red <- merge(fread(file.path(REPO, "human/pm/gt_eval_v9/rediportal_v9/scoring_mismatch76/positives.tsv")), REDI,
             by = "variant")[, .(
  dataset = "Human teratoma", truth_set, definition, t_pct = NA_integer_,
  positives_total = n_pos, positives_AG = n_pos_forward, positives_TC = n_pos_reverse,
  negatives_germline = germ[truth == "REDIportal", n_germline_called], negatives_other = n_neg_universe,
  negatives_other_type = "NTC candidate sites not in REDIportal")]

aor <- merge(fread(file.path(REPO, "mouse/aorta_gse255661/truth/eval_v9/excluded/truth_counts.tsv")), GATE,
             by = "gate")[, .(
  dataset = "Mouse aorta (held out)", truth_set, definition, t_pct = t_of(t),
  positives_total = positive, positives_AG = positive_AG, positives_TC = positive_TC,
  negatives_germline = 0L, negatives_other = informative_negative, negatives_other_type = "informative negatives")]

tab <- rbind(ter, red, aor)
tab[, negatives_total := negatives_germline + negatives_other]
stopifnot(tab[, all(positives_total == positives_AG + positives_TC)], !anyNA(tab[, !"t_pct"]))
tab[, `:=`(dataset = factor(dataset, c("Human teratoma", "Mouse aorta (held out)")),
           truth_set = factor(truth_set, c(GATE$truth_set, REDI$truth_set)))]
setorder(tab, dataset, truth_set, t_pct)
setcolorder(tab, c("dataset", "truth_set", "definition", "t_pct", "positives_total", "positives_AG", "positives_TC",
                   "negatives_total", "negatives_germline", "negatives_other", "negatives_other_type"))
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(tab, OUT, sep = "\t", na = "NA", quote = FALSE)
print(tab)
