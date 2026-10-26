#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_nomm_cluster_retention_source_data.R   (2026-10-18; split from the _plot.R script 2026-10-25)
#
# R2 major 3 (does the nearby-mismatch pass discard hyper-edited clusters?). Pool = truth candidates (all 12-type
# gate sites, PGP1 SNV positions removed).
#   a_bins  KO-depleted A-to-G sites binned by the number of other A-to-G pool sites within 200 nt on the same
#           strand (0 / 1 / 2 / 3 / > 3): sites and % kept by "remove if near other MM types" (w1 = 200), gate x t
#   b_kept  both passes, alone and together, on the same KO-responsive candidate sites (best w2 per t): A-to-G vs
#           the other 11 mismatch types kept, % of no filter (ADAR KO cannot act on the other types)
#
# INPUT : human/pm/gt_eval_v8/truth/support/window_cluster_retention.tsv (human/scripts/gt_v8_window_cluster_retention.R)
#         human/pm/gt_eval_v8/truth/support/window_truth_compare.tsv
# OUTPUT: source_data/pm_gt/GTDEF_v8_nomm_cluster_retention.{a_bins,b_kept}.tsv
#         (drawn by GTDEF_v8_nomm_cluster_retention_plot.R: response figure 8)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_v8_nomm_cluster_retention"

F_CLU <- file.path(V8_SUPPORT, "window_cluster_retention.tsv")
F_CMP <- file.path(V8_SUPPORT, "window_truth_compare.tsv")
stopifnot(file.exists(F_CLU), file.exists(F_CMP))
GATE_LAB <- c(expected = "Power-gated", depth = "Depth-gated")
T_LAB <- c(t25 = "\u2265 25% editing\nreduction in KO", t50 = "\u2265 50% editing\nreduction in KO",
           t75 = "\u2265 75% editing\nreduction in KO")
CL_LEV <- c("0", "1", "2", "3", ">3")

cr <- fread(F_CLU)[cl != "all"]
stopifnot(setequal(cr$cl, CL_LEV), cr[gate == "expected" & t == "t50", sum(n)] == 1686L)
cr[, `:=`(clf = factor(cl, levels = CL_LEV), g = factor(GATE_LAB[gate], levels = GATE_LAB),
          tl = factor(T_LAB[t], levels = T_LAB))]
write_sd(cr[, .(g, tl, clf, n, pct = kept_pct)], SRC_DIR, STEM, "a_bins")

PASS_LAB <- c(none = "No local\ncontext filters", mismatch = "+ Remove if near\nother MM types",
              clustered = "+ Clustered", both = "+ Remove if near\nother MM types\n+ Clustered")
CLS_LAB <- c(ag = "KO-responsive A-to-G sites", other = "Other 11 mismatch types,\ndecreased mismatch rate in KO")
cm <- fread(F_CMP)
stopifnot(setequal(cm$state, names(PASS_LAB)), nrow(cm) == 24L)
cm <- merge(cm, cm[state == "none", .(gate, t, ag0 = n_ag, oth0 = n_depleted - n_ag)], by = c("gate", "t"))
m <- rbind(cm[, .(gate, t, state, w, cls = CLS_LAB[["ag"]], n = n_ag, pct = 100 * n_ag / ag0)],
           cm[, .(gate, t, state, w, cls = CLS_LAB[["other"]], n = n_depleted - n_ag,
                  pct = 100 * (n_depleted - n_ag) / oth0)])
m[, `:=`(g = factor(GATE_LAB[gate], levels = GATE_LAB), st = factor(PASS_LAB[state], levels = PASS_LAB),
         cls = factor(cls, levels = CLS_LAB), tl = factor(T_LAB[t], levels = T_LAB))]
write_sd(m, SRC_DIR, STEM, "b_kept")
