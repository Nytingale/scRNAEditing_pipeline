#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_sj_indel_source_data.R   (2026-10-23; split from the _plot.R script 2026-10-25)
#
# R2 major 3 (false-positive SNVs near exon-exon junctions and microindels). KO ground truth sites (>= 50% reduction,
# both gates; unique sites, independent of any caller): positives and informative negatives at each distance
# 0-10 nt and > 10 nt to the nearest splice junction (STAR SJ.out.tab, >= 5 unique reads or annotated) and to the
# nearest indel in the merged NTC reads (>= 3 reads).
#
# INPUT : human/pm/gt_eval_v8/indel_sj/truth_distance_nt.tsv (human/scripts/gt_v8_truth_sj_indel_nt.R)
# OUTPUT: source_data/pm_gt/GTDEF_v8_sj_indel.distances.tsv (drawn by GTDEF_v8_sj_indel_plot.R: response figure 13)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_v8_sj_indel"

D <- normalizePath(file.path(dirname(dirname(V8_SUPPORT)), "indel_sj"))
F_NT <- file.path(D, "truth_distance_nt.tsv")
stopifnot(file.exists(F_NT))
D_LEV <- c(as.character(0:10), "> 10")
GATE_LAB <- c(expected = "Power-gated", depth = "Depth-gated")

n <- fread(F_NT)
n[, `:=`(d = factor(d, D_LEV), g = factor(GATE_LAB[gate], GATE_LAB),
         k = factor(c(positive = "Positive", informative_negative = "Informative negative")[class],
                    names(cmap_truth_class)))]
n[, rng := factor(fifelse(d == D_LEV[length(D_LEV)], "far", "near"), c("near", "far"))]
write_sd(n, SRC_DIR, STEM, "distances")
