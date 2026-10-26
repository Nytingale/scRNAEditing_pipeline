#!/usr/bin/env Rscript
# ============================================================================
# LUNGTME_depmap_exADAR_source_data.R   (2026-10-04; split from LUNGTME_depmap_plot_exADAR.R 2026-10-25)
# Reviewer 3, Major 10. DepMap 26Q1 lung cell lines; ISG scores exclude only ADAR (dsRNA sensors stay in the
# gene sets). Dependent = DepMap gene dependency probability >= 0.5.
#   a_prob   ADAR and ADARB1 dependency probability per lung line
#   b_rho    per-gene Spearman's rho (expression vs ADAR Chronos) in dependent and non-dependent lines, for each
#            ISG set and for all expressed genes (gzipped)
#   b_tests  paired Wilcoxon signed-rank test, dependent vs non-dependent rho of the same genes, per facet
#   c_genes  dsRNA sensors and the N_TOP ISGs with the most negative rho in dependent lines
#   c_ref    cut-off of the 5% lowest rho values over all expressed genes in dependent lines
# INPUT : human/dl/depmap_results_exADAR/ (human/dl/scripts/depmap_adar_isg_exADAR.sh)
# OUTPUT: source_data/dl/LUNGTME_panel_depmap_adar_isg_exADAR.<part>.tsv (drawn by LUNGTME_depmap_exADAR_plot.R:
#         response figure 19)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "LUNGTME_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
lungtme_set_dataset("dl")
STEM <- "LUNGTME_panel_depmap_adar_isg_exADAR"

RES <- file.path(REPO_DIR, "human", "dl", "depmap_results_exADAR")
IN_LINE <- file.path(RES, "depmap_lung_per_line.tsv.gz")
IN_GENE <- file.path(RES, "lung_isg_gene_vs_chronos_by_dependency.tsv")
IN_GW   <- file.path(RES, "lung_genomewide_expr_vs_adar_chronos_by_dependency.tsv.gz")
PROB_CUT <- 0.5
N_TOP <- 20
SET_LAB <- c(HALLMARK_INTERFERON_ALPHA_RESPONSE = "Hallmark\nIFN-alpha",
             HALLMARK_INTERFERON_GAMMA_RESPONSE = "Hallmark\nIFN-gamma",
             HALLMARK_IL6_JAK_STAT3_SIGNALING = "Hallmark\nIL6-JAK-STAT3",
             HCNetlas_423_ISGs = "HCNetlas\nISGs")
DEP_LEV <- c("Dependent", "Not dependent")

status_of <- function(dep) factor(fifelse(dep, "Dependent", "Not dependent"), DEP_LEV)
pl <- fread(IN_LINE)
pg <- fread(IN_GENE)
gw <- fread(IN_GW)
stopifnot(!"ADAR" %in% pg$gene)

# a
la <- melt(pl, id.vars = "ModelID", measure.vars = c("ADAR_prob", "ADARB1_prob"),
           variable.name = "gene", value.name = "prob")[is.finite(prob)]
la[, gene := factor(sub("_prob", "", gene), c("ADAR", "ADARB1"))]
la[, status := status_of(prob >= PROB_CUT)]
write_sd(la, SRC_DIR, STEM, "a_prob")

# b
BG_LAB <- "All expressed\ngenes"
lb <- rbind(
  rbindlist(lapply(names(SET_LAB), function(s)
    pg[get(paste0("in__", s)) == TRUE & is.finite(rho), .(isg_set = SET_LAB[[s]], gene, ADAR_dependent, rho)])),
  gw[is.finite(rho), .(isg_set = BG_LAB, gene, ADAR_dependent, rho)])
n_genes <- lb[ADAR_dependent == TRUE, .(n_genes = uniqueN(gene)), by = isg_set]
facet_lev <- c(SET_LAB, BG_LAB)
strip_of <- function(s) factor(sprintf("%s\n(%s genes)", s, prettyNum(n_genes[match(s, isg_set), n_genes], big.mark = ",")),
                               sprintf("%s\n(%s genes)", facet_lev, prettyNum(n_genes[match(facet_lev, isg_set), n_genes], big.mark = ",")))
lb[, `:=`(strip = strip_of(isg_set), status = status_of(ADAR_dependent), bg = isg_set == BG_LAB)]
pw <- dcast(lb, strip + gene ~ ADAR_dependent, value.var = "rho")[!is.na(`TRUE`) & !is.na(`FALSE`)]
p_b_test <- pw[, .(n_genes = .N, median_diff = median(`TRUE` - `FALSE`),
                   p = wilcox.test(`TRUE`, `FALSE`, paired = TRUE)$p.value), by = strip]
p_b_test[, lab := ifelse(p < 1e-3, sprintf("p = %.1e", p), sprintf("p = %.2g", p))]
print(p_b_test)
write_sd(lb, SRC_DIR, STEM, "b_rho", gzip = TRUE)
write_sd(p_b_test, SRC_DIR, STEM, "b_tests")

# c
dep <- pg[ADAR_dependent == TRUE & is.finite(rho)]
top <- dep[sensor == FALSE][order(rho)][seq_len(min(N_TOP, .N)), gene]
sens <- dep[sensor == TRUE][order(rho), gene]
grp_lab <- c(sensor = "dsRNA sensors", isg = sprintf("ISGs with the most negative rho (top %d)", N_TOP))
lc <- pg[gene %in% c(sens, top) & is.finite(rho), .(gene, sensor, ADAR_dependent, rho)]
lc[, group := factor(fifelse(sensor, grp_lab[["sensor"]], grp_lab[["isg"]]), grp_lab)]
lc[, gene := factor(gene, c(sens, top))]
lc[, status := status_of(ADAR_dependent)]
q05 <- gw[ADAR_dependent == TRUE, quantile(rho, 0.05)]
lab_q <- data.table(group = factor(grp_lab[["isg"]], grp_lab), gene = top[length(top)], rho = q05,
                    lab = "5% lowest rho values across all genes (dependent lines)")
write_sd(lc, SRC_DIR, STEM, "c_genes")
write_sd(lab_q, SRC_DIR, STEM, "c_ref")
