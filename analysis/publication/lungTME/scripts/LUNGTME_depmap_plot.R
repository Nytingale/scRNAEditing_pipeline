#!/usr/bin/env Rscript
# ============================================================================
# LUNGTME_depmap_plot.R   (2026-10-04; Reviewer 3, Major 10)
# DepMap 26Q1 lung cell lines: ADAR dependency vs interferon state, four gene sets of the Down
# syndrome analysis (Hallmark IFN-alpha, IFN-gamma, IL6-JAK-STAT3; HCNetlas 423 ISGs).
# Dependent = DepMap gene dependency probability >= 0.5.
#   a  ADAR and ADARB1 dependency probability, all lung lines; dependent lines in red
#   b  ADAR Chronos vs per-line fgsea NES, ADAR-dependent lines only, one facet per gene set
#      (lines with NA NES dropped per set; n in the strip)
#   c  ADAR Chronos vs expression of the dsRNA sensors / effectors, dependent / not dependent
# Input : human/dl/depmap_results/ (human/scripts/depmap_adar_isg.R, depmap_fgsea_lines.R)
# Output: figures/dl/LUNGTME_panel_depmap_adar_isg.svg (+ source data)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "LUNGTME_common.R"))
lungtme_set_dataset("dl")
SCRIPT <- "publication/lungTME/scripts/LUNGTME_depmap_plot.R"

IN_LINE <- file.path(DEPMAP_RES, "depmap_lung_per_line.tsv.gz")
IN_ISG  <- file.path(DEPMAP_RES, "depmap_lung_isg_scores.tsv.gz")
IN_DEP  <- file.path(DEPMAP_RES, "lung_isg_vs_chronos_by_dependency.tsv")
IN_EXPR <- file.path(DEPMAP_RES, "lung_expr_vs_chronos_by_dependency.tsv")
IN_GSEA <- file.path(DEPMAP_RES, "lung_fgsea_per_line.tsv")
IN_GSUM <- file.path(DEPMAP_RES, "lung_fgsea_summary.tsv")
PROB_CUT <- 0.5
SET_LAB <- c(HALLMARK_INTERFERON_ALPHA_RESPONSE = "Hallmark\nIFN-alpha response",
             HALLMARK_INTERFERON_GAMMA_RESPONSE = "Hallmark\nIFN-gamma response",
             HALLMARK_IL6_JAK_STAT3_SIGNALING = "Hallmark\nIL6-JAK-STAT3",
             HCNetlas_423_ISGs = "HCNetlas\n423 ISGs")
GENES_C <- c("EIF2AK2", "IFIH1", "RIGI", "MAVS", "ZBP1")
cmap_dep <- c(Dependent = "#CB181D", `Not dependent` = "grey55")
LAB_SIZE <- 3
ADAR_EFFECT_LAB <- "Growth effect of ADAR KO"

relabel <- function(d) d[, isg_set := factor(SET_LAB[as.character(isg_set)], SET_LAB)]
status_of <- function(dep) factor(fifelse(dep, "Dependent", "Not dependent"), names(cmap_dep))
pl <- fread(IN_LINE)
lg <- relabel(fread(IN_ISG))
by_dep <- relabel(fread(IN_DEP))
bde <- fread(IN_EXPR)
gs <- relabel(fread(IN_GSEA))
gsum <- relabel(fread(IN_GSUM))
stopifnot(!anyNA(lg$isg_set), !anyNA(by_dep$isg_set), !anyNA(gs$isg_set), !anyNA(gsum$isg_set))
fmt_p <- function(p) ifelse(p < 1e-3, sprintf("p = %.1e", p), sprintf("p = %.2g", p))
fmt_lab <- function(rho, p) sprintf("Spearman's rho = %.2f\n%s", rho, fmt_p(p))
jit <- position_jitter(width = 0.2, height = 0, seed = 1)
base <- theme_manuscript(base_size = 14, facets = TRUE)

scatter_layers <- function(lab) list(
  geom_label(data = lab, aes(x = -Inf, y = Inf, label = fmt_lab(rho, p)), hjust = 0, vjust = 1,
             size = LAB_SIZE, lineheight = 0.9, fill = NA, label.size = 0, label.padding = unit(0.5, "lines")),
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.4)))
)

# a
la <- melt(pl, id.vars = "ModelID", measure.vars = c("ADAR_prob", "ADARB1_prob"),
           variable.name = "gene", value.name = "prob")[is.finite(prob)]
la[, gene := factor(sub("_prob", "", gene), c("ADAR", "ADARB1"))]
la[, status := status_of(prob >= PROB_CUT)]
lab_a <- la[, .(lab = sprintf("%d / %d\ndependent", sum(prob >= PROB_CUT), .N)), by = gene]
p_a <- ggplot(la, aes(gene, prob)) +
  geom_hline(yintercept = PROB_CUT, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = status), position = jit, size = 1.6, alpha = 0.7) +
  geom_boxplot(fill = NA, outlier.shape = NA, colour = "black", width = 0.5) +
  geom_text(data = lab_a, aes(y = 1.13, label = lab), size = LAB_SIZE, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 1.22), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.03, 0))) +
  scale_colour_manual(values = cmap_dep, name = NULL) +
  labs(x = NULL, y = "Dependency probability") + base + grid_y

# b
lb <- gs[ADAR_dependent == TRUE & is.finite(ADAR_chronos) & is.finite(NES)]
lab_b <- lb[, .(rho = cor(NES, ADAR_chronos, method = "spearman"),
                p = suppressWarnings(cor.test(NES, ADAR_chronos, method = "spearman")$p.value), n = .N), by = isg_set]
chk <- gsum[, .(isg_set, rho_sum = rho_NES_effect_dep)][lab_b, on = "isg_set"]
stopifnot(all(abs(chk$rho_sum - chk$rho) < 1e-6))
strip_b <- setNames(sprintf("%s\n(n = %d)", lab_b$isg_set, lab_b$n), as.character(lab_b$isg_set))
lb[, strip := factor(strip_b[as.character(isg_set)], strip_b[levels(isg_set)[levels(isg_set) %in% names(strip_b)]])]
lab_b[, strip := factor(strip_b[as.character(isg_set)], levels(lb$strip))]
p_b <- ggplot(lb, aes(NES, ADAR_chronos)) +
  geom_point(colour = cmap_dep[["Dependent"]], size = 1.8, alpha = 0.75) +
  scatter_layers(lab_b) +
  facet_wrap(~strip, nrow = 1, scales = "free_x") +
  labs(x = "Per-line enrichment (fgsea NES)", y = ADAR_EFFECT_LAB,
       subtitle = "ADAR-dependent lung lines (dependency probability >= 0.5)") +
  base + grid_both + theme(aspect.ratio = 1)

# c
lc <- melt(pl[is.finite(ADAR_chronos) & !is.na(ADAR_dependent)], id.vars = c("ModelID", "ADAR_chronos", "ADAR_dependent"),
           measure.vars = paste0(GENES_C, "_expr"), variable.name = "gene", value.name = "expr")
lc[, gene := factor(sub("_expr$", "", gene), GENES_C)]
lc[, status := status_of(ADAR_dependent)]
lab_c <- bde[variable %in% GENES_C, .(gene = factor(variable, GENES_C), status = status_of(ADAR_dependent), rho, p)]
chk <- lc[, .(rho_plot = cor(expr, ADAR_chronos, method = "spearman")), by = .(gene, status)][lab_c, on = c("gene", "status")]
stopifnot(all(abs(chk$rho_plot - chk$rho) < 1e-6))
n_c <- lc[, .(n = uniqueN(ModelID)), by = status]
lc[, status_lab := factor(sprintf("%s\n(n = %d)", status, n_c[match(lc$status, status), n]))]
lab_c[, status_lab := factor(sprintf("%s\n(n = %d)", status, n_c[match(lab_c$status, status), n]), levels(lc$status_lab))]
p_c <- ggplot(lc, aes(expr, ADAR_chronos)) +
  geom_point(aes(colour = status), size = 1.2, alpha = 0.7, show.legend = FALSE) +
  scatter_layers(lab_c) +
  scale_colour_manual(values = cmap_dep) +
  facet_grid(status_lab ~ gene, scales = "free") +
  labs(x = "Expression, log2(TPM + 1)", y = ADAR_EFFECT_LAB) +
  base + grid_both + theme(aspect.ratio = 1, strip.text.y = element_text(angle = 270))

p <- patchwork::wrap_plots(patchwork::wrap_plots(p_a + theme(legend.position = "none"), p_b, widths = c(0.28, 1)),
                           p_c, ncol = 1, heights = c(1, 1.25)) +
  patchwork::plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 18, face = "bold"))
save_figure(p, "LUNGTME_panel_depmap_adar_isg", width = 16, height = 12,
            sources = c(IN_LINE, IN_GSEA, IN_GSUM, IN_EXPR),
            data = merge(pl, dcast(gs, ModelID ~ isg_set, value.var = "NES"), by = "ModelID"),
            script = SCRIPT)
