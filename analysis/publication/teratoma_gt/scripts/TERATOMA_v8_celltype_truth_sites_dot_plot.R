#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_celltype_truth_sites_dot_plot.R   (2026-10-21)
# Supplementary: positives of each v8 truth set in the per-unit recall denominator of each of the 23 NTC cell
# types (power-gated: cell-type coverage x merged NTC editing level >= 3; depth-gated and REDIportal: >= 10
# stranded reads), one point per truth; cell types ordered by their minimum over the nine truths (ringed).
# Dashed line = the cell-type minimum (MIN_SITES under every truth) for the cell-type tier.
#
# Input : V8_CT_COVERED (plotting/gt_v8_paths.R; headline: unit_cov/celltype_covered_positives_pwcov.tsv, power-gated
#         positives also need >= 10 reads; gt_v8_celltype_keep_pw.R gate_cov)
# Output: figures/pm_gt/TERATOMA_v8_dot_celltype_truth_sites_min<MIN_SITES>.svg
# Usage : Rscript TERATOMA_v8_celltype_truth_sites_dot_plot.R [MIN_SITES (V8_CT_MIN_SITES)]
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_celltype_truth_sites_dot_plot.R"
args <- commandArgs(trailingOnly = TRUE)
MIN_SITES <- if (length(args)) as.integer(args[1]) else V8_CT_MIN_SITES

SRC <- V8_CT_COVERED
stopifnot(file.exists(SRC))
BS <- 16

w <- fread(SRC)
w[, min_all := do.call(pmin, .SD), .SDcols = V8_TRUTH_IDS]
w[, kept := min_all >= MIN_SITES]
lev <- w[order(min_all, unit), unit]
d <- melt(w, id.vars = c("unit", "min_all", "kept"), measure.vars = V8_TRUTH_IDS, variable.name = "truth",
          value.name = "n_sites")
d[, `:=`(Truth = factor(V8_TRUTH_LABELS[as.character(truth)], levels = V8_TRUTH_LABELS[V8_TRUTH_IDS]),
         unit = factor(unit, levels = lev), x = log10(n_sites + 1))]
d[, family := factor(fifelse(startsWith(as.character(truth), "expected"), "power",
                     fifelse(startsWith(as.character(truth), "depth"), "depth", "redi")), levels = c("power", "depth", "redi"))]
mins <- unique(d[, .(unit, x = log10(min_all + 1), kept)])
s <- sort(w$min_all)
gap <- c(max(s[s < MIN_SITES]), min(s[s >= MIN_SITES]))
message(sprintf("%d of %d cell types with >= %d sites under every truth; largest excluded %d, smallest kept %d",
                sum(w$kept), nrow(w), MIN_SITES, gap[1], gap[2]))

cmap <- setNames(NINE_TRUTHS$colour, V8_TRUTH_LABELS[V8_TRUTH_IDS])
smap <- setNames(rep(c(16, 17, 15), each = 3), V8_TRUTH_LABELS[V8_TRUTH_IDS])
brk <- sort(unique(c(0, 1, 10, 100, 1000, 10000, 1e5, MIN_SITES)))
xr <- c(0, max(d$x))
kept_rows <- data.table(unit = factor(w[kept == TRUE, unit], levels = lev))
p <- ggplot(d, aes(x, unit)) +
  geom_tile(data = kept_rows, aes(x = mean(xr), y = unit), width = diff(xr) * 1.2, height = 1,
            inherit.aes = FALSE, fill = "grey90", colour = NA) +
  geom_vline(xintercept = log10(MIN_SITES + 1), linetype = 2, linewidth = 0.7, colour = "black") +
  geom_point(data = mins, aes(x, unit), inherit.aes = FALSE, shape = 21, size = 5.2, stroke = 0.9,
             colour = "black", fill = NA) +
  geom_point(aes(colour = Truth, shape = Truth), size = 2.6, alpha = 0.9) +
  scale_y_discrete(limits = lev) +
  scale_x_continuous(breaks = log10(brk + 1), labels = format(brk, big.mark = ",", scientific = FALSE, trim = TRUE)) +
  coord_cartesian(xlim = xr, expand = TRUE) +
  scale_colour_manual(values = cmap, name = NULL) + scale_shape_manual(values = smap, name = NULL) +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  labs(x = "Truth sites observable in the cell type", y = NULL,
       title = sprintf("Truth sites per cell type (%d NTC cell types)", nrow(w)),
       subtitle = sprintf("Ringed: minimum over the nine truth sets\nShaded: the %d cell types with \u2265 %d under every truth set",
                          sum(w$kept), MIN_SITES)) +
  theme_manuscript(base_size = BS) + grid_both +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = BS - 3),
        legend.position = "bottom", legend.text = element_text(size = BS - 3), legend.key.size = unit(0.9, "lines"))

save_figure(p, sprintf("TERATOMA_v8_dot_celltype_truth_sites_min%d", MIN_SITES), width = 9, height = 11,
            sources = SRC, data = d[, .(unit, truth, Truth, n_sites, min_all, kept)], grid = "both", script = SCRIPT)
