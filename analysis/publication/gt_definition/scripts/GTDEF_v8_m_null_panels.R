# ============================================================================
# GTDEF_v8_m_null_panels.R   (2026-10-18)
# Panel builders for the empirical check of the m = 3 gate, shared by GTDEF_gtdefinition_plot.R (Figure 2 e-f,
# GTDEF_v8_informative_ko) and GTDEF_v8_m_null_check_plot.R (standalone copy and source data).
#   m_null_zero_panel(z)    report/pm_gt/GTDEF_v8_m_null_zero.tsv
#   m_null_excess_panel(w)  report/pm_gt/GTDEF_v8_m_null_excess.tsv
# ============================================================================
M_NULL_COL = c("Non-AG/TC sites (10 mismatch types)" = "grey40", "A-to-G candidate sites" = "#D55E00")
M_NULL_XMAX = 8; M_NULL_YMIN = 1e-3

m_null_zero_panel = function(zf, base_size) {
  zf = as.data.frame(zf)
  zf$grp = factor(zf$grp, levels = names(M_NULL_COL))
  zf$lo = pmax(zf$lo, M_NULL_YMIN)
  M_NULL_KO = c(ADAR1 = "ADAR1 KO", ADAR2 = "ADAR2 KO")
  zf$ko = factor(M_NULL_KO[zf$arm], levels = M_NULL_KO)
  e_lab = "e\u207B\u1D50: chance of 0 variant reads if unchanged"
  v_lab = "m = 3 (power-gated)"
  lv = c(names(M_NULL_COL), e_lab, v_lab)
  col = c(M_NULL_COL, setNames(c("black", "grey45"), c(e_lab, v_lab)))
  lt = setNames(c("solid", "solid", "dashed", "dotted"), lv)
  curve = data.frame(m = seq(0.05, M_NULL_XMAX, 0.05)); curve$fit = exp(-curve$m); curve$grp = e_lab
  vl = data.frame(m = c(3, 3), fit = c(M_NULL_YMIN, 1), grp = v_lab)
  both = function(d) do.call(rbind, lapply(M_NULL_KO, function(k) transform(d, ko = k)))
  ln = rbind(zf[, c("m", "fit", "grp")] |> transform(ko = as.character(zf$ko)), both(curve), both(vl))
  ln$grp = factor(ln$grp, levels = lv); ln$ko = factor(ln$ko, levels = M_NULL_KO)
  ggplot(zf, aes(m, fit)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = grp), alpha = 0.2, show.legend = FALSE) +
    geom_line(data = ln, aes(colour = grp, linetype = grp, linewidth = grp)) +
    scale_colour_manual(values = col, name = NULL) +
    scale_linetype_manual(values = lt, name = NULL) +
    scale_linewidth_manual(values = setNames(c(1, 1, 0.6, 0.5), lv), name = NULL) +
    scale_fill_manual(values = M_NULL_COL, guide = "none") +
    facet_wrap(~ ko, nrow = 1) +
    scale_y_log10(breaks = c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001),
                  labels = function(x) paste0(100 * x, "%")) +
    scale_x_continuous(breaks = 0:M_NULL_XMAX) +
    coord_cartesian(xlim = c(0, M_NULL_XMAX), ylim = c(M_NULL_YMIN, 1)) +
    labs(x = "m = expected edited KO reads\n(KO coverage \u00D7 NTC allele fraction)",
         y = "Sites with 0 variant reads\nin the KO (%)", title = "Sites with no variant reads in the KO: non-AG/TC vs A-to-G") +
    guides(colour = guide_legend(ncol = 1, keywidth = unit(2.2, "lines"))) +
    theme_manuscript(base_size = base_size, facets = TRUE) + grid_both +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
}

m_null_excess_panel = function(w, base_size) {
  w = as.data.frame(w)
  w$lab = sprintf("%.0f%%", 100 * w$excess)
  ggplot(w, aes(m_mid, 100 * excess)) +
    geom_vline(xintercept = 3, linetype = 2, colour = "grey45") +
    geom_hline(yintercept = 50, linetype = 3, colour = "grey55") +
    geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 0, linewidth = 0.5, colour = "#01665E") +
    geom_line(colour = "#01665E", linewidth = 0.6) +
    geom_point(size = 2.8, colour = "#01665E") +
    geom_text(aes(label = lab), vjust = -1.1, hjust = 1.15, size = 3.8) +
    annotate("text", x = 3.1, y = 5, label = "m = 3", hjust = 0, size = 4.5, colour = "grey30") +
    scale_x_continuous(breaks = 0:M_NULL_XMAX, limits = c(0, M_NULL_XMAX)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
    labs(x = "m = expected edited KO reads",
         y = "A-to-G zero-read KO arms\nbeyond chance (1 \u2212 non-AG/TC / A-to-G rate)",
         title = "Zero KO reads that reflect real loss") +
    theme_manuscript(base_size = base_size) + grid_both +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
