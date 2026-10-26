#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_artifact_null_plot.R   (2026-10-08)
#
# What are the non-A>G / non-T>C sites that pass the v6 ground-truth rules, and what do they say
# about chance among the A>G / T>C positives? The rules applied to all 12 mismatch types
# (gt_spectrum12_eval_v6.R); positive = KO depletion >= t in a covered KO arm.
#   a  positives per 1-Mb bin along the genome, per group (A>G, T>C, other 10 types) x t
#   b  sequence context of the positives (first match): satellite; simple repeat or low complexity;
#      within 2 Mb of a centromere (acen) or heterochromatin (gvar) band; Alu; other repeat; unique.
#      The first three are the artifact-prone mask.
# Source data only (not plotted):
#   chance table: observed A>G / T>C positives (outside the mask) vs the number expected by chance:
#      sum over strata s of n_type(s) x r0(s), r0 = pass rate of other-type sites outside the mask.
#      Strata: expected KO alt reads per arm (not covered, <3, 3-5, 5-10, 10-20, >=20) x
#      NTC edit (<0.15, 0.15-0.3, >=0.3); a cell with < 30 null sites backs off to the per-arm cells
#      pooled over NTC edit, then to the most-covered-arm bin. 95% CI: 1,000 bootstrap draws of the
#      null sites (seed 20261008).
# One figure per cover gate (expected-gated = frozen v6; depth-gated).
#
# INPUT : human/pm/gt_eval_v6/spectrum12/{eval,eval_depth}/truth12_sites.tsv.gz
#         human/ref_annots/{rmsk,cytoBand}.txt.gz (UCSC hg38; download_provenance_repeats.txt)
# OUTPUT: figures/<DATASET>/GTDEF_panel_artifact_null{,_depth}.svg (+ source data)
#         source_data/<DATASET>/GTDEF_artifact_null_context{,_depth}.tsv   counts per group x t x context
#         source_data/<DATASET>/GTDEF_artifact_null_chance{,_depth}.tsv    observed vs expected by chance
#         source_data/<DATASET>/GTDEF_artifact_null_sites{,_depth}.tsv.gz  per gate site: context,
#           stratum, r0 and type pass rate at t50, chance probability min(1, r0 / r_type)
# Usage : Rscript GTDEF_artifact_null_plot.R [DATASET]
# ============================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "GTDEF_common.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
argv = commandArgs(trailingOnly = TRUE)
gtdef_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT = caller_script()

REF = file.path(REPO_DIR, "human", "ref_annots")
RMSK = file.path(REF, "rmsk.txt.gz"); CYTO = file.path(REF, "cytoBand.txt.gz"); FAI = file.path(REF, "hg38.analysisSet.fa.fai")
GATES = c(expected = "eval", depth = "eval_depth")
DS = c(t25 = 0.25, t50 = 0.50, t75 = 0.75)
CEN_PAD = 2e6
MIN_NULL = 30
N_BOOT = 1000
LAM_BR = c(0, 3, 5, 10, 20, Inf); LAM_LAB = c("<3", "3-5", "5-10", "10-20", "\u226520")
EDIT_BR = c(0, 0.15, 0.3, Inf)
MASK = c("Satellite", "Simple repeat / low complexity", "Pericentromeric")
CTX = c(MASK, "Alu", "Other repeat", "Unique")
CTX_PLOT = c("Alu", MASK, "Other repeat", "Unique")
cmap_ctx = setNames(c("#332288", "#E69F00", "#56B4E9", "#CC79A7", "grey60", "grey88"), CTX)
MASK_COL = "#E41A1C"
MASK_LAB = paste0("Artifact-prone bin: \u2265 half of its sites are in a RepeatMasker satellite, simple repeat or\n",
                  "low-complexity region, or within 2 Mb of a centromere (acen) or heterochromatin (gvar) band")
cmap_grp = c("A>G" = "#D55E00", "T>C" = "#0072B2", "Other" = "grey35")
grp_of = function(ty) factor(fifelse(ty == "AG", "A>G", fifelse(ty == "TC", "T>C", "Other")), names(cmap_grp))
CHR = paste0("chr", c(1:22, "X", "Y"))

# ---- annotation: context of every position --------------------------------------------------
rm = fread(cmd = paste("gzip -cd", shQuote(RMSK)), select = c(6, 7, 8, 12, 13), header = FALSE,
           col.names = c("chr", "start", "end", "cls", "fam"))[chr %in% CHR]
rm[, ctx := fifelse(cls == "Satellite", 1L, fifelse(cls %in% c("Simple_repeat", "Low_complexity"), 2L,
              fifelse(fam == "Alu", 4L, 5L)))]
rm[, start := start + 1L]
setkey(rm, chr, start, end)
cy = fread(cmd = paste("gzip -cd", shQuote(CYTO)), header = FALSE, col.names = c("chr", "start", "end", "band", "stain"))
cen = cy[stain %in% c("acen", "gvar") & chr %in% CHR, .(start = max(0, min(start) - CEN_PAD) + 1, end = max(end) + CEN_PAD), by = chr]
setkey(cen, chr, start, end)
context_of = function(chr, pos) {
  q = unique(data.table(chr = chr, start = pos, end = pos)); setkey(q, chr, start, end)
  r = foverlaps(q, rm, nomatch = NULL)[, .(rep = min(ctx)), by = .(chr, pos = i.start)]
  c3 = foverlaps(q, cen, nomatch = NULL)[, .(chr, pos = i.start, peri = TRUE)]
  q = merge(q[, .(chr, pos = start)], r, by = c("chr", "pos"), all.x = TRUE)
  q = merge(q, unique(c3), by = c("chr", "pos"), all.x = TRUE)
  q[, k := fifelse(!is.na(rep) & rep <= 2L, rep, fifelse(!is.na(peri), 3L, fifelse(!is.na(rep), rep, 6L)))]
  q[, .(chr, pos, context = factor(CTX[k], CTX))]
}
fai = fread(FAI, header = FALSE, select = 1:2, col.names = c("chr", "len"))[chr %in% CHR]
fai[, offset := c(0, cumsum(as.numeric(len))[-.N])]

lam_bin = function(cov, lam, covered) factor(fifelse(!covered, "nc", as.character(cut(lam, LAM_BR, LAM_LAB, right = FALSE))),
                                            c("nc", LAM_LAB))
for (g in names(GATES)) {
  SRC = file.path(GT_V6_DIR, "spectrum12", GATES[[g]], "truth12_sites.tsv.gz")
  x = fread(SRC)
  x[, grp := grp_of(type)]
  x = merge(x, context_of(x$chr, x$pos), by = c("chr", "pos"))
  x[, masked := context %in% MASK]

  # covered arms and d_max parity with gt_spectrum12_eval_v6.R
  x[, `:=`(lam1 = adar1_cov * ntc_edit, lam2 = adarb1_cov * ntc_edit)]
  if (g == "expected") x[, `:=`(c1 = lam1 >= 3, c2 = lam2 >= 3)] else x[, `:=`(c1 = adar1_cov >= 10, c2 = adarb1_cov >= 10)]
  x[, d_chk := pmax(fifelse(c1, 1 - adar1_edit / ntc_edit, -Inf), fifelse(c2, 1 - adarb1_edit / ntc_edit, -Inf))]
  n_bad = x[abs(d_chk - d_max) > 1e-6, .N]
  message(g, ": ", nrow(x), " gate sites; d_max parity mismatches: ", n_bad)
  stopifnot(n_bad == 0, x[, all(c1 | c2)])
  x[, `:=`(s1 = lam_bin(adar1_cov, lam1, c1), s2 = lam_bin(adarb1_cov, lam2, c2),
           se = cut(ntc_edit, EDIT_BR, right = FALSE, labels = c("<0.15", "0.15-0.3", "\u22650.3")))]
  x[, lmax := factor(as.character(cut(pmax(fifelse(c1, lam1, -1), fifelse(c2, lam2, -1)), LAM_BR, LAM_LAB, right = FALSE)), LAM_LAB)]
  for (tn in names(DS)) set(x, j = paste0("pos_", tn), value = x$d_max >= DS[[tn]])

  # ---- b: context counts ----------------------------------------------------------------------
  ctx = rbindlist(lapply(names(DS), function(tn) x[get(paste0("pos_", tn)) == TRUE, .N, by = .(grp, context)][, t := tn]))
  ctx = rbind(ctx, x[, .N, by = .(grp, context)][, t := "eligible"])
  ctx = merge(CJ(grp = levels(x$grp), context = CTX, t = c("eligible", names(DS))), ctx, by = c("grp", "context", "t"), all.x = TRUE)
  ctx[is.na(N), N := 0L][, pct := 100 * N / sum(N), by = .(grp, t)]
  ctx[, `:=`(grp = factor(grp, names(cmap_grp)), context = factor(context, CTX))]
  fwrite(ctx, file.path(SRC_DIR, sprintf("GTDEF_artifact_null_context%s.tsv", if (g == "depth") "_depth" else "")), sep = "\t")
  print(dcast(ctx[t != "eligible"], t + grp ~ context, value.var = "N"))
  bt = ctx[t %in% names(DS)][, t_lab := factor(paste0("\u2265 ", 100 * DS[t], "%"), paste0("\u2265 ", 100 * DS, "%"))]
  bn = bt[, .(n = sum(N), n_mask = sum(N[context %in% MASK])), by = .(grp, t_lab)]
  bt[, context := factor(context, CTX_PLOT)]
  pb = ggplot(bt, aes(grp, pct, fill = context)) +
    geom_col(width = 0.75, colour = "grey20", linewidth = 0.15) +
    geom_text(data = bn, aes(grp, 100, label = paste0(n_mask, "\nof ", n)), inherit.aes = FALSE, vjust = -0.15,
              size = 5.2, lineheight = 0.9, colour = "grey20") +
    facet_wrap(~ t_lab, nrow = 1) +
    scale_fill_manual(values = cmap_ctx, breaks = CTX_PLOT, name = NULL) +
    scale_y_continuous(limits = c(0, 122), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Sites passing the\nground truth criteria (%)", title = "Editing reduction in KO",
         subtitle = "Numbers: sites in satellite, simple repeat / low complexity\nor pericentromeric context, of all") +
    theme_manuscript(base_size = 22, facets = TRUE) + grid_y +
    theme(legend.position = "right", legend.key.size = unit(0.8, "cm"), legend.text = element_text(size = 20),
          plot.title = element_text(size = 22), plot.subtitle = element_text(size = 18)) +
    guides(fill = guide_legend(ncol = 1))

  # ---- a: genome map per t ---------------------------------------------------------------------
  gm = rbindlist(lapply(names(DS), function(tn)
    x[get(paste0("pos_", tn)) == TRUE, .(n = .N, n_mask = sum(masked)), by = .(grp, chr, bin = pos %/% 1e6)][, t := tn]))
  gm = merge(gm, fai[, .(chr, offset)], by = "chr")
  gm[, `:=`(xg = (offset + (bin + 0.5) * 1e6) / 1e6, Mask = factor(fifelse(n_mask >= n / 2, MASK_LAB, "Outside mask"),
                                                              c(MASK_LAB, "Outside mask")),
            row = factor(paste0(grp, "\n\u2265 ", 100 * DS[t], "%"),
                         paste0(rep(names(cmap_grp), each = 3), "\n\u2265 ", 100 * DS, "%")))]
  lab = gm[grp == "Other"][order(-n), .SD[1:3], by = t][, lab := sprintf("%s: %d\u2013%d Mb, %d sites", chr, bin, bin + 1L, n)]
  lab[, `:=`(right = xg > 0.75 * max(gm$xg), side = n == max(n) & xg > 0.75 * max(gm$xg)), by = t]
  chr_mid = fai[chr %in% unique(gm$chr), .(chr, mid = (offset + len / 2) / 1e6, lo = offset / 1e6, hi = (offset + len) / 1e6)]
  chr_mid[, `:=`(shade = seq_len(.N) %% 2 == 0, short = sub("^chr", "", chr))]
  pa = ggplot(gm) +
    geom_rect(data = chr_mid[shade == TRUE], aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf), fill = "grey95", inherit.aes = FALSE) +
    geom_segment(aes(x = xg, xend = xg, y = 0, yend = n, colour = Mask), linewidth = 0.45) +
    geom_point(aes(xg, n, colour = Mask), size = 1.3) +
    geom_text(data = lab, aes(xg, n, label = lab, hjust = fifelse(right, 1.04, -0.06), vjust = fifelse(side, 0.5, -0.35)),
              size = 7, colour = "grey15") +
    facet_grid(row ~ ., scales = "free_y") +
    scale_colour_manual(values = setNames(c(MASK_COL, "grey45"), c(MASK_LAB, "Outside mask")), name = NULL) +
    scale_x_continuous(breaks = chr_mid$mid, labels = chr_mid$short, expand = expansion(mult = c(0.005, 0.005))) +
    scale_y_continuous(breaks = function(l) { b = pretty(c(0, l[2] / 1.25), n = 2); b[b <= l[2]] }, expand = expansion(mult = c(0, 0.45))) +
    labs(x = "Chromosome", y = "Sites passing the ground truth criteria per 1-Mb bin",
         subtitle = "Rows: mismatch group and editing reduction in KO") +
    theme_manuscript(base_size = 22, facets = TRUE) + grid_y +
    theme(legend.position = "top", panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
          axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          strip.text.y = element_text(angle = 0, size = 18), legend.text = element_text(size = 18),
          plot.subtitle = element_text(size = 18)) +
    guides(colour = guide_legend(override.aes = list(size = 4, linewidth = 1.2)))

  # ---- chance rate from the other types outside the mask (source data only) ----------------
  null = x[grp == "Other" & masked == FALSE]
  test = x[grp %in% c("A>G", "T>C") & masked == FALSE]
  rates = function(d, tn) {
    d = d[, .(s1, s2, se, lmax, k = get(paste0("pos_", tn)))]
    list(full = d[, .(n0 = .N, k0 = sum(k)), by = .(s1, s2, se)],
         arm = d[, .(n0 = .N, k0 = sum(k)), by = .(s1, s2)],
         top = d[, .(n0 = .N, k0 = sum(k)), by = .(lmax)])
  }
  r0_of = function(R, d) {
    f = R$full[d[, .(s1, s2, se)], on = .(s1, s2, se)]; a = R$arm[d[, .(s1, s2)], on = .(s1, s2)]; tp = R$top[d[, .(lmax)], on = .(lmax)]
    fifelse(!is.na(f$n0) & f$n0 >= MIN_NULL, f$k0 / f$n0,
            fifelse(!is.na(a$n0) & a$n0 >= MIN_NULL, a$k0 / a$n0, tp$k0 / tp$n0))
  }
  set.seed(20261008)
  boot_idx = replicate(N_BOOT, sample.int(nrow(null), replace = TRUE), simplify = FALSE)
  ch = rbindlist(lapply(names(DS), function(tn) rbindlist(lapply(c("A>G", "T>C"), function(gp) {
    d = test[grp == gp]
    e = sum(r0_of(rates(null, tn), d), na.rm = TRUE)
    eb = vapply(boot_idx, function(i) sum(r0_of(rates(null[i], tn), d), na.rm = TRUE), numeric(1))
    obs = sum(d[[paste0("pos_", tn)]])
    data.table(t = tn, grp = gp, n_gate = nrow(d), observed = obs, expected_chance = e,
               lo = quantile(eb, 0.025), hi = quantile(eb, 0.975), chance_pct = 100 * e / obs,
               chance_pct_lo = 100 * quantile(eb, 0.025) / obs, chance_pct_hi = 100 * quantile(eb, 0.975) / obs,
               null_n = nrow(null), null_pass = sum(null[[paste0("pos_", tn)]]))
  }))))
  print(ch, digits = 3)
  fwrite(ch, file.path(SRC_DIR, sprintf("GTDEF_artifact_null_chance%s.tsv", if (g == "depth") "_depth" else "")), sep = "\t")

  # per-site chance probability at t50 (for weighting or filtering the truth later)
  R0 = rates(null, "t50")
  x[, r0_t50 := r0_of(R0, x)]
  for (gp in c("A>G", "T>C")) {
    Rg = rates(test[grp == gp], "t50"); i = which(x$grp == gp)
    x[i, r_type_t50 := r0_of(Rg, x[i])]
  }
  x[, chance_prob_t50 := fifelse(grp %in% c("A>G", "T>C") & r_type_t50 > 0, pmin(1, r0_t50 / r_type_t50), NA_real_)]
  fwrite(x[, .(chr, pos, strand, type, ntc_edit, ntc_cov, d_max, context, masked, stratum_adar1 = s1, stratum_adarb1 = s2,
               stratum_edit = se, pos_t25, pos_t50, pos_t75, r0_t50, r_type_t50, chance_prob_t50)],
         file.path(SRC_DIR, sprintf("GTDEF_artifact_null_sites%s.tsv.gz", if (g == "depth") "_depth" else "")), sep = "\t")

  p = pa / (pb + plot_spacer() + plot_layout(widths = c(1.6, 0.4))) + plot_layout(heights = c(2.2, 1)) +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 28, face = "bold"))
  save_figure(p, paste0("GTDEF_panel_artifact_null", if (g == "depth") "_depth" else ""), width = 18, height = 18,
              sources = c(SRC, RMSK, CYTO), data = gm, grid = "y", script = SCRIPT)
}
message("Done")
