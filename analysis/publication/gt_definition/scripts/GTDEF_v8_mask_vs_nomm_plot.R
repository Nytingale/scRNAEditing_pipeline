#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_v8_mask_vs_nomm_plot.R   (2026-10-13)
#
# Why v8 drops the artifact-prone region mask: once the nearby-mismatch filter is applied to the
# truth, the mask removes few further sites and does not change the benchmark.
# All four cleanings remove PGP1 SNV positions, multiallelic and dbSNP common sites
# (gt_truth_v7_variants.R *_pgp1; v8 = nomm_pgp1).
#   a  positives (t = 25 / 50 / 75%) and informative negatives per cleaning, both gates
#   b  A>G share of the 12-type depleted gate sites per cleaning (purity of what remains)
#   c  12-type depleted gate sites (t = 50%) by fate: removed by mask only / filter only / both / kept,
#      A>G vs other types; label = share of the masked sites that the filter also removes
#   d  configuration rank with vs without the mask (v8 vs mask + filter), call state mismatch
#   e  PGP1 SNV positions kept vs removed under the nearby-mismatch filter (v7_nomm vs v8)
# Supersedes TERATOMA_v7_mask_vs_nomm{,_ag_share,_spectrum12{,_depth}}, GTDEF_panel_nomm_removed{,_depth},
# GTDEF_bar_truth_spectrum12_masked{,_depth}, TERATOMA_v7_pgp1.
#
# INPUT : human/pm/gt_eval_v8/truth/support/{truth_variant_counts,spectrum12_cleaning}.tsv,
#         truth12_site_flags.tsv.gz; gt_eval_v8/compare_config_rank.tsv;
#         gt_eval_v8/mask_compare/compare_config_rank.tsv (mask + filter comparison variant, per-unit recall)
# OUTPUT: figures/<DATASET>/GTDEF_v8_mask_vs_nomm.svg (+ source data)
# Usage : Rscript GTDEF_v8_mask_vs_nomm_plot.R [DATASET]
# ============================================================================
source(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "GTDEF_common.R"))
argv <- commandArgs(trailingOnly = TRUE)
gtdef_set_dataset(if (length(argv) >= 1) argv[1] else "pm_gt")
SCRIPT <- caller_script()

F_CNT  <- file.path(V8_SUPPORT, "truth_variant_counts.tsv")
F_SPEC <- file.path(V8_SUPPORT, "spectrum12_cleaning.tsv")
F_T12  <- file.path(V8_SUPPORT, "truth12_site_flags.tsv.gz")
F_CMP7 <- file.path(V8_DIR, "mask_compare", "compare_config_rank.tsv")
SRC <- c(F_CNT, F_SPEC, F_T12, V8_CMP_RANK, F_CMP7)
stopifnot(all(file.exists(SRC)))
BS <- 13
CLEAN_LAB <- c(none_pgp1 = "None", mask_pgp1 = "Region mask", nomm_pgp1 = "Removed if near other\nmismatch types",
               mask_nomm_pgp1 = "Both")
cmap_clean <- setNames(c("grey60", "#E69F00", "#56B4E9", "#009E73"), CLEAN_LAB)
GATE_LAB <- c(expected = "Power-gated", depth = "Depth-gated")
T_LAB <- c(t25 = "\u2265 25%", t50 = "\u2265 50%", t75 = "\u2265 75%")
T_FACET <- setNames(paste(T_LAB, "editing reduction in KO"), names(T_LAB))

# ---- a: truth sizes ----------------------------------------------------------------------------
cnt <- fread(F_CNT)[cleaning %in% names(CLEAN_LAB)]
stopifnot(cnt[cleaning == "nomm_pgp1" & gate == "expected" & t == "t50", positive] == 1462L)
ca <- rbind(cnt[, .(gate, cleaning, what = T_LAB[t], n = positive)],
            unique(cnt[, .(gate, cleaning, what = "Informative\nnegatives", n = informative_negative)]))
ca[, `:=`(cl = factor(CLEAN_LAB[cleaning], levels = CLEAN_LAB), g = factor(GATE_LAB[gate], levels = GATE_LAB),
          what = factor(what, levels = c(T_LAB, "Informative\nnegatives")))]
p_a <- ggplot(ca, aes(what, n, fill = cl)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = formatC(n, format = "d", big.mark = ",")), position = position_dodge(width = 0.85),
            vjust = -0.3, size = 2.7, angle = 0) +
  facet_grid(g ~ ., scales = "free_y") +
  scale_fill_manual(values = cmap_clean, name = "Cleaning") +
  scale_y_continuous(labels = shorten_axis, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Editing reduction in KO", y = "Ground-truth sites", title = "Ground-truth sizes per cleaning (PGP1 SNVs removed)") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y + theme(legend.position = "bottom")

# ---- b: A>G share of 12-type depleted gate sites -------------------------------------------------
sp <- fread(F_SPEC)[cleaning %in% names(CLEAN_LAB)]
sb <- sp[, .(ag_pct = 100 * sum(N[type == "AG"]) / sum(N), n = sum(N)), by = .(gate, t, cleaning)]
sb[, `:=`(cl = factor(CLEAN_LAB[cleaning], levels = CLEAN_LAB), g = factor(GATE_LAB[gate], levels = GATE_LAB),
          tl = factor(T_LAB[t], levels = T_LAB))]
p_b <- ggplot(sb, aes(tl, ag_pct, fill = cl)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f", ag_pct)), position = position_dodge(width = 0.85), vjust = -0.3, size = 2.9) +
  facet_grid(g ~ ., scales = "free_y") +
  scale_fill_manual(values = cmap_clean, guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Editing reduction in KO (all 12 mismatch types)", y = "A-to-G share of sites (%)",
       title = "Purity: A-to-G share across all 12 mismatch types") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y

# ---- c: fate of the 12-type depleted gate sites (t = 50%) ------------------------------------
t12 <- fread(F_T12)[pgp1 == FALSE & d_max >= 0.5]
FATE <- c("Kept", "Region mask only", "Near other mismatch types only", "Both")
t12[, fate := factor(fifelse(!masked & nomm_pass_pgp1, FATE[1], fifelse(masked & nomm_pass_pgp1, FATE[2],
                     fifelse(!masked & !nomm_pass_pgp1, FATE[3], FATE[4]))), levels = FATE)]
t12[, ty := factor(fifelse(type == "AG", "A-to-G", "Other 11 types"), levels = c("A-to-G", "Other 11 types"))]
cc <- t12[, .N, by = .(gate, ty, fate)]
cc <- merge(CJ(gate = names(GATE_LAB), ty = levels(t12$ty), fate = FATE), cc, all.x = TRUE)[is.na(N), N := 0L]
cc[, `:=`(g = factor(GATE_LAB[gate], levels = GATE_LAB), fate = factor(fate, levels = rev(FATE)),
          ty = factor(ty, levels = c("A-to-G", "Other 11 types")))]
cmap_fate <- setNames(c("grey80", "#E69F00", "#56B4E9", "#009E73"), FATE)
mc <- cc[, .(masked = sum(N[fate %in% FATE[c(2, 4)]]), both = N[fate == FATE[4]],
             total = sum(N)), by = .(gate, g, ty)]
mc[, `:=`(pct_filter = 100 * both / masked,
          lab = sprintf("%d of %d masked\nare also near other\nmismatch types (%.0f%%)", both, masked, 100 * both / masked))]
print(mc)
p_c <- ggplot(cc, aes(ty, N, fill = fate)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = ifelse(N > 0, formatC(N, format = "d", big.mark = ","), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  geom_text(data = mc, aes(ty, total, label = lab), inherit.aes = FALSE, vjust = -0.15, size = 3, lineheight = 0.9) +
  facet_grid(g ~ ., scales = "free_y") +
  scale_fill_manual(values = cmap_fate, breaks = FATE, name = "Removed by") +
  scale_y_continuous(labels = shorten_axis, expand = expansion(mult = c(0, 0.3))) +
  labs(x = NULL, y = "Sites with \u2265 50% editing reduction in KO\n(all 12 mismatch types)", title = "What each cleaning removes") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y + theme(legend.position = "bottom", legend.direction = "vertical")

# ---- d: does the mask change the ranking? -------------------------------------------------------
r8 <- fread(V8_CMP_RANK)[state == V8_STATE, .(gate, t, config, rank_v8 = overall_rank)]
r7 <- fread(F_CMP7)[variant == "v7_mask_nomm_pgp1" & state == V8_STATE, .(gate, t, config, rank_mask = overall_rank)]
rr <- merge(r8, r7, by = c("gate", "t", "config"))
stopifnot(rr[, .N, by = .(gate, t)][, all(N == 58L)])
rs <- rr[, .(rho = cor(rank_v8, rank_mask, method = "spearman"),
             top5 = length(intersect(config[rank_v8 <= 5], config[rank_mask <= 5])),
             n_top5 = sum(rank_v8 <= 5), same_top1 = all(config[rank_v8 == 1] %in% config[rank_mask == 1]),
             max_shift = max(abs(rank_v8 - rank_mask))), by = .(gate, t)]
rs[, lab := sprintf("%s\ntop 5 shared: %d / %d", spearman_lab(rho, 3), top5, n_top5)]
rr[, `:=`(g = factor(GATE_LAB[gate], levels = GATE_LAB), tl = factor(T_FACET[t], levels = T_FACET))]
rs[, `:=`(g = factor(GATE_LAB[gate], levels = GATE_LAB), tl = factor(T_FACET[t], levels = T_FACET))]
p_d <- ggplot(rr, aes(rank_mask, rank_v8)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 1.4, alpha = 0.7, colour = unname(cmap_clean[3])) +
  geom_text(data = rs, aes(x = 1, y = 57, label = lab), hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
  facet_grid(g ~ tl, labeller = labeller(tl = label_wrap_gen(22))) + coord_equal() +
  scale_x_continuous(breaks = c(1, 20, 40, 58)) + scale_y_continuous(breaks = c(1, 20, 40, 58)) +
  labs(x = "Configuration rank, with region mask", y = "Configuration rank, without region mask",
       title = "Benchmark ranking with vs without the mask") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_both

# ---- e: PGP1 SNV positions kept vs removed, nearby-mismatch filter (v8 vs v8 without PGP1 removal) ----
pe <- fread(F_CNT)[cleaning %in% c("nomm", "nomm_pgp1")]
pe <- rbind(pe[, .(gate, cleaning, what = T_LAB[t], n = positive)],
            unique(pe[, .(gate, cleaning, what = "Informative\nnegatives", n = informative_negative)]))
pe[, `:=`(pg = factor(fifelse(cleaning == "nomm", "PGP1 SNVs kept", "PGP1 SNVs removed"),
                      levels = c("PGP1 SNVs kept", "PGP1 SNVs removed")),
          g = factor(GATE_LAB[gate], levels = GATE_LAB), what = factor(what, levels = c(T_LAB, "Informative\nnegatives")))]
p_e <- ggplot(pe, aes(what, n, fill = pg)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = formatC(n, format = "d", big.mark = ",")), position = position_dodge(width = 0.85),
            vjust = -0.3, size = 2.7) +
  facet_grid(g ~ ., scales = "free_y") +
  scale_fill_manual(values = c("PGP1 SNVs kept" = "grey75", "PGP1 SNVs removed" = unname(cmap_clean[3])), name = NULL) +
  scale_y_continuous(labels = shorten_axis, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Editing reduction in KO", y = "Ground-truth sites", title = "PGP1 donor SNVs") +
  theme_manuscript(base_size = BS, facets = TRUE) + grid_y +
  theme(legend.position = "bottom", legend.direction = "vertical", axis.text.x = element_text(size = 9))

fwrite(pe, file.path(REPORT_DIR, "GTDEF_v8_pgp1_counts.tsv"), sep = "\t")
fwrite(ca, file.path(REPORT_DIR, "GTDEF_v8_mask_vs_nomm_counts.tsv"), sep = "\t")
fwrite(sb, file.path(REPORT_DIR, "GTDEF_v8_mask_vs_nomm_ag_share.tsv"), sep = "\t")
fwrite(cc, file.path(REPORT_DIR, "GTDEF_v8_mask_vs_nomm_fate.tsv"), sep = "\t")
fwrite(mc, file.path(REPORT_DIR, "GTDEF_v8_mask_vs_nomm_mask_covered.tsv"), sep = "\t")
fwrite(rs, file.path(REPORT_DIR, "GTDEF_v8_mask_vs_nomm_rank.tsv"), sep = "\t")
fig <- cowplot::plot_grid(
  cowplot::plot_grid(p_a, p_b, p_c, nrow = 1, rel_widths = c(1.35, 1, 0.75), labels = c("a", "b", "c"),
                     label_size = 18, align = "h", axis = "t"),
  cowplot::plot_grid(p_d, p_e, nrow = 1, rel_widths = c(1.5, 1), labels = c("d", "e"), label_size = 18),
  ncol = 1, rel_heights = c(1, 0.85))
save_figure(fig, "GTDEF_v8_mask_vs_nomm", width = 17, height = 13, sources = SRC, data = cc, script = SCRIPT)
print(sb[t == "t50"][order(gate, cleaning)])
print(rs)
