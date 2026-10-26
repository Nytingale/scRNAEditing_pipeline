#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_mismatch_composition_plot.R
# Mismatch-type composition of every reported configuration's NTC call set in the three
# local-mismatch-context states (HANDOFF 2026-09-29 UPDATE-24). Mismatch types are after
# stranding; counts are unique sites, pooled over the units of a tier.
#   none     = step 4, mismatch = step 5 (near other MM types excluded), both = step 6 (+ clustered)
#
# Input : human/pm/gt_eval_v6/prewindow/vartype/vartype_counts_all_tiers.tsv.gz (gt_vartype_counts_v6.R)
# Outputs (figures/pm_gt/)
#   TERATOMA_mismatch_class_by_window.svg     % A-to-G / T-to-C / other transitions / transversions per configuration
#   TERATOMA_mismatch_retention_by_type.svg   % of step-4 sites of each mismatch type kept by the window filters
#   TERATOMA_mismatch_composition_arm.svg     per-configuration composition, arm tier
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_mismatch_composition_plot.R"

IN <- Sys.getenv("VARTYPE_COUNTS", file.path(REPO_DIR, "human", "pm", "gt_eval_v6", "prewindow", "vartype",
                                            "vartype_counts_all_tiers.tsv.gz"))
TIER_LAB <- c(arm = "Arm", sample = "Sample", celltype = "Cell type")
UNREPORTED_BAM <- "sample"
UNREPORTED_CALLER_BAM <- list(scallele = c("recal", "recalq20"), cellsnp = c("recal", "recalq20"))
TRANSITIONS <- c("CT", "GA")
CLASS_LAB <- c(AG = "A-to-G", TC = "T-to-C", ts = "Other transitions (C-to-T, G-to-A)",
               tv = "Transversions (8 types)")
cmap_class <- setNames(c("#D55E00", "#F0B37E", "#7F7F7F", "#C8C8C8"), CLASS_LAB)
WIN_LAB <- label_windowFilter

cnt <- fread(IN)
cnt <- cnt[!(bam_type %in% UNREPORTED_BAM |
             mapply(function(c, b) b %in% UNREPORTED_CALLER_BAM[[c]], caller, bam_type) &
             caller %in% names(UNREPORTED_CALLER_BAM))]
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")
cnt[, config := do.call(paste, c(.SD, sep = "__")), .SDcols = COMBO]
pool <- cnt[, .(n = sum(n)), by = c("window_filter", "tier", "config", COMBO, "varType")]
full <- CJ(window_filter = names(WIN_LAB), tier = names(TIER_LAB), config = unique(pool$config),
           varType = all_vartypes())
full[, (COMBO) := tstrsplit(config, "__", fixed = TRUE)]
pool <- merge(full, pool, by = c("window_filter", "tier", "config", COMBO, "varType"), all.x = TRUE)
pool[is.na(n), n := 0]
pool[, `:=`(Tier = factor(TIER_LAB[tier], TIER_LAB),
            Window = factor(WIN_LAB[window_filter], WIN_LAB),
            method = std_caller(caller))]

# ---- (a) share of each mismatch class per configuration --------------------------
pool[, cls := fcase(varType == "AG", "AG", varType == "TC", "TC", varType %in% TRANSITIONS, "ts", default = "tv")]
cl <- pool[, .(n = sum(n)), by = .(window_filter, Window, tier, Tier, config, caller, method, strand_mode, cls)]
cl[, total := sum(n), by = .(window_filter, tier, config)]
cl <- cl[total > 0][, pct := 100 * n / total]
cl[, Class := factor(CLASS_LAB[cls], CLASS_LAB)]
fwrite(cl, file.path(SRC_DIR, "TERATOMA_mismatch_class_by_window.tsv"), sep = "\t")

STEM <- "TERATOMA_mismatch_class_by_window"
p <- ggplot(cl, aes(Window, pct)) +
  geom_line(aes(group = config), colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = method), position = position_jitter(width = 0.08, height = 0, seed = 1),
             size = 1.1, alpha = 0.8) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = NA, colour = "grey30") +
  facet_grid(Class ~ Tier, scales = "free_y", labeller = label_wrap_gen(width = 22)) +
  scale_colour_manual_cvd(values = cmap_method, name = NULL) +
  scale_x_discrete(labels = function(x) gsub("Near other MM types excluded", "Near other MM\ntypes excluded", x)) +
  labs(x = "Local mismatch context", y = "% of called sites (after stranding)",
       title = "Mismatch-type composition of the call set, every reported configuration") +
  theme_manuscript() + theme(legend.position = "bottom", axis.text.x = element_text(size = 7))
save_figure(p, STEM, width = 9, height = 10, sources = IN, data = cl, grid = "y", script = SCRIPT)

# ---- (b) retention of each mismatch type by the window filters ------------------
ret <- dcast(pool, tier + Tier + config + caller + method + strand_mode + varType ~ window_filter, value.var = "n")
ret <- melt(ret, measure.vars = c("mismatch", "both"), variable.name = "window_filter", value.name = "n_kept")
ret <- ret[none >= 20][, pct_kept := 100 * n_kept / none]
ret[, `:=`(Window = factor(WIN_LAB[as.character(window_filter)], WIN_LAB[c("mismatch", "both")]),
           varType = factor(varType, c("AG", "TC", TRANSITIONS, setdiff(all_vartypes(), c("AG", "TC", TRANSITIONS)))),
           Class = factor(CLASS_LAB[fcase(varType == "AG", "AG", varType == "TC", "TC",
                                           varType %in% TRANSITIONS, "ts", default = "tv")], CLASS_LAB))]
fwrite(ret, file.path(SRC_DIR, "TERATOMA_mismatch_retention_by_type.tsv"), sep = "\t")

STEM <- "TERATOMA_mismatch_retention_by_type"
p <- ggplot(ret, aes(varType, pmax(pct_kept, 0.01), fill = Class)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  facet_grid(Window ~ Tier, labeller = label_wrap_gen(width = 16)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c("<=0.01", "0.1", "1", "10", "100")) +
  scale_fill_manual(values = cmap_class, name = NULL) +
  labs(x = "Mismatch type (after stranding)", y = "% of sites without window filters that are kept",
       title = "Sites kept by the local mismatch-context filters, by mismatch type") +
  theme_manuscript() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = 2))
save_figure(p, STEM, width = 9, height = 7, sources = IN, data = ret, grid = "y", script = SCRIPT)

# ---- (c) per-configuration composition, arm tier --------------------------------
SHORT_BAM <- c(sampleRmDup = "Dedup (coord.)", sampleUmiDup = "Dedup (UMI)", recal = "Recal",
               recalq20 = "Recal (Q20)", recaldefault = "Recal (default)")
SHORT_STRAND <- c(annotation = "Annotation", hybrid = "Hybrid", caller = "Caller's own")
SHORT_SNP <- c(off = "dbSNP", on = "dbSNP + blacklist")
ca <- cl[tier == "arm"]
ca[, c("bam_type", "snp_subtract") := tstrsplit(config, "__", fixed = TRUE)[c(2, 4)]]
ca[, label := paste(relabel(caller, label_caller), relabel(bam_type, SHORT_BAM),
                    relabel(strand_mode, SHORT_STRAND), relabel(snp_subtract, SHORT_SNP), sep = " | ")]
ord <- ca[window_filter == "none" & cls == "AG"][order(pct), label]
ca[, label := factor(label, ord)]
STEM <- "TERATOMA_mismatch_composition_arm"
p <- ggplot(ca, aes(pct, label, fill = Class)) +
  geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
  geom_vline(xintercept = c(50, 70, 90), linetype = 3, colour = "black", linewidth = 0.3) +
  facet_grid(~ Window) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  scale_fill_manual(values = cmap_class, name = NULL) +
  labs(x = "% of called sites (after stranding)", y = NULL,
       title = "Arm-tier call-set composition per configuration (ordered by A-to-G % without window filters)") +
  theme_manuscript() + theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5),
                             panel.spacing.x = unit(0.8, "lines")) +
  guides(fill = guide_legend(nrow = 2))
save_figure(p, STEM, width = 11, height = 9.5, sources = IN, data = ca, grid = "none", script = SCRIPT)

# ---- summary for the text -------------------------------------------------------
print(cl[cls == "AG", .(max_AG = round(max(pct), 1), median_AG = round(median(pct), 1),
                        n_ge70 = sum(pct >= 70), n = .N), by = .(Tier, window_filter)][order(Tier, window_filter)])
print(ret[, .(median_pct_kept = round(median(pct_kept), 1)), by = .(Tier, window_filter, Class)][order(Tier, window_filter, Class)])
