#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_linkage_ased_thresholds.R   (2026-10-18)
# Review R3.9: is single-cell data too sparse for allele-specific RNA editing (ASED)? A site x unit pair is testable
# at threshold k when one PGP1 heterozygous SNV of the site has >= k linked molecules (cell barcode + UMI covering
# both the site and the SNV) on each SNV allele. k = 5 is the bulk criterion of Zhou et al. 2018 (Nucleic Acids Res;
# >= 5 reads per allele, i.e. >= 10 linked molecules); k = 3, 2, 1 relax it down to observing both alleles once
# (the per-replicate minimum of Park et al. 2017, Genome Biol). Units: merged NTC arm, NTC teratomas (4), NTC cell
# types (23 kept), single NTC cells. Sites: v8 KO positives (union) and merged NTC arm candidate A-to-G sites.
# Input : human/pm/percell_gt_linkage_v8/T<1-4>_NTC_ctrl.allele_cell.tsv.gz
#         (human/pm/scripts/percell_gt_linkage_allele_v8.sh -> percell_gt_linkage_allele.py)
# Output: figures/pm_gt/TERATOMA_v8_linkage_ased_thresholds.svg (+ source data);
#         human/pm/percell_gt_linkage_v8/ased_threshold_summary.tsv
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
SCRIPT <- "publication/teratoma_gt/scripts/TERATOMA_v8_linkage_ased_thresholds.R"

PM <- file.path(REPO_DIR, "human", "pm")
LDIR <- file.path(PM, "percell_gt_linkage_v8")
VCF <- file.path(PM, "pgp1_snps", "PGP1_hu43860C_hg38_snvs.vcf.gz")
CELLS_F <- file.path(PM, "seurat_analysis", "ground_truth_cells_paper.tsv")
SITES_F <- file.path(LDIR, "linkage_sites.tsv")
AL_F <- file.path(LDIR, paste0("T", 1:4, "_NTC_ctrl.allele_cell.tsv.gz"))
stopifnot(all(file.exists(c(VCF, CELLS_F, SITES_F, AL_F))))
K <- c(1L, 2L, 3L, 5L)

v <- fread(cmd = paste("zcat", shQuote(VCF), "| grep -v '^##'"), select = 1:2, col.names = c("chrom", "pos"))
s <- fread(SITES_F)
s[, c("chrom", "pos") := tstrsplit(site, ":", fixed = TRUE)][, pos := as.integer(pos)]
s <- s[!v, on = .(chrom, pos)]
SET_LAB <- c(ko = "KO-based positives", cand = "All candidate A-to-G sites")
sets <- list(ko = s[ko_positive == 1, site], cand = s[candidate == 1, site])

cells <- fread(CELLS_F)[arm == "NTC_ctrl", .(sample = paste0(sample, "_NTC_ctrl"), barcode, celltype)]
CT_KEEP <- sub("^NTC_ctrl__(.*)\\.bam$", "\\1",
               list.files(file.path(PM, "02_gt_bams_celltype_umi_v6"), "^NTC_ctrl__.*\\.bam$"))
stopifnot(length(CT_KEEP) == 23)
N_CELLS <- nrow(cells)

al <- rbindlist(lapply(AL_F, fread))[site %in% s$site]
al[cells, celltype := i.celltype, on = .(sample, barcode)]
stopifnot(al[is.na(celltype), .N] == 0)
al[, n := unedited + edited]

LEVELS <- c(arm = "Merged NTC arm", tera = "NTC teratoma", ct = "NTC cell type", cell = "Single NTC cell")
units <- c(arm = 1L, tera = 4L, ct = length(CT_KEEP), cell = N_CELLS)
by_unit <- list(arm = character(0), tera = "sample", ct = "celltype", cell = c("sample", "barcode"))
# per unit x site x SNV: molecules on the less covered allele; the site must also be edited in the unit (>= 1 edited
# linked molecule and editing >= 1%, Park et al.'s per-replicate minimum), otherwise there is no editing to compare
minal <- rbindlist(lapply(setNames(nm = names(LEVELS)), function(l) {
  d <- if (l == "ct") al[celltype %in% CT_KEEP] else al
  d <- d[, .(n = sum(n), e = sum(edited)), by = c(by_unit[[l]], "site", "snv", "snv_allele")]
  d <- dcast(d, as.formula(paste(paste(c(by_unit[[l]], "site", "snv"), collapse = " + "), "~ snv_allele")),
             value.var = c("n", "e"), fill = 0)
  for (cc in setdiff(c("n_ref", "n_alt", "e_ref", "e_alt"), names(d))) d[, (cc) := 0]
  d[, unit := if (length(by_unit[[l]])) do.call(paste, .SD) else "arm", .SDcols = by_unit[[l]]]
  d[, edited := (e_ref + e_alt) >= 1 & (e_ref + e_alt) / (n_ref + n_alt) >= 0.01]
  d[, .(min_allele = max(pmin(n_ref, n_alt)), min_allele_edited = max(fifelse(edited, pmin(n_ref, n_alt), 0)),
        linked = max(n_ref + n_alt), linked_edited = any(edited)), by = .(unit, site)]
}), idcol = "level")

summ <- rbindlist(lapply(setNames(nm = names(sets)), function(k) rbindlist(lapply(setNames(nm = names(LEVELS)), function(l) {
  a <- minal[level == l & site %in% sets[[k]]]
  data.table(k = K, n_pairs = as.numeric(length(sets[[k]])) * units[[l]],
             n_hit = vapply(K, function(x) sum(a$min_allele_edited >= x), numeric(1)),
             n_sites = vapply(K, function(x) uniqueN(a[min_allele_edited >= x, site]), numeric(1)),
             n_hit_any_editing = vapply(K, function(x) sum(a$min_allele >= x), numeric(1)),
             n_linked = nrow(a), n_linked_edited = sum(a$linked_edited))
}), idcol = "level")), idcol = "set")
summ[, pct := 100 * n_hit / n_pairs]
fwrite(summ, file.path(LDIR, "ased_threshold_summary.tsv"), sep = "\t")
print(summ[, .(set, level, k, n_pairs, n_linked, n_hit_any_editing, n_hit, n_sites, pct = signif(pct, 3))])

# --- plot: pairs testable by threshold ----------------------------------------
FLOOR <- 1e-6
pd <- copy(summ)[, `:=`(Set = factor(SET_LAB[set], SET_LAB), Level = factor(LEVELS[level], LEVELS),
                        K = factor(sprintf("\u2265 %d\n(\u2265 %d in total)", k, 2 * k),
                                   sprintf("\u2265 %d\n(\u2265 %d in total)", K, 2 * K)),
                        y = pmax(pct, FLOOR), lab = format(n_hit, big.mark = ",", trim = TRUE))]
cmap_lev <- setNames(c("#08306B", "#2171B5", "#6BAED6", "#D94801"), LEVELS)
p <- ggplot(pd, aes(K, y, colour = Level, group = Level)) +
  geom_hline(yintercept = FLOOR, colour = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(data = pd[level == "cell"], aes(label = lab), vjust = -0.9, size = 3.4, show.legend = FALSE) +
  facet_wrap(~Set, nrow = 1) +
  scale_colour_manual(values = cmap_lev, name = NULL) +
  scale_y_log10(limits = c(FLOOR, 100), breaks = 10^(-6:2),
                labels = function(x) ifelse(x <= FLOOR, "0", sub("\\.?0+$", "", formatC(x, format = "f", digits = 6)))) +
  labs(x = "Linked molecules required on each allele of a heterozygous SNV",
       y = "Edited site \u00d7 unit pairs testable (%)") +
  theme_manuscript(base_size = 13, axis_text_size = 11) + grid_y +
  theme(legend.position = "bottom")
save_figure(p, "TERATOMA_v8_linkage_ased_thresholds", width = 12, height = 5.5,
            sources = c(SITES_F, AL_F, VCF, CELLS_F),
            data = pd[, .(Set, Level, k, n_pairs, n_hit, n_sites, pct)], grid = "y", script = SCRIPT)
message(sprintf("%d NTC cells", N_CELLS))
