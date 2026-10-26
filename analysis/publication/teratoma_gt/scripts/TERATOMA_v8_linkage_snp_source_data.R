#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_linkage_snp_source_data.R   (2026-10-17; split from TERATOMA_v8_linkage_snp_plot.R 2026-10-25)
# Review R3.9 (linked allele-specific markers): can a candidate editing site be linked to a heterozygous donor SNP
# on the same reads in single cells? (a) share of sites with a PGP1 heterozygous SNV within one read length;
# (b) share of site x unit pairs with >= 1 and >= 10 linking reads (reads with an aligned base at the site and at a
# heterozygous SNV) at each pseudobulk tier and in single cells, NTC arm; (c) cells needed for 10 linking reads;
# (d) site x unit pairs testable by bulk allele-specific editing criteria (ased_threshold_summary.tsv,
# TERATOMA_v8_linkage_ased_thresholds.R); ko_steps: KO-based positives with an edited linked molecule and with
# both SNV alleles covered (response figure 18).
# Sites: v8 KO positives (union) and the merged NTC arm GIREMI candidate A-to-G sites; sites at PGP1 variant
# positions removed.
# The inputs are read-level tables from the per-cell-type NTC BAMs, which are not deposited; the source data are.
# Input : human/pm/percell_gt_linkage_v8/<T>_<arm>.{span,readlen}.tsv
#         (human/pm/scripts/percell_gt_linkage_v8.sh -> percell_gt_linkage.py),
#         linkage_cell_summary.tsv (human/pm/scripts/percell_gt_linkage_cell_summary.R), ased_threshold_summary.tsv
# Output: source_data/pm_gt/TERATOMA_v8_linkage_snp.{a_sites,b_pairs,c_cells,c_refs,d_ased,ko_steps}.tsv
#         (drawn by TERATOMA_v8_linkage_snp_plot.R); human/pm/percell_gt_linkage_v8/linkage_summary.tsv
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
teratoma_set_dataset("pm_gt")
STEM <- "TERATOMA_v8_linkage_snp"

PM <- file.path(REPO_DIR, "human", "pm")
LDIR <- file.path(PM, "percell_gt_linkage_v8")
VCF <- file.path(PM, "pgp1_snps", "PGP1_hu43860C_hg38_snvs.vcf.gz")
CELLS_F <- file.path(PM, "seurat_analysis", "ground_truth_cells_paper.tsv")
SITES_F <- file.path(LDIR, "linkage_sites.tsv")
TERAS <- paste0("T", 1:4)
SPAN_F <- file.path(LDIR, paste0(TERAS, "_NTC_ctrl.span.tsv.gz"))
RLEN_F <- list.files(LDIR, "\\.readlen\\.tsv$", full.names = TRUE)
stopifnot(all(file.exists(c(VCF, CELLS_F, SITES_F, SPAN_F))), length(RLEN_F) == 20)

# --- sites, SNVs, read length ------------------------------------------------
v <- fread(cmd = paste("zcat", shQuote(VCF), "| grep -v '^##'"), select = c(1, 2, 4, 5, 7, 10),
           col.names = c("chrom", "pos", "ref", "alt", "filter", "gt"))
v[, gt := sub(":.*", "", gt)]
het <- v[filter == "PASS" & gt == "0/1" & nchar(ref) == 1 & nchar(alt) == 1, .(chrom, snp = pos)]
s <- fread(SITES_F)
s[, c("chrom", "pos") := tstrsplit(site, ":", fixed = TRUE)][, pos := as.integer(pos)]
s <- s[!v[, .(chrom, pos)], on = .(chrom, pos)]
rl <- rbindlist(lapply(RLEN_F, fread))[, .(n = sum(n)), by = read_length]
READ_LEN <- rl[which.max(n), read_length]
setkey(het, chrom, snp)
s[, snp := pos]
near <- het[s, on = .(chrom, snp), roll = "nearest", .(site, d = abs(x.snp - i.pos))]
s[near, d_snp := i.d, on = "site"]
SET_LAB <- c(ko = "KO-based positives", cand = "All candidate A-to-G sites")
sets <- list(ko = s[ko_positive == 1, site], cand = s[candidate == 1, site])

# --- linking reads, NTC arm --------------------------------------------------
cells <- fread(CELLS_F)[arm == "NTC_ctrl", .(sample = paste0(sample, "_NTC_ctrl"), barcode, celltype)]
CT_KEEP <- sub("^NTC_ctrl__(.*)\\.bam$", "\\1",
               list.files(file.path(PM, "02_gt_bams_celltype_umi_v6"), "^NTC_ctrl__.*\\.bam$"))
stopifnot(length(CT_KEEP) == 23)
sp <- rbindlist(lapply(SPAN_F, function(f) fread(cmd = paste("zcat", shQuote(f)))))
sp <- sp[site %in% s$site]
sp[cells, celltype := i.celltype, on = .(sample, barcode)]
N_CELLS <- nrow(cells)
stopifnot(sp[is.na(celltype), .N] == 0)

LEVELS <- c(arm = "Merged NTC arm", tera = "NTC teratoma", ct = "NTC cell type", cell = "Single NTC cell")
units <- c(arm = 1L, tera = 4L, ct = length(CT_KEEP), cell = N_CELLS)
agg <- list(arm  = sp[, .(n_reads = sum(n_reads), n_umis = sum(n_umis)), by = site],
            tera = sp[, .(n_reads = sum(n_reads), n_umis = sum(n_umis)), by = .(sample, site)],
            ct   = sp[celltype %in% CT_KEEP, .(n_reads = sum(n_reads), n_umis = sum(n_umis)), by = .(celltype, site)],
            cell = sp[, .(site, n_reads, n_umis)])
summ <- rbindlist(lapply(setNames(nm = names(sets)), function(k) rbindlist(lapply(setNames(nm = names(LEVELS)), function(l) {
  a <- agg[[l]][site %in% sets[[k]]]
  n_pairs <- as.numeric(length(sets[[k]])) * units[[l]]
  rbindlist(lapply(c(reads = "n_reads", umis = "n_umis"), function(col) data.table(
    min_reads = c(1L, 10L), n_pairs = n_pairs,
    n_hit = c(sum(a[[col]] >= 1), sum(a[[col]] >= 10)))), idcol = "count")
}), idcol = "level")), idcol = "set")
summ[, pct := 100 * n_hit / n_pairs]
summ[, `:=`(Set = factor(SET_LAB[set], SET_LAB), Level = factor(LEVELS[level], LEVELS))]
fwrite(summ, file.path(LDIR, "linkage_summary.tsv"), sep = "\t")
print(summ[, .(set, level, count, min_reads, n_pairs, n_hit, pct = signif(pct, 3))])

# --- a: sites with a heterozygous SNV within one read length / with a linking read -----
pa_d <- rbindlist(lapply(names(sets), function(k) {
  x <- s[site %in% sets[[k]]]
  hit <- c(sum(x$d_snp <= READ_LEN - 1, na.rm = TRUE), sum(x$site %in% agg$arm[n_reads >= 1, site]))
  data.table(Set = SET_LAB[[k]], n = nrow(x), n_hit = hit, pct = 100 * hit / nrow(x),
             what = c(sprintf("Heterozygous SNV\nwithin %d nt", READ_LEN - 1), "\u2265 1 linking read\n(merged NTC arm)"))
}))
pa_d[, `:=`(Set = factor(Set, SET_LAB), what = factor(what, unique(what)))]
write_sd(pa_d, SRC_DIR, STEM, "a_sites")

# --- b: site x unit pairs with linked molecules, pseudobulks and single cells at four linkage levels --
# (UMIs; human/pm/scripts/percell_gt_linkage_cell_summary.R)
CELL_SUM_F <- file.path(LDIR, "linkage_cell_summary.tsv")
CELL_LEV <- c("Same read" = "Same\nread", "Same molecule" = "Same\nmolecule",
              "Same cell and gene" = "SNV in\nsame gene",
              "Same cell and gene, monoallelic" = "SNV in\nsame gene,\nmonoallelic")
B_LEV <- c("Merged\nNTC arm", "NTC\nteratoma", "NTC\ncell type", CELL_LEV)
cs <- fread(CELL_SUM_F, integer64 = "double")[level %in% names(CELL_LEV)]
pb_d <- rbind(
  summ[count == "umis" & level != "cell", .(set, x = sub(" ", "\n", LEVELS[level]), min_reads, n_pairs, n_hit)],
  rbind(cs[, .(set, x = CELL_LEV[level], min_reads = 1L, n_hit = pairs_ge1)],
        cs[, .(set, x = CELL_LEV[level], min_reads = 10L, n_hit = pairs_ge10)])[
    , n_pairs := as.numeric(lengths(sets)[set]) * N_CELLS])
pb_d[, `:=`(pct = 100 * n_hit / n_pairs, Set = factor(SET_LAB[set], SET_LAB), x = factor(x, B_LEV),
            Threshold = factor(sprintf("\u2265 %d", min_reads), c("\u2265 1", "\u2265 10")),
            lab = format(n_hit, big.mark = ",", trim = TRUE))]
stopifnot(pb_d[, all(n_hit > 0)], nrow(pb_d) == 2 * 2 * length(B_LEV))
write_sd(pb_d, SRC_DIR, STEM, "b_pairs")

# --- c: NTC cells needed for 10 linked molecules ------------------------------
pc_d <- rbindlist(lapply(names(sets), function(k) {
  a <- agg$arm[site %in% sets[[k]] & n_umis >= 1]
  data.table(Set = SET_LAB[[k]], site = a$site, n_umis_arm = a$n_umis, n_sites = length(sets[[k]]),
             cells_needed = 10 * N_CELLS / a$n_umis)
}))
pc_d[, Set := factor(Set, SET_LAB)]
setorder(pc_d, Set, cells_needed)
pc_d[, pct_sites := 100 * seq_len(.N) / n_sites, by = Set]
write_sd(pc_d, SRC_DIR, STEM, "c_cells")
ct_size <- cells[celltype %in% CT_KEEP, .N, by = celltype][, median(N)]
write_sd(data.table(ref = c("median_celltype", "all_cells"), x = c(ct_size, N_CELLS),
                    lab = c(sprintf("median NTC\ncell type (%s)", ct_size),
                            sprintf("all NTC cells\n(%s)", format(N_CELLS, big.mark = ",")))),
         SRC_DIR, STEM, "c_refs")

# --- d: site x unit pairs testable by bulk allele-specific editing criteria (TERATOMA_v8_linkage_ased_thresholds.R) --
ASED_F <- file.path(LDIR, "ased_threshold_summary.tsv")
K_LAB <- function(k) sprintf("\u2265 %d\n(\u2265 %d in total)", k, 2 * k)
pd_d <- fread(ASED_F, integer64 = "double")
pd_d[, `:=`(Set = factor(SET_LAB[set], SET_LAB), Level = factor(LEVELS[level], LEVELS),
            K = factor(K_LAB(k), K_LAB(sort(unique(k)))), lab = format(n_hit, big.mark = ",", trim = TRUE))]
write_sd(pd_d, SRC_DIR, STEM, "d_ased")

# --- ko_steps: KO-based positives, site x unit pairs with an edited linked molecule, and among them both SNV
# alleles covered by >= 1 and >= 3 linked molecules each; per_unit = pairs per unit at that level
STEPS <- c("Edited and linked to a SNV", "Both SNV alleles covered (\u2265 1 molecule each)",
           "Both SNV alleles covered (\u2265 3 molecules each)")
fk <- pd_d[set == "ko"]
fn <- rbind(fk[k == 1, .(level, n_pairs, step = 1L, n = n_linked_edited)], fk[k == 1, .(level, n_pairs, step = 2L, n = n_hit)],
            fk[k == 3, .(level, n_pairs, step = 3L, n = n_hit)])
fn[, `:=`(Level = factor(LEVELS[level], LEVELS), Step = factor(STEPS[step], STEPS), pct = 100 * n / n_pairs,
          per_unit = n / units[level])]
write_sd(fn, SRC_DIR, STEM, "ko_steps")
message(sprintf("read length %d; %d NTC cells; median kept NTC cell type %s cells", READ_LEN, N_CELLS, ct_size))
