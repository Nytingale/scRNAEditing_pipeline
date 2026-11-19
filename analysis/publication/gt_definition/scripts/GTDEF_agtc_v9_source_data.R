#!/usr/bin/env Rscript
# ============================================================================
# GTDEF_agtc_v9_source_data.R
#
# v9 version of GTDEF_agtc_truth_plot.R for the headline teratoma truth (umi_tools --random-seed 8787). Teratoma only:
# this figure supports a ground-truth construction decision, and the aorta is held out. Mismatch types are in read orientation (10x 3'), so read-orientation T-to-C on
# one strand is the same genomic change as A-to-G on the other. All 12-type gate sites before cleaning:
#   pass   per type, % of gate sites with editing reduction in KO >= t (t = 25-75%)
#   gene   per A-to-G / T-to-C / other and KO-reduction bin, where the gene bodies lie relative to the read strand
#          (gene on the read strand only / on the opposite strand only / on both / no gene)
#
# INPUT : human/pm/gt_eval_v9/build/variants/truth12_site_flags.tsv.gz              (teratoma, GRCh38)
#         annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData
# OUTPUT: source_data/pm_gt/GTDEF_agtc_v9.{pass,gene}.tsv (drawn by GTDEF_agtc_v9_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "GTDEF_common.R"))
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
suppressPackageStartupMessages(library(GenomicRanges))
gtdef_set_dataset("pm_gt")
STEM <- "GTDEF_agtc_v9"

IN <- list(
  Teratoma = list(sites = file.path(REPO_DIR, "human/pm/gt_eval_v9/build/variants/truth12_site_flags.tsv.gz"),
                  annot = file.path(REPO_DIR, "annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData")))
TS <- seq(0.25, 0.75, by = 0.05)
BINS <- c(-Inf, 0, 0.25, 0.5, 0.75, Inf)
BIN_LAB <- c("< 0", "0-25", "25-50", "50-75", "\u2265 75")
GENE_LEV <- c("Gene on the opposite strand", "Gene on the read strand", "Genes on both strands", "No gene")

read_ds <- function(ds) {
  x <- fread(IN[[ds]]$sites, select = c("gate", "chr", "pos", "strand", "type", "d_max"))
  stopifnot(all(x$type %chin% all_vartypes()), all(x$strand %chin% c("+", "-")))
  x <- x[!is.na(d_max)]
  e <- new.env()
  load(IN[[ds]]$annot, envir = e)
  gr <- GRanges(x$chr, IRanges(x$pos, width = 1))
  gp <- overlapsAny(gr, e$plus.annot[e$plus.annot$geneType != "intergenic"])
  gm <- overlapsAny(gr, e$minus.annot[e$minus.annot$geneType != "intergenic"])
  stopifnot(any(gp), any(gm))
  read_gene <- fifelse(x$strand == "+", gp, gm)
  opp_gene <- fifelse(x$strand == "+", gm, gp)
  x[, gene_class := factor(fifelse(opp_gene & !read_gene, GENE_LEV[1], fifelse(read_gene & !opp_gene, GENE_LEV[2],
                                   fifelse(read_gene & opp_gene, GENE_LEV[3], GENE_LEV[4]))), GENE_LEV)]
  x[, `:=`(dataset = ds, grp = fifelse(type %in% c("AG", "TC"), type, "other"))]
  x
}
x <- rbindlist(lapply(names(IN), read_ds))
x[, gate := fifelse(gate == "expected", "Power-gated", "Depth-gated")]

pass <- rbindlist(lapply(TS, function(t) x[, .(t = t, n_gate = .N, n_pass = sum(d_max >= t)), by = .(dataset, gate, type)]))
pass[, pct := 100 * n_pass / n_gate]
write_sd(pass, SRC_DIR, STEM, "pass")

x[, bin := factor(as.character(cut(d_max, BINS, labels = BIN_LAB, right = FALSE)), BIN_LAB)]
gene <- x[, .(n = .N), by = .(dataset, gate, grp, bin, gene_class)]
gene <- merge(CJ(dataset = unique(x$dataset), gate = unique(x$gate), grp = c("AG", "TC", "other"), bin = BIN_LAB,
                 gene_class = GENE_LEV), gene, by = c("dataset", "gate", "grp", "bin", "gene_class"), all.x = TRUE)
gene[is.na(n), n := 0L]
gene[, `:=`(n_bin = sum(n), pct = if (sum(n) > 0) 100 * n / sum(n) else NA_real_), by = .(dataset, gate, grp, bin)]
write_sd(gene, SRC_DIR, STEM, "gene")

chk <- dcast(pass[, .(dataset, gate, t, type, pct)], dataset + gate + t ~ type, value.var = "pct")
others <- setdiff(all_vartypes(), c("AG", "TC"))
chk[, max_other := do.call(pmax, .SD), .SDcols = others]
print(chk[, .(n_t = .N, agtc_above_all_others = sum(pmin(AG, TC) > max_other),
              min_ratio = round(min(pmin(AG, TC) / max_other), 1)), by = .(dataset, gate)])
print(dcast(gene[bin %in% c("50-75", "\u2265 75"), .(n = sum(n)), by = .(dataset, gate, grp, gene_class)][
  , pct := round(100 * n / sum(n), 1), by = .(dataset, gate, grp)], dataset + gate + grp ~ gene_class, value.var = "pct"))
message("Done")
