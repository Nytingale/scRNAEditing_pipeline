#!/usr/bin/env Rscript
# =============================================================================
# ding_3prime_bias_rediportal.R   (2026-10-18; Reviewer 3, major 2: editing sites missed because of 3' bias)
# Caller-independent: coverage and edited reads at every REDIportal A-to-I site in Ding et al. human PBMC data,
# 10x Chromium 3' v3 (PBMC1, 4 libraries) vs full-length Smart-seq2 (768 cells, PBMC1 + PBMC2), from the pileups of
# the WASP reference-bias analysis (all reads, before WASP filtering). Smart-seq2 counts are binomially thinned to the
# 10x total read number (474 M vs 723 M reads), so the platforms are compared at equal sequencing depth.
#   covered : coverage >= 10;   edited : coverage >= 10 and >= 2 edited reads
# Sites by RefSeq region (REDIportal Func.refGene) and by genomic distance to the nearest RefSeq transcript 3' end
# (TES) on the site's strand.
# Inputs : human/jzlxv3/wasp_refbias/pileup_all.tsv, human/jzlss2/wasp_refbias/pileup_all.tsv (+ readloss.tsv),
#          ref_annots/TABLE1_hg38_v2.txt.gz, ref_annots/hg38.ncbiRefSeq.gtf
# Output : human/jzlxv3/three_prime_bias/{three_prime_bias_sites.tsv.gz, three_prime_bias_region.tsv,
#          three_prime_bias_tes.tsv}
# =============================================================================
suppressPackageStartupMessages(library(data.table))
H <- file.path(Sys.getenv("SCRNA_DIR", "."), "human")
OUT <- file.path(H, "jzlxv3", "three_prime_bias")
dir.create(OUT, showWarnings = FALSE)
MIN_COV <- 10L
MIN_ALT <- 2L
set.seed(20261018)

reads <- function(d) fread(file.path(H, d, "wasp_refbias", "readloss.tsv"))[metric == "reads_all", value]
f_thin <- reads("jzlxv3") / reads("jzlss2")
stopifnot(f_thin > 0, f_thin < 1)
pu <- function(d) fread(file.path(H, d, "wasp_refbias", "pileup_all.tsv"), select = c("chr", "pos", "coverage", "variant_cov"))
x <- pu("jzlxv3")
s <- pu("jzlss2")
s[, alt_thin := rbinom(.N, variant_cov, f_thin)]
s[, cov_thin := alt_thin + rbinom(.N, coverage - variant_cov, f_thin)]

# TABLE1 has 25 header names but 22 data columns (the three GeneDetail columns are absent): Func.refGene is column 14
rp <- fread(cmd = paste("zcat", file.path(H, "ref_annots", "TABLE1_hg38_v2.txt.gz"), "| cut -f1,2,5,14"))
setnames(rp, c("chr", "pos", "strand", "func"))
stopifnot(rp[, mean(func %chin% c("intergenic", "intronic", "UTR3", "exonic", "UTR5", "ncRNA_intronic", "ncRNA_exonic",
                                  "downstream", "upstream", "upstream;downstream", "splicing", "exonic;splicing",
                                  "ncRNA_exonic;splicing", "ncRNA_splicing", "UTR5;UTR3"))] > 0.99)
d <- merge(rp, x[, .(chr, pos, cov_10x = coverage, alt_10x = variant_cov)], by = c("chr", "pos"), all.x = TRUE)
d <- merge(d, s[, .(chr, pos, cov_ss2 = cov_thin, alt_ss2 = alt_thin)], by = c("chr", "pos"), all.x = TRUE)
for (v in c("cov_10x", "alt_10x", "cov_ss2", "alt_ss2")) set(d, which(is.na(d[[v]])), v, 0L)
d[, region := fcase(func == "UTR3", "3' UTR", func == "exonic", "Coding exon", func == "UTR5", "5' UTR",
                    func == "intronic", "Intron", grepl("^ncRNA", func), "Noncoding RNA",
                    func %chin% c("downstream", "upstream", "upstream;downstream"), "Up- or downstream",
                    default = "Intergenic and other")]

gtf <- fread(cmd = paste("awk -F'\\t' '$3==\"transcript\"'", file.path(H, "ref_annots", "hg38.ncbiRefSeq.gtf")),
             header = FALSE, select = c(1, 4, 5, 7))
setnames(gtf, c("chr", "start", "end", "strand"))
tes_dt <- unique(gtf[, .(chr, strand, tes = fifelse(strand == "+", end, start))])[, tes_pos := tes]
setkey(tes_dt, chr, strand, tes)
d[, tes := pos]
d[, d_tes := abs(tes_dt[d, on = .(chr, strand, tes), roll = "nearest", tes_pos] - pos)]
d[, tes := NULL]
print(d[, .N, by = .(strand, has_tes = !is.na(d_tes))])
d[, `:=`(covered_10x = cov_10x >= MIN_COV, covered_ss2 = cov_ss2 >= MIN_COV,
         edited_10x = cov_10x >= MIN_COV & alt_10x >= MIN_ALT, edited_ss2 = cov_ss2 >= MIN_COV & alt_ss2 >= MIN_ALT)]
fwrite(d, file.path(OUT, "three_prime_bias_sites.tsv.gz"), sep = "\t")

cls <- function(a, b) fcase(a & b, "Both", a, "10x 3' v3 only", b, "Smart-seq2 only", default = "Neither")
reg <- rbind(
  d[, .(what = "covered", class = cls(covered_10x, covered_ss2), region)],
  d[, .(what = "edited", class = cls(edited_10x, edited_ss2), region)])[class != "Neither", .N, by = .(what, region, class)]
reg <- rbind(reg, reg[, .(region = "All", N = sum(N)), by = .(what, class)])
reg[, pct := 100 * N / sum(N), by = .(what, region)]
fwrite(reg, file.path(OUT, "three_prime_bias_region.tsv"), sep = "\t")
print(dcast(reg, what + region ~ class, value.var = "N"))

BRK <- c(-1, 500, 1000, 2000, 5000, 10000, 20000, 50000, Inf)
LAB <- c("0-0.5", "0.5-1", "1-2", "2-5", "5-10", "10-20", "20-50", ">50")
d[, tes_bin := factor(LAB[findInterval(d_tes, BRK, left.open = TRUE)], LAB)]
EXONIC <- c("3' UTR", "Coding exon", "5' UTR")
tb <- rbind(d[!is.na(tes_bin)][, sites := "All"], d[!is.na(tes_bin) & region %chin% EXONIC][, sites := "Exonic"])[
  , .(n_sites = .N, covered_10x = sum(covered_10x), covered_ss2 = sum(covered_ss2), edited_10x = sum(edited_10x),
      edited_ss2 = sum(edited_ss2)), by = .(sites, tes_bin)][order(sites, tes_bin)]
fwrite(tb, file.path(OUT, "three_prime_bias_tes.tsv"), sep = "\t")
print(tb)
message(sprintf("Smart-seq2 thinned by %.3f; %d REDIportal sites", f_thin, nrow(d)))
