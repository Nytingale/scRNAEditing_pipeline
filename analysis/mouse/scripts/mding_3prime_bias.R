#!/usr/bin/env Rscript
# =============================================================================
# mding_3prime_bias.R   (2026-11-07; Reviewer 3, major 2; mouse version of human/scripts/ding_3prime_bias_rediportal_pu.R
# (DEPTH = subsampled) + ding_3prime_bias_genebody.R)
# Ding mouse cortex, per mouse (MOUSE = Cortex1 | Cortex2) and Smart-seq2 seed (SEED): REDIportal mm39 sites covered /
# edited on the 10x 3' UMI-aware MarkDuplicates vs the Smart-seq2 coordinate-deduplicated BAM, both subsampled to the
# same number of reads (mding_dedup_subsample.sh, so no thinning; 2026-11-14, was 10x umi_tools: x10us / ss2ds).
#   covered : coverage >= 10;   edited : coverage >= 10 and >= 3 edited reads (REDIportal Ed base; unstranded pileups)
#   region  : REDIportal Func.refGene; distance to the nearest RefSeq transcript 3' end on the site's strand
#   gene body: relative distance from the gene 3' end = distance / gene length (gene = union of its RefSeq transcripts
#             on the site's strand; sites in > 1 gene on their strand or outside RefSeq genes left out)
# Input : mouse/mjzl10x/three_prime_bias/pileup/<M>_{x10ms20261024,ss2dsm<SEED>}_<chr>.tsv.gz
# Output: mouse/mjzl10x/three_prime_bias/rp_<M>_x10ms20261024_ss2dsm<SEED>_subsampled_{sites.tsv.gz, region.tsv, tes.tsv,
#         genebody_sites.tsv.gz, genebody_bins.tsv, genebody_tests.tsv}
# Usage : MOUSE=Cortex1 SEED=20261024 Rscript mding_3prime_bias.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
MO <- file.path(Sys.getenv("SCRNA_DIR", "."), "mouse")
OUT <- file.path(MO, "mjzl10x", "three_prime_bias"); PU <- file.path(OUT, "pileup")
MIN_COV <- 10L; MIN_ALT <- 3L; N_BIN <- 10L
M <- Sys.getenv("MOUSE"); SEED <- Sys.getenv("SEED")
stopifnot(M %in% c("Cortex1", "Cortex2"), SEED %in% c("20261024", "20261025", "20261026"))
X10 <- "x10ms20261024"; SS2 <- paste0("ss2dsm", SEED)
PREFIX <- paste("rp", M, X10, SS2, "subsampled", sep = "_")

n_10x <- fread(file.path(MO, "mjzl10x/results", M, "02_dedup/subsample_umimd", paste0(M, "_umiMd_sub20261024.count")))$reads_primary
n_ss2 <- fread(file.path(MO, "mjzlss2/results", M, "02_dedup/subsample_umimd", sprintf("%s_rmDup_sub%s.count", M, SEED)))$reads_primary
stopifnot(n_10x == n_ss2)
message(sprintf("%s: %d reads per platform", PREFIX, n_10x))

# 25 header names but 26 data fields: Func.refGene is data column 14
rp <- fread(file.path(MO, "ref_annots", "REDIportal_mm39.tsv"), header = FALSE, skip = 1, select = c(1:3, 5, 14))
setnames(rp, c("chr", "pos", "strand", "ed", "func"))
CHRS <- intersect(paste0("chr", c(1:19, "X")), unique(rp$chr))

read_pu <- function(plat, cc) {
  p <- fread(file.path(PU, sprintf("%s_%s_%s.tsv.gz", M, plat, cc)), colClasses = list(character = "ref"))
  for (b in c("A", "C", "G", "T")) set(p, j = b, value = p[[b]] + p[[tolower(b)]])
  m <- merge(rp[chr == cc, .(chr, pos, ed)], p[, .(chr, pos, coverage, A, C, G, T)], by = c("chr", "pos"))
  m[, alt := fcase(ed == "G", G, ed == "C", C, ed == "A", A, ed == "T", T)]
  m[, .(chr, pos, ed, cov = coverage, alt)]
}
pu <- rbindlist(lapply(CHRS, function(cc) {
  x <- read_pu(X10, cc); s <- read_pu(SS2, cc)
  merge(x[, .(chr, pos, ed, cov_10x = cov, alt_10x = alt)], s[, .(chr, pos, ed, cov_ss2 = cov, alt_ss2 = alt)],
        by = c("chr", "pos", "ed"), all = TRUE)
}))
d <- merge(rp[chr %chin% CHRS], pu, by = c("chr", "pos", "ed"), all.x = TRUE)
for (v in c("cov_10x", "alt_10x", "cov_ss2", "alt_ss2")) set(d, which(is.na(d[[v]])), v, 0L)
d[, region := fcase(func == "UTR3", "3' UTR", func == "exonic", "Coding exon", func == "UTR5", "5' UTR",
                    func == "intronic", "Intron", grepl("^ncRNA", func), "Noncoding RNA",
                    func %chin% c("downstream", "upstream", "upstream;downstream"), "Up- or downstream",
                    default = "Intergenic and other")]

gtf <- fread(cmd = paste("awk -F'\\t' '$3==\"transcript\"'", file.path(MO, "ref_annots", "mm39.ncbiRefSeq.gtf")),
             header = FALSE, select = c(1, 4, 5, 7, 9), sep = "\t", quote = "")
setnames(gtf, c("chr", "start", "end", "strand", "attr"))
tes_dt <- unique(gtf[, .(chr, strand, tes = fifelse(strand == "+", end, start))])[, tes_pos := tes]
setkey(tes_dt, chr, strand, tes)
d[, tes := pos]
d[, d_tes := abs(tes_dt[d, on = .(chr, strand, tes), roll = "nearest", tes_pos] - pos)]
d[, tes := NULL]
d[, `:=`(covered_10x = cov_10x >= MIN_COV, covered_ss2 = cov_ss2 >= MIN_COV,
         edited_10x = cov_10x >= MIN_COV & alt_10x >= MIN_ALT, edited_ss2 = cov_ss2 >= MIN_COV & alt_ss2 >= MIN_ALT)]
fwrite(d[covered_10x | covered_ss2], file.path(OUT, paste0(PREFIX, "_sites.tsv.gz")), sep = "\t")

cls <- function(a, b) fcase(a & b, "Both", a, "10x 3' only", b, "Smart-seq2 only", default = "Neither")
reg <- rbind(
  d[, .(what = "covered", class = cls(covered_10x, covered_ss2), region)],
  d[, .(what = "edited", class = cls(edited_10x, edited_ss2), region)])[class != "Neither", .N, by = .(what, region, class)]
reg <- rbind(reg, reg[, .(region = "All", N = sum(N)), by = .(what, class)])
reg[, pct := 100 * N / sum(N), by = .(what, region)]
reg[, `:=`(mouse = M, x10 = X10, ss2 = SS2, reads = n_10x)]
fwrite(reg, file.path(OUT, paste0(PREFIX, "_region.tsv")), sep = "\t")
print(dcast(reg, what + region ~ class, value.var = "pct", fun.aggregate = function(z) if (length(z)) round(z[1], 1) else NA_real_))

BRK <- c(-1, 500, 1000, 2000, 5000, 10000, 20000, 50000, Inf)
LAB <- c("0-0.5", "0.5-1", "1-2", "2-5", "5-10", "10-20", "20-50", ">50")
d[, tes_bin := factor(LAB[findInterval(d_tes, BRK, left.open = TRUE)], LAB)]
EXONIC <- c("3' UTR", "Coding exon", "5' UTR")
tb <- rbind(d[!is.na(tes_bin)][, sites := "All"], d[!is.na(tes_bin) & region %chin% EXONIC][, sites := "Exonic"])[
  , .(n_sites = .N, covered_10x = sum(covered_10x), covered_ss2 = sum(covered_ss2), edited_10x = sum(edited_10x),
      edited_ss2 = sum(edited_ss2)), by = .(sites, tes_bin)][order(sites, tes_bin)]
tb[, `:=`(mouse = M, x10 = X10, ss2 = SS2, reads = n_10x)]
fwrite(tb, file.path(OUT, paste0(PREFIX, "_tes.tsv")), sep = "\t")

# ---- gene body ----
g0 <- d[(edited_10x | edited_ss2) & strand %chin% c("+", "-")]
gtf[, gene := sub('^gene_id "([^"]+)".*', "\\1", attr)]
genes <- gtf[chr %chin% unique(g0$chr), .(start = min(start), end = max(end)), by = .(chr, strand, gene)]
genes[, len := as.numeric(end - start + 1L)]
setkey(genes, chr, strand, start, end)
g0[, `:=`(s = pos, e = pos)]
ov <- foverlaps(g0[, .(chr, strand, pos, s, e)], genes, by.x = c("chr", "strand", "s", "e"), nomatch = NULL)
ov[, n_gene := .N, by = .(chr, strand, pos)]
message(sprintf("%d sites inside more than one gene on their strand, left out", uniqueN(ov[n_gene > 1], by = c("chr", "strand", "pos"))))
ov <- ov[n_gene == 1L]
ov[, rel_3p := fifelse(strand == "+", end - pos, pos - start) / len]
g0 <- merge(g0, ov[, .(chr, strand, pos, gene, gene_len = len, rel_3p)], by = c("chr", "strand", "pos"), all.x = TRUE)
g0[, c("s", "e") := NULL]
message(sprintf("%s: %d edited sites, %d inside a RefSeq gene on their strand", PREFIX, nrow(g0), g0[!is.na(rel_3p), .N]))
fwrite(g0, file.path(OUT, paste0(PREFIX, "_genebody_sites.tsv.gz")), sep = "\t")

REG <- c("Coding exon", "5' UTR", "3' UTR", "Noncoding RNA", "Intron")
g <- g0[!is.na(rel_3p) & region %chin% REG]
long <- rbind(g[edited_ss2 == TRUE, .(region, rel_3p, gene_len, platform = "Smart-seq2")],
              g[edited_10x == TRUE, .(region, rel_3p, gene_len, platform = "10x 3'")])
long <- rbind(long, copy(long)[, region := "All genic sites"])
long[, bin := pmin(N_BIN, floor(rel_3p * N_BIN) + 1L)]
bins <- long[, .N, by = .(region, platform, bin)][CJ(region = unique(long$region), platform = unique(long$platform),
                                                      bin = seq_len(N_BIN)), on = .(region, platform, bin)]
bins[is.na(N), N := 0L]
bins[, pct := 100 * N / sum(N), by = .(region, platform)]
bins[, `:=`(bin_lo = (bin - 1) / N_BIN, bin_hi = bin / N_BIN)]
fwrite(bins, file.path(OUT, paste0(PREFIX, "_genebody_bins.tsv")), sep = "\t")
tests <- long[, {
  a <- rel_3p[platform == "10x 3'"]; b <- rel_3p[platform == "Smart-seq2"]
  .(n_10x = length(a), n_ss2 = length(b), median_10x = median(a), median_ss2 = median(b),
    frac_3p_half_10x = mean(a < 0.5), frac_3p_half_ss2 = mean(b < 0.5),
    median_gene_len_10x = median(gene_len[platform == "10x 3'"]),
    median_gene_len_ss2 = median(gene_len[platform == "Smart-seq2"]),
    wilcox_p = if (length(a) && length(b)) wilcox.test(a, b)$p.value else NA_real_)
}, by = region]
fwrite(tests, file.path(OUT, paste0(PREFIX, "_genebody_tests.tsv")), sep = "\t")
print(tests)
