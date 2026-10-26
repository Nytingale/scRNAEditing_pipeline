#!/usr/bin/env Rscript
# =============================================================================
# ding_3prime_bias_genebody.R   (2026-10-23; Reviewer 3, major 2: 3' bias normalised by gene length)
# Position of each edited REDIportal site along its gene body, so that genes of different lengths are comparable:
#   relative distance from the 3' end = genomic distance from the gene 3' end / gene length
#   (0 = gene 3' end, 1 = gene 5' end; gene = union of the gene's RefSeq transcripts on the site's strand)
# Sites inside more than one gene on their strand (the read's gene is unknown) and sites outside RefSeq genes are
# left out (user 2026-10-25; previously such sites were given the shortest gene).
# Input : human/jzlxv3/three_prime_bias/rp_<X10>_<SS2>_<DEPTH>_sites.tsv.gz (ding_3prime_bias_rediportal_pu.R)
# Output: human/jzlxv3/three_prime_bias/rp_<X10>_<SS2>_<DEPTH>_genebody_{sites.tsv.gz,bins.tsv,tests.tsv}
# Usage : X10=x10u SS2=ss2p1 DEPTH=sequenced Rscript ding_3prime_bias_genebody.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
H <- file.path(Sys.getenv("SCRNA_DIR", "."), "human")
OUT <- file.path(H, "jzlxv3", "three_prime_bias")
PREFIX <- paste("rp", Sys.getenv("X10", "x10u"), Sys.getenv("SS2", "ss2p1"), Sys.getenv("DEPTH", "sequenced"), sep = "_")
N_BIN <- 10L

d <- fread(file.path(OUT, paste0(PREFIX, "_sites.tsv.gz")))[edited_10x | edited_ss2]
d <- d[strand %chin% c("+", "-")]

gtf <- fread(cmd = paste("awk -F'\\t' '$3==\"transcript\"'", file.path(H, "ref_annots", "hg38.ncbiRefSeq.gtf")),
             header = FALSE, select = c(1, 4, 5, 7, 9), sep = "\t", quote = "")
setnames(gtf, c("chr", "start", "end", "strand", "attr"))
gtf[, gene := sub('^gene_id "([^"]+)".*', "\\1", attr)]
genes <- gtf[chr %chin% unique(d$chr), .(start = min(start), end = max(end)), by = .(chr, strand, gene)]
genes[, len := as.numeric(end - start + 1L)]
setkey(genes, chr, strand, start, end)

d[, `:=`(s = pos, e = pos)]
ov <- foverlaps(d[, .(chr, strand, pos, s, e)], genes, by.x = c("chr", "strand", "s", "e"), nomatch = NULL)
ov[, n_gene := .N, by = .(chr, strand, pos)]
message(sprintf("%d sites inside more than one gene on their strand, left out", uniqueN(ov[n_gene > 1], by = c("chr", "strand", "pos"))))
ov <- ov[n_gene == 1L]
ov[, rel_3p := fifelse(strand == "+", end - pos, pos - start) / len]
d <- merge(d, ov[, .(chr, strand, pos, gene, gene_len = len, rel_3p)], by = c("chr", "strand", "pos"), all.x = TRUE)
d[, c("s", "e") := NULL]
message(sprintf("%s: %d edited sites, %d inside a RefSeq gene on their strand", PREFIX, nrow(d), d[!is.na(rel_3p), .N]))
fwrite(d, file.path(OUT, paste0(PREFIX, "_genebody_sites.tsv.gz")), sep = "\t")

REG <- c("Coding exon", "5' UTR", "3' UTR", "Noncoding RNA", "Intron")
g <- d[!is.na(rel_3p) & region %chin% REG]
long <- rbind(g[edited_ss2 == TRUE, .(region, rel_3p, gene_len, platform = "Smart-seq2")],
              g[edited_10x == TRUE, .(region, rel_3p, gene_len, platform = "10x 3' v3")])
long <- rbind(long, copy(long)[, region := "All genic sites"])
long[, bin := pmin(N_BIN, floor(rel_3p * N_BIN) + 1L)]
bins <- long[, .N, by = .(region, platform, bin)][CJ(region = unique(long$region), platform = unique(long$platform),
                                                      bin = seq_len(N_BIN)), on = .(region, platform, bin)]
bins[is.na(N), N := 0L]
bins[, pct := 100 * N / sum(N), by = .(region, platform)]
bins[, `:=`(bin_lo = (bin - 1) / N_BIN, bin_hi = bin / N_BIN)]
fwrite(bins, file.path(OUT, paste0(PREFIX, "_genebody_bins.tsv")), sep = "\t")

tests <- long[, {
  a <- rel_3p[platform == "10x 3' v3"]; b <- rel_3p[platform == "Smart-seq2"]
  .(n_10x = length(a), n_ss2 = length(b), median_10x = median(a), median_ss2 = median(b),
    frac_3p_half_10x = mean(a < 0.5), frac_3p_half_ss2 = mean(b < 0.5),
    median_gene_len_10x = median(gene_len[platform == "10x 3' v3"]),
    median_gene_len_ss2 = median(gene_len[platform == "Smart-seq2"]),
    wilcox_p = wilcox.test(a, b)$p.value)
}, by = region]
fwrite(tests, file.path(OUT, paste0(PREFIX, "_genebody_tests.tsv")), sep = "\t")
print(tests)
