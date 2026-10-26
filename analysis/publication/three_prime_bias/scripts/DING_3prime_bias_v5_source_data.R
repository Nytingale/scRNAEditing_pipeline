#!/usr/bin/env Rscript
# =============================================================================
# DING_3prime_bias_v5_source_data.R   (2026-10-24; split from human/scripts/ding_3prime_bias_fig_v5.R 2026-10-25)
# Source data for response figure 16 (Reviewer 3, major 2: editing sites missed because of 3' bias).
# Ding et al. human PBMC1, 10x Chromium 3' v3 vs full-length Smart-seq2, the pipeline's dedup BAMs (10x umiDup,
# Smart-seq2 rmDup) subsampled to the same number of unique reads (215,322,584; Smart-seq2 as whole proper pairs,
# seeds 20261024-26; 10x one read dropped, the same reads for every seed). REDIportal sites; covered = >= 10 reads,
# edited = covered and >= 3 variant reads. Smart-seq2 values = mean over seeds, ranges = min-max.
# Gene-based regions of ding_3prime_bias_genebody.R (sites inside exactly one RefSeq gene on their strand).
#   edited_by_region  edited sites per platform and region
#   dist3p_bins       edited sites by absolute distance from the gene 3' end (nt; log10 bins of 0.5) as % of each
#                     platform's edited sites in the region
#   dist3p_tests      median distance per platform and Wilcoxon test, per seed and region
#   overlap_by_region edited sites found by Smart-seq2 only, both platforms or 10x only, per region of the site
#                     annotation (rp_*_sites.tsv.gz; all edited sites, genic or not), mean / min / max over seeds
# The linear-bin and other bin-width layouts are in human/scripts/_superseded/pre_split_2026-10-25/.
# Input : human/jzlxv3/three_prime_bias/rp_x10us20261024_ss2p1ds<seed>_subsampled_{genebody_sites,sites}.tsv.gz
# Output: publication/three_prime_bias/source_data/ding_3prime_bias_v5.<part>.tsv
#         (drawn by DING_3prime_bias_v5_plot.R)
# =============================================================================
B <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "..", "..", ".."))
source(file.path(B, "plotting", "plotting_functions.R"))
source(file.path(B, "plotting", "source_data_io.R"))
SRC_DIR <- file.path(B, "publication/three_prime_bias/source_data")
STEM <- "ding_3prime_bias_v5"
IN <- file.path(B, "human/jzlxv3/three_prime_bias")

PLAT <- names(cmap_3prime_platform)
SEEDS <- 20261024:20261026
GB_REG <- c("All genic sites", "5' UTR", "Coding exon", "3' UTR", "Intron", "Noncoding RNA")
LOG_W <- 0.5
stem <- function(s) file.path(IN, sprintf("rp_x10us20261024_ss2p1ds%d_subsampled", s))

gs <- rbindlist(lapply(SEEDS, function(s) fread(paste0(stem(s), "_genebody_sites.tsv.gz"))[!is.na(rel_3p)][, seed := s]))
gs[, dist := round(rel_3p * gene_len)]
gs <- rbind(gs[region %chin% GB_REG], copy(gs)[, region := "All genic sites"])
long <- rbind(gs[edited_ss2 == TRUE, .(seed, region, dist, platform = PLAT[1])],
              gs[edited_10x == TRUE, .(seed, region, dist, platform = PLAT[2])])
long[, `:=`(region = factor(region, GB_REG), platform = factor(platform, PLAT),
            bin = as.integer(floor(log10(dist + 1) / LOG_W)))]
c_bin <- long[, .N, by = .(seed, region, platform, bin)]
c_bin <- c_bin[CJ(seed = SEEDS, region = factor(GB_REG, GB_REG), platform = factor(PLAT, PLAT),
                  bin = seq(0L, max(long$bin))), on = .(seed, region, platform, bin)][is.na(N), N := 0L]
c_bin[, x := 10^((bin + 0.5) * LOG_W) - 1]
c_bin[, pct := 100 * N / sum(N), by = .(seed, region, platform)]
c_dt <- c_bin[, .(pct = mean(pct), pct_min = min(pct), pct_max = max(pct)), by = .(region, platform, bin, x)]
c_te <- long[, {
  a <- dist[platform == PLAT[2]]; b <- dist[platform == PLAT[1]]
  .(n_10x = length(a), n_ss2 = length(b), median_10x = median(a), median_ss2 = median(b),
    wilcox_p = wilcox.test(a, b)$p.value)
}, by = .(seed, region)]
cnt <- long[, .N, by = .(seed, region, platform)][, .(N = mean(N), N_min = min(N), N_max = max(N)),
                                                    by = .(region, platform)]
OV_REG <- c("Coding exon", "5' UTR", "3' UTR", "Noncoding RNA", "Intron", "Up- or downstream",
            "Intergenic and other", "All regions")
st <- rbindlist(lapply(SEEDS, function(s) fread(paste0(stem(s), "_sites.tsv.gz"))[edited_10x | edited_ss2][, seed := s]))
st <- rbind(st, copy(st)[, region := "All regions"])
stopifnot(all(st$region %chin% OV_REG))
ov <- st[, .(ss2_only = sum(edited_ss2 & !edited_10x), shared = sum(edited_ss2 & edited_10x),
             x10_only = sum(edited_10x & !edited_ss2)), by = .(seed, region)]
ov <- melt(ov, id.vars = c("seed", "region"), variable.name = "set", value.name = "N")[
  , .(N = mean(N), N_min = min(N), N_max = max(N)), by = .(region, set)][, region := factor(region, OV_REG)]
write_sd(cnt, SRC_DIR, STEM, "edited_by_region")
write_sd(ov[order(region, set)], SRC_DIR, STEM, "overlap_by_region")
write_sd(c_dt, SRC_DIR, STEM, "dist3p_bins")
write_sd(c_te, SRC_DIR, STEM, "dist3p_tests")
print(dcast(cnt, region ~ platform, value.var = "N"))
