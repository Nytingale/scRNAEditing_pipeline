#!/usr/bin/env Rscript
# Summary of single-cell linkage levels (R3.9), NTC arm: site x cell pairs with >= 1 / >= 10 site UMIs linked to a
# heterozygous PGP1 SNV by (1) the same read, (2) the same molecule (CB + UB), (3) the same cell and gene (any SNV
# UMI of the gene in the cell), (4) as 3 with monoallelic evidence (every SNV of the gene with >= 2 UMIs in the cell
# shows one allele; >= 1 such SNV). Also per site: linked site UMIs pooled over cells.
# Input : human/pm/percell_gt_linkage_v8/T*_NTC_ctrl.{span,cell}.tsv.gz
# Output: human/pm/percell_gt_linkage_v8/linkage_cell_summary.tsv, linkage_cell_pooled.tsv
suppressMessages(library(data.table))
LDIR <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/percell_gt_linkage_v8")
PM <- dirname(LDIR)
S <- paste0("T", 1:4, "_NTC_ctrl")
s <- fread(file.path(LDIR, "linkage_sites.tsv"))
# all KO-positive rows; candidate rows only with a SNV UMI of the gene in the cell (others are unlinked at every level)
ko_f <- tempfile()
writeLines(s[ko_positive == 1, site], ko_f)
cl <- rbindlist(lapply(S, function(x) fread(cmd = sprintf(
  "zcat %s | awk -F'\\t' 'NR==FNR{k[$1]=1;next} FNR==1 || $7>=1 || ($3 in k)' %s -",
  file.path(LDIR, paste0(x, ".cell.tsv.gz")), ko_f))))
sp <- rbindlist(lapply(S, function(x) fread(cmd = paste("zcat", file.path(LDIR, paste0(x, ".span.tsv.gz"))))))
N_CELLS <- sum(sapply(S, function(x) length(readLines(file.path(PM, "gt_barcodes", paste0(x, ".tsv"))))))
cl[sp, read_linked_umis := i.n_umis, on = .(sample, barcode, site)]
cl[is.na(read_linked_umis), read_linked_umis := 0L]
cl[, `:=`(gene_umis = fifelse(gene_snv_umis >= 1, site_umis, 0L),
          mono_umis = fifelse(gene_snvs_ge2 >= 1 & gene_snvs_mono == gene_snvs_ge2, site_umis, 0L))]
LEV <- c(read_linked_umis = "Same read", mol_linked_umis = "Same molecule", gene_umis = "Same cell and gene",
         mono_umis = "Same cell and gene, monoallelic", site_umis = "Site covered (any)")
sets <- list(ko = s[ko_positive == 1, site], cand = s[candidate == 1, site])
n_sites_set <- fread(file.path(LDIR, "linkage_summary.tsv"), integer64 = "double")[level == "arm" & count == "reads" & min_reads == 1,
                                                            setNames(n_pairs, set)]
res <- rbindlist(lapply(names(sets), function(k) {
  x <- cl[site %in% sets[[k]]]
  n_sites <- n_sites_set[[k]]
  cols <- if (k == "ko") names(LEV) else setdiff(names(LEV), "site_umis")
  rbindlist(lapply(cols, function(col) data.table(
    set = k, level = LEV[[col]], n_sites = n_sites, n_pairs = as.numeric(n_sites) * N_CELLS,
    pairs_ge1 = sum(x[[col]] >= 1), pairs_ge10 = sum(x[[col]] >= 10),
    pooled_sites_ge10 = x[, sum(get(col)), by = site][V1 >= 10, .N],
    pooled_median = x[, sum(get(col)), by = site][V1 >= 1, as.numeric(median(V1))])))
}))
res[, `:=`(pct_ge1 = 100 * pairs_ge1 / n_pairs, pct_ge10 = 100 * pairs_ge10 / n_pairs)]
fwrite(res, file.path(LDIR, "linkage_cell_summary.tsv"), sep = "\t")
options(width = 200)
print(res)
cat("NTC cells:", N_CELLS, "\n")
print(cl[site %in% sets$ko, .(pairs = .N, any_gene_snv = mean(gene_snv_umis >= 1), ge2 = mean(gene_snvs_ge2 >= 1),
                               mono_given_ge2 = sum(gene_snvs_ge2 >= 1 & gene_snvs_mono == gene_snvs_ge2) /
                                 sum(gene_snvs_ge2 >= 1))])
