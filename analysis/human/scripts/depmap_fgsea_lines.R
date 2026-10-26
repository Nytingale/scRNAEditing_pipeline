#!/usr/bin/env Rscript
# ============================================================================
# depmap_fgsea_lines.R   (2026-10-04; Reviewer 3, Major 10)
# Per-cell-line gene set enrichment of the Down syndrome gene sets in DepMap 26Q1 lung lines.
# For each line, genes are ranked by expression z-score across the lung lines (log2(TPM + 1),
# default profile per model) and fgsea gives one NES per gene set. Dependent lines
# (ADAR dependency probability >= 0.5) are then compared with non-dependent lines.
# Gene sets exclude ADAR genes and the dsRNA sensors tested individually (as depmap_adar_isg.R).
# Run with the perturbseq env (fgsea).
# Input : human/dl/depmap_data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
#         human/dl/depmap_results/depmap_lung_per_line.tsv.gz (depmap_adar_isg.R)
#         ISG_FILES (as depmap_adar_isg.R)
# Output: human/dl/depmap_results/lung_fgsea_per_line.tsv, lung_fgsea_summary.tsv
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

B <- Sys.getenv("SCRNA_DIR", ".")
DM <- file.path(B, "human", "dl", "depmap_data")
OUT <- file.path(B, "human", "dl", "depmap_results")
TESTED <- c("ADAR", "ADARB1", "ADARB2", "IFIH1", "EIF2AK2", "RIGI", "DDX58", "MAVS", "ZBP1")
PADJ_CUT <- 0.05
NPROC <- as.integer(Sys.getenv("NSLOTS", "1"))
set.seed(1)

read_geneset <- function(entry) {
  nm <- NULL
  if (grepl("=", entry)) { nm <- sub("=.*$", "", entry); entry <- sub("^[^=]*=", "", entry) }
  l <- suppressWarnings(readLines(entry))
  l <- l[nzchar(trimws(l)) & !startsWith(l, "#")]
  if (grepl("\\.gmt$", entry)) {
    s <- strsplit(l, "\t")
    return(setNames(lapply(s, function(x) unique(x[-(1:2)])), vapply(s, `[`, "", 1)))
  }
  g <- trimws(sub("[\t,].*$", "", l))
  g <- g[!tolower(g) %in% c("gene", "genes", "symbol", "gene_symbol")]
  setNames(list(unique(g)), if (is.null(nm)) tools::file_path_sans_ext(basename(entry)) else nm)
}
gsets <- do.call(c, lapply(strsplit(Sys.getenv("ISG_FILES"), ";")[[1]], read_geneset))
ALIAS <- c(DDX58 = "RIGI", MB21D1 = "CGAS", FAM46A = "TENT5A", METTL7B = "TMT1B", PHB = "PHB1")
gsets <- lapply(gsets, function(g) unique(ifelse(g %in% names(ALIAS), ALIAS[g], g)))
gsets <- lapply(gsets, setdiff, TESTED)

pl <- fread(file.path(OUT, "depmap_lung_per_line.tsv.gz"))[!is.na(ADAR_dependent)]
ex <- fread(file.path(DM, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"))
ex <- ex[IsDefaultEntryForModel %in% c("Yes", "True", "TRUE", TRUE) & ModelID %in% pl$ModelID]
meta_cols <- intersect(names(ex), c("V1", "", "SequencingID", "ModelConditionID", "ModelID",
                                     "IsDefaultEntryForMC", "IsDefaultEntryForModel"))
X <- as.matrix(ex[, setdiff(names(ex), meta_cols), with = FALSE])
colnames(X) <- sub(" \\(.*\\)$", "", colnames(X))
X <- X[, !duplicated(colnames(X))]
rownames(X) <- ex$ModelID
rm(ex); invisible(gc())
X <- X[, apply(X, 2, sd) > 0]
Z <- scale(X)
gsets <- lapply(gsets, intersect, colnames(Z))
message(nrow(Z), " lung lines x ", ncol(Z), " genes; set sizes: ",
        paste(sprintf("%s %d", names(gsets), lengths(gsets)), collapse = ", "))

res <- rbindlist(lapply(rownames(Z), function(m) {
  r <- Z[m, ]
  f <- suppressWarnings(fgsea(gsets, r, minSize = 10, maxSize = 1000, nproc = NPROC, nPermSimple = 10000, eps = 0))
  data.table(ModelID = m, isg_set = f$pathway, size = f$size, ES = f$ES, NES = f$NES, pval = f$pval)
}))
message("lines x sets with NA NES: ", res[is.na(NES), .N])
res[, padj_lines := p.adjust(pval, "BH"), by = isg_set]
res <- merge(res, pl[, .(ModelID, ADAR_chronos, ADAR_prob, ADAR_dependent)], by = "ModelID")
res[, enriched := NES > 0 & padj_lines < PADJ_CUT]
fwrite(res, file.path(OUT, "lung_fgsea_per_line.tsv"), sep = "\t")

spear <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  list(rho = unname(ct$estimate), p = ct$p.value)
}
summ <- res[!is.na(NES), {
  d <- ADAR_dependent
  ft <- fisher.test(table(factor(d, c(TRUE, FALSE)), factor(enriched, c(TRUE, FALSE))))
  sd_ <- spear(NES[d], ADAR_chronos[d])
  sn_ <- spear(NES[!d], ADAR_chronos[!d])
  list(n_dep = sum(d), n_not = sum(!d),
       median_NES_dep = median(NES[d]), median_NES_not = median(NES[!d]),
       p_wilcox_NES = wilcox.test(NES[d], NES[!d])$p.value,
       n_enriched_dep = sum(enriched & d), n_enriched_not = sum(enriched & !d),
       pct_enriched_dep = 100 * mean(enriched[d]), pct_enriched_not = 100 * mean(enriched[!d]),
       odds_ratio = unname(ft$estimate), p_fisher = ft$p.value,
       rho_NES_effect_dep = sd_$rho, p_NES_effect_dep = sd_$p,
       rho_NES_effect_not = sn_$rho, p_NES_effect_not = sn_$p)
}, by = isg_set]
fwrite(summ, file.path(OUT, "lung_fgsea_summary.tsv"), sep = "\t")
print(summ)
