#!/usr/bin/env Rscript
# ============================================================================
# depmap_adar_isg.R   (2026-10-04; Reviewer 3, Major 10)
# DepMap 26Q1 context for the lung TME editing results, lung cell lines only (OncotreeLineage == "Lung").
# Default expression entry per model; Chronos from CRISPRGeneEffect.
#   1. ADAR / ADARB1 Chronos vs each ISG score (Spearman)
#   2. ISG-score tertiles: median ADAR Chronos, fraction dependent (dependency probability >= 0.5,
#      the DepMap recommendation); ISG vs Chronos within dependent / non-dependent lines and
#      over a scan of gene-effect and probability cut-offs
#   3. Chronos vs expression of ADAR, ADARB1 and the dsRNA sensors / effectors
#   4. sensor modifier: lm(ADAR Chronos ~ ISG * sensor expression), per gene set
#   5. genome-wide Spearman of ADAR Chronos with expression
#   6. each gene of the ISG sets and each sensor vs ADAR Chronos, dependent / non-dependent lines
# ISG score: mean per-gene z (across lung lines) of each gene set of the Down syndrome analysis,
# excluding ISG_EXCLUDE (comma-separated; default: the genes tested individually, TESTED).
# DEPMAP_OUT overrides the output directory.
# Input : DepMap 26Q1 CSVs in human/dl/depmap_data/{Model,CRISPRGeneEffect,CRISPRGeneDependency,
#                               OmicsExpressionTPMLogp1HumanProteinCodingGenes}.csv
#         ISG_FILES: ';'-separated, each a .gmt or "name=path" to a one-gene-per-line file
# Output: human/dl/depmap_results/
# ============================================================================
suppressPackageStartupMessages(library(data.table))

B <- Sys.getenv("SCRNA_DIR", ".")
DM <- file.path(B, "human", "dl", "depmap_data")
OUT <- Sys.getenv("DEPMAP_OUT", file.path(B, "human", "dl", "depmap_results"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DEP_CUT <- -0.5
PROB_CUT <- 0.5
TESTED <- c("ADAR", "ADARB1", "ADARB2", "IFIH1", "EIF2AK2", "RIGI", "DDX58", "MAVS", "ZBP1")
SENSORS <- c("IFIH1", "EIF2AK2", "RIGI", "MAVS", "ZBP1")
EXCLUDE <- if (nzchar(Sys.getenv("ISG_EXCLUDE"))) strsplit(Sys.getenv("ISG_EXCLUDE"), ",")[[1]] else TESTED
message("output: ", OUT, "; excluded from ISG scores: ", paste(EXCLUDE, collapse = ","))

read_geneset <- function(entry) {
  nm <- NULL
  if (grepl("=", entry)) { nm <- sub("=.*$", "", entry); entry <- sub("^[^=]*=", "", entry) }
  l <- readLines(entry)
  l <- l[nzchar(trimws(l)) & !startsWith(l, "#")]
  if (grepl("\\.gmt$", entry)) {
    s <- strsplit(l, "\t")
    return(setNames(lapply(s, function(x) unique(x[-(1:2)])), vapply(s, `[`, "", 1)))
  }
  g <- trimws(sub("[\t,].*$", "", l))
  g <- g[!tolower(g) %in% c("gene", "genes", "symbol", "gene_symbol")]
  setNames(list(unique(g)), if (is.null(nm)) tools::file_path_sans_ext(basename(entry)) else nm)
}
isg_files <- strsplit(Sys.getenv("ISG_FILES"), ";")[[1]]
if (!length(isg_files)) stop("set ISG_FILES")
gsets <- do.call(c, lapply(isg_files, read_geneset))
# older symbols in the gene sets -> DepMap (HGNC current) symbols
ALIAS <- c(DDX58 = "RIGI", MB21D1 = "CGAS", FAM46A = "TENT5A", METTL7B = "TMT1B", PHB = "PHB1")
gsets <- lapply(gsets, function(g) unique(ifelse(g %in% names(ALIAS), ALIAS[g], g)))
message("gene set(s): ", paste(sprintf("%s (%d)", names(gsets), lengths(gsets)), collapse = ", "))

strip_sym <- function(x) sub(" \\(.*\\)$", "", x)
model <- fread(file.path(DM, "Model.csv"),
               select = c("ModelID", "CellLineName", "OncotreeLineage", "OncotreePrimaryDisease",
                          "OncotreeSubtype", "OncotreeCode"))[OncotreeLineage == "Lung"]

chr <- fread(file.path(DM, "CRISPRGeneEffect.csv"))
setnames(chr, 1, "ModelID")
setnames(chr, -1, strip_sym(names(chr)[-1]))
chr <- chr[ModelID %in% model$ModelID, .(ModelID, ADAR_chronos = ADAR, ADARB1_chronos = ADARB1)]
dep_cols <- names(fread(file.path(DM, "CRISPRGeneDependency.csv"), nrows = 0))
dep_cols <- c(dep_cols[1], dep_cols[strip_sym(dep_cols) %in% c("ADAR", "ADARB1")])
prob <- fread(file.path(DM, "CRISPRGeneDependency.csv"), select = dep_cols)
setnames(prob, c("ModelID", paste0(strip_sym(dep_cols[-1]), "_prob")))
chr <- merge(chr, prob, by = "ModelID", all.x = TRUE)
chr[, ADAR_dependent := ADAR_prob >= PROB_CUT]

ex <- fread(file.path(DM, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"))
ex <- ex[IsDefaultEntryForModel %in% c("Yes", "True", "TRUE", TRUE) & ModelID %in% chr$ModelID]
stopifnot(!anyDuplicated(ex$ModelID))
meta_cols <- intersect(names(ex), c("V1", "", "SequencingID", "ModelConditionID", "ModelID",
                                     "IsDefaultEntryForMC", "IsDefaultEntryForModel"))
X <- as.matrix(ex[, setdiff(names(ex), meta_cols), with = FALSE])
colnames(X) <- strip_sym(colnames(X))
X <- X[, !duplicated(colnames(X))]
rownames(X) <- ex$ModelID
rm(ex); invisible(gc())
message("lung lines with Chronos and expression: ", nrow(X), " x ", ncol(X), " genes")

Z <- scale(X)
isg <- rbindlist(lapply(names(gsets), function(s) {
  g <- setdiff(intersect(gsets[[s]], colnames(Z)), EXCLUDE)
  g <- g[is.finite(colSums(Z[, g, drop = FALSE]))]
  data.table(set = s, n_in_set = length(gsets[[s]]), n_used = length(g),
             missing = paste(setdiff(gsets[[s]], colnames(X)), collapse = ","),
             genes = paste(g, collapse = ","))
}))
fwrite(isg, file.path(OUT, "isg_geneset_used.tsv"), sep = "\t")
print(isg[, .(set, n_in_set, n_used)])
score <- sapply(names(gsets), function(s) rowMeans(Z[, strsplit(isg[set == s, genes], ",")[[1]], drop = FALSE]))
score <- as.data.table(score, keep.rownames = "ModelID")
gz <- rbindlist(lapply(names(gsets), function(s) {
  g <- strsplit(isg[set == s, genes], ",")[[1]]
  data.table(isg_set = s, ModelID = rep(rownames(Z), length(g)), gene = rep(g, each = nrow(Z)), z = as.vector(Z[, g]))
}))
fwrite(gz, file.path(OUT, "depmap_lung_isg_gene_z.tsv.gz"), sep = "\t")

expr_genes <- intersect(c(TESTED, SENSORS), colnames(X))
pl <- data.table(ModelID = rownames(X), X[, expr_genes, drop = FALSE])
setnames(pl, expr_genes, paste0(expr_genes, "_expr"))
pl <- merge(merge(model, chr, by = "ModelID"), pl, by = "ModelID")
ID_COLS <- c("ModelID", "ADAR_chronos", "ADARB1_chronos", "ADAR_prob", "ADARB1_prob", "ADAR_dependent")
long <- melt(merge(pl[, ..ID_COLS], score, by = "ModelID"),
             id.vars = ID_COLS, variable.name = "isg_set", value.name = "isg")
long[, isg_set := factor(isg_set, names(gsets))]
long[, isg_tertile := cut(isg, quantile(isg, 0:3 / 3), labels = c("Low", "Mid", "High"), include.lowest = TRUE),
     by = isg_set]
fwrite(pl, file.path(OUT, "depmap_lung_per_line.tsv.gz"), sep = "\t")
fwrite(long, file.path(OUT, "depmap_lung_isg_scores.tsv.gz"), sep = "\t")

spear <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(list(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  list(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

# 1. dependency vs ISG score
cor_isg <- long[, rbindlist(lapply(c("ADAR_chronos", "ADARB1_chronos"), function(v)
  c(list(dependency = v), spear(isg, get(v))))), by = isg_set]
fwrite(cor_isg, file.path(OUT, "lung_chronos_vs_isg.tsv"), sep = "\t")

# 2. ISG tertiles; dependent = dependency probability >= PROB_CUT
ter <- long[is.finite(ADAR_chronos) & !is.na(ADAR_dependent),
            .(n = .N, median_chronos = median(ADAR_chronos),
              n_dependent = sum(ADAR_dependent), pct_dependent = 100 * mean(ADAR_dependent),
              n_effect_lt_cut = sum(ADAR_chronos < DEP_CUT)), keyby = .(isg_set, isg_tertile)]
hl <- long[isg_tertile %in% c("Low", "High") & is.finite(ADAR_chronos) & !is.na(ADAR_dependent),
           .(p_high_vs_low_fisher = fisher.test(table(isg_tertile == "High", ADAR_dependent))$p.value,
             p_high_vs_low_wilcox = wilcox.test(ADAR_chronos[isg_tertile == "High"], ADAR_chronos[isg_tertile == "Low"])$p.value),
           by = isg_set]
ter <- merge(ter, hl, by = "isg_set")
fwrite(ter, file.path(OUT, "lung_isg_tertiles.tsv"), sep = "\t")

# 2b. ISG score vs gene effect within dependent / non-dependent lines, and over a range of cut-offs
cls <- pl[is.finite(ADAR_chronos) & !is.na(ADAR_dependent),
          .N, by = .(prob_dependent = ADAR_dependent, effect_lt_cut = ADAR_chronos < DEP_CUT)]
fwrite(cls, file.path(OUT, "lung_dependency_class.tsv"), sep = "\t")
dl <- long[is.finite(ADAR_chronos) & !is.na(ADAR_dependent)]
by_dep <- dl[, c(list(variable = "isg"), spear(isg, ADAR_chronos)), by = .(isg_set, ADAR_dependent)]
by_dep_expr <- rbindlist(lapply(c(TRUE, FALSE), function(k) {
  d <- pl[ADAR_dependent == k & is.finite(ADAR_chronos)]
  rbindlist(lapply(expr_genes, function(g) c(list(ADAR_dependent = k, variable = g), spear(d[[paste0(g, "_expr")]], d$ADAR_chronos))))
}))
fwrite(by_dep, file.path(OUT, "lung_isg_vs_chronos_by_dependency.tsv"), sep = "\t")
fwrite(by_dep_expr, file.path(OUT, "lung_expr_vs_chronos_by_dependency.tsv"), sep = "\t")
scan <- rbind(
  rbindlist(lapply(seq(-0.3, -1.0, by = -0.1), function(k)
    dl[ADAR_chronos < k, c(list(measure = "gene_effect", cut = k), spear(isg, ADAR_chronos)), by = isg_set])),
  rbindlist(lapply(seq(0.3, 0.9, by = 0.1), function(k)
    dl[ADAR_prob >= k, c(list(measure = "probability", cut = k), spear(isg, ADAR_chronos)), by = isg_set])))
fwrite(scan, file.path(OUT, "lung_isg_vs_chronos_cut_scan.tsv"), sep = "\t")

# 3. dependency vs expression of single genes
cor_expr <- rbindlist(lapply(c("ADAR_chronos", "ADARB1_chronos"), function(v)
  rbindlist(lapply(expr_genes, function(g) c(list(dependency = v, gene = g), spear(pl[[paste0(g, "_expr")]], pl[[v]]))))))
cor_expr[, fdr := p.adjust(p, "BH"), by = dependency]
fwrite(cor_expr, file.path(OUT, "lung_chronos_vs_expression.tsv"), sep = "\t")

# 4. sensors as modifiers of the ISG-dependency link
sens <- pl[, c("ModelID", paste0(intersect(SENSORS, expr_genes), "_expr")), with = FALSE]
modif <- rbindlist(lapply(levels(long$isg_set), function(s) rbindlist(lapply(intersect(SENSORS, expr_genes), function(g) {
  d <- merge(long[isg_set == s], sens, by = "ModelID")[, sensor := as.numeric(scale(get(paste0(g, "_expr"))))]
  m <- summary(lm(ADAR_chronos ~ isg * sensor, data = d))$coefficients
  data.table(isg_set = s, sensor = g, term = c("isg", "sensor", "isg:sensor"),
             beta = m[c("isg", "sensor", "isg:sensor"), 1], p = m[c("isg", "sensor", "isg:sensor"), 4])
}))))
fwrite(modif, file.path(OUT, "lung_sensor_modifier_lm.tsv"), sep = "\t")

# 5. genome-wide expression correlates of ADAR dependency
d <- pl[is.finite(ADAR_chronos)]
M <- X[d$ModelID, , drop = FALSE]
M <- M[, apply(M, 2, sd) > 0 & colMeans(M > 1) >= 0.1, drop = FALSE]
rho <- cor(apply(M, 2, rank), rank(d$ADAR_chronos))[, 1]
n <- nrow(M)
tt <- rho * sqrt((n - 2) / (1 - rho^2))
gw <- data.table(gene = colnames(M), n = n, rho = rho, p = 2 * pt(-abs(tt), n - 2))
gw[, fdr := p.adjust(p, "BH")]
for (s in names(gsets)) set(gw, j = paste0("in__", s), value = gw$gene %in% gsets[[s]])
setorder(gw, rho)
fwrite(gw, file.path(OUT, "lung_genomewide_expr_vs_adar_chronos.tsv.gz"), sep = "\t")
enr <- rbindlist(lapply(names(gsets), function(s) {
  w <- suppressWarnings(wilcox.test(gw$rho[gw[[paste0("in__", s)]]], gw$rho[!gw[[paste0("in__", s)]]]))
  data.table(isg_set = s, n_genes = sum(gw[[paste0("in__", s)]]), median_rho_set = median(gw$rho[gw[[paste0("in__", s)]]]),
             median_rho_other = median(gw$rho[!gw[[paste0("in__", s)]]]), p_wilcox = w$p.value)
}))
fwrite(enr, file.path(OUT, "lung_genomewide_set_shift.tsv"), sep = "\t")

# 6. single genes of the ISG sets (and the sensors) vs ADAR Chronos, within dependent / non-dependent lines
ig <- setdiff(intersect(unique(c(unlist(gsets), SENSORS)), colnames(X)), EXCLUDE)
pg <- rbindlist(lapply(c(TRUE, FALSE), function(k) {
  d <- pl[ADAR_dependent == k & is.finite(ADAR_chronos)]
  M <- X[d$ModelID, ig, drop = FALSE]
  rbindlist(lapply(ig, function(g) c(list(ADAR_dependent = k, gene = g), spear(M[, g], d$ADAR_chronos))))
}))
pg[, fdr := p.adjust(p, "BH"), by = ADAR_dependent]
pg[, sensor := gene %in% SENSORS]
for (s in names(gsets)) set(pg, j = paste0("in__", s), value = pg$gene %in% gsets[[s]])
setorder(pg, -ADAR_dependent, rho)
fwrite(pg, file.path(OUT, "lung_isg_gene_vs_chronos_by_dependency.tsv"), sep = "\t")
# background: every expressed gene, within dependent / non-dependent lines; shift of the ISG genes
gwd <- rbindlist(lapply(c(TRUE, FALSE), function(k) {
  d <- pl[ADAR_dependent == k & is.finite(ADAR_chronos)]
  M <- X[d$ModelID, , drop = FALSE]
  M <- M[, apply(M, 2, sd) > 0 & colMeans(M > 1) >= 0.1, drop = FALSE]
  r <- cor(apply(M, 2, rank), rank(d$ADAR_chronos))[, 1]
  data.table(ADAR_dependent = k, gene = colnames(M), n = nrow(M), rho = r)
}))
gwd[, `:=`(isg = gene %in% ig, sensor = gene %in% SENSORS)]
fwrite(gwd, file.path(OUT, "lung_genomewide_expr_vs_adar_chronos_by_dependency.tsv.gz"), sep = "\t")
bg <- gwd[, .(n_genes = .N, n_isg = sum(isg), median_rho_isg = median(rho[isg]), median_rho_other = median(rho[!isg]),
              pct_neg_isg = 100 * mean(rho[isg] < 0), pct_neg_other = 100 * mean(rho[!isg] < 0),
              p_wilcox = wilcox.test(rho[isg], rho[!isg])$p.value,
              isg_pct_rank_median = 100 * mean(rho[!isg] < median(rho[isg]))), by = ADAR_dependent]
fwrite(bg, file.path(OUT, "lung_isg_vs_background_by_dependency.tsv"), sep = "\t")

message("\n== Chronos vs ISG"); print(cor_isg)
message("\n== ISG tertiles"); print(ter)
message("\n== dependency classes"); print(cls)
message("\n== ISG vs Chronos by dependency"); print(by_dep)
message("\n== expression vs Chronos by dependency"); print(by_dep_expr)
message("\n== cut-off scan"); print(scan)
message("\n== Chronos vs expression"); print(cor_expr)
message("\n== sensor modifiers (isg:sensor)"); print(modif[term == "isg:sensor"])
message("\n== genome-wide rho shift of each set"); print(enr)
message("\n== genome-wide top 20 (FDR < 0.1: ", gw[fdr < 0.1, .N], ")"); print(head(gw[, .(gene, rho, fdr)], 20))
message("\n== ISG / sensor genes vs Chronos, dependent lines, top 25 (FDR < 0.1: ", pg[ADAR_dependent == TRUE & fdr < 0.1, .N], ")")
print(head(pg[ADAR_dependent == TRUE, .(gene, sensor, n, rho, p, fdr)], 25))
print(pg[ADAR_dependent == TRUE & sensor == TRUE, .(gene, n, rho, p, fdr)])
message("\n== ISG genes vs all expressed genes, rho with ADAR Chronos"); print(bg, digits = 3)
