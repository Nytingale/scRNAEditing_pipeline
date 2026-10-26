#!/usr/bin/env Rscript
# assign_guides_paper.R — reproduce the authors' guide calling EXACTLY
# (scScreens_ADARs Figure6/F06_02_ADARKO_GuideCalling.R):
#   * total guide UMI >= 5
#   * assign cell to its single max-UMI guide; ties -> multiplet (dropped)
#   * drop guide ADARB1_1 (inadequate cutting efficiency)
#   * collapse guide_id -> target by prefix (ADAR1/ADARB1/ADARB2/AAVS1/NTC)
# Then join our annotation and report per-arm CELL + UMI depth.

suppressPackageStartupMessages({ library(Matrix); library(data.table) })
geo   <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/geo_guides")
outdir<- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm/seurat_analysis")
samples <- c(T1="GSM7923327", T2="GSM7923328", T3="GSM7923329", T4="GSM7923330")
MIN_GUIDE <- 5
fread_gz <- function(p, ...) fread(cmd=paste("zcat", shQuote(p)), ...)
strip_bc <- function(x) sub("-[0-9]+$", "", x)

per_cell <- list()
for (s in names(samples)) {
  g <- samples[[s]]
  bc <- fread_gz(file.path(geo, sprintf("%s_%s_barcodes.tsv.gz", g, s)), header=FALSE)$V1
  ft <- fread_gz(file.path(geo, sprintf("%s_%s_features.tsv.gz", g, s)), header=FALSE)
  m  <- readMM(gzfile(file.path(geo, sprintf("%s_%s_matrix.mtx.gz", g, s))))
  gidx <- which(ft$V3 == "CRISPR Guide Capture"); eidx <- which(ft$V3 == "Gene Expression")
  gex_umi <- Matrix::colSums(m[eidx, , drop=FALSE])
  gm <- as.matrix(m[gidx, , drop=FALSE]); rownames(gm) <- ft$V1[gidx]  # 20 x cells
  tot <- colSums(gm); top <- apply(gm, 2, max)
  gid <- rep(NA_character_, ncol(gm))
  for (it in seq_len(ncol(gm))) {
    if (tot[it] >= MIN_GUIDE) {
      w <- which(gm[, it] == top[it])
      gid[it] <- if (length(w) == 1) rownames(gm)[w] else "multiplet"
    }
  }
  per_cell[[s]] <- data.table(sample=s, barcode=strip_bc(bc), gex_umi=gex_umi,
                              tot_gUMI=tot, top_gUMI=top, guide_id=gid)
}
pc <- rbindlist(per_cell)
# drop ADARB1_1 guide (paper) + multiplets/unassigned
pc <- pc[!is.na(guide_id) & guide_id != "multiplet" & guide_id != "ADARB1_1"]
pc[, target := tstrsplit(guide_id, split="-|_", keep=1L)[[1]]]
pc[, arm := fcase(target=="ADAR1","ADAR1_KO", target=="AAVS1","AAVS1_ctrl",
                  target=="NTC","NTC_ctrl", target=="ADARB1","ADARB1_KO",
                  target=="ADARB2","ADARB2_KO", default=target)]

ann <- rbindlist(lapply(names(samples), function(s) {
  f <- file.path(outdir, sprintf("barcodes_celltype_scored_%s.tsv", s))
  if (file.exists(f)) fread(f)[, .(sample=s, barcode=strip_bc(barcode), celltype)] }), fill=TRUE)
j <- merge(pc, ann, by=c("sample","barcode"))   # annotated (human) cells only
fwrite(j[, .(sample, barcode, celltype, guide_id, target, arm, tot_gUMI, top_gUMI, gex_umi)],
       file.path(outdir, "ground_truth_cells_paper.tsv"), sep="\t")

fmtM <- function(x) sprintf("%.2fM", x/1e6)
armf <- c("ADAR1_KO","ADARB1_KO","ADARB2_KO","AAVS1_ctrl","NTC_ctrl")
cat("\n=== PAPER thresholds (gUMI>=5, max-guide, ADARB1_1 dropped) ===\n")
cat(sprintf("Assigned singlet cells (all): %d ; annotated: %d\n", nrow(pc), nrow(j)))
cat("\n--- MERGED across T1-T4 (annotated) ---\n")
print(j[arm %in% armf, .(cells=.N, total_UMI=fmtM(sum(gex_umi)),
        med_UMI=as.integer(median(gex_umi))), by=arm][order(match(arm,armf))])
cat("\n--- PER SAMPLE (cells) ---\n")
print(dcast(j[arm %in% armf, .N, by=.(arm,sample)], arm ~ sample, value.var="N", fill=0))
cat("\n--- ADAR1_KO vs AAVS1_ctrl vs NTC_ctrl per cell type (cells; merged) ---\n")
ct <- dcast(j[arm %in% c("ADAR1_KO","AAVS1_ctrl","NTC_ctrl"), .N, by=.(celltype,arm)],
            celltype ~ arm, value.var="N", fill=0); setorder(ct, -ADAR1_KO)
print(head(ct, 15))
cat("\nWrote ground_truth_cells_paper.tsv\n")
