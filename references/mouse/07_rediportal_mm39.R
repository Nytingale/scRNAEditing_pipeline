# Mouse REDIportal (raerdata EH8233, mm10) -> mm39 via UCSC mm10ToMm39 chain.
# Keeps sites that lift 1:1 to a single base. Writes:
#   REDIportal_mm39.tsv  full annotation (Func.wgEncodeGencodeBasicVM16 = region)
#   REDIportal_mm39.vcf  minimal VCF for the pipeline's --rediportal
# Needs mm39.fa.fai (01). Requires R with GenomicRanges, rtracklayer, data.table.
#
# Usage: REFDIR=<dir> Rscript 07_rediportal_mm39.R
suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(data.table)
})
setwd(Sys.getenv("REFDIR", "."))
# ExperimentHub EH8233 (raerdata 1.0.0); the fetch URL takes the file id 8292, not the EH number.
# md5 ba646db389e78eecfa7af41476faabc0
if (!file.exists("rediportal_full_mm10.rda"))
  download.file("https://experimenthub.bioconductor.org/fetch/8292", "rediportal_full_mm10.rda", mode = "wb")
if (!file.exists("mm10ToMm39.over.chain.gz"))
  download.file("https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz",
                "mm10ToMm39.over.chain.gz", mode = "wb")

e = new.env(); load("rediportal_full_mm10.rda", envir = e)
gr = e$rediportal_full_mm10
mcols(gr)$mm10_pos = start(gr)
chain_txt = "mm10ToMm39.over.chain"
if (!file.exists(chain_txt)) writeLines(readLines(gzfile("mm10ToMm39.over.chain.gz")), chain_txt)
chain = import.chain(chain_txt)
lo = liftOver(gr, chain)
ok = elementNROWS(lo) == 1
gr39 = unlist(lo[ok])
gr39 = gr39[width(gr39) == 1]
fai = fread("mm39.fa.fai", header = FALSE, select = 1)$V1
gr39 = gr39[as.character(seqnames(gr39)) %in% fai]
cat(sprintf("mm10 sites %d; lifted 1:1 %d (%.2f%%)\n", length(gr), length(gr39),
            100 * length(gr39) / length(gr)))

dt = as.data.table(gr39)
setnames(dt, c("seqnames", "start"), c("chr", "pos"))
dt[, c("end", "width") := NULL]
dt = dt[order(chr, pos)]
dt = unique(dt, by = c("chr", "pos", "strand"))
fwrite(dt, "REDIportal_mm39.tsv", sep = "\t")

vcf = unique(dt[, .(`#CHROM` = chr, POS = pos, ID = ".", REF = Ref, ALT = Ed,
                    QUAL = ".", FILTER = ".", INFO = paste0("STRAND=", strand))],
             by = c("#CHROM", "POS"))
con = file("REDIportal_mm39.vcf", "w")
writeLines(c("##fileformat=VCFv4.2",
             "##source=REDIportal mouse (raerdata EH8233, mm10) lifted to mm39 with UCSC mm10ToMm39.over.chain"),
           con)
close(con)
fwrite(vcf, "REDIportal_mm39.vcf", sep = "\t", append = TRUE, col.names = TRUE)
cat(sprintf("VCF positions %d\n", nrow(vcf)))
print(dt[, .N, by = Func.wgEncodeGencodeBasicVM16][order(-N)])
