#!/usr/bin/env Rscript
# REDIportal overlap summariser (human datasets).
# For every filtered pre100 site table produced by downstream_filter.R, count how
# many discovered sites fall on a known REDIportal editing position vs how many do
# not (novel), plus the fractions — emitted as one tidy table for figure plotting.
#
# Overlap is POSITION-level (chr:pos): "site is a known editing position". We
# report both ALL discovered varTypes and the A>G-only subset (canonical A-to-I).
# REDIportal is loaded ONCE; every pre100 file staged in the CWD is scored.
options(scipen = 999)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt = parse_args(OptionParser(option_list = list(
  make_option("--rediportal", type = "character",
              help = "REDIportal VCF(.gz): #CHROM POS ... (position-level known editing sites)"),
  make_option("--inputs_glob", type = "character", default = "*_pre100.tsv.gz",
              help = "Glob for pre100 site tables staged in CWD [default %default]"),
  make_option("--output", type = "character", default = "rediportal_overlap.tsv",
              help = "Output tidy TSV [default %default]")
)))

# ── Load REDIportal positions once ──────────────────────────────────────────
if (is.null(opt$rediportal) || !file.exists(opt$rediportal)) {
  stop("REDIportal file not found: ", opt$rediportal)
}
cat("Loading REDIportal:", opt$rediportal, "\n")
redi = fread(opt$rediportal, skip = "#CHROM", header = TRUE, select = 1:2,
             col.names = c("chr", "pos"))
redi[, chr := as.character(chr)]
redi[, pos := as.integer(pos)]
redi = unique(redi)
setkey(redi, chr, pos)
cat("  REDIportal unique positions:", nrow(redi), "\n")

# ── Score every pre100 table ────────────────────────────────────────────────
files = Sys.glob(opt$inputs_glob)
cat("Found", length(files), "pre100 files\n")

# Parse combo identity from filename:
#   <sample>_<caller>_<bam_type>_<strand_mode>_snp<off|on>_pre100.tsv.gz
# sample id may contain underscores (e.g. T1_ADAR1_KO), so anchor on the known
# tail tokens: caller, bam_type, strand_mode, snpX.
CALLERS = c("bcftools", "reditools3", "gatk", "jacusa2", "scallele")
parse_combo = function(fn) {
  base = sub("_pre100\\.tsv\\.gz$", "", basename(fn))
  toks = strsplit(base, "_", fixed = TRUE)[[1]]
  # snp token is last: snpoff/snpon
  snp_tok = toks[length(toks)]
  snp_sub = sub("^snp", "", snp_tok)
  strand_mode = toks[length(toks) - 1]
  bam_type    = toks[length(toks) - 2]
  # caller is the token before bam_type
  caller      = toks[length(toks) - 3]
  sample      = paste(toks[1:(length(toks) - 4)], collapse = "_")
  list(sample = sample, caller = caller, bam_type = bam_type,
       strand_mode = strand_mode, snp_sub = snp_sub)
}

summarise_one = function(fn) {
  meta = parse_combo(fn)
  dt = tryCatch(fread(fn, sep = "\t"), error = function(e) NULL)
  # Fallbacks for identity if file is empty / unreadable
  if (is.null(dt) || nrow(dt) == 0 || !all(c("chr", "pos") %in% names(dt))) {
    rows = rbindlist(lapply(c("all", "AG"), function(st) data.table(
      sample = meta$sample, caller = meta$caller, bam_type = meta$bam_type,
      strand_mode = meta$strand_mode, snp_sub = meta$snp_sub, set = st,
      n_total = 0L, n_rediportal = 0L, n_novel = 0L,
      frac_rediportal = NA_real_, frac_novel = NA_real_)))
    return(rows)
  }
  dt[, chr := as.character(chr)]
  dt[, pos := as.integer(pos)]
  # Sample id from the table itself is authoritative when present
  if ("sample" %in% names(dt) && nrow(dt) > 0) meta$sample = as.character(dt$sample[1])
  dt[, in_redi := FALSE]
  dt[redi, in_redi := TRUE, on = .(chr, pos)]

  mk = function(sub_dt, st) {
    n_total = nrow(sub_dt)
    n_redi  = sum(sub_dt$in_redi)
    n_novel = n_total - n_redi
    data.table(sample = meta$sample, caller = meta$caller, bam_type = meta$bam_type,
               strand_mode = meta$strand_mode, snp_sub = meta$snp_sub, set = st,
               n_total = n_total, n_rediportal = n_redi, n_novel = n_novel,
               frac_rediportal = if (n_total > 0) n_redi / n_total else NA_real_,
               frac_novel      = if (n_total > 0) n_novel / n_total else NA_real_)
  }
  ag = if ("varType" %in% names(dt)) dt[varType == "AG"] else dt[0]
  rbind(mk(dt, "all"), mk(ag, "AG"))
}

out = rbindlist(lapply(files, summarise_one))
setorder(out, sample, caller, bam_type, strand_mode, snp_sub, set)
fwrite(out, opt$output, sep = "\t")
cat("Wrote", nrow(out), "rows to", opt$output, "\n")
