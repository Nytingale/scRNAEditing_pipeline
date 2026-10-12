#!/usr/bin/env Rscript
options(scipen = 999)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
  library(GenomicRanges)
  library(janitor)
  library(Rcpp)
  library(stringr)
})

# ── Argument parsing ─────────────────────────────────────────────────────────
option_list = list(
  make_option("--input",          type = "character", help = "Input variants file (gzipped or plain text)"),
  make_option("--caller",         type = "character", help = "Caller: bcftools, reditools3, gatk, jacusa2, scallele, cellsnp"),
  make_option("--sample",         type = "character", help = "Sample identifier"),
  make_option("--bam_label",      type = "character", default = "", help = "BAM processing label"),
  make_option("--genome_fasta",   type = "character", help = "Reference genome FASTA"),
  make_option("--gtf",            type = "character", default = NULL, help = "Gene annotation GTF (unused, kept for interface)"),
  make_option("--snp_rdata",      type = "character", default = NULL, help = "Known SNPs as GRanges RData (snp_gr object)"),
  make_option("--snp_bed",        type = "character", default = NULL, help = "Known SNP BED file (legacy fallback)"),
  make_option("--snp_blacklist",  type = "character", default = NULL, help = "Per-sample genetic-variant blacklist (TSV: chr<TAB>pos, no header). Sites removed in addition to known SNPs. Populated by callers that genotype (scAllele PASS set; optionally GATK/BCFtools hom-ALT)."),
  make_option("--strand_annots",  type = "character", default = NULL, help = "RData with plus.annot/minus.annot GRanges"),
  make_option("--alu_annots",     type = "character", default = NULL, help = "RData with alu.plus.annot/alu.minus.annot GRanges"),
  make_option("--main_chroms",    type = "character", default = NULL, help = "File with one chromosome name per line"),
  make_option("--windows",        type = "character", default = "0,20,40,60,80,100,120,140,160,180,200",
              help = "Comma-separated window sizes"),
  make_option("--mismatch_window_filter", type = "integer", default = 1,
              help = paste("1 = pass 1 searches w1 over --windows (remove sites with >= 2 mismatch types within w1 bp);",
                           "0 = w1 fixed at 0, so only the clustered-sites pass (w2) is searched [default 1]")),
  make_option("--min_cov",        type = "integer",   default = 10, help = "Minimum total coverage"),
  make_option("--frac_weight",    type = "double",    default = 1.0, help = "Weight for AG fraction in scoring"),
  make_option("--output_filtered",type = "character", help = "Output filtered sites (gzipped TSV)"),
  make_option("--output_summary", type = "character", help = "Output summary statistics (plain text)"),
  make_option("--output_stepwise",type = "character", default = NULL, help = "Output stepwise mismatch counts (TSV)"),
  make_option("--output_pre100",  type = "character", default = NULL, help = "Output sites before 100%% editing removal (TSV)"),
  make_option("--save_intermediates", type = "integer", default = 0,
              help = "Save full site tables from this step onwards (1-7, 0=off). Step 1=multiMM, 2=SNP, 3=strand, 4=cov/edit, 5=noMM(w1), 6=clustered(w2), 7=no100pct"),
  make_option("--strand_mode", type = "character", default = "annotation",
              help = paste("Strand assignment mode: 'annotation' (gene annotation only),",
                           "'hybrid' (read-strand I16/SB first, annotation fallback), or",
                           "'caller' (keep the caller's OWN strand call, no annotation override)")),
  make_option("--strand_ambiguity_threshold", type = "double", default = 0.3,
              help = "For hybrid mode: if abs(fwd-rev)/(fwd+rev) < this, site is ambiguous and falls back to annotation"),
  make_option("--min_alt",        type = "integer",   default = 3, help = "Minimum alt allele read count (used when --alt_rule count)"),
  make_option("--min_edit_level", type = "double",    default = 0.01, help = "Minimum editing level (cov_alt / total)"),
  # ── alt-read support rule ─────────────────────────────────────────────────
  # `count` is the historical behaviour: cov_alt >= min_alt.
  #
  # >>> WHY THIS IS PARAMETERISED. <<<
  # `min_alt = 3` is the single most consequential threshold in the pipeline and
  # it does not look like one. Because it implies editing >= 3/coverage, the
  # EFFECTIVE editing floor is set by depth, not by --min_edit_level:
  #     coverage  16 (25th pct) -> 18.8% editing required
  #     coverage  36 (median)   ->  8.3%
  #     coverage 142 (90th pct) ->  2.1%
  # Meanwhile --min_edit_level 0.01 binds on only 2.58% of sites (it needs >300
  # reads to bind at all) and --min_cov 10 excludes 13.8%. Measured on the merged
  # teratoma NTC arm.
  # The cost: at well-covered sites (>=50 reads) the median editing level is just
  # 5.4%, and 76.9% of edited sites are below 10%. So at median depth the 3-read
  # rule discards ~71% of the edited sites that are visible where depth allows.
  #
  # `binom` is the principled alternative: keep a site when its edited-read count
  # is significantly above the sequencing error rate, which is what the fixed
  # count is groping at. At 36 reads, 2 edited reads is already highly
  # significant against a 0.1% error rate, so low-level editing is recovered
  # while single-read noise is still rejected (min_alt_floor).
  make_option("--alt_rule", type = "character", default = "count",
              help = "Alt-read support rule: 'count' (cov_alt >= min_alt) or 'binom' (one-sided binomial vs --seq_error_rate)"),
  make_option("--seq_error_rate", type = "double", default = 0.001,
              help = "For --alt_rule binom: per-base sequencing error rate under the null"),
  make_option("--alt_binom_alpha", type = "double", default = 0.05,
              help = "For --alt_rule binom: significance threshold for the alt-read count"),
  make_option("--min_alt_floor", type = "integer", default = 2,
              help = "For --alt_rule binom: absolute minimum alt reads, so a single read is never accepted"),
  make_option("--jacusa2_min_score", type = "double", default = 1.25, help = "Minimum JACUSA2 quality score"),
  make_option("--remove_100pct_editing", type = "logical", default = TRUE,
              help = "Remove sites with 100%% editing level [default TRUE]"),
  make_option("--output_metrics", type = "character", default = NULL,
              help = "Output single-row run metrics TSV for cross-run collation"),
  make_option("--snp_subtract_label", type = "character", default = "off",
              help = "Benchmark-axis label recorded in outputs: 'off' (dbSNP only), 'on' (dbSNP + per-sample blacklist), or 'NA' (caller cannot genotype). [default off]")
)
opt = parse_args(OptionParser(option_list = option_list))

# ── Constants ────────────────────────────────────────────────────────────────
BASES    = c("A", "C", "G", "T")
VARTYPES = as.character(outer(BASES, BASES, paste0))
VARTYPES = VARTYPES[!substr(VARTYPES, 1, 1) == substr(VARTYPES, 2, 2)]

COMP_MAP = c(A = "T", T = "A", C = "G", G = "C")
complement_vartype = function(vt) {
  paste0(COMP_MAP[substr(vt, 1, 1)], COMP_MAP[substr(vt, 2, 2)])
}

# REDItools3 `analyze` writes Reference, BaseCount[A,C,G,T] and AllSubs in the
# forward genomic frame even when Strand is "-". Genomic T>C on minus is RNA A>G.
# Swap counts A<->T, C<->G and complement ref + varType so the row is in the RNA
# frame of the current strand. Unstranded (2/*) stay mapped to "+" (plus frame).
reditools3_complement_frame = function(df, idx) {
  if (!length(idx)) return(df)
  tmp = df$counts_A[idx]
  df$counts_A[idx] = df$counts_T[idx]
  df$counts_T[idx] = tmp
  tmp = df$counts_C[idx]
  df$counts_C[idx] = df$counts_G[idx]
  df$counts_G[idx] = tmp
  ref = df$ref[idx]
  hit = !is.na(ref) & ref %in% names(COMP_MAP)
  df$ref[idx[hit]] = unname(COMP_MAP[ref[hit]])
  vt = df$varType[idx]
  hit = !is.na(vt) & nchar(vt) == 2L
  df$varType[idx[hit]] = complement_vartype(vt[hit])
  df
}

# Keep counts and varType in the RNA frame of the assigned strand.
apply_strand_flip = function(df, needs_flip) {
  if (!any(needs_flip)) return(df)
  if (identical(opt$caller, "reditools3")) {
    df = reditools3_complement_frame(df, which(needs_flip))
  } else {
    df$varType[needs_flip] = complement_vartype(df$varType[needs_flip])
  }
  df
}

numThreads = min(parallel::detectCores(), 8)
windows    = as.numeric(unlist(strsplit(opt$windows, ",")))
w1_windows = if (opt$mismatch_window_filter == 1) windows else 0

summary_lines = c()
log_msg = function(...) {
  msg = paste0(...)
  cat(msg, "\n")
  summary_lines <<- c(summary_lines, msg)
}

write_empty_stepwise = function() {
  if (!is.null(opt$output_stepwise)) {
    empty_df = data.frame(varType = character(), n = integer(), percent = numeric(),
                          totalSites = integer(), setType = character(),
                          sample = character(), caller = character(),
                          bam_type = character(), strand_mode = character(),
                          step_order = integer())
    fwrite(empty_df, opt$output_stepwise, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(opt$output_pre100)) {
    fwrite(data.table(), opt$output_pre100, sep = "\t")
  }
  if (!is.null(opt$output_metrics)) {
    empty_metrics = data.frame(
      sample = opt$sample, caller = opt$caller, bam_type = opt$bam_label,
      strand_mode = opt$strand_mode, snp_subtract = opt$snp_subtract_label,
      best_w1 = NA_integer_, best_w2 = NA_integer_,
      sites_step1_multiMM = 0L, sites_step2_noSNP = 0L,
      sites_step3_stranded = 0L, sites_step4_covEdit = 0L,
      sites_step5_noMM = 0L, sites_step6_clustered = 0L,
      sites_step7_no100pct = 0L,
      ag_n_step1 = 0L, ag_pct_step1 = 0, ag_n_step2 = 0L, ag_pct_step2 = 0,
      ag_n_step3 = 0L, ag_pct_step3 = 0, ag_n_step4 = 0L, ag_pct_step4 = 0,
      ag_n_step5 = 0L, ag_pct_step5 = 0, ag_n_step6 = 0L, ag_pct_step6 = 0,
      ag_n_step7 = 0L, ag_pct_step7 = 0,
      n_strand_from_read = NA_integer_, n_strand_from_annot = NA_integer_,
      stringsAsFactors = FALSE)
    fwrite(empty_metrics, opt$output_metrics, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

out_prefix = tools::file_path_sans_ext(tools::file_path_sans_ext(opt$output_filtered))
save_intermediate = function(df, step_num, step_label) {
  if (opt$save_intermediates > 0 && step_num >= opt$save_intermediates && nrow(df) > 0) {
    fname = paste0(out_prefix, "_step", step_num, "_", step_label, ".tsv.gz")
    fwrite(df, fname, sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("  [intermediate] Saved ", nrow(df), " sites to ", fname)
  }
}

collate_vartypes_tbl = function(df, setTypeName) {
  if (nrow(df) == 0) {
    return(data.frame(varType = VARTYPES, n = 0L, percent = 0,
                      totalSites = 0L, setType = setTypeName,
                      stringsAsFactors = FALSE))
  }
  tbl = as.data.frame(janitor::tabyl(df, varType))
  missing = VARTYPES[!VARTYPES %in% tbl$varType]
  if (length(missing) > 0) {
    tbl = rbind(tbl, data.frame(varType = missing, n = 0L, percent = 0))
  }
  tbl = tbl[order(tbl$varType), ]
  tbl$totalSites = sum(tbl$n)
  tbl$setType    = setTypeName
  tbl
}

# ── Compile C++ sitesInWindow function ───────────────────────────────────────
cpp_code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void filReliable(int window, CharacterVector chr, NumericVector pos, CharacterVector site, CharacterVector type,
                 IntegerVector sites_in_wind_n, IntegerVector types_in_wind_n) {
  for (int i=0; i<chr.size(); i++){
    int start=pos[i]-window;
    int end=pos[i]+window;
    CharacterVector listsite;
    CharacterVector listtype;
    int j0;
    int jn;
    if (i>=window*4){j0=i-window*4;} else {j0=0;}
    if (i+window*4<=chr.size()){jn=i+window*4;} else {jn=chr.size();}
    for (int j=j0; j<jn;j++){
      if ((chr[j]==chr[i]) && (pos[j]>=start) && (pos[j]<=end)){
        listsite.push_back(site[j]);
        listtype.push_back(type[j]);
      }
    }
    CharacterVector listsite_uniq=unique(listsite);
    CharacterVector listtype_uniq=unique(listtype);
    sites_in_wind_n[i]=listsite_uniq.size();
    types_in_wind_n[i]=listtype_uniq.size();
  }
}
'
cpp_tmpfile = tempfile(fileext = ".cpp")
writeLines(cpp_code, cpp_tmpfile)
Rcpp::sourceCpp(cpp_tmpfile)

# ── Helper: window filtering wrappers ────────────────────────────────────────
filterForReliableSites = function(sample.df, window_size) {
  sample.df$sites_in_wind_n      = NA_integer_
  sample.df$uniqueTypesInWindow   = NA_integer_
  sample.df$pos = as.numeric(sample.df$pos)
  sample.df = sample.df[order(sample.df$chr, sample.df$pos), ]
  filReliable(window_size, sample.df$chr, sample.df$pos,
              sample.df$site, sample.df$varType,
              sample.df$sites_in_wind_n, sample.df$uniqueTypesInWindow)
  return(sample.df)
}

applyFilterForReliableSites = function(df, window_size) {
  result_plus  = NULL
  result_minus = NULL
  df_plus  = df[df$strand == "+", ]
  df_minus = df[df$strand == "-", ]
  if (nrow(df_plus)  > 0) result_plus  = filterForReliableSites(df_plus,  window_size)
  if (nrow(df_minus) > 0) result_minus = filterForReliableSites(df_minus, window_size)
  rbind(result_plus, result_minus)
}

# ── Helper: gene-type ranking for strand annotation ─────────────────────────
processAnnotatedBED = function(df) {
  colnames(df) = c("chr", "pos", "gene_name", "strand", "gene_type")
  df$gene_type[is.na(df$gene_type)] = "NA"
  type_names = c("NA", "gene", "transcript", "intergenic", "intron",
                 "5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR")
  ranks = c(0, 0, 0, 1, 2, rep(3, 6))
  map = setNames(ranks, type_names)
  df$ranktype = map[df$gene_type]
  df$ranktype[is.na(df$ranktype)] = 0
  df = as.data.table(df)
  df = df[, .(genetype  = paste(unique(gene_type), collapse = ","),
              genename  = paste(unique(gene_name), collapse = ","),
              ranktype  = max(ranktype)),
          by = .(chr, pos)]
  df$genetype = stringr::str_replace_na(df$genetype)
  df$genename = stringr::str_replace_na(df$genename)
  as.data.frame(df)
}

# ── STEP 0: Read input ──────────────────────────────────────────────────────
log_msg("=== Downstream filter: ", opt$caller, " | sample: ", opt$sample,
        " | strand_mode: ", opt$strand_mode, " ===")
log_msg("Input: ", opt$input)

if (!file.exists(opt$input)) {
  log_msg("ERROR: Input file not found: ", opt$input)
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

file_info = file.info(opt$input)
if (file_info$size == 0) {
  log_msg("Input file is empty.")
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

# bcftools rows below min_cov are dropped first thing in STEP 1; dropping them
# while streaming keeps deep arms (Replogle NTC, >100 M rows) out of RAM.
df = tryCatch(
  if (opt$caller == "bcftools") {
    fread(cmd = sprintf("zcat -f %s | awk -F'\\t' '$5+0 >= %s'",
                        shQuote(opt$input), opt$min_cov),
          data.table = FALSE, sep = "\t")
  } else {
    fread(opt$input, data.table = FALSE, sep = "\t")
  },
  error = function(e) {
    log_msg("ERROR reading input: ", e$message)
    data.frame()
  }
)
if (nrow(df) == 0) {
  log_msg("No rows in input.")
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

# ── STEP 1: Parse caller-specific input ─────────────────────────────────────
log_msg("STEP 1: Parsing ", opt$caller, " output (", nrow(df), " rows)")

if (opt$caller == "bcftools") {
  # bcftools query -f output: chr, pos, ref, alt, dp, ad, I16
  # AD field is comma-separated: ref_count,alt_count[,other...]
  # I16: 16 comma-separated floats; first 4 are ref_fwd,ref_rev,alt_fwd,alt_rev
  colnames(df) = c("chr", "pos", "ref", "alt", "cov", "ad_field", "i16_field")
  df$site   = paste0(df$chr, ":", df$pos)
  df$strand = "+"
  df$cov    = as.numeric(df$cov)
  df        = df[df$cov >= opt$min_cov, ]
  df$strand_orig = df$strand
  df$alt    = gsub(",<*>", "", df$alt, fixed = TRUE)
  df$varType = paste0(df$ref, df$alt)

  # Parse AD: first value = ref count, second = alt count
  ad_split = str_split_fixed(df$ad_field, ",", 3)
  df$cov_ref  = as.numeric(ad_split[, 1])
  df$cov_alt  = as.numeric(ad_split[, 2])
  df$ad_field = NULL

  # Parse I16 for read-strand evidence (only used in hybrid mode)
  if (opt$strand_mode == "hybrid" && nrow(df) > 0) {
    i16_split = str_split_fixed(df$i16_field, ",", 16)
    df$read_ref_fwd = suppressWarnings(as.numeric(i16_split[, 1]))
    df$read_ref_rev = suppressWarnings(as.numeric(i16_split[, 2]))
    df$read_alt_fwd = suppressWarnings(as.numeric(i16_split[, 3]))
    df$read_alt_rev = suppressWarnings(as.numeric(i16_split[, 4]))
    # Missing I16 (reported as ".") becomes NA -> treat as 0
    df$read_ref_fwd[is.na(df$read_ref_fwd)] = 0
    df$read_ref_rev[is.na(df$read_ref_rev)] = 0
    df$read_alt_fwd[is.na(df$read_alt_fwd)] = 0
    df$read_alt_rev[is.na(df$read_alt_rev)] = 0
  }
  df$i16_field = NULL
  # Base data.frame rejects recycling length-1 into 0 rows (tiny cell-type
  # BAMs often lose every site at min_cov). Empty frames fall through to the
  # shared nrow==0 exit below and write empty outputs with status 0.
  if (nrow(df) > 0) df$sample = opt$sample

} else if (opt$caller == "reditools3") {
  # REDItools3 `analyze` native output = 14 columns WITH a header row:
  #   Region Position Reference Strand Coverage MeanQ BaseCount[A,C,G,T] AllSubs
  #   Frequency gCoverage gMeanQ gBaseCount[A,C,G,T] gAllSubs gFrequency
  # BaseCount is a single bracketed, space-padded field, e.g. "[5, 0, 8, 0]" =
  # per-base A,C,G,T read counts in the SAME FORWARD GENOMIC frame as Reference
  # (verified: char1 of AllSubs == Reference for every single-sub row), including
  # when Strand is "-". Split it into counts_A/C/G/T, derive genomic varType from
  # ref + the dominant non-ref base (robust to multi-allelic AllSubs like "AG AT";
  # mirrors the JACUSA2 branch), then convert minus-strand rows into the RNA frame.
  # A legacy pre-split 10-col schema (chr,pos,ref,strand,cov,A,C,G,T,varType) is
  # still accepted for back-compat.
  if (ncol(df) >= 14) {
    df = df[, 1:8]
    colnames(df) = c("chr", "pos", "ref", "strand", "cov",
                     "meanq", "basecount", "allsubs")
    bc = gsub("[][[:space:]]", "", df$basecount)   # "[5, 0, 8, 0]" -> "5,0,8,0"
    bc = str_split_fixed(bc, ",", 4)
    df$counts_A = suppressWarnings(as.numeric(bc[, 1]))
    df$counts_C = suppressWarnings(as.numeric(bc[, 2]))
    df$counts_G = suppressWarnings(as.numeric(bc[, 3]))
    df$counts_T = suppressWarnings(as.numeric(bc[, 4]))
    df$meanq = NULL; df$basecount = NULL; df$allsubs = NULL
  } else {
    colnames(df) = c("chr", "pos", "ref", "strand", "cov",
                     "counts_A", "counts_C", "counts_G", "counts_T", "varType")
  }
  df$cov = suppressWarnings(as.numeric(df$cov))
  # varType = ref + dominant non-ref base (from per-base counts)
  cnt_mat = as.matrix(df[, c("counts_A", "counts_C", "counts_G", "counts_T")])
  colnames(cnt_mat) = BASES
  ref_idx = match(df$ref, BASES)
  alt_mat = cnt_mat
  ok_rows = which(!is.na(ref_idx))
  alt_mat[cbind(ok_rows, ref_idx[ok_rows])] = 0        # zero out the ref base
  alt_idx  = max.col(alt_mat, ties.method = "first")
  alt_base = BASES[alt_idx]
  alt_base[rowSums(alt_mat, na.rm = TRUE) == 0] = NA_character_   # no alt reads
  df$varType = ifelse(is.na(ref_idx) | is.na(alt_base),
                      NA_character_, paste0(df$ref, alt_base))
  df$strand = gsub("0", "-", df$strand, fixed = TRUE)
  df$strand = gsub("1", "+", df$strand, fixed = TRUE)
  # "U" (unstranded, value 2 or *) counts are in the reference/plus frame,
  # so we treat them as "+" from the start — simplifies downstream flip logic
  df$strand = gsub("2", "+", df$strand, fixed = TRUE)
  df$strand = gsub("*", "+", df$strand, fixed = TRUE)
  minus = which(df$strand == "-")
  if (length(minus) > 0L) {
    df = reditools3_complement_frame(df, minus)
    log_msg("  REDItools3: RNA-frame complement for ", length(minus),
            " minus-strand rows (ref, counts A<->T C<->G, varType)")
  }
  df$strand_orig = df$strand
  df$site   = paste0(df$chr, ":", df$pos)
  if (nrow(df) > 0) df$sample = opt$sample

} else if (opt$caller == "gatk") {
  # GATK output: chr, pos, ref, alt, cov_ref, cov_alt, cov, sb_field
  # SB = ref_fwd,ref_rev,alt_fwd,alt_rev (from StrandBiasBySample)
  colnames(df) = c("chr", "pos", "ref", "alt", "cov_ref", "cov_alt", "cov", "sb_field")
  df$cov_ref    = as.numeric(df$cov_ref)
  df$cov_alt    = as.numeric(df$cov_alt)
  df$cov        = as.numeric(df$cov)
  df$site       = paste0(df$chr, ":", df$pos)
  df$strand     = "+"
  df$strand_orig = df$strand
  df$varType    = paste0(df$ref, df$alt)

  # Parse SB for read-strand evidence (only used in hybrid mode)
  if (opt$strand_mode == "hybrid" && nrow(df) > 0) {
    sb_split = str_split_fixed(df$sb_field, ",", 4)
    df$read_ref_fwd = suppressWarnings(as.numeric(sb_split[, 1]))
    df$read_ref_rev = suppressWarnings(as.numeric(sb_split[, 2]))
    df$read_alt_fwd = suppressWarnings(as.numeric(sb_split[, 3]))
    df$read_alt_rev = suppressWarnings(as.numeric(sb_split[, 4]))
    df$read_ref_fwd[is.na(df$read_ref_fwd)] = 0
    df$read_ref_rev[is.na(df$read_ref_rev)] = 0
    df$read_alt_fwd[is.na(df$read_alt_fwd)] = 0
    df$read_alt_rev[is.na(df$read_alt_rev)] = 0
  }
  df$sb_field = NULL
  if (nrow(df) > 0) df$sample = opt$sample

} else if (opt$caller == "jacusa2") {
  # Header-only jacusa2 (tiny cell-type BAMs): default fread treats the ## version
  # banner as the sole column name and the #contig header as the only "data" row,
  # so bases11 is missing and str_split_fixed crashes with
  # "replacement has 0 items, need 4". Re-read skipping the banner; 0 call rows
  # then fall through to the shared empty handler.
  if (!"bases11" %in% names(df)) {
    log_msg("  jacusa2: re-reading with skip=1 (##-banner / header-only mishandle)")
    df = tryCatch(
      fread(opt$input, data.table = FALSE, sep = "\t", skip = 1L),
      error = function(e) data.frame()
    )
  }
  if (nrow(df) == 0L || !"bases11" %in% names(df)) {
    df = data.frame()
  } else {
  colnames(df)[1] = "chr"
  df$site   = paste0(df$chr, ":", df$end)
  df[, c("cts_A", "cts_C", "cts_G", "cts_T")] = str_split_fixed(df$bases11, ",", 4)
  df$cts_A = as.numeric(df$cts_A); df$cts_C = as.numeric(df$cts_C)
  df$cts_G = as.numeric(df$cts_G); df$cts_T = as.numeric(df$cts_T)
  df$cov = df$cts_A + df$cts_C + df$cts_G + df$cts_T
  df = df[df$cov >= opt$min_cov, ]
  if (nrow(df) == 0) {
    # fall through to the shared empty handler after the caller switch
  } else {
  df$sample = opt$sample
  df$strand = gsub(".", "+", df$strand, fixed = TRUE)
  df$strand_orig = df$strand

  counts_mat = as.matrix(df[, c("cts_A", "cts_C", "cts_G", "cts_T")])
  colnames(counts_mat) = BASES
  rows    = seq_len(nrow(df))
  ref_idx = match(df$ref, colnames(counts_mat))
  alt_mat = counts_mat
  alt_mat[cbind(rows, ref_idx)] = 0
  n_alts  = rowSums(alt_mat > 0)
  alt_idx = max.col(alt_mat, ties.method = "first")

  df$alt     = ifelse(n_alts == 1, colnames(counts_mat)[alt_idx], NA_character_)
  df$cov_alt = ifelse(n_alts == 1, alt_mat[cbind(rows, alt_idx)], NA_integer_)
  df$cov_ref = counts_mat[cbind(rows, ref_idx)]

  df$varType = paste0(df$ref, df$alt)

  # JACUSA2 score filter
  if ("score" %in% colnames(df)) {
    df = df[!is.na(df$score) & df$score >= opt$jacusa2_min_score, ]
  }

  # Rename end -> pos and drop extra columns
  df$pos = df$end
  df = df[, intersect(c("site","chr","pos","ref","alt","strand","strand_orig",
                         "cov","cov_ref","cov_alt","varType","sample"),
                       colnames(df))]
  } # end nrow>0 after min_cov
  } # end has bases11 / call rows

} else if (opt$caller == "scallele") {
  # scAllele VCF parsed to: chr, pos, ref, alt, cov_ref (RC), cov_alt (AC), cov (DP)
  # SNVs only; no I16/SB — hybrid strand_mode falls back to annotation.
  colnames(df) = c("chr", "pos", "ref", "alt", "cov_ref", "cov_alt", "cov")
  df$cov_ref     = as.numeric(df$cov_ref)
  df$cov_alt     = as.numeric(df$cov_alt)
  df$cov         = as.numeric(df$cov)
  df$site        = paste0(df$chr, ":", df$pos)
  df$strand      = "+"
  df$strand_orig = df$strand
  df$varType     = paste0(df$ref, df$alt)
  if (nrow(df) > 0) df$sample = opt$sample

} else if (opt$caller == "cellsnp") {
  # CELLSNP_LITE TSV: chr, pos, ref, alt, cov_ref (DP - AD), cov_alt (AD), cov (DP + OTH).
  # Genomic frame, no strand or read-strand counts: annotation stranding only.
  colnames(df) = c("chr", "pos", "ref", "alt", "cov_ref", "cov_alt", "cov")
  df$cov_ref     = as.numeric(df$cov_ref)
  df$cov_alt     = as.numeric(df$cov_alt)
  df$cov         = as.numeric(df$cov)
  df             = df[df$cov >= opt$min_cov, ]
  df$site        = paste0(df$chr, ":", df$pos)
  df$strand      = rep("+", nrow(df))
  df$strand_orig = df$strand
  df$varType     = paste0(df$ref, df$alt)
  if (nrow(df) > 0) df$sample = opt$sample

} else {
  stop("Unknown caller: ", opt$caller)
}

if (nrow(df) == 0) {
  log_msg("No rows after caller-specific parsing.")
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}
log_msg("  After parsing: ", nrow(df), " sites")

# ── STEP 2: Filter to main chromosomes (mandatory) ─────────────────────────
log_msg("STEP 2: Filtering to main chromosomes")
if (!is.null(opt$main_chroms) && file.exists(opt$main_chroms)) {
  main_chrs = readLines(opt$main_chroms)
  main_chrs = trimws(main_chrs[nchar(trimws(main_chrs)) > 0])
  log_msg("  Using ", length(main_chrs), " chromosomes from --main_chroms file")
} else {
  main_chrs = unique(df$chr)
  log_msg("  No --main_chroms file, using all ", length(main_chrs), " chromosomes in data")
}
df = df[df$chr %in% main_chrs, ]
log_msg("  After chr filter: ", nrow(df), " sites")
if (nrow(df) == 0) {
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

# ── STEP 3: Remove multi-mismatch sites ─────────────────────────────────────
log_msg("STEP 3: Removing multi-mismatch sites")
n_before = nrow(df)
df = df[df$varType %in% VARTYPES, ]
log_msg("  Removed ", n_before - nrow(df), " multi-mismatch sites; ", nrow(df), " remain")

tbl.raw = collate_vartypes_tbl(df, "After_multiMM_removal")
save_intermediate(df, 1, "multiMM_removed")

if (nrow(df) == 0) {
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

# ── STEP 4: Remove known SNPs ───────────────────────────────────────────────
log_msg("STEP 4: Removing known SNPs")

snp_gr = NULL
if (!is.null(opt$snp_rdata) && file.exists(opt$snp_rdata)) {
  e_snp = new.env()
  load(opt$snp_rdata, envir = e_snp)
  if ("snp_gr" %in% ls(e_snp)) {
    snp_gr = e_snp$snp_gr
    log_msg("  Loaded SNP GRanges: ", length(snp_gr), " positions from RData")
  } else {
    for (nm in ls(e_snp)) {
      obj = e_snp[[nm]]
      if (is.data.frame(obj) && "site" %in% colnames(obj)) {
        log_msg("  Loaded legacy SNP data.frame '", nm, "' with ", nrow(obj), " rows, converting to GRanges")
        parts = data.table::tstrsplit(obj$site, ":", fixed = TRUE)
        snp_gr = GRanges(seqnames = parts[[1]], ranges = IRanges(start = as.integer(parts[[2]]), end = as.integer(parts[[2]])))
        snp_gr = sort(unique(snp_gr))
        break
      }
    }
  }
  rm(e_snp); gc()
} else if (!is.null(opt$snp_bed) && file.exists(opt$snp_bed)) {
  snp_df = fread(opt$snp_bed, header = FALSE, select = c(1, 3))
  colnames(snp_df) = c("chr", "pos")
  snp_gr = GRanges(seqnames = snp_df$chr, ranges = IRanges(start = snp_df$pos, end = snp_df$pos))
  snp_gr = sort(unique(snp_gr))
  log_msg("  Loaded SNP BED: ", length(snp_gr), " positions")
  rm(snp_df); gc()
}

if (!is.null(snp_gr)) {
  n_before = nrow(df)
  df_gr = GRanges(seqnames = df$chr, ranges = IRanges(start = as.integer(df$pos), end = as.integer(df$pos)))
  hits = queryHits(findOverlaps(df_gr, snp_gr, type = "equal"))
  if (length(hits) > 0) {
    df = df[-hits, ]
  }
  log_msg("  Removed ", n_before - nrow(df), " SNP sites; ", nrow(df), " remain")
  rm(snp_gr, df_gr); gc()
} else {
  log_msg("  No SNP file provided or found, skipping")
}
# ── STEP 4b: Remove per-sample genetic-variant blacklist ────────────────────
# Sample-specific SNP set produced by callers that genotype/classify. For
# scAllele this is its PASS set (GLM trained without editing -> PASS = genetic;
# A-to-I is always LowQual, so this subtraction is editing-safe). Supplements
# static dbSNP, catching somatic/rare variants absent from any database.
if (!is.null(opt$snp_blacklist) && file.exists(opt$snp_blacklist) &&
    file.info(opt$snp_blacklist)$size > 0 && nrow(df) > 0) {
  bl = tryCatch(fread(opt$snp_blacklist, header = FALSE, sep = "\t"),
                error = function(e) NULL)
  if (!is.null(bl) && ncol(bl) >= 2 && nrow(bl) > 0) {
    setnames(bl, 1:2, c("chr", "pos"))
    bl_gr = GRanges(seqnames = as.character(bl$chr),
                    ranges = IRanges(start = as.integer(bl$pos),
                                     end = as.integer(bl$pos)))
    bl_gr = sort(unique(bl_gr))
    n_before = nrow(df)
    df_gr = GRanges(seqnames = df$chr,
                    ranges = IRanges(start = as.integer(df$pos),
                                     end = as.integer(df$pos)))
    hits = queryHits(findOverlaps(df_gr, bl_gr, type = "equal"))
    if (length(hits) > 0) df = df[-hits, ]
    log_msg("STEP 4b: Removed ", n_before - nrow(df),
            " per-sample blacklist sites (", length(bl_gr),
            " in list); ", nrow(df), " remain")
    rm(bl_gr, df_gr); gc()
  }
} else {
  log_msg("STEP 4b: No per-sample SNP blacklist provided, skipping")
}

tbl.noSNP = collate_vartypes_tbl(df, "SNPs_removed")
save_intermediate(df, 2, "SNPs_removed")

if (nrow(df) == 0) {
  log_msg("No sites remain after SNP removal.")
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

# ── STEP 5: Strand assignment ────────────────────────────────────────────────
log_msg("STEP 5: Strand assignment (mode: ", opt$strand_mode, ")")

# --- 5a: Read-based strand inference (hybrid mode only) ----------------------
read_strand_sites = NULL
if (opt$strand_mode == "hybrid" && opt$caller %in% c("bcftools", "gatk") &&
    "read_alt_fwd" %in% colnames(df)) {
  log_msg("  Hybrid: inferring strand from read-level counts (I16/SB)")
  thresh = opt$strand_ambiguity_threshold

  # Use ALT allele strand counts: alt_fwd vs alt_rev
  total_alt = df$read_alt_fwd + df$read_alt_rev
  # Strand imbalance ratio: (fwd - rev) / (fwd + rev), ranges from -1 to +1
  # Positive = mostly forward-strand reads, Negative = mostly reverse-strand reads
  imbalance = ifelse(total_alt > 0,
                     (df$read_alt_fwd - df$read_alt_rev) / total_alt,
                     0)

  # Classify: clear forward, clear reverse, or ambiguous
  df$read_strand = ifelse(imbalance > thresh, "+",
                   ifelse(imbalance < -thresh, "-", "U"))

  # Sites with total_alt == 0 are ambiguous by definition
  df$read_strand[total_alt == 0] = "U"

  n_plus  = sum(df$read_strand == "+")
  n_minus = sum(df$read_strand == "-")
  n_ambig = sum(df$read_strand == "U")
  log_msg("  Read-strand: +strand=", n_plus, ", -strand=", n_minus,
          ", ambiguous=", n_ambig, " (threshold=", thresh, ")")

  # For non-ambiguous sites, record them; they skip annotation lookup
  read_strand_sites = df$site[df$read_strand != "U"]
}

# --- 5b: Annotation-based strand assignment ----------------------------------
# >>> 'caller' MODE BYPASSES THIS ENTIRELY. <<<
# The point of caller-native stranding is to report what a user actually gets
# from the tool, so the annotation must NOT override it. df$strand is already
# the caller's own call (df$strand_orig, set in the per-caller parse above), so
# this branch only has to label it and leave it alone.
#
# NO SITES ARE DROPPED, matching the original misc_scripts/reditools_ByRedi.R:
# that script maps REDItools strand 2/* to "U" and then never filters on strand
# — the label is informational. Our parse maps 2/* to "+" instead, which is
# equivalent in effect, because REDItools reports counts in the reference/plus
# frame when it cannot determine a strand. Note this is CONSERVATIVE, not
# generous: a genuine minus-strand A>G event at an unstranded site appears as
# T>C in the plus frame and is therefore excluded from the A>G set.
#
# This also applies to UNSTRANDED libraries on purpose. A caller run on an
# unstranded library returns no strand, so this mode degenerates to "everything
# plus" and will misassign mismatch types — which is exactly the result a user
# would get, and is itself the argument for annotation-based stranding as the
# default. Gating it away would hide that finding.
#
# No further complement here: strand == strand_orig by construction. Minus-strand
# REDItools3 rows were already converted to the RNA frame at parse.
if (identical(opt$strand_mode, "caller")) {
  log_msg("  Caller-native stranding: keeping ", opt$caller, "'s own strand call, ",
          "no annotation override, no sites dropped")
  if (nrow(df) > 0) {
    df$strand = df$strand_orig
    df$strand_source = "caller"
  }
  n_plus  = sum(df$strand == "+", na.rm = TRUE)
  n_minus = sum(df$strand == "-", na.rm = TRUE)
  log_msg("  Caller strand: +strand=", n_plus, ", -strand=", n_minus,
          " (unstranded calls were mapped to '+' at parse time)")
} else if (!is.null(opt$strand_annots) && file.exists(opt$strand_annots)) {
  load(opt$strand_annots)  # loads plus.annot, minus.annot

  # In hybrid mode, only look up annotation for ambiguous sites
  if (opt$strand_mode == "hybrid" && !is.null(read_strand_sites)) {
    annot_df = df[!df$site %in% read_strand_sites, ]
    log_msg("  Annotation lookup for ", nrow(annot_df), " sites (ambiguous by read-strand)")
  } else {
    annot_df = df
    log_msg("  Annotation lookup for all ", nrow(annot_df), " sites")
  }

  if (nrow(annot_df) > 0) {
    # ── Annotation overlap, CHUNKED BY SITE ─────────────────────────────────
    # >>> PEAK MEMORY MUST NOT SCALE WITH TOTAL SITE COUNT. <<<
    # findOverlaps returns one hit PER (site, annotation feature) pair, so
    # `sub.gr[queryHits(hits)]` expands the GRanges from one row per site to one
    # row per overlap — and gene/transcript/exon features overlap heavily, so the
    # expansion is large. `as.data.frame()` then `unique()` each take a further
    # full copy of that expanded table, and the whole thing is done twice (plus
    # and minus). Peak was therefore several copies of a table many times larger
    # than the input.
    #
    # MEASURED: the Replogle merged NTC arm with reditools3 reached 31,238,406
    # sites after SNP removal and was OOM-killed here at MEM_LIMIT_G 128.00 on
    # the 4th and final attempt of the memory ladder, which flipped
    # errorStrategy to `finish` and DRAINED the whole run (09-15 18:57). A second
    # task died the same way at 24,673,054 sites. Raising the ladder only buys
    # headroom; this removes the scaling.
    #
    # WHY CHUNKING IS EXACT, NOT AN APPROXIMATION. processAnnotatedBED()
    # aggregates strictly `by = .(chr, pos)`, and findOverlaps is run against the
    # FULL annotation for every chunk, so a site's complete feature set is always
    # within its own chunk. Per-chunk voting therefore gives byte-identical rows
    # to voting over everything at once, and `unique()` per chunk equals
    # `unique()` overall because duplicate (site, gene, strand, type) rows can
    # only arise within one site. Results are consumed by
    # merge(by = c("chr","pos")), so row order does not matter either.
    # The one requirement is that a site never spans two chunks — asserted below.
    CHUNK_SITES = suppressWarnings(as.integer(Sys.getenv("DOWNSTREAM_ANNOT_CHUNK", "2000000")))
    if (is.na(CHUNK_SITES) || CHUNK_SITES < 1L) CHUNK_SITES = 2000000L

    # Contiguous row blocks keep a site intact only if each site is one row.
    # That is the normal case; if it ever is not, fall back to a single block,
    # which is the original behaviour and cannot split a site.
    if (anyDuplicated(annot_df$site)) {
      log_msg("  NOTE: annot_df has repeated sites; using one block so no site is split")
      CHUNK_SITES = nrow(annot_df)
    }

    # One strand's overlap -> per-site vote. Returns the same columns as
    # processAnnotatedBED so the callers below are unchanged.
    vote_overlaps = function(gr, annot) {
      hits = findOverlaps(gr, annot)
      if (length(hits) == 0L) {
        return(data.frame(chr = character(0), pos = numeric(0),
                          genetype = character(0), genename = character(0),
                          ranktype = numeric(0), stringsAsFactors = FALSE))
      }
      matched = gr[queryHits(hits)]
      mcols(matched) = cbind.data.frame(mcols(matched),
                                        mcols(annot[subjectHits(hits)]))
      d = as.data.frame(matched)[, c("seqnames", "end", "id", "strand", "geneType")]
      colnames(d) = c("chr", "pos", "gene_name", "strand", "gene_type")
      # character, not the seqnames factor: chunks would otherwise carry
      # different level sets into the rbind.
      d$chr = as.character(d$chr)
      processAnnotatedBED(unique(d))
    }

    n_sites = nrow(annot_df)
    block_starts = seq(1L, n_sites, by = CHUNK_SITES)
    plus_parts  = vector("list", length(block_starts))
    minus_parts = vector("list", length(block_starts))
    if (length(block_starts) > 1L)
      log_msg("  Chunking annotation lookup: ", n_sites, " sites in ",
              length(block_starts), " block(s) of <=", CHUNK_SITES)

    for (k in seq_along(block_starts)) {
      rows = block_starts[k]:min(block_starts[k] + CHUNK_SITES - 1L, n_sites)
      pos_k = as.numeric(annot_df$pos[rows])
      sub.gr = GRanges(annot_df$chr[rows], IRanges(pos_k, pos_k))
      plus_parts[[k]]  = vote_overlaps(sub.gr, plus.annot)
      minus_parts[[k]] = vote_overlaps(sub.gr, minus.annot)
      rm(sub.gr, pos_k); invisible(gc(FALSE))
      if (length(block_starts) > 1L)
        log_msg("    block ", k, "/", length(block_starts), ": ", length(rows),
                " sites -> +", nrow(plus_parts[[k]]), " / -", nrow(minus_parts[[k]]),
                " voted")
    }

    gr.plus.voted  = as.data.frame(rbindlist(plus_parts,  use.names = TRUE, fill = TRUE))
    gr.minus.voted = as.data.frame(rbindlist(minus_parts, use.names = TRUE, fill = TRUE))
    rm(plus_parts, minus_parts); invisible(gc(FALSE))

    # Merge plus and minus
    u.merge = merge(gr.plus.voted, gr.minus.voted,
                    by = c("chr", "pos"), all = TRUE, suffixes = c(".p", ".m"))
    u.merge$site = paste0(u.merge$chr, ":", u.merge$pos)

    # Assign strand by rank comparison
    u.merge$ranktype.p[is.na(u.merge$ranktype.p)] = -1
    u.merge$ranktype.m[is.na(u.merge$ranktype.m)] = -1
    u.merge$genename.p[is.na(u.merge$genename.p)] = "NA"
    u.merge$genename.m[is.na(u.merge$genename.m)] = "NA"
    u.merge$genetype.p[is.na(u.merge$genetype.p)] = "NA"
    u.merge$genetype.m[is.na(u.merge$genetype.m)] = "NA"

    u.merge$new_strand = ifelse(u.merge$ranktype.m > u.merge$ranktype.p, "-",
                          ifelse(u.merge$ranktype.m < u.merge$ranktype.p, "+", "U"))

    # Break ties using gene names
    u_rows = which(u.merge$new_strand == "U")
    if (length(u_rows) > 0) {
      has_plus_gene  = !grepl("^NA$", u.merge$genename.p[u_rows]) & u.merge$genename.m[u_rows] == "NA"
      has_minus_gene = !grepl("^NA$", u.merge$genename.m[u_rows]) & u.merge$genename.p[u_rows] == "NA"
      u.merge$new_strand[u_rows[has_plus_gene]]  = "+"
      u.merge$new_strand[u_rows[has_minus_gene]] = "-"
    }

    # Try Alu annotations for remaining "U" sites
    still_u = u.merge$site[u.merge$new_strand == "U"]
    if (length(still_u) > 0 && !is.null(opt$alu_annots) && file.exists(opt$alu_annots)) {
      load(opt$alu_annots)  # loads alu.plus.annot, alu.minus.annot
      alu_sites = annot_df[annot_df$site %in% still_u, ]
      if (nrow(alu_sites) > 0) {
        alu.bed = data.frame(chrom = alu_sites$chr,
                              start = as.numeric(alu_sites$pos),
                              end   = as.numeric(alu_sites$pos))
        alu.gr = with(alu.bed, GRanges(chrom, IRanges(start, end)))

        alu.plus.hits  = findOverlaps(alu.gr, alu.plus.annot)
        alu.minus.hits = findOverlaps(alu.gr, alu.minus.annot)

        plus_sites  = unique(alu_sites$site[queryHits(alu.plus.hits)])
        minus_sites = unique(alu_sites$site[queryHits(alu.minus.hits)])

        alu_plus_only  = setdiff(plus_sites, minus_sites)
        alu_minus_only = setdiff(minus_sites, plus_sites)

        u.merge$new_strand[u.merge$site %in% alu_plus_only]  = "+"
        u.merge$new_strand[u.merge$site %in% alu_minus_only] = "-"
      }
      rm(alu.plus.annot, alu.minus.annot); gc()
    }

    annot_stranded = u.merge[u.merge$new_strand != "U", c("site", "new_strand")]
    n_annot_unstranded = sum(u.merge$new_strand == "U")
    log_msg("  Annotation: stranded=", nrow(annot_stranded),
            ", unstranded=", n_annot_unstranded)
  } else {
    annot_stranded = data.frame(site = character(), new_strand = character())
  }

  # --- 5c: Combine read-strand and annotation-strand assignments ---------------
  if (opt$strand_mode == "hybrid" && "read_strand" %in% colnames(df)) {
    # Read-strand assignments for non-ambiguous sites.
    # Tiny BAMs can leave every remaining site ambiguous (all U) — assigning a
    # length-1 "read" into a 0-row frame crashes base R the same way empty
    # df$sample did. Mirror the annot_stranded empty guard below.
    read_assigned = df[df$read_strand != "U", c("site", "read_strand")]
    colnames(read_assigned)[2] = "new_strand"
    if (nrow(read_assigned) > 0) {
      read_assigned$strand_source = "read"
    } else {
      read_assigned$strand_source = character()
    }

    # Annotation assignments for ambiguous sites
    if (nrow(annot_stranded) > 0) {
      annot_stranded$strand_source = "annotation"
    } else {
      annot_stranded$strand_source = character()
    }

    # Combine: read-strand takes priority (non-ambiguous sites never went to annotation)
    all_stranded = rbind(read_assigned, annot_stranded)
    n_read  = sum(all_stranded$strand_source == "read")
    n_annot = sum(all_stranded$strand_source == "annotation")
    log_msg("  Hybrid total: ", n_read, " from reads, ", n_annot, " from annotation")
  } else {
    all_stranded = annot_stranded
    if (nrow(all_stranded) > 0) all_stranded$strand_source = "annotation"
  }

  # Apply strand assignments and drop unresolved sites
  n_unstranded = nrow(df) - nrow(df[df$site %in% all_stranded$site, ])
  log_msg("  Unstranded sites dropped: ", n_unstranded)

  df = df[df$site %in% all_stranded$site, ]
  if (nrow(df) > 0) {
    df$strand = all_stranded$new_strand[match(df$site, all_stranded$site)]
    df$strand_source = all_stranded$strand_source[match(df$site, all_stranded$site)]

    # Flip into the newly assigned strand's RNA frame
    df = apply_strand_flip(df, df$strand != df$strand_orig)
  }
  rm(plus.annot, minus.annot); gc()
  if (exists("u.merge")) rm(u.merge)
  if (exists("all_stranded")) rm(all_stranded)
} else {
  log_msg("  No strand annotations provided, keeping original strand")
  # In hybrid mode without annotations, use read-strand only
  if (opt$strand_mode == "hybrid" && "read_strand" %in% colnames(df)) {
    log_msg("  Using read-strand only (no annotation file)")
    df = df[df$read_strand != "U", ]
    if (nrow(df) > 0) {
      df$strand = df$read_strand
      df$strand_source = "read"
      df = apply_strand_flip(df, df$strand != df$strand_orig)
    }
  }
}

log_msg("  After strand assignment: ", nrow(df), " sites")
strand_label = if (identical(opt$strand_mode, "hybrid")) {
  "Stranded_hybrid"
} else if (identical(opt$strand_mode, "caller")) {
  "Stranded_by_caller"
} else {
  "Stranded_by_annotations"
}
tbl.stranded = collate_vartypes_tbl(df, strand_label)
save_intermediate(df, 3, "stranded")

if (nrow(df) == 0) {
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

# RNA-frame cov_alt/cov_ref for reditools3 (counts match assigned strand after parse + any flip)
if (opt$caller == "reditools3" && nrow(df) > 0) {
  var_base = substr(df$varType, 2, 2)
  ref_base = substr(df$varType, 1, 1)
  df$cov_alt = ifelse(var_base == "A", df$counts_A,
               ifelse(var_base == "C", df$counts_C,
               ifelse(var_base == "G", df$counts_G,
               ifelse(var_base == "T", df$counts_T, NA_real_))))
  df$cov_ref = ifelse(ref_base == "A", df$counts_A,
               ifelse(ref_base == "C", df$counts_C,
               ifelse(ref_base == "G", df$counts_G,
               ifelse(ref_base == "T", df$counts_T, NA_real_))))
}

# ── STEP 6: Coverage and editing level filters ──────────────────────────────
log_msg("STEP 6: Coverage and editing level filters")
df$cov_alt = as.numeric(df$cov_alt)
df$cov_ref = as.numeric(df$cov_ref)

n_before = nrow(df)
if (identical(opt$alt_rule, "binom")) {
  # One-sided binomial: is cov_alt more edited reads than sequencing error alone
  # would produce at this depth?  P(X >= cov_alt | n, p = seq_error_rate)
  # `lower.tail = FALSE` at cov_alt - 1 gives P(X >= cov_alt).
  n_tot = df$cov_ref + df$cov_alt
  df$alt_pvalue = pbinom(df$cov_alt - 1, size = n_tot,
                         prob = opt$seq_error_rate, lower.tail = FALSE)
  keep = !is.na(df$cov_alt) & df$cov_alt >= opt$min_alt_floor &
         !is.na(df$alt_pvalue) & df$alt_pvalue <= opt$alt_binom_alpha
  df = df[keep, ]
  log_msg("  Removed ", n_before - nrow(df), " sites failing the binomial alt rule",
          " (error=", opt$seq_error_rate, ", alpha=", opt$alt_binom_alpha,
          ", floor=", opt$min_alt_floor, ")")
} else {
  df = df[!is.na(df$cov_alt) & df$cov_alt >= opt$min_alt, ]
  log_msg("  Removed ", n_before - nrow(df), " sites with cov_alt < ", opt$min_alt,
          "  [NB this implies editing >= ", opt$min_alt, "/coverage, which is the",
          " binding sensitivity floor at typical depth — see --alt_rule]")
}

# Editing level: consistent across all callers — cov_alt / (cov_ref + cov_alt)
# This ignores reads supporting other bases (e.g., C or T at an A>G site).
df$edit_level = df$cov_alt / (df$cov_ref + df$cov_alt)
n_before = nrow(df)
df = df[!is.na(df$edit_level) & df$edit_level >= opt$min_edit_level, ]
log_msg("  Removed ", n_before - nrow(df), " sites with edit_level < ", opt$min_edit_level)
log_msg("  After coverage/edit filters: ", nrow(df), " sites")
tbl.editCov = collate_vartypes_tbl(df, "Cov_and_edit_filtered")
save_intermediate(df, 4, "cov_edit_filtered")

if (nrow(df) == 0) {
  log_msg("No sites remain after coverage/edit filters.")
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

# ── STEP 7: Two-pass adaptive window filtering ──────────────────────────────
log_msg("STEP 7: Two-pass adaptive window filtering")
log_msg("  Window sizes: ", paste(windows, collapse = ","))
log_msg("  Nearby-mismatch filter (pass 1): ", if (opt$mismatch_window_filter == 1) "on, w1 searched" else "off, w1 = 0")

# Pre-compute sitesInWindow for each window size, per chromosome
work_cols = c("site", "chr", "pos", "varType", "strand", "sample")
df_work = df[, intersect(work_cols, colnames(df))]
df_work$pos = as.numeric(df_work$pos)
split_chroms = split(df_work, df_work$chr)

siw_cache = list()
for (w in w1_windows) {
  cat("  Computing sitesInWindow for window =", w, "...\n")
  if (w == 0) {
    tmp = df_work
    tmp$sites_in_wind_n    = 1L
    tmp$uniqueTypesInWindow = 1L
  } else {
    tmp = rbindlist(mclapply(split_chroms, function(subdf) {
      applyFilterForReliableSites(subdf, w)
    }, mc.cores = numThreads))
    tmp = as.data.frame(tmp)
  }
  siw_cache[[as.character(w)]] = tmp
}

# Test all w1 x w2 combinations
cat("  Testing window combinations...\n")
combination_results = list()
idx = 1
for (w1 in w1_windows) {
  df_w1 = siw_cache[[as.character(w1)]]
  if (w1 != 0) {
    nomm = df_w1[df_w1$uniqueTypesInWindow < 2, ]
  } else {
    nomm = df_w1
  }
  if (nrow(nomm) == 0) {
    for (w2 in windows) {
      combination_results[[idx]] = data.frame(n = 0L, percent = 0,
                                               window_sizes = paste0(w1, ",", w2),
                                               stringsAsFactors = FALSE)
      idx = idx + 1
    }
    next
  }

  nomm_split = split(nomm[, intersect(c("site","chr","pos","varType","strand","sample"), colnames(nomm))],
                      nomm$chr)

  for (w2 in windows) {
    if (w2 != 0) {
      min2 = rbindlist(mclapply(nomm_split, function(subdf) {
        applyFilterForReliableSites(subdf, w2)
      }, mc.cores = max(numThreads %/% 4, 1)))
      min2 = as.data.frame(min2)
      min2 = min2[min2$sites_in_wind_n >= 2, ]
    } else {
      min2 = nomm
    }
    if (nrow(min2) > 0) {
      tbl = as.data.frame(janitor::tabyl(min2, varType))
      if ("AG" %in% tbl$varType) {
        ag_row = tbl[tbl$varType == "AG", ]
        combination_results[[idx]] = data.frame(n = ag_row$n,
                                                 percent = ag_row$percent,
                                                 window_sizes = paste0(w1, ",", w2),
                                                 stringsAsFactors = FALSE)
      } else {
        combination_results[[idx]] = data.frame(n = 0L, percent = 0,
                                                 window_sizes = paste0(w1, ",", w2),
                                                 stringsAsFactors = FALSE)
      }
    } else {
      combination_results[[idx]] = data.frame(n = 0L, percent = 0,
                                               window_sizes = paste0(w1, ",", w2),
                                               stringsAsFactors = FALSE)
    }
    idx = idx + 1
  }
}
combination_df = rbindlist(combination_results) |> as.data.frame()
combination_df$sample = opt$sample

# Score combinations
n_weight = 1 - opt$frac_weight
if (length(unique(combination_df$percent)) == 1) {
  combination_df$scaled_n       = scale(combination_df$n)
  combination_df$scaled_percent = combination_df$scaled_n
} else {
  combination_df$scaled_n       = scale(combination_df$n)
  combination_df$scaled_percent = scale(combination_df$percent)
}
combination_df$scaled_n[is.na(combination_df$scaled_n)]             = 0
combination_df$scaled_percent[is.na(combination_df$scaled_percent)] = 0
combination_df$weighted_score = n_weight * combination_df$scaled_n +
                                 opt$frac_weight * combination_df$scaled_percent

best_row = combination_df[which.max(combination_df$weighted_score), ]
best_w1  = as.numeric(strsplit(best_row$window_sizes, ",")[[1]][1])
best_w2  = as.numeric(strsplit(best_row$window_sizes, ",")[[1]][2])
log_msg("  Best window combination: w1=", best_w1, ", w2=", best_w2,
        " (AG n=", best_row$n, ", AG%=", round(best_row$percent * 100, 2), "%)")

# Apply best window combination
df_w1_best = siw_cache[[as.character(best_w1)]]
if (best_w1 != 0) {
  nomm_best = df_w1_best[df_w1_best$uniqueTypesInWindow < 2, ]
} else {
  nomm_best = df_w1_best
}
log_msg("  After pass 1 (noMM, w1=", best_w1, "): ", nrow(nomm_best), " sites")
tbl.noMM = collate_vartypes_tbl(nomm_best, "Removed_variable_regions")
save_intermediate(nomm_best, 5, "variable_regions_removed")

if (nrow(nomm_best) == 0) {
  log_msg("No sites remain after pass 1.")
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

if (best_w2 != 0) {
  nomm_split = split(nomm_best[, intersect(c("site","chr","pos","varType","strand","sample"), colnames(nomm_best))],
                      nomm_best$chr)
  min2_best = rbindlist(mclapply(nomm_split, function(subdf) {
    applyFilterForReliableSites(subdf, best_w2)
  }, mc.cores = numThreads)) |> as.data.frame()
  min2_best = min2_best[min2_best$sites_in_wind_n >= 2, ]
} else {
  min2_best = nomm_best
}
log_msg("  After pass 2 (min2sites, w2=", best_w2, "): ", nrow(min2_best), " sites")
tbl.min2sites = collate_vartypes_tbl(min2_best, "Clustered_sites")
save_intermediate(min2_best, 6, "clustered")

# Merge back to full df to get all columns
if (nrow(min2_best) > 0) {
  keep_sites = unique(min2_best$site)
  df = df[df$site %in% keep_sites, ]
} else {
  df = df[0, ]
}

rm(siw_cache, split_chroms, df_work); gc()

if (nrow(df) == 0) {
  log_msg("No sites remain after window filtering.")
  fwrite(data.table(), opt$output_filtered, sep = "\t")
  writeLines(summary_lines, opt$output_summary)
  write_empty_stepwise()
  quit(status = 0)
}

# ── STEP 8: Remove 100% edited sites ────────────────────────────────────────
# Save pre-100% output (this is what the manuscript reports as the final result)
if (!is.null(opt$output_pre100) && nrow(df) > 0) {
  fwrite(df, opt$output_pre100, sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("  Pre-100% removal sites written to: ", opt$output_pre100)
}

log_msg("STEP 8: Removing 100% edited sites")
n_before = nrow(df)
if (opt$remove_100pct_editing) {
  df = df[df$edit_level < 1, ]
  log_msg("  Removed ", n_before - nrow(df), " fully edited sites; ", nrow(df), " remain")
} else {
  log_msg("  Skipped (--remove_100pct_editing FALSE); ", nrow(df), " sites retained")
}
tbl.no100 = collate_vartypes_tbl(df, "No_100pct_editing")

# ── STEP 9: Output ──────────────────────────────────────────────────────────
log_msg("STEP 9: Writing output")

# Clean up temporary read-strand columns
for (col in c("read_strand", "read_ref_fwd", "read_ref_rev",
              "read_alt_fwd", "read_alt_rev", "i16_field", "sb_field")) {
  if (col %in% colnames(df)) df[[col]] = NULL
}

log_msg("  Final sites: ", nrow(df))

if (nrow(df) > 0) {
  ag_final = sum(df$varType == "AG")
  ag_pct   = round(ag_final / nrow(df) * 100, 2)
  log_msg("  Final A>G count: ", ag_final, " (", ag_pct, "%)")
}

# Write filtered output
fwrite(df, opt$output_filtered, sep = "\t", quote = FALSE, row.names = FALSE)

# Write summary
summary_lines = c(summary_lines, "",
  "=== Mismatch type tables at each step ===")

write_tbl = function(tbl, label) {
  summary_lines <<- c(summary_lines, paste0("\n--- ", label, " ---"))
  for (i in seq_len(nrow(tbl))) {
    summary_lines <<- c(summary_lines,
      paste0(tbl$varType[i], "\t", tbl$n[i], "\t", round(tbl$percent * 100, 2)[i], "%"))
  }
  summary_lines <<- c(summary_lines, paste0("Total: ", tbl$totalSites[1]))
}

write_tbl(tbl.raw,       "After multi-mismatch removal")
write_tbl(tbl.noSNP,     "After SNP removal")
write_tbl(tbl.stranded,  "After strand assignment")
write_tbl(tbl.editCov,   "After coverage/edit filters")
write_tbl(tbl.noMM,      "After variable region removal (w1)")
write_tbl(tbl.min2sites, "After cluster requirement (w2)")
write_tbl(tbl.no100,     "After 100% editing removal")

summary_lines = c(summary_lines, "",
  paste0("Best_window_combination\t", best_w1, ",", best_w2),
  paste0("Best_AG_n\t", best_row$n),
  paste0("Best_AG_percent\t", round(best_row$percent * 100, 4)))

writeLines(summary_lines, opt$output_summary)

# ── Write stepwise mismatch count table (machine-readable) ──────────────────
if (!is.null(opt$output_stepwise)) {
  all_tbls = rbind(tbl.raw, tbl.noSNP,
                   tbl.stranded, tbl.editCov,
                   tbl.noMM, tbl.min2sites, tbl.no100)
  all_tbls$sample      = opt$sample
  all_tbls$caller      = opt$caller
  all_tbls$bam_type    = opt$bam_label
  all_tbls$strand_mode = opt$strand_mode
  all_tbls$step_order = match(all_tbls$setType,
    c("After_multiMM_removal", "SNPs_removed", strand_label,
      "Cov_and_edit_filtered", "Removed_variable_regions",
      "Clustered_sites", "No_100pct_editing"))
  fwrite(all_tbls, opt$output_stepwise, sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("  Stepwise mismatch table written to: ", opt$output_stepwise)
}

# ── Write single-row run metrics (machine-readable, for cross-run collation) ─
if (!is.null(opt$output_metrics)) {
  extract_ag = function(tbl) {
    row = tbl[tbl$varType == "AG", ]
    if (nrow(row) == 0) return(list(n = 0L, pct = 0))
    list(n = row$n, pct = row$percent)
  }

  ag_raw       = extract_ag(tbl.raw)
  ag_noSNP     = extract_ag(tbl.noSNP)
  ag_stranded  = extract_ag(tbl.stranded)
  ag_editCov   = extract_ag(tbl.editCov)
  ag_noMM      = extract_ag(tbl.noMM)
  ag_min2      = extract_ag(tbl.min2sites)
  ag_no100     = extract_ag(tbl.no100)

  # Count strand sources if available.
  # 'caller' is counted explicitly: without it, caller-native stranding reported
  # read=0 and annotation=0, which reads as "no strand was assigned" when in fact
  # every site was stranded by the caller itself.
  n_strand_read   = if ("strand_source" %in% colnames(df)) sum(df$strand_source == "read", na.rm = TRUE) else NA_integer_
  n_strand_annot  = if ("strand_source" %in% colnames(df)) sum(df$strand_source == "annotation", na.rm = TRUE) else NA_integer_
  n_strand_caller = if ("strand_source" %in% colnames(df)) sum(df$strand_source == "caller", na.rm = TRUE) else NA_integer_

  metrics = data.frame(
    sample       = opt$sample,
    caller       = opt$caller,
    bam_type     = opt$bam_label,
    strand_mode  = opt$strand_mode,
    snp_subtract = opt$snp_subtract_label,
    min_cov      = opt$min_cov,
    min_alt      = opt$min_alt,
    min_edit_level = opt$min_edit_level,
    jacusa2_min_score = opt$jacusa2_min_score,
    remove_100pct = opt$remove_100pct_editing,
    best_w1      = best_w1,
    best_w2      = best_w2,
    sites_step1_multiMM   = tbl.raw$totalSites[1],
    sites_step2_noSNP     = tbl.noSNP$totalSites[1],
    sites_step3_stranded  = tbl.stranded$totalSites[1],
    sites_step4_covEdit   = tbl.editCov$totalSites[1],
    sites_step5_noMM      = tbl.noMM$totalSites[1],
    sites_step6_clustered = tbl.min2sites$totalSites[1],
    sites_step7_no100pct  = tbl.no100$totalSites[1],
    ag_n_step1   = ag_raw$n,       ag_pct_step1   = ag_raw$pct,
    ag_n_step2   = ag_noSNP$n,     ag_pct_step2   = ag_noSNP$pct,
    ag_n_step3   = ag_stranded$n,  ag_pct_step3   = ag_stranded$pct,
    ag_n_step4   = ag_editCov$n,   ag_pct_step4   = ag_editCov$pct,
    ag_n_step5   = ag_noMM$n,      ag_pct_step5   = ag_noMM$pct,
    ag_n_step6   = ag_min2$n,      ag_pct_step6   = ag_min2$pct,
    ag_n_step7   = ag_no100$n,     ag_pct_step7   = ag_no100$pct,
    n_strand_from_read   = n_strand_read,
    n_strand_from_annot  = n_strand_annot,
    n_strand_from_caller = n_strand_caller,
    stringsAsFactors = FALSE
  )
  fwrite(metrics, opt$output_metrics, sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("  Run metrics written to: ", opt$output_metrics)
}

log_msg("Done.")
