#!/usr/bin/env Rscript
# =============================================================================
# gt_score_germline_v8.R   (2026-10-13)
# Germline false-positive rate (user 2026-10-13), a truth-independent companion to the v8 F1.
#   unit         as gt_score_v7.R: arm = NTC_ctrl, sample = T1-T4_NTC_ctrl, celltype = NTC_ctrl__<ct>
#                (GT_UNIT_KEEP: the 11 kept cell types)
#   covered      germline sites (gt_v8_germline_sites.R) with >= GT_UNIT_MIN_COV (10) reads on their read
#                orientation in the unit (gt_v8_unit_cov.sh, OUT = gt_eval_v8/germline/unit_cov)
#   germline FP  a call with RNA-frame varType AG and strand = the site's read orientation at a covered site
#   unit rate    germline FP / covered sites (NA if no covered site); also for GT 1/1 and non-1/1 sites
#   germ_frac    germline FP / A-to-G sites called in the unit (share of the discovered A-to-G sites that are
#                germline SNVs; NA if none called); tier value = mean of the unit values
#   tier rate    mean of the unit rates
# Call sets: the scored state (default mismatch = results_gt{,_persample,_celltype}_v8/*_step5_nomm), all
# before the 100%-editing removal (downstream_filter.R step 7), so homozygous sites are not removed by it.
# Usage : Rscript gt_score_germline_v8.R <tier> [results_dir] [out_dir] [call_regex]
# Output: <out_dir>/germline_unit_scores.tsv, germline_scores.tsv (config x tier)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
TIER <- args[1]
PM <- file.path(Sys.getenv("SCRNA_DIR", "."), "human/pm")
G <- file.path(PM, "gt_eval_v8", "germline")
UNIT_RE <- list(arm = "^NTC_ctrl$", sample = "^T[0-9]+_NTC_ctrl$", celltype = "^NTC_ctrl__")
stopifnot(TIER %in% names(UNIT_RE))
RES <- if (length(args) >= 2) args[2] else file.path(PM, c(arm = "results_gt_v8", sample = "results_gt_persample_v8", celltype = "results_gt_celltype_v8")[[TIER]])
OUT <- if (length(args) >= 3) args[3] else file.path(G, Sys.getenv("SCN", "scoring"), "mismatch", TIER)
CALL_RE <- if (length(args) >= 4) args[4] else "_step5_nomm\\.tsv\\.gz$"
MIN_COV <- as.integer(Sys.getenv("GT_UNIT_MIN_COV", "10"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
COMBO <- c("caller", "bam_type", "strand_mode", "snp_subtract")

gs <- fread(file.path(G, "germline_sites.tsv.gz"), select = c("site", "read_strand", "gt"))
gs[, hom := gt == "1/1"]
ucov <- fread(file.path(G, "unit_cov", "ko_truth_unit_cov.tsv.gz"), select = c("unit", "site", "coverage"))

units <- grep(UNIT_RE[[TIER]], list.dirs(RES, recursive = FALSE, full.names = FALSE), value = TRUE)
if (nzchar(Sys.getenv("GT_UNIT_KEEP")) && TIER == "celltype") {
  keep <- readLines(Sys.getenv("GT_UNIT_KEEP")); stopifnot(all(keep %in% units)); units <- intersect(units, keep)
}
stopifnot(length(units) > 0, all(units %chin% ucov$unit))
cov <- merge(ucov[unit %chin% units & coverage >= MIN_COV, .(unit, site)], gs, by = "site")
den <- cov[, .(n_cov = .N, n_cov_hom = sum(hom), n_cov_het = sum(!hom)), by = unit]

files <- rbindlist(lapply(units, function(u) {
  f <- list.files(file.path(RES, u, "04_filtered"), CALL_RE, recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NULL)
  p <- tstrsplit(sub(paste0("^", file.path(RES, u, "04_filtered"), "/"), "", f), "/", fixed = TRUE)
  data.table(unit = u, caller = p[[1]], bam_type = p[[2]], strand_mode = p[[3]], snp_subtract = sub("^snp_", "", p[[4]]), file = f)
}))
files[, config := do.call(paste, c(.SD, sep = "__")), .SDcols = COMBO]
message("tier ", TIER, ": ", length(units), " units, ", nrow(files), " call files, ", uniqueN(files$config), " configurations")

n_ag <- integer(nrow(files))
hits <- rbindlist(lapply(seq_len(nrow(files)), function(i) {
  if (file.size(files$file[i]) == 0) return(NULL)
  x <- fread(files$file[i], select = c("site", "varType", "strand"), showProgress = FALSE)
  x <- unique(x[varType == "AG", .(site, strand)])
  n_ag[i] <<- nrow(x)
  x <- x[site %chin% gs$site]
  if (!nrow(x)) return(NULL)
  x[, `:=`(unit = files$unit[i], config = files$config[i])]
}))
files[, n_ag := n_ag]
hits <- merge(hits, cov, by = c("unit", "site"))[strand == read_strand]
fp <- hits[, .(n_fp = uniqueN(site), n_fp_hom = uniqueN(site[hom]), n_fp_het = uniqueN(site[!hom]),
               n_fp_gt01 = uniqueN(site[gt == "0/1"]), n_fp_gtdot1 = uniqueN(site[gt == "./1"]),
               n_fp_gt1 = uniqueN(site[gt == "1"])), by = .(unit, config)]
stopifnot(all(gs$gt %in% c("1/1", "0/1", "./1", "1")))

us <- CJ(unit = units, config = unique(files$config))
us <- merge(us, unique(files[, c("unit", "config", COMBO, "n_ag"), with = FALSE]), by = c("unit", "config"))
us <- merge(merge(us, den, by = "unit", all.x = TRUE), fp, by = c("unit", "config"), all.x = TRUE)
for (cn in c("n_cov", "n_cov_hom", "n_cov_het", "n_fp", "n_fp_hom", "n_fp_het", "n_fp_gt01", "n_fp_gtdot1", "n_fp_gt1"))
  set(us, which(is.na(us[[cn]])), cn, 0L)
rt <- function(a, b) fifelse(b > 0, a / b, NA_real_)
us[, `:=`(fp_rate = rt(n_fp, n_cov), fp_rate_hom = rt(n_fp_hom, n_cov_hom), fp_rate_het = rt(n_fp_het, n_cov_het),
          germ_frac = rt(n_fp, n_ag), germ_frac_hom = rt(n_fp_hom, n_ag), germ_frac_het = rt(n_fp_het, n_ag),
          germ_frac_gt01 = rt(n_fp_gt01, n_ag), germ_frac_gtdot1 = rt(n_fp_gtdot1, n_ag), germ_frac_gt1 = rt(n_fp_gt1, n_ag))]
fwrite(us, file.path(OUT, "germline_unit_scores.tsv"), sep = "\t")
mn <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sm <- us[, .(tier = TIER, n_units = .N, fp_rate = mn(fp_rate), fp_rate_hom = mn(fp_rate_hom), fp_rate_het = mn(fp_rate_het),
             germ_frac = mn(germ_frac), germ_frac_hom = mn(germ_frac_hom), germ_frac_het = mn(germ_frac_het),
             germ_frac_gt01 = mn(germ_frac_gt01), germ_frac_gtdot1 = mn(germ_frac_gtdot1), germ_frac_gt1 = mn(germ_frac_gt1),
             n_fp = sum(n_fp), n_cov = sum(n_cov), n_ag = sum(n_ag)), by = c("config", COMBO)]
setorder(sm, fp_rate)
fwrite(sm, file.path(OUT, "germline_scores.tsv"), sep = "\t")
print(head(sm[, .(config, fp_rate = round(fp_rate, 4), fp_rate_het = round(fp_rate_het, 4), fp_rate_hom = round(fp_rate_hom, 4))], 10))
print(tail(sm[, .(config, fp_rate = round(fp_rate, 4))], 5))
