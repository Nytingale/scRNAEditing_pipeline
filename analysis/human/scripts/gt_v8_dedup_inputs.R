#!/usr/bin/env Rscript
# =============================================================================
# gt_v8_dedup_inputs.R   (2026-10-13)
#
# Inputs for the v8 BAM-processing figures (TERATOMA_v8_dedup_compare_plot.R,
# GTDEF_v8_bam_processing_plot.R).
#   discovery  every site called in NTC per configuration, all 12 varTypes, v8 call state
#              (mismatch = nearby-mismatch pass only: results_gt{,_persample,_celltype}_v8/
#              <NTC unit>/04_filtered/**/*_step5_nomm.tsv.gz). Same columns as
#              gt_eval_v6/dedup_compare/discovery (gt_eval_sites.R --all_vartypes true), which
#              used the final pre100 call sets (both window filters).
#              varType is in the RNA frame of the call's strand for every caller (REDItools3
#              after the RNA-frame parse fix; gt_eval_v7/variants/reditools3_frame_check.tsv).
#   read counts / site depth  truth-independent read counts and samtools depth at the v6
#              truth positions (gt_eval_v6/dedup_compare), restricted here to the v8
#              expected-gated t = 50% positives (all 1,462 are among the v6 positions).
# Outputs: gt_eval_v8/dedup_compare/
#          discovery/{arm,sample,celltype}/gt_eval_discovery.tsv.gz
#          dedup_read_counts.tsv, dedup_truth_site_depth.tsv.gz
# =============================================================================
suppressPackageStartupMessages(library(data.table))
B = Sys.getenv("SCRNA_DIR", ".")
PM = file.path(B, "human/pm")
PUB = c(arm = "results_gt_v8", sample = "results_gt_persample_v8", celltype = "results_gt_celltype_v8")
DD6 = file.path(PM, "gt_eval_v6/dedup_compare")
OUT = file.path(PM, "gt_eval_v8/dedup_compare")
msg = function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")

for (tier in c("arm", "sample", "celltype")) {
  root = file.path(PM, PUB[[tier]])
  units = grep("NTC_ctrl", list.dirs(root, recursive = FALSE, full.names = FALSE), value = TRUE)
  files = unlist(lapply(units, function(u)
    list.files(file.path(root, u, "04_filtered"), pattern = "_step5_nomm\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)))
  d = rbindlist(lapply(files, function(f) {
    p = strsplit(sub(paste0("^", root, "/"), "", f), "/", fixed = TRUE)[[1]]
    # <unit>/04_filtered/<caller>/<bam_type>/<strand_mode>/snp_<on|off>/<file>
    x = fread(f, select = c("chr", "pos", "strand", "varType", "site"), showProgress = FALSE)
    if (!nrow(x)) return(NULL)
    x[, `:=`(caller = p[3], bam_type = p[4], strand_mode = p[5], snp_subtract = sub("^snp_", "", p[6]), arm = p[1])]
  }))
  d[, combo := paste(caller, bam_type, strand_mode, snp_subtract, sep = "__")]
  setcolorder(d, c("combo", "caller", "bam_type", "strand_mode", "snp_subtract", "arm", "chr", "pos", "strand", "varType", "site"))
  dir.create(file.path(OUT, "discovery", tier), recursive = TRUE, showWarnings = FALSE)
  fwrite(d, file.path(OUT, "discovery", tier, "gt_eval_discovery.tsv.gz"), sep = "\t")
  msg(tier, ": ", length(units), " NTC units, ", length(files), " files, ", nrow(d), " calls, ",
      uniqueN(d$combo), " configurations")
}

file.copy(file.path(DD6, "dedup_read_counts.tsv"), file.path(OUT, "dedup_read_counts.tsv"), overwrite = TRUE)
v8 = fread(file.path(PM, "gt_eval_v8/truth/gt_v7_truth_sites.tsv.gz"))[gate == "expected" & class_t50 == "positive"]
dep = fread(file.path(DD6, "dedup_truth_site_depth.tsv.gz"))
dep = dep[paste0(chr, ":", pos) %chin% v8$site]
stopifnot(dep[, uniqueN(paste(chr, pos))] == nrow(v8))
fwrite(dep, file.path(OUT, "dedup_truth_site_depth.tsv.gz"), sep = "\t")
msg("depth: ", nrow(v8), " v8 expected t50 positives x ", uniqueN(dep$unit), " units x ", uniqueN(dep$bam_type), " BAM types")
