#!/usr/bin/env Rscript
# ============================================================================
# TERATOMA_v8_caller_resources_source_data.R   (2026-10-15; split from the _plot.R script 2026-10-25)
# Source data for response figure 12:
# Compute cost of each caller in the v8 teratoma runs (review R1 minor 9): per calling task (one BAM type of
# one unit), wall time, CPU hours (realtime x %CPU) and peak resident memory, from the Nextflow traces, and
# storage = bytes of the files each task wrote to its work directory (human/pm/scripts/caller_storage_v8.sh;
# scAllele = per-chromosome outputs + merged output, as both are on disk at the end of a run).
# scAllele = SCALLELE_CHROM (per-chromosome scatter) summed per unit x BAM type + SCALLELE_GATHER; GATK =
# HaplotypeCaller only (SplitNCigarReads / BQSR are shared preprocessing). Failed attempts are counted.
# Inputs : human/pm/results_gt_{v8,persample_v8,celltype_v8}/pipeline_trace.txt,
#          human/pm/gt_eval_v8/caller_storage/storage_{arm,sample,celltype}.tsv
# Output : source_data/pm_gt/TERATOMA_v8_caller_resources.{tasks,summary}.tsv
#          (drawn by TERATOMA_v8_caller_resources_plot.R)
# ============================================================================
source(file.path(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "TERATOMA_common.R"))
teratoma_set_dataset("pm_gt")
source(file.path(REPO_DIR, "plotting", "source_data_io.R"))
STEM <- "TERATOMA_v8_caller_resources"

PM <- file.path(REPO_DIR, "human", "pm")
SRC <- c(arm = "results_gt_v8", sample = "results_gt_persample_v8", celltype = "results_gt_celltype_v8")
SRC <- setNames(file.path(PM, SRC, "pipeline_trace.txt"), names(SRC))
STORE <- setNames(file.path(PM, "gt_eval_v8", "caller_storage", sprintf("storage_%s.tsv", names(SRC))), names(SRC))
stopifnot(all(file.exists(SRC)), all(file.exists(STORE)))
PROC <- c(BCFTOOLS_MPILEUP = "bcftools", CELLSNP_LITE = "cellsnp", GATK_HAPLOTYPECALLER = "gatk",
          JACUSA2_CALL1 = "jacusa2", REDITOOLS3 = "reditools3", SCALLELE_CHROM = "scallele", SCALLELE_GATHER = "scallele")

dur_h <- function(x) vapply(x, function(s) {
  if (is.na(s) || s %in% c("-", "")) return(NA_real_)
  u <- c(d = 24, h = 1, m = 1 / 60, s = 1 / 3600, ms = 1 / 3.6e6)
  p <- regmatches(s, gregexpr("[0-9.]+(ms|d|h|m|s)", s))[[1]]
  sum(as.numeric(sub("[a-z]+$", "", p)) * u[sub("^[0-9.]+", "", p)])
}, numeric(1), USE.NAMES = FALSE)
mem_gb <- function(x) {
  v <- as.numeric(sub(" .*", "", x)); u <- sub("^[0-9.]+ ", "", x)
  v * c(B = 1 / 1024^3, KB = 1 / 1024^2, MB = 1 / 1024, GB = 1, TB = 1024)[u]
}

tr <- rbindlist(lapply(names(SRC), function(t) fread(SRC[[t]], sep = "\t", quote = "")[process %in% names(PROC)][, tier := t]))
tr[, `:=`(caller = PROC[process], wall_h = dur_h(realtime), rss_gb = mem_gb(peak_rss),
          cpu_frac = as.numeric(sub("%", "", `%cpu`)) / 100)]
tr[, cpu_h := wall_h * cpu_frac]
tr[, unit_bam := fifelse(process == "SCALLELE_CHROM", sub(":chr[^:]*$", "", tag), tag)]
fails <- tr[status != "COMPLETED", .N, by = .(tier, caller, status)]
cat("failed / non-completed attempts:\n"); print(fails)
ok <- tr[status == "COMPLETED"]
st <- rbindlist(lapply(STORE, fread, na.strings = "NA"))
stopifnot(!anyNA(st$bytes))
ok <- merge(ok, st[, .(tier, hash, bytes)], by = c("tier", "hash"), all.x = TRUE)
stopifnot(!anyNA(ok$bytes))
task <- ok[, .(wall_h = if (caller[1] == "scallele") max(wall_h[process == "SCALLELE_CHROM"]) + sum(wall_h[process == "SCALLELE_GATHER"])
                        else sum(wall_h),
               cpu_h = sum(cpu_h, na.rm = TRUE), rss_gb = max(rss_gb, na.rm = TRUE),
               storage_mb = sum(bytes) / 1024^2), by = .(tier, caller, unit_bam)]
summ <- task[, .(n_tasks = .N, median_wall_h = round(median(wall_h), 3), median_cpu_h = round(median(cpu_h), 3),
                 total_cpu_h = round(sum(cpu_h), 1), median_rss_gb = round(median(rss_gb), 2), max_rss_gb = round(max(rss_gb), 2),
                 median_storage_mb = round(median(storage_mb), 1), total_storage_gb = round(sum(storage_mb) / 1024, 2)),
             by = .(tier, caller)][order(tier, -total_cpu_h)]
summ <- merge(summ, fails[, .(failed_attempts = sum(N)), by = .(tier, caller)], by = c("tier", "caller"), all.x = TRUE)
summ[is.na(failed_attempts), failed_attempts := 0L]
print(summ)

units <- unique(task[, .(tier, unit = sub(":[^:]*$", "", unit_bam))])
units[, `:=`(arm = sub("^T[1-4]_", "", sub("__.*$", "", unit)), ct = sub("^.*__", "", unit))]
n_arm <- units[tier == "arm", uniqueN(arm)]
n_ter <- units[tier == "sample", uniqueN(sub("_.*$", "", unit))]
ct_per_arm <- units[tier == "celltype", .N, by = arm]$N
stopifnot(units[tier == "sample", uniqueN(arm)] == n_arm, units[tier == "celltype", uniqueN(arm)] == n_arm)
TIER_LAB <- c(arm = sprintf("Arm\n(%d arms)", n_arm),
              sample = sprintf("Sample\n(%d arms \u00d7\n%d teratomas)", n_arm, n_ter),
              celltype = sprintf("Cell type\n(%d arms \u00d7\n%d\u2013%d cell types)", n_arm, min(ct_per_arm), max(ct_per_arm)))
task[, wall_min := 60 * wall_h]
METRIC_LAB <- c(wall_min = "Run time (min)", rss_gb = "Peak memory (GB)", storage_mb = "Storage (MB)")
long <- melt(task, id.vars = c("tier", "caller", "unit_bam"), measure.vars = names(METRIC_LAB))
long[, `:=`(Tier = factor(TIER_LAB[tier], TIER_LAB), method = caller_factor(caller, rev = TRUE),
            Metric = factor(METRIC_LAB[variable], METRIC_LAB))]

write_sd(long, SRC_DIR, STEM, "tasks")
write_sd(summ, SRC_DIR, STEM, "summary")
