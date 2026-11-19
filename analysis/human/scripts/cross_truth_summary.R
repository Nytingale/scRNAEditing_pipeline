#!/usr/bin/env Rscript
# cross_truth_summary.R
# Cross-truth scoring summary (human/pm/scripts/cross_truth_score.sh): every call set against every truth.
# Truths: v9 (unseeded umi_tools BAMs), seed1, seed2 (umi_tools --random-seed 1 / 2), umimd (UMI-aware MarkDuplicates),
#         agree (sites labelled the same by all four; cross_truth_agreement.R).
# umi_tools call sources for the 4 BCFtools umi_tools configurations: pub (unseeded, the published run), seed1, seed2.
# A truth and a call source are "matched" when built from the same BAMs (v9-pub, seed1-seed1, seed2-seed2).
# F1 per configuration = mean over the 6 truths (gate x t) of the cross-tier mean (as umimd_compare.R, harmonic76.R).
# Configurations and the sanity reference: human/pm/gt_eval_v9/configs76/configs76_seeded.tsv (gt_v10_configs76.R).
# Outputs (human/pm/gt_eval_v9/cross_truth/):
#   cross_truth_f1.tsv        ds, truth, calls, config, f1, precision, recall (76 configurations + seed BCFtools sets)
#   cross_truth_bcftools.tsv  BCFtools umi_tools / umimd / coordinate F1 per truth x umi_tools source, matched flag
#   cross_truth_bias.tsv      shared-read bias: matched minus crossed umi_tools F1 (same truth), per BCFtools config
#   cross_truth_top.tsv       top 10 of the 76 per truth (umi_tools BCFtools from each source)
#   cross_truth_spearman.tsv  rank correlation of the 76 between truths
#   cross_truth_summary.log   printed summary
suppressPackageStartupMessages(library(data.table))
B <- Sys.getenv("SCRNA_DIR", ".")
setwd(B)
O <- "human/pm/gt_eval_v9/cross_truth"
TRUTHS <- c("v9", "seed1", "seed2", "umimd", "agree", "coord")
CALLS <- c("pub", "ch", "umimd", "seed1", "seed2")
TIERS <- list(teratoma = c("arm", "sample", "celltype"), aorta = c("arm", "sample"))
MATCH <- c(v9 = "pub", seed1 = "seed1", seed2 = "seed2")
cols <- c("gate", "t", "caller", "bam_type", "strand_mode", "snp_subtract", "precision", "recall", "f1")
C76 <- fread("human/pm/gt_eval_v9/configs76/configs76_seeded.tsv")   # human/scripts/gt_v10_configs76.R
c60 <- grep("__sampleUmiMd__", C76$config, value = TRUE, invert = TRUE)
stopifnot(nrow(C76) == 76, length(c60) == 60)
ch_cfg <- grep("^cellsnp__.*__hybrid__", c60, value = TRUE)

s <- rbindlist(lapply(names(TIERS), function(d) rbindlist(lapply(TRUTHS, function(tr) rbindlist(lapply(CALLS, function(cl)
  rbindlist(lapply(TIERS[[d]], function(ti) {
    f <- file.path(O, "scores", d, tr, cl, ti, "gt_v7_scores.tsv")
    if (!file.exists(f)) stop("missing ", f)
    unique(fread(f)[, ..cols])[, `:=`(ds = d, truth = tr, calls = cl, tier = ti)]
  }))))))))
s[calls == "umimd", bam_type := "sampleUmiMd"]
s[, config := paste(caller, bam_type, strand_mode, snp_subtract, sep = "__")]
s <- s[(calls == "pub" & config %in% setdiff(c60, ch_cfg)) | (calls == "ch" & config %in% ch_cfg) | calls == "umimd" |
       (calls %like% "^seed" & caller == "bcftools" & bam_type == "sampleUmiDup")]
stopifnot(s[, .N, by = .(ds, truth, config, calls, gate, t)][, all(N == length(TIERS[[ds[1]]])), by = ds]$V1)
bt <- s[, .(f1 = mean(f1), precision = mean(precision, na.rm = TRUE), recall = mean(recall)),
        by = .(ds, truth, calls, config, caller, bam_type, strand_mode, snp_subtract, gate, t)]
f1 <- bt[, .(f1 = mean(f1), precision = mean(precision), recall = mean(recall)),
         by = .(ds, truth, calls, config, caller, bam_type, strand_mode, snp_subtract)]
stopifnot(f1[calls %in% c("pub", "ch", "umimd"), .N, by = .(ds, truth)][, all(N == 76)],
          f1[calls %like% "^seed", .N, by = .(ds, truth, calls)][, all(N == 4)])
fwrite(f1, file.path(O, "cross_truth_f1.tsv"), sep = "\t")

sink(file.path(O, "cross_truth_summary.log"), split = TRUE)
options(width = 220)

# sanity: v9 truth x published / umimd call sets reproduces the existing headline scores
ref <- melt(C76[, .(config, teratoma = f1_te, aorta = f1_ao)],
            id.vars = "config", variable.name = "ds", value.name = "f1_ref", variable.factor = FALSE)
chk <- merge(f1[truth == "v9" & calls %in% c("pub", "ch", "umimd")], ref, by = c("ds", "config"))
cat("Sanity (v9 truth, published + umimd call sets vs configs76_seeded.tsv):", nrow(chk), "rows; max |dF1| =",
    signif(chk[, max(abs(f1 - f1_ref))], 3), "\n")

# BCFtools: umi_tools (per source), umimd, coordinate under each truth
bc <- f1[caller == "bcftools" & bam_type %in% c("sampleUmiDup", "sampleRmDup", "sampleUmiMd")]
bc[, dedup := fcase(bam_type == "sampleRmDup", "coordinate", bam_type == "sampleUmiMd", "umimd",
                    calls == "pub", "umi_tools_unseeded", default = paste0("umi_tools_", calls))]
bc[, cfg := paste(strand_mode, snp_subtract, sep = "/")]
bc[, matched := (dedup %like% "^umi_tools" & MATCH[truth] == calls) | (truth == "umimd" & dedup == "umimd") |
     (truth == "coord" & dedup == "coordinate")]
bc[is.na(matched), matched := FALSE]
fwrite(bc, file.path(O, "cross_truth_bcftools.tsv"), sep = "\t")
for (d in names(TIERS)) {
  cat("\n=====", d, ": BCFtools F1 by deduplication (rows) and truth (columns); * = call BAMs built the truth\n")
  w <- dcast(bc[ds == d], cfg + dedup ~ truth, value.var = "f1")
  m <- dcast(bc[ds == d], cfg + dedup ~ truth, value.var = "matched")
  for (tr in TRUTHS) set(w, j = tr, value = sprintf("%.4f%s", w[[tr]], fifelse(m[[tr]], "*", " ")))
  print(w[, c("cfg", "dedup", TRUTHS), with = FALSE])
}

# shared-read bias: same truth, umi_tools calls from the truth's own BAMs minus calls from the other seeds' BAMs
ut <- bc[dedup %like% "^umi_tools" & truth %in% names(MATCH)]
bias <- ut[, .(f1_matched = f1[calls == MATCH[truth[1]]], f1_crossed = mean(f1[calls != MATCH[truth[1]]])),
           by = .(ds, truth, cfg)][, bias := f1_matched - f1_crossed]
fwrite(bias, file.path(O, "cross_truth_bias.tsv"), sep = "\t")
cat("\nShared-read bias (umi_tools F1, matched minus crossed calls, same truth):\n")
print(bias[, .(mean_bias = mean(bias), min = min(bias), max = max(bias), n = .N), by = ds])
gap <- dcast(bc[dedup %in% c("umimd", "umi_tools_unseeded", "coordinate")], ds + truth + cfg ~ dedup, value.var = "f1")
gap[, `:=`(umimd_minus_umi_tools = umimd - umi_tools_unseeded, umimd_minus_coordinate = umimd - coordinate)]
cat("\nBCFtools umimd minus published umi_tools, by truth (mean over 4 configurations):\n")
print(dcast(gap[, .(d = round(mean(umimd_minus_umi_tools), 4)), by = .(ds, truth)], ds ~ truth, value.var = "d"))
cat("\nBCFtools coordinate minus umimd, by truth (mean over 4 configurations):\n")
print(dcast(gap[, .(d = round(mean(-umimd_minus_coordinate), 4)), by = .(ds, truth)], ds ~ truth, value.var = "d"))
cg <- gap[, .(coord_truth = -umimd_minus_coordinate[truth == "coord"],
              other_truths = mean(-umimd_minus_coordinate[truth %in% c("v9", "seed1", "seed2", "umimd")])), by = .(ds, cfg)]
cat("\nCoordinate shared-read bias: (coordinate - umimd) under the coordinate truth minus under the four other truths:\n")
print(cg[, .(mean_bias = round(mean(coord_truth - other_truths), 5), min = round(min(coord_truth - other_truths), 5),
             max = round(max(coord_truth - other_truths), 5)), by = ds])

# rankings of the 76 per truth, umi_tools BCFtools from each source
r76 <- rbindlist(lapply(c("pub", "seed1", "seed2"), function(src) {
  x <- f1[calls %in% c("pub", "ch", "umimd")]
  if (src != "pub") {
    x <- x[!(caller == "bcftools" & bam_type == "sampleUmiDup")]
    x <- rbind(x, f1[calls == src])
  }
  x[, umi_src := src]
}))
r76[, rk := frank(-f1, ties.method = "min"), by = .(ds, truth, umi_src)]
stopifnot(r76[, .N, by = .(ds, truth, umi_src)][, all(N == 76)])
top <- r76[rk <= 10][order(ds, truth, umi_src, rk)]
fwrite(top, file.path(O, "cross_truth_top.tsv"), sep = "\t")
cat("\n#1 of the 76 per truth and umi_tools source:\n")
print(dcast(r76[rk == 1, .(ds, truth, umi_src, top = paste0(config, " ", sprintf("%.4f", f1)))], ds + truth ~ umi_src,
            value.var = "top", fun.aggregate = function(v) paste(v, collapse = "; ")))
cat("\nRank of each BCFtools annotation/SNP-on deduplication among the 76:\n")
print(dcast(r76[caller == "bcftools" & strand_mode == "annotation" & snp_subtract == "on" &
                bam_type %in% c("sampleUmiDup", "sampleRmDup", "sampleUmiMd")],
            ds + umi_src + bam_type ~ truth, value.var = "rk"))

# rank agreement of the 76 between truths (published umi_tools source)
sp <- rbindlist(lapply(names(TIERS), function(d) {
  w <- dcast(r76[ds == d & umi_src == "pub"], config ~ truth, value.var = "f1")
  rbindlist(combn(TRUTHS, 2, function(p) data.table(ds = d, a = p[1], b = p[2],
    spearman = cor(w[[p[1]]], w[[p[2]]], method = "spearman")), simplify = FALSE))
}))
fwrite(sp, file.path(O, "cross_truth_spearman.tsv"), sep = "\t")
cat("\nSpearman of the 76 F1 between truths:\n")
print(sp[, .(ds, a, b, spearman = round(spearman, 4))])

# umimd vs published umi_tools across the 16 re-run configurations, per truth
w16 <- merge(f1[calls == "umimd", .(ds, truth, caller, strand_mode, snp_subtract, f1_umimd = f1)],
             f1[calls %in% c("pub", "ch") & bam_type == "sampleUmiDup", .(ds, truth, caller, strand_mode, snp_subtract, f1_ut = f1)],
             by = c("ds", "truth", "caller", "strand_mode", "snp_subtract"))
cat("\numimd F1 > published umi_tools F1 (of 16 configurations), per truth:\n")
print(dcast(w16[, .(wins = sum(f1_umimd > f1_ut), n = .N), by = .(ds, truth)][, v := paste0(wins, "/", n)],
            ds ~ truth, value.var = "v"))
sink()
