#!/bin/bash
# =============================================================================
# run_gt_editlevel_by_bamtype_v6.sh  [submit|scan|beds|pile|collect|correlate]
#
# v6 coordinate vs UMI de-duplication: A-to-G editing level and KO/NTC
# depletion ratio at a caller-independent site universe, quantified on both
# dedup BAMs (v5 equivalent: run_gt_editlevel_by_bamtype.sh, which also had
# Unprocessed; Unprocessed is not reported from v6 on).
#
# Each tier is quantified on its own BAMs:
#   arm       merged_<arm>.bam                  (1 group: merged)
#   sample    T{1-4}_<arm>.bam                  (groups: teratoma)
#   celltype  <arm>__<celltype>.bam             (groups: cell type)
# under 02_gt_bams{,_celltype}_{rmdup,umi}_v6/, for NTC_ctrl, ADAR1_KO, ADARB1_KO.
#
# Site universe (one list, all tiers; no caller involved): the v6 truth
# candidate floor (merged NTC cov >= 5, edited >= 1, editing >= 1%, BQ/MQ 30)
# applied to the UMI merged NTC (gt_eval_v6/candidates/candidates_all.bed) and,
# by the `scan` stage here, to the coordinate merged NTC; union of the two,
# minus dbSNP156 common. Taking the union keeps the selection symmetric
# between the two BAMs being compared.
#
# Pileups: gt_site_pileup.sh, forward library, BQ/MQ 30 (as v5 and the v6 truth);
# only rows with cov >= 10 are kept (MIN_COV of correlate).
#
# correlate (per tier):
#   edit level   per (site, group) both BAMs cov >= 10; pooled = mean over groups
#   depletion    per BAM: per (site, group) NTC and KO cov >= 10; ratio =
#                mean KO edit / mean NTC edit over those groups (NTC edit > 0)
#
# Outputs: human/pm/gt_eval_v6/dedup_compare/editlevel/   Logs: logs/
# Do NOT run samtools on the login node; `submit` only calls qsub.
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 2
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00

set -eo pipefail
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:/usr/local/bin:/usr/bin:/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"

B="${SCRNA_DIR:-.}"
PM="$B/human/pm"
SELF="$PM/scripts/run_gt_editlevel_by_bamtype_v6.sh"
V6="$PM/gt_eval_v6"
OUT="$V6/dedup_compare/editlevel"
FASTA="$B/human/ref_annots/hg38.analysisSet.fa"
PERL="$B/misc_scripts/parse_pileup_stranded.pl"
REGIONS="$V6/scan/regions.tsv"
UMI_CBED="$V6/candidates/candidates_all.bed"
SNP_BED="$V6/snp/dbSNP156_common.bed"
UNIVERSE="$OUT/beds/universe_sites.bed"
UNITS="$OUT/_units.tsv"
MIN_COV=10

ARMS=(NTC_ctrl ADAR1_KO ADARB1_KO)
BAMTYPES=(coord_dedup umi_dedup)

mkdir -p "$OUT"/{beds,pileups,frags,scan_coord}
cd "$B"
STAGE="${1:-submit}"

bam_dir() {  # tier bam_type
  local sfx; [ "$2" = "coord_dedup" ] && sfx=rmdup || sfx=umi
  if [ "$1" = "celltype" ]; then echo "$PM/02_gt_bams_celltype_${sfx}_v6"; else echo "$PM/02_gt_bams_${sfx}_v6"; fi
}

# ------------------------------------------------------------------ submit ---
if [ "$STAGE" = "submit" ]; then
  for f in "$REGIONS" "$UMI_CBED" "$SNP_BED"; do [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
  # unit table: tier group arm bam_type bam
  : > "$UNITS"
  for bt in "${BAMTYPES[@]}"; do
    for arm in "${ARMS[@]}"; do
      echo -e "arm\tmerged\t$arm\t$bt\t$(bam_dir arm "$bt")/merged_${arm}.bam" >> "$UNITS"
      for tt in T1 T2 T3 T4; do
        echo -e "sample\t$tt\t$arm\t$bt\t$(bam_dir sample "$bt")/${tt}_${arm}.bam" >> "$UNITS"
      done
      for f in "$(bam_dir celltype "$bt")/${arm}__"*.bam; do
        ct=$(basename "$f" .bam); ct="${ct#${arm}__}"
        echo -e "celltype\t$ct\t$arm\t$bt\t$f" >> "$UNITS"
      done
    done
  done
  miss=0
  while IFS=$'\t' read -r _ _ _ _ bam; do [ -s "$bam" ] && [ -s "$bam.bai" ] || { echo "ERROR: missing $bam(.bai)" >&2; miss=1; }; done < "$UNITS"
  [ "$miss" -eq 0 ] || exit 1
  n=$(wc -l < "$UNITS")
  echo "units: $n"; cut -f1,4 "$UNITS" | sort | uniq -c
  nreg=$(wc -l < "$REGIONS")
  j0=$(qsub -terse -N el6_scan -q slice -pe threaded 1 -l h_vmem=8G -t 1-"$nreg" "$SELF" scan | cut -d. -f1)
  j1=$(qsub -terse -N el6_beds -q slice -pe threaded 1 -l h_vmem=32G -hold_jid "$j0" "$SELF" beds | cut -d. -f1)
  j2=$(qsub -terse -N el6_pile -q slice -hold_jid "$j1" -t 1-"$n" -tc 60 "$SELF" pile | cut -d. -f1)
  j3=$(qsub -terse -N el6_coll -q slice -pe threaded 1 -l h_vmem=48G -hold_jid "$j2" "$SELF" collect | cut -d. -f1)
  j4=$(qsub -terse -N el6_corr -q slice -pe threaded 1 -l h_vmem=64G -hold_jid "$j3" "$SELF" correlate | cut -d. -f1)
  echo "scan $j0 (1-$nreg) | beds $j1 | pile $j2 (1-$n) | collect $j3 | correlate $j4"
  exit 0
fi

# -------------------------------------------------------------------- scan ---
# The v6 candidate floor on the coordinate merged NTC (same chunks, flags and
# floor as run_gt_truth_rebuild_v6.sh scan; candidate rows only, no histogram).
if [ "$STAGE" = "scan" ]; then
  read -r CHUNK CHR S E REG < <(sed -n "${SGE_TASK_ID:?}p" "$REGIONS")
  BAM="$PM/02_gt_bams_rmdup_v6/merged_NTC_ctrl.bam"
  CAND="$OUT/scan_coord/${CHUNK}.cand.bed"
  rm -f "$OUT/scan_coord/_ok_${CHUNK}"
  echo "=== $(date) === coord scan ${CHUNK} ${REG} $(samtools --version | head -1)"
  samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$REG" -f "$FASTA" "$BAM" \
    | perl "$PERL" forward \
    | awk -F'\t' -v OFS='\t' 'NR > 1 && $4 >= 5 && $5 >= 1 && $16 + 0 >= 0.01 {
        s = (toupper($3) == "A") ? "+" : "-"; print $1, $2 - 1, $2, $1 ":" $2, ".", s }' > "$CAND"
  echo "candidates: $(wc -l < "$CAND")"
  touch "$OUT/scan_coord/_ok_${CHUNK}"
  exit 0
fi

# -------------------------------------------------------------------- beds ---
if [ "$STAGE" = "beds" ]; then
  nreg=$(wc -l < "$REGIONS"); nok=$(ls "$OUT"/scan_coord/_ok_* 2>/dev/null | wc -l)
  [ "$nok" -eq "$nreg" ] || { echo "ERROR: $nok/$nreg coord scan chunks OK" >&2; exit 1; }
  COORD_CBED="$OUT/beds/candidates_coord.bed"
  cut -f1 "$REGIONS" | while read -r c; do cat "$OUT/scan_coord/${c}.cand.bed"; done \
    | sort -k1,1 -k2,2n -S 4G > "$COORD_CBED"
  cat "$UMI_CBED" "$COORD_CBED" | sort -k1,1 -k2,2n -u -S 4G \
    | bedtools intersect -v -a - -b "$SNP_BED" > "$UNIVERSE"
  # provenance: which scan each universe site passed
  awk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { u[$4] = 1; next } FILENAME == ARGV[2] { c[$4] = 1; next }
    FNR == 1 { print "site", "in_umi_scan", "in_coord_scan" } { print $4, ($4 in u) ? 1 : 0, ($4 in c) ? 1 : 0 }' \
    "$UMI_CBED" "$COORD_CBED" "$UNIVERSE" | gzip -c > "$OUT/beds/universe_sites_source.tsv.gz"
  {
    printf "metric\tn\n"
    printf "umi_candidates\t%s\n"   "$(wc -l < "$UMI_CBED")"
    printf "coord_candidates\t%s\n" "$(wc -l < "$COORD_CBED")"
    printf "universe_noSNP\t%s\n"   "$(wc -l < "$UNIVERSE")"
    zcat "$OUT/beds/universe_sites_source.tsv.gz" | awk -F'\t' 'NR > 1 { n[$2 $3]++ }
      END { printf "universe_both\t%d\nuniverse_umi_only\t%d\nuniverse_coord_only\t%d\n", n["11"], n["10"], n["01"] }'
  } | tee "$OUT/beds/universe_stats.tsv"
  exit 0
fi

# -------------------------------------------------------------------- pile ---
if [ "$STAGE" = "pile" ]; then
  IFS=$'\t' read -r TIER GRP ARM BT BAM < <(sed -n "${SGE_TASK_ID:?}p" "$UNITS")
  BED="$UNIVERSE"
  tag="${TIER}__${GRP}__${ARM}__${BT}"
  OUTF="$OUT/pileups/pileup_${tag}.tsv"
  FRAG="$OUT/frags/frag_${tag}.tsv.gz"
  [ -s "$BED" ] || { echo "ERROR: missing $BED" >&2; exit 1; }
  echo "=== $(date) === task $SGE_TASK_ID: $tag sites=$(wc -l < "$BED") bam=$BAM"
  bash "$B/human/scripts/gt_site_pileup.sh" \
    --bam "$BAM" --bed "$BED" --fasta "$FASTA" --out "$OUTF" \
    --libtype forward --min-bq 30 --min-mq 30
  MIN_COV="$MIN_COV" TIER="$TIER" GRP="$GRP" ARM="$ARM" BT="$BT" OUTF="$OUTF" FRAG="$FRAG" Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))
p <- fread(Sys.getenv("OUTF"), sep = "\t", showProgress = FALSE)
out <- p[coverage >= as.integer(Sys.getenv("MIN_COV")), .(tier = Sys.getenv("TIER"), group = Sys.getenv("GRP"),
                         arm = Sys.getenv("ARM"), bam_type = Sys.getenv("BT"),
                         site = paste0(chr, ":", pos), cov = as.integer(coverage),
                         alt = as.integer(variant_cov), edit = as.numeric(editing_level))]
fwrite(out, sub("\\.gz$", "", Sys.getenv("FRAG")), sep = "\t")
cat("frag rows:", nrow(out), "\n")
RS
  # fwrite of a zero-row table straight to .gz leaves a truncated gzip stream
  gzip -f "${FRAG%.gz}"
  rm -f "$OUTF"
  echo "=== $(date) === done"
  exit 0
fi

# ----------------------------------------------------------------- collect ---
if [ "$STAGE" = "collect" ]; then
  n_exp=$(wc -l < "$UNITS")
  n_got=$(ls "$OUT"/frags/frag_*.tsv.gz 2>/dev/null | wc -l)
  [ "$n_got" -eq "$n_exp" ] || { echo "ERROR: $n_got/$n_exp frags present" >&2; exit 1; }
  OUT="$OUT" Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))
out <- Sys.getenv("OUT")
files <- list.files(file.path(out, "frags"), pattern = "^frag_.*\\.tsv\\.gz$", full.names = TRUE)
fr <- lapply(files, function(f) fread(cmd = paste("gzip -cd", shQuote(f)),
                                      colClasses = list(character = c("group"))))
empty <- vapply(fr, nrow, 1L) == 0L
if (any(empty)) cat("no sites at the coverage floor:", paste(basename(files[empty]), collapse = ", "), "\n")
d <- rbindlist(fr[!empty])
fwrite(d, file.path(out, "editlevel_long.tsv.gz"), sep = "\t")
cat("editlevel_long.tsv.gz rows=", nrow(d), " sites=", uniqueN(d$site), "\n")
print(d[, .(n_units = uniqueN(group), n_rows = .N), by = .(tier, arm, bam_type)][order(tier, arm, bam_type)])
RS
  exit 0
fi

# --------------------------------------------------------------- correlate ---
if [ "$STAGE" = "correlate" ]; then
  MIN_COV="$MIN_COV" OUT="$OUT" Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))
out <- Sys.getenv("OUT")
MIN_COV <- as.integer(Sys.getenv("MIN_COV"))
d <- fread(cmd = paste("gzip -cd", shQuote(file.path(out, "editlevel_long.tsv.gz"))),
           colClasses = list(character = "group"), showProgress = FALSE)
w <- dcast(d, tier + group + arm + site ~ bam_type, value.var = c("cov", "edit"))
w <- w[cov_coord_dedup >= MIN_COV & cov_umi_dedup >= MIN_COV]

cor_row <- function(a, b) data.table(n_sites = length(a),
  pearson = cor(a, b), spearman = cor(a, b, method = "spearman"),
  mae = mean(abs(a - b)), median_abs_delta = median(abs(a - b)))

by_grp <- w[, if (.N >= 20L) cor_row(edit_coord_dedup, edit_umi_dedup), by = .(tier, arm, group)]
fwrite(by_grp, file.path(out, "editlevel_correlation_by_group.tsv"), sep = "\t")

pts <- w[, .(n_groups = .N, edit_coord = mean(edit_coord_dedup), edit_umi = mean(edit_umi_dedup)),
         by = .(tier, arm, site)]
fwrite(pts, file.path(out, "editlevel_pooled_points.tsv.gz"), sep = "\t")
pool <- pts[, cor_row(edit_coord, edit_umi), by = .(tier, arm)]
fwrite(pool, file.path(out, "editlevel_correlation_pooled.tsv"), sep = "\t")
cat("\n=== edit level, coordinate vs UMI (pooled over groups) ===\n"); print(pool)

# depletion ratio per BAM type, NTC and KO matched within group
dep <- rbindlist(lapply(c("ADAR1_KO", "ADARB1_KO"), function(ko) {
  n <- d[arm == "NTC_ctrl" & cov >= MIN_COV, .(tier, group, site, bam_type, edit_ntc = edit)]
  k <- d[arm == ko & cov >= MIN_COV, .(tier, group, site, bam_type, edit_ko = edit)]
  m <- merge(n, k, by = c("tier", "group", "site", "bam_type"))
  m <- m[, .(n_groups = .N, edit_ntc = mean(edit_ntc), edit_ko = mean(edit_ko)),
         by = .(tier, site, bam_type)][edit_ntc > 0]
  m[, `:=`(ko_arm = ko, ratio = edit_ko / edit_ntc)]
}))
dw <- dcast(dep, tier + ko_arm + site ~ bam_type, value.var = c("ratio", "edit_ntc", "edit_ko", "n_groups"))
dw <- dw[!is.na(ratio_coord_dedup) & !is.na(ratio_umi_dedup)]
fwrite(dw, file.path(out, "depletion_ratio_points.tsv.gz"), sep = "\t")
conc <- dw[, .(n_sites = .N,
               agree_50 = mean((ratio_coord_dedup <= 0.5) == (ratio_umi_dedup <= 0.5)),
               both_dep_50 = mean(ratio_coord_dedup <= 0.5 & ratio_umi_dedup <= 0.5),
               only_coord_50 = mean(ratio_coord_dedup <= 0.5 & ratio_umi_dedup > 0.5),
               only_umi_50 = mean(ratio_coord_dedup > 0.5 & ratio_umi_dedup <= 0.5),
               spearman_ratio = cor(ratio_coord_dedup, ratio_umi_dedup, method = "spearman")),
           by = .(tier, ko_arm)]
fwrite(conc, file.path(out, "depletion_concordance.tsv"), sep = "\t")
cat("\n=== KO/NTC depletion ratio, coordinate vs UMI ===\n"); print(conc)
RS
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
