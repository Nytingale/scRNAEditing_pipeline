#!/bin/bash
# =============================================================================
# run_gt_editlevel_by_bamtype_v9.sh  [submit|scan|beds|pile|collect|correlate]   (2026-11-16)
#
# run_gt_editlevel_by_bamtype_v6.sh on the seed-8787 umi_tools BAMs, with UMI-aware MarkDuplicates as a third
# de-duplication (Response Figure 15): A-to-G editing level and KO/NTC depletion ratio at a caller-independent site
# universe, each pair of the three BAM types compared.
#
# BAM types (tier BAMs as in v6):
#   coord_dedup  02_gt_bams{,_celltype}_rmdup_v6/
#   umi_dedup    02_gt_bams{,_celltype}_umi_v6/     (seed 8787)
#   umimd_dedup  02_gt_bams_umimd/, 02_gt_bams_celltype_umimd/ (cell types: NTC only, so no UMI-aware
#                MarkDuplicates depletion ratio at the cell-type tier)
# Arms NTC_ctrl, ADAR1_KO, ADARB1_KO.
#
# Site universe: the v6 truth candidate floor (merged NTC cov >= 5, edited >= 1, editing >= 1%, BQ/MQ 30) on each
# merged NTC BAM (umi_tools: gt_eval_v6/candidates/candidates_all.bed, rebuilt from the seed-8787 BAM; coordinate
# and UMI-aware MarkDuplicates: `scan` here), union, minus dbSNP156 common.
# Pileups: gt_site_pileup.sh, forward library, BQ/MQ 30; rows with cov >= 10 kept.
# correlate: per pair of BAM types and tier, as v6 (edit level: both cov >= 10, pooled = mean over groups;
#            depletion: per BAM NTC and KO cov >= 10 in the group, ratio = mean KO edit / mean NTC edit).
#
# Outputs: human/pm/gt_eval_v9/dedup_compare/editlevel/   Logs: logs/
# Do NOT run samtools on the login node; `submit` only calls qsub.
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 2
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00

set -eo pipefail
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true

B="${SCRNA_DIR:-.}"
PM="$B/human/pm"
SELF="$PM/scripts/run_gt_editlevel_by_bamtype_v9.sh"
V6="$PM/gt_eval_v6"
OUT="$PM/gt_eval_v9/dedup_compare/editlevel"
FASTA="$B/human/ref_annots/hg38.analysisSet.fa"
PERL="$B/misc_scripts/parse_pileup_stranded.pl"
REGIONS="$V6/scan/regions.tsv"
UMI_CBED="$V6/candidates/candidates_all.bed"
SNP_BED="$V6/snp/dbSNP156_common.bed"
UNIVERSE="$OUT/beds/universe_sites.bed"
UNITS="$OUT/_units.tsv"
MIN_COV=10
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript

ARMS=(NTC_ctrl ADAR1_KO ADARB1_KO)
BAMTYPES=(coord_dedup umi_dedup umimd_dedup)
SCANS=(coord_dedup umimd_dedup)

mkdir -p "$OUT"/{beds,pileups,frags,scan_coord_dedup,scan_umimd_dedup}
cd "$B"
STAGE="${1:-submit}"

bam_dir() {  # tier bam_type
  case "$2" in
    coord_dedup) [ "$1" = celltype ] && echo "$PM/02_gt_bams_celltype_rmdup_v6" || echo "$PM/02_gt_bams_rmdup_v6" ;;
    umi_dedup)   [ "$1" = celltype ] && echo "$PM/02_gt_bams_celltype_umi_v6"   || echo "$PM/02_gt_bams_umi_v6" ;;
    umimd_dedup) [ "$1" = celltype ] && echo "$PM/02_gt_bams_celltype_umimd"    || echo "$PM/02_gt_bams_umimd" ;;
  esac
}

# ------------------------------------------------------------------ submit ---
if [ "$STAGE" = "submit" ]; then
  for f in "$REGIONS" "$UMI_CBED" "$SNP_BED"; do [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done
  : > "$UNITS"
  for bt in "${BAMTYPES[@]}"; do
    for arm in "${ARMS[@]}"; do
      echo -e "arm\tmerged\t$arm\t$bt\t$(bam_dir arm "$bt")/merged_${arm}.bam" >> "$UNITS"
      for tt in T1 T2 T3 T4; do
        echo -e "sample\t$tt\t$arm\t$bt\t$(bam_dir sample "$bt")/${tt}_${arm}.bam" >> "$UNITS"
      done
      for f in "$(bam_dir celltype "$bt")/${arm}__"*.bam; do
        [ -e "$f" ] || continue
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
  j0=$(qsub -terse -N el9_scan -pe threaded 1 -l h_vmem=8G -t 1-"$((2 * nreg))" "$SELF" scan | cut -d. -f1)
  j1=$(qsub -terse -N el9_beds -pe threaded 1 -l h_vmem=32G -hold_jid "$j0" "$SELF" beds | cut -d. -f1)
  j2=$(qsub -terse -N el9_pile -hold_jid "$j1" -t 1-"$n" -tc 60 "$SELF" pile | cut -d. -f1)
  j3=$(qsub -terse -N el9_coll -pe threaded 1 -l h_vmem=48G -hold_jid "$j2" "$SELF" collect | cut -d. -f1)
  j4=$(qsub -terse -N el9_corr -pe threaded 1 -l h_vmem=64G -hold_jid "$j3" "$SELF" correlate | cut -d. -f1)
  echo "scan $j0 (1-$((2 * nreg))) | beds $j1 | pile $j2 (1-$n) | collect $j3 | correlate $j4"
  exit 0
fi

# -------------------------------------------------------------------- scan ---
# tasks 1..nreg: coordinate merged NTC; nreg+1..2*nreg: UMI-aware MarkDuplicates merged NTC
if [ "$STAGE" = "scan" ]; then
  nreg=$(wc -l < "$REGIONS")
  i=$(( (${SGE_TASK_ID:?} - 1) / nreg )); r=$(( (SGE_TASK_ID - 1) % nreg + 1 ))
  bt="${SCANS[$i]}"
  read -r CHUNK CHR S E REG < <(sed -n "${r}p" "$REGIONS")
  BAM="$(bam_dir arm "$bt")/merged_NTC_ctrl.bam"
  D="$OUT/scan_${bt}"
  CAND="$D/${CHUNK}.cand.bed"
  rm -f "$D/_ok_${CHUNK}"
  echo "=== $(date) === $bt scan ${CHUNK} ${REG} $(samtools --version | head -1)"
  samtools mpileup -B -q 30 -Q 30 -d 1000000 -r "$REG" -f "$FASTA" "$BAM" \
    | perl "$PERL" forward \
    | awk -F'\t' -v OFS='\t' 'NR > 1 && $4 >= 5 && $5 >= 1 && $16 + 0 >= 0.01 {
        s = (toupper($3) == "A") ? "+" : "-"; print $1, $2 - 1, $2, $1 ":" $2, ".", s }' > "$CAND"
  echo "candidates: $(wc -l < "$CAND")"
  touch "$D/_ok_${CHUNK}"
  exit 0
fi

# -------------------------------------------------------------------- beds ---
if [ "$STAGE" = "beds" ]; then
  nreg=$(wc -l < "$REGIONS")
  for bt in "${SCANS[@]}"; do
    nok=$(ls "$OUT/scan_${bt}"/_ok_* 2>/dev/null | wc -l)
    [ "$nok" -eq "$nreg" ] || { echo "ERROR: $bt $nok/$nreg scan chunks OK" >&2; exit 1; }
    cut -f1 "$REGIONS" | while read -r c; do cat "$OUT/scan_${bt}/${c}.cand.bed"; done \
      | sort -k1,1 -k2,2n -S 4G > "$OUT/beds/candidates_${bt}.bed"
  done
  cp "$UMI_CBED" "$OUT/beds/candidates_umi_dedup.bed"
  cat "$OUT"/beds/candidates_{umi_dedup,coord_dedup,umimd_dedup}.bed | sort -k1,1 -k2,2n -u -S 4G \
    | bedtools intersect -v -a - -b "$SNP_BED" > "$UNIVERSE"
  awk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] { u[$4] = 1; next } FILENAME == ARGV[2] { c[$4] = 1; next }
    FILENAME == ARGV[3] { m[$4] = 1; next }
    FNR == 1 { print "site", "in_umi_scan", "in_coord_scan", "in_umimd_scan" }
    { print $4, ($4 in u) ? 1 : 0, ($4 in c) ? 1 : 0, ($4 in m) ? 1 : 0 }' \
    "$OUT"/beds/candidates_{umi_dedup,coord_dedup,umimd_dedup}.bed "$UNIVERSE" \
    | gzip -c > "$OUT/beds/universe_sites_source.tsv.gz"
  {
    printf "metric\tn\n"
    for bt in umi_dedup coord_dedup umimd_dedup; do printf "%s_candidates\t%s\n" "$bt" "$(wc -l < "$OUT/beds/candidates_${bt}.bed")"; done
    printf "universe_noSNP\t%s\n" "$(wc -l < "$UNIVERSE")"
    zcat "$OUT/beds/universe_sites_source.tsv.gz" | awk -F'\t' 'NR > 1 { n[$2 $3 $4]++ }
      END { printf "universe_all3\t%d\n", n["111"] }'
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
  MIN_COV="$MIN_COV" TIER="$TIER" GRP="$GRP" ARM="$ARM" BT="$BT" OUTF="$OUTF" FRAG="$FRAG" $RS - <<'RS'
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
  OUT="$OUT" $RS - <<'RS'
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
  MIN_COV="$MIN_COV" OUT="$OUT" $RS - <<'RS'
suppressPackageStartupMessages(library(data.table))
out <- Sys.getenv("OUT")
MIN_COV <- as.integer(Sys.getenv("MIN_COV"))
d <- fread(cmd = paste("gzip -cd", shQuote(file.path(out, "editlevel_long.tsv.gz"))),
           colClasses = list(character = "group"), showProgress = FALSE)
PAIRS <- list(c("coord_dedup", "umi_dedup"), c("coord_dedup", "umimd_dedup"), c("umi_dedup", "umimd_dedup"))
pair_id <- function(p) paste(p, collapse = "_vs_")

cor_row <- function(a, b) data.table(n_sites = length(a),
  pearson = cor(a, b), spearman = cor(a, b, method = "spearman"),
  mae = mean(abs(a - b)), median_abs_delta = median(abs(a - b)))

w_all <- dcast(d, tier + group + arm + site ~ bam_type, value.var = c("cov", "edit"))
by_grp <- list(); pts <- list(); pool <- list()
for (p in PAIRS) {
  cx <- paste0("cov_", p); ex <- paste0("edit_", p)
  w <- w_all[get(cx[1]) >= MIN_COV & get(cx[2]) >= MIN_COV, c("tier", "group", "arm", "site", ex), with = FALSE]
  setnames(w, ex, c("edit_x", "edit_y"))
  by_grp[[pair_id(p)]] <- w[, if (.N >= 20L) cor_row(edit_x, edit_y), by = .(tier, arm, group)]
  pt <- w[, .(n_groups = .N, edit_x = mean(edit_x), edit_y = mean(edit_y)), by = .(tier, arm, site)]
  pts[[pair_id(p)]] <- pt
  pool[[pair_id(p)]] <- pt[, cor_row(edit_x, edit_y), by = .(tier, arm)]
}
fwrite(rbindlist(by_grp, idcol = "pair"), file.path(out, "editlevel_correlation_by_group.tsv"), sep = "\t")
fwrite(rbindlist(pts, idcol = "pair"), file.path(out, "editlevel_pooled_points.tsv.gz"), sep = "\t")
pool <- rbindlist(pool, idcol = "pair")
fwrite(pool, file.path(out, "editlevel_correlation_pooled.tsv"), sep = "\t")
cat("\n=== edit level, pairs of BAM types (pooled over groups) ===\n"); print(pool)

dep <- rbindlist(lapply(c("ADAR1_KO", "ADARB1_KO"), function(ko) {
  n <- d[arm == "NTC_ctrl" & cov >= MIN_COV, .(tier, group, site, bam_type, edit_ntc = edit)]
  k <- d[arm == ko & cov >= MIN_COV, .(tier, group, site, bam_type, edit_ko = edit)]
  m <- merge(n, k, by = c("tier", "group", "site", "bam_type"))
  m <- m[, .(n_groups = .N, edit_ntc = mean(edit_ntc), edit_ko = mean(edit_ko)),
         by = .(tier, site, bam_type)][edit_ntc > 0]
  m[, `:=`(ko_arm = ko, ratio = edit_ko / edit_ntc)]
}))
dw_all <- dcast(dep, tier + ko_arm + site ~ bam_type, value.var = "ratio")
dpts <- list(); conc <- list()
for (p in PAIRS) {
  dw <- dw_all[!is.na(get(p[1])) & !is.na(get(p[2])), c("tier", "ko_arm", "site", p), with = FALSE]
  setnames(dw, p, c("ratio_x", "ratio_y"))
  dpts[[pair_id(p)]] <- dw
  conc[[pair_id(p)]] <- dw[, .(n_sites = .N,
               agree_50 = mean((ratio_x <= 0.5) == (ratio_y <= 0.5)),
               both_dep_50 = mean(ratio_x <= 0.5 & ratio_y <= 0.5),
               only_x_50 = mean(ratio_x <= 0.5 & ratio_y > 0.5),
               only_y_50 = mean(ratio_x > 0.5 & ratio_y <= 0.5),
               spearman_ratio = cor(ratio_x, ratio_y, method = "spearman")), by = .(tier, ko_arm)]
}
fwrite(rbindlist(dpts, idcol = "pair"), file.path(out, "depletion_ratio_points.tsv.gz"), sep = "\t")
conc <- rbindlist(conc, idcol = "pair")
fwrite(conc, file.path(out, "depletion_concordance.tsv"), sep = "\t")
cat("\n=== KO/NTC depletion ratio, pairs of BAM types ===\n"); print(conc)
RS
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
