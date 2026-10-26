#!/bin/bash
# =============================================================================
# run_dedup_compare_v6.sh  [submit|depth|collect|discover TIER]
#
# Inputs for the UMI vs coordinate de-duplication figure on the v6 data
# (publication/teratoma_gt/scripts/TERATOMA_dedup_compare_plot.R).
#
#   depth     per teratoma x arm (20 units, array): primary-read counts and
#             per-site samtools depth (BQ/MQ 30) at the 1,692 v6 ground-truth
#             sites, for three BAMs:
#               unprocessed  results_gt_persample_v6/<U>/01_merged/<U>_merged.bam
#               coord_dedup  02_gt_bams_rmdup_v6/<U>.bam
#               umi_dedup    02_gt_bams_umi_v6/<U>.bam
#             Unprocessed is recorded (it is the denominator for reads removed)
#             but is not a reported configuration.
#   collect   concatenates the per-unit fragments.
#   discover  all-varType pre100 discovery on the NTC units only (the scored
#             arm): TIER = sample (T1-T4 NTC, results_gt_persample_v6) or
#             arm (merged NTC, results_gt_v6) or celltype (NTC_ctrl__<celltype>,
#             results_gt_celltype_v6). Feeds the mismatch spectrum.
#             Cell type is not in `submit`:
#               qsub -N dd_disc_c -q slice -pe threaded 1 -l h_vmem=64G \
#                 human/pm/scripts/run_dedup_compare_v6.sh discover celltype
#
# Outputs: human/pm/gt_eval_v6/dedup_compare/   Logs: logs/
# Do NOT run samtools on the login node; `submit` only calls qsub.
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -j y
#$ -pe threaded 4
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00

set -eo pipefail
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:/usr/local/bin:/usr/bin:/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"; source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing 2>/dev/null || true
export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${PATH}"

B="${SCRNA_DIR:-.}"
PM="$B/human/pm"
SELF="$PM/scripts/run_dedup_compare_v6.sh"
OUT="$PM/gt_eval_v6/dedup_compare"
BED="$PM/gt_eval_v6/truth/v6_truth_TRUE.bed"
PS="$PM/results_gt_persample_v6"
UNITS="$OUT/units.txt"

mkdir -p "$OUT"/{depth,discovery}
cd "$B"
STAGE="${1:-submit}"

if [ "$STAGE" = "submit" ]; then
  [ -s "$BED" ] || { echo "ERROR: missing $BED" >&2; exit 1; }
  : > "$UNITS"
  for a in NTC_ctrl AAVS1_ctrl ADAR1_KO ADARB1_KO ADARB2_KO; do
    for t in T1 T2 T3 T4; do echo "${t}_${a}" >> "$UNITS"; done
  done
  n=$(wc -l < "$UNITS")
  jd=$(qsub -terse -N dd_depth -q slice -t 1-"$n" "$SELF" depth | cut -d. -f1)
  jc=$(qsub -terse -N dd_collect -q slice -pe threaded 1 -l h_vmem=4G -hold_jid "$jd" "$SELF" collect)
  js=$(qsub -terse -N dd_disc_s -q slice -pe threaded 1 -l h_vmem=32G "$SELF" discover sample)
  ja=$(qsub -terse -N dd_disc_a -q slice -pe threaded 1 -l h_vmem=32G "$SELF" discover arm)
  echo "depth $jd (1-$n) | collect $jc | discover sample $js | discover arm $ja"
  exit 0
fi

if [ "$STAGE" = "depth" ]; then
  U=$(sed -n "${SGE_TASK_ID:?}p" "$UNITS")
  declare -A BAM=(
    [unprocessed]="$PS/$U/01_merged/${U}_merged.bam"
    [coord_dedup]="$PM/02_gt_bams_rmdup_v6/$U.bam"
    [umi_dedup]="$PM/02_gt_bams_umi_v6/$U.bam"
  )
  for bt in unprocessed coord_dedup umi_dedup; do
    [ -s "${BAM[$bt]}" ] || { echo "ERROR: missing ${BAM[$bt]}" >&2; exit 1; }
  done
  READS="$OUT/depth/_reads_${U}.tsv"; SITES="$OUT/depth/_sites_${U}.tsv"
  echo -e "unit\tbam_type\tn_reads" > "$READS"
  echo -e "unit\tbam_type\tchr\tpos\tdepth" > "$SITES"
  for bt in unprocessed coord_dedup umi_dedup; do
    echo "=== $(date) === $U $bt"
    # primary, mapped, not QC-fail, not duplicate-flagged
    n=$(samtools view -c -@ 4 -F 0xF04 "${BAM[$bt]}")
    echo -e "$U\t$bt\t$n" >> "$READS"
    samtools depth -a -q 30 -Q 30 -b "$BED" "${BAM[$bt]}" \
      | awk -v u="$U" -v bt="$bt" 'BEGIN{OFS="\t"} {print u, bt, $1, $2, $3}' >> "$SITES"
  done
  nexp=$(( $(wc -l < "$BED") * 3 + 1 ))
  [ "$(wc -l < "$SITES")" -eq "$nexp" ] || { echo "ERROR: $SITES has $(wc -l < "$SITES") lines, expected $nexp" >&2; exit 1; }
  touch "$OUT/depth/_ok_${U}"
  cat "$READS"
  exit 0
fi

if [ "$STAGE" = "collect" ]; then
  miss=0
  while read -r U; do [ -f "$OUT/depth/_ok_${U}" ] || { echo "ERROR: incomplete $U" >&2; miss=1; }; done < "$UNITS"
  [ "$miss" -eq 0 ] || exit 1
  first=1
  while read -r U; do
    if [ "$first" -eq 1 ]; then
      cat "$OUT/depth/_reads_${U}.tsv" > "$OUT/dedup_read_counts.tsv"
      cat "$OUT/depth/_sites_${U}.tsv" > "$OUT/dedup_truth_site_depth.tsv"
      first=0
    else
      tail -n +2 "$OUT/depth/_reads_${U}.tsv" >> "$OUT/dedup_read_counts.tsv"
      tail -n +2 "$OUT/depth/_sites_${U}.tsv" >> "$OUT/dedup_truth_site_depth.tsv"
    fi
  done < "$UNITS"
  gzip -f "$OUT/dedup_truth_site_depth.tsv"
  column -t "$OUT/dedup_read_counts.tsv"
  exit 0
fi

if [ "$STAGE" = "discover" ]; then
  t="${2:?discover needs TIER (sample|arm|celltype)}"
  case "$t" in
    sample) res="$PS"; arms="T1_NTC_ctrl,T2_NTC_ctrl,T3_NTC_ctrl,T4_NTC_ctrl" ;;
    arm)    res="$PM/results_gt_v6"; arms="NTC_ctrl" ;;
    celltype)
      res="$PM/results_gt_celltype_v6"
      arms=$(ls -d "$res"/NTC_ctrl__*/04_filtered | awk -F/ '{print $(NF-1)}' | paste -sd,)
      [ -n "$arms" ] || { echo "ERROR: no NTC_ctrl__* units under $res" >&2; exit 1; }
      echo "celltype NTC units: $(echo "$arms" | tr ',' '\n' | wc -l)" ;;
    *) echo "ERROR: unknown tier $t" >&2; exit 1 ;;
  esac
  Rscript "$B/human/scripts/gt_eval_sites.R" --results_dir "$res" \
    --outdir "$OUT/discovery/$t" --arms "$arms" --table pre100 --all_vartypes true
  exit 0
fi

echo "ERROR: unknown stage '$STAGE'" >&2; exit 1
