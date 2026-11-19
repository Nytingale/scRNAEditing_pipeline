#!/bin/bash
# molecule_check_v9.sh   STAGE=manifest|sites|count|apply|score|rank|submit|submit2   (2026-11-17)
# v9 / 76-configuration copy of molecule_check_v8.sh (Reviewer 2, major 1; response paragraph 366): >= 3 variant
# molecules (distinct cell barcode + UMI) instead of >= 3 variant reads, on the step-5 tables of the 76 ranked
# configurations (gt_eval_v9/filter_sweep76/results/centre: the reruns on seed-8787 umi_tools and UMI-aware
# MarkDuplicates BAMs, plus links to the unchanged configurations), NTC units of the three tiers (cell types: the 8 of
# celltype_keep_pwcov.txt). Each unit x BAM type's BAM is traced from a calling task of filter_sweep76 (sampleUmiDup,
# sampleUmiMd) or filter_sweep2 (the unchanged BAM types). Scored on the v9 truth as scoring_pwcov76 (GT_FP_RULE=ag).
#   manifest  MW/manifest.tsv (id, tier, caller, prefix, rel_dir, calls, work_dir)
#   sites     MW/bams.tsv, MW/sites/<tier>/<unit>__<bam_type>.tsv.gz (called positions; arm also the v9 truth positions)
#   count / apply / score / rank as molecule_check_v8.sh -> MW/{counts,results,summary,scoring_pwcov}/,
#             gt_v9_molecule_compare.R -> MW/compare_rank.tsv, positives_molecules.tsv
#   On the Slurm login node (cluster/slurm on PATH): STAGE=submit bash human/pm/scripts/molecule_check_v9.sh
#$ -S /bin/bash
#$ -N molchk9
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=8G
#$ -l h_rt=12:00:00
set -euo pipefail
export LC_ALL=C
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
V9=$PM/gt_eval_v9
MW=$V9/molecules
CEN=$V9/filter_sweep76/results/centre
KEEP=$V8/unit_cov/celltype_keep_pwcov.txt
SELF=$PM/scripts/molecule_check_v9.sh
ENV=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin
RS=$ENV/Rscript
FA=$PM/../ref_annots/hg38.analysisSet.fa
export PATH=$ENV:$PATH
export MOL_DIR=$MW MOL_MANIFEST=$MW/manifest.tsv
STAGE=${STAGE:-${1:-}}
TIERS=(arm sample celltype)
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
declare -A TIER_OF=([results_gt_v8]=arm [results_gt_persample_v8]=sample [results_gt_celltype_v8]=celltype)

case "$STAGE" in
manifest)
  mkdir -p "$MW"
  # work dirs per tier x unit x bam type (any caller): the calling task whose input links to the BAM
  awk -F'\t' -v OFS='\t' 'NR>1{split($7, a, "/"); print $2, a[2], a[5], $3}' \
    "$V9/filter_sweep76/manifest.tsv" "$V8/filter_sweep2/manifest.tsv" > "$MW/workdirs.tsv"
  { printf "id\ttier\tcaller\tprefix\trel_dir\tcalls\n"
    for t in results_gt_v8 results_gt_persample_v8 results_gt_celltype_v8; do
      find -L "$CEN/$t" -name '*_step5_nomm.tsv.gz' | sort | while read -r f; do
        rel=$(dirname "${f#$CEN/}"); IFS=/ read -r _ unit _ caller bt _ <<< "$rel"
        [ "$bt" = sample ] && continue
        case $unit in NTC_ctrl|T[1-4]_NTC_ctrl) ;; NTC_ctrl__*) grep -qx "$unit" "$KEEP" || continue ;; *) continue ;; esac
        printf "%s\t%s\t%s\t%s\t%s\n" "${TIER_OF[$t]}" "$caller" "$(basename "$f" _step5_nomm.tsv.gz)" "$rel" "$f"
      done
    done | awk -v OFS='\t' '{print NR, $0}'; } > "$MW/manifest.tsv"
  awk -F'\t' 'NR>1{split($5, a, "/"); print $2, a[5]}' "$MW/manifest.tsv" | sort | uniq -c
  awk -F'\t' 'NR>1{split($5, a, "/"); c[$2" "a[2]]++} END{for (k in c) n[c[k]]++; for (v in n) print n[v], "units with", v, "call sets"}' "$MW/manifest.tsv"
  ;;

sites)
  mkdir -p "$MW/sites"
  printf "tier\tunit\tbam_type\tbam\tmin_bq\tsites\tcounts\n" > "$MW/bams.tsv"
  awk -F'\t' -v OFS='\t' 'NR>1{split($5, a, "/"); print $2, a[2], a[5]}' "$MW/manifest.tsv" | sort -u |
  while IFS=$'\t' read -r tier unit bt; do
    bam=""
    while IFS=$'\t' read -r _ _ _ wd; do
      i=$(grep -oP -- '--input \K\S+' "$wd/.command.sh") || continue
      t=$(readlink -f "$wd/$i") || continue
      b=$(ls "$(dirname "$t")"/*.bam 2>/dev/null | head -1) || true
      [ -n "$b" ] && { bam=$(readlink -f "$b"); break; }
    done < <(awk -F'\t' -v t="$tier" -v u="$unit" -v b="$bt" '$1==t && $2==u && $3==b' "$MW/workdirs.tsv")
    [ -n "$bam" ] && [ -s "$bam.bai" -o -s "${bam%.bam}.bai" ] || { echo "no BAM for $tier $unit $bt" >&2; exit 1; }
    case $bt in recalq20) q=20 ;; recaldefault) q=10 ;; *) q=30 ;; esac
    mkdir -p "$MW/sites/$tier" "$MW/counts/$tier"
    s=$MW/sites/$tier/${unit}__$bt.tsv.gz
    {
      awk -F'\t' -v t="$tier" -v u="$unit" -v b="$bt" 'NR>1{split($5, a, "/"); if ($2==t && a[2]==u && a[5]==b) print $6}' \
        "$MW/manifest.tsv" | while read -r f; do zcat "$f" | awk -F'\t' 'NR>1{print $2"\t"$3}'; done
      if [ "$tier" = arm ]; then zcat "$V9/truth/gt_v9_truth_sites.tsv.gz" | awk -F'\t' 'NR>1{print $3"\t"$4}'; fi
    } | sort -u -k1,1 -k2,2n | gzip > "$s"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tier" "$unit" "$bt" "$bam" "$q" "$s" "$MW/counts/$tier/${unit}__$bt.tsv.gz" \
      >> "$MW/bams.tsv"
    echo "$tier $unit $bt: $(zcat "$s" | wc -l) sites, BQ >= $q, $bam"
  done
  ;;

count)
  IFS=$'\t' read -r tier unit bt bam q s c < <(sed -n "$(( SGE_TASK_ID + 1 ))p" "$MW/bams.tsv")
  [ -s "$c" ] && { echo "skip (done): $c"; exit 0; }
  mkdir -p "$(dirname "$c")"
  case $bt in recalq20|recaldefault) mq=20 ;; *) mq=30 ;; esac
  echo "$tier $unit $bt MAPQ>=$mq BQ>=$q threads=${NSLOTS:-1}"
  $ENV/python $B/human/scripts/molecule_counts_v8.py --bam "$bam" --fasta "$FA" --sites "$s" --min_bq "$q" \
    --min_mapq "$mq" --threads "${NSLOTS:-1}" --out "$c.tmp.gz"
  mv "$c.tmp.gz" "$c"
  ;;

apply)
  cd "$B"; $RS human/scripts/molecule_filter_v8.R "${TIERS[$(( SGE_TASK_ID - 1 ))]}"
  ;;

score)
  tier=${TIERS[$(( SGE_TASK_ID - 1 ))]}; cd "$B"
  GT_UNIT_KEEP=$KEEP GT_UNIT_RULE=gate_cov GT_FP_RULE=ag GT_V7_TRUTH=$V9/truth/gt_v9_truth_sites.tsv.gz \
  GT_UNIT_COV=$V9/unit_cov/ko_truth_unit_cov.tsv.gz \
    $RS human/scripts/gt_score_v7.R "$tier" "$MW/results/${PUB[$tier]}" "$MW/scoring_pwcov/mismatch/$tier" '_step5_nomm\.tsv\.gz$'
  ;;

rank)
  cd "$B"
  $RS human/scripts/gt_rank_v7.R "$MW/scoring_pwcov/mismatch" > "$MW/scoring_pwcov/mismatch/rank.log" 2>&1
  $RS human/scripts/gt_v9_molecule_compare.R
  ;;

submit)
  cd "$B"; mkdir -p "$MW"
  js=$(qsub -terse -N molchk9_sites -l h_vmem=4G -l h_rt=4:00:00 -v STAGE=sites "$SELF")
  jn=$(qsub -terse -N molchk9_sub2 -l h_vmem=2G -l h_rt=0:30:00 -hold_jid "$js" -v STAGE=submit2 "$SELF")
  echo "sites $js | submit2 $jn (count, apply, score, rank)"
  ;;

submit2)
  # after the sites stage: the arm BAMs are large, so the arm rows get more threads
  cd "$B"; n=$(( $(wc -l < "$MW/bams.tsv") - 1 ))
  arm=$(awk -F'\t' 'NR>1 && $1=="arm"' "$MW/bams.tsv" | wc -l)
  jc=$(qsub -terse -N molchk9_cnt_arm -t 1-$arm -pe threaded 8 -l h_vmem=3G -l h_rt=24:00:00 -v STAGE=count "$SELF" | cut -d. -f1)
  jd=$(qsub -terse -N molchk9_cnt -t $(( arm + 1 ))-$n -tc 100 -pe threaded 4 -l h_vmem=3G -l h_rt=12:00:00 \
         -v STAGE=count "$SELF" | cut -d. -f1)
  ja=$(qsub -terse -N molchk9_apply -t 1-3 -l h_vmem=24G -l h_rt=6:00:00 -hold_jid "$jc,$jd" -v STAGE=apply "$SELF" | cut -d. -f1)
  jsc=$(qsub -terse -N molchk9_score -t 1-3 -l h_vmem=48G -l h_rt=12:00:00 -hold_jid "$ja" -v STAGE=score "$SELF" | cut -d. -f1)
  jr=$(qsub -terse -N molchk9_rank -l h_vmem=16G -l h_rt=2:00:00 -hold_jid "$jsc" -v STAGE=rank "$SELF")
  echo "count arm $jc, others $jd | apply $ja | score $jsc | rank $jr"
  ;;

*) echo "unknown STAGE '$STAGE'" >&2; exit 1 ;;
esac
