#!/bin/bash
# baqoff_v9.sh   (2026-11-17) v9 / 76-configuration copy of baqoff_v8.sh (response paragraph 382, BAQ on vs off):
#   new BAQ-off pileups only for the BAM types that changed since the v8 run (sampleUmiDup: seed-8787 umi_tools;
#   sampleUmiMd: new), rows of gt_eval_v9/filter_sweep76/manifest.tsv; the v8 BAQ-off call sets of the unchanged
#   BAM types (sampleRmDup, recal, recalq20; gt_eval_v8/baq_compare/current) are linked in (STAGE=link_old);
#   all 20 BCFtools configurations are scored on the v9 truth (GT_FP_RULE=ag, as scoring_pwcov76) and compared with
#   scoring_pwcov76 (gt_v9_baq_compare.R). Submit from the Slurm login node: STAGE=submit_pile bash baqoff_v9.sh
# baqoff_v8.sh   STAGE=manifest|filter|score|compare|submit   (Reviewer 2, major 3; user 2026-10-22)
# BCFtools with BAQ off (-B) through the current downstream filters, so that BAQ is the only difference from the
# published call sets (gt_v8_baq_compare.R used the 2026-10-14 filter state). BAQ-off pileups (user 2026-11-06,
# TSO-clipped tree): baq_compare/pile_nobaq/<tier>/<work hash>/ = the BCFTOOLS_MPILEUP task's own .command.sh with
# `nobaq` (same BAMs and options). The 16 reported BCFtools configurations, every NTC unit of
# gt_eval_v8/filter_sweep2/manifest.tsv.
#   piles     O/piles.tsv (tier, BCFTOOLS_MPILEUP work dir) of those configurations
#   pile      array over O/piles.tsv -> PILE/<tier>/<hash>/; then STAGE=after_pile (held): manifest + submit
#   manifest  B/manifest.tsv  (filter_sweep2 manifest columns + BAQ-off pileup dir)
#   filter    array: filter_sweep2 union + apply (centre) with the task's variant and blacklist inputs swapped
#             -> B/results/centre/<rel_dir>/
#   score     array 1-4: gt_score_v7.R (arm / sample / celltype; headline scoring_pwcov settings) and
#             gt_score_rediportal_v8.R -> B/scoring_pwcov/<tier>/, B/rediportal/
#   compare   gt_v8_singlecell_modes_compare.R baqoff
#   STAGE=manifest bash baqoff_v8.sh; STAGE=submit bash baqoff_v8.sh
#$ -S /bin/bash
#$ -N baqoff9
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=12:00:00
set -uo pipefail
export LC_ALL=C
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
V9=$PM/gt_eval_v9
W=$V9/filter_sweep76
O=$V9/baq_compare
PILE=$O/pile_nobaq
OLD=$V8/baq_compare/current
SELF=$PM/scripts/baqoff_v9.sh
SWEEP_R=$B/human/scripts/downstream_filter_sweep.R
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
STAGE=${STAGE:-${1:-}}
TIERS=(arm sample celltype)
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
row() { sed -n "$(( SGE_TASK_ID + 1 ))p" "${M:-$O/manifest.tsv}"; }
sweep_cmd() {  # <work_dir> <extra args>: the task's own command with the sweep script (as filter_sweep2_v8.sh)
  awk -v r="$SWEEP_R" -v x="$2" '
    /^Rscript .*downstream_filter(_v8|_scallele)?\.R/ { sub(/[^ ]*downstream_filter(_v8|_scallele)?\.R/, r) }
    { lines[++n] = $0; if ($0 ~ /[^[:space:]]/) last = n }
    END { sub(/[ \\]*$/, " \\", lines[last]); for (i = 1; i <= last; i++) print lines[i]; print "    " x }
  ' "$1/.command.sh"
}

case "$STAGE" in
piles)
  set -e; mkdir -p "$O"
  awk -F'\t' 'NR>1 && $5=="bcftools" && $7 !~ /\/bcftools\/sample\//{print $2"\t"$3}' "$W/manifest.tsv" |
    while IFS=$'\t' read -r t wd; do
      i=$(grep -oP -- '--input \K\S+' "$wd/.command.sh"); printf "%s\t%s\n" "$t" "$(dirname "$(readlink -f "$wd/$i")")"
    done | sort -u > "$O/piles.tsv"
  cut -f1 "$O/piles.tsv" | uniq -c
  ;;

pile)
  export PATH=$B/nextflow/bin:${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
  IFS=$'\t' read -r tier wd < <(sed -n "${SGE_TASK_ID}p" "$O/piles.tsv")
  d=$PILE/$tier/$(basename "$(dirname "$wd")")_$(basename "$wd")
  [ -e "$d/done" ] && { echo "skip (done): $d"; exit 0; }
  rm -rf "$d"; mkdir -p "$d"
  for e in "$wd"/*; do
    n=$(basename "$e")
    case "$n" in *_bcftools_variants.txt.gz|*_bcftools_snpblacklist.txt.gz) continue ;; esac
    ln -s "$(readlink -f "$e")" "$d/$n"
  done
  sed 's/ baq$/ nobaq/' "$wd/.command.sh" > "$d/run.sh"
  grep -q ' nobaq$' "$d/run.sh" || { echo "no baq argument in $wd/.command.sh" >&2; exit 2; }
  ( cd "$d" && bash run.sh > run.log 2>&1 ) || { echo "pileup failed: $d" >&2; exit 1; }
  find "$d" -maxdepth 1 -type l -delete
  touch "$d/done"
  ;;

submit_pile)
  STAGE=piles bash "$SELF"
  n=$(wc -l < "$O/piles.tsv")
  jp=$(qsub -terse -N baqoff9_pile -t 1-$n -tc 30 -pe threaded 12 -l h_vmem=4G -l h_rt=24:00:00 -v STAGE=pile "$SELF" | cut -d. -f1)
  ja=$(qsub -terse -N baqoff9_after -l h_vmem=4G -l h_rt=1:00:00 -hold_jid "$jp" -v STAGE=after_pile "$SELF")
  echo "pile $jp (1-$n) | after_pile $ja"
  ;;

link_old)
  # v8 BAQ-off call sets of the BAM types that did not change (same BAMs and downstream filter as published)
  n=0
  for d in "$OLD"/results/centre/results_gt*/*/04_filtered/bcftools/{sampleRmDup,recal,recalq20}; do
    [ -d "$d" ] || continue
    t=$O/results/centre/${d#$OLD/results/centre/}
    mkdir -p "$(dirname "$t")"; [ -e "$t" ] || { ln -s "$d" "$t"; n=$(( n + 1 )); }
  done
  echo "linked $n v8 BAQ-off directories"
  ;;

after_pile)
  n=$(wc -l < "$O/piles.tsv"); k=$(ls "$PILE"/*/*/done 2>/dev/null | wc -l)
  [ "$k" -eq "$n" ] || { echo "only $k of $n BAQ-off pileups done" >&2; exit 1; }
  STAGE=manifest bash "$SELF" && STAGE=link_old bash "$SELF" && STAGE=submit bash "$SELF"
  ;;

manifest)
  set -e; mkdir -p "$O"
  { head -1 "$W/manifest.tsv" | awk -v OFS='\t' '{print $0, "pile"}'
    awk -F'\t' 'NR>1 && $5=="bcftools" && $7 !~ /\/bcftools\/sample\//' "$W/manifest.tsv" |
    while IFS=$'\t' read -r id tier wd pub caller prefix rel mem; do
      i=$(grep -oP -- '--input \K\S+' "$wd/.command.sh")
      p=$(ls "$PILE/$tier"/*/"$i" 2>/dev/null)
      [ "$(echo "$p" | grep -c .)" = 1 ] || { echo "BAQ-off pileup for $id ($tier $i): '$p'" >&2; exit 1; }
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$id" "$tier" "$wd" "$pub" "$caller" "$prefix" "$rel" "$mem" "$(dirname "$p")"
    done; } > "$O/manifest.tsv"
  cut -f2 "$O/manifest.tsv" | tail -n +2 | sort | uniq -c
  ;;

filter)
  export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
  mkdir -p "$O/status" "$O/runs"
  IFS=$'\t' read -r id tier wd pub caller prefix rel mem pile < <(row)
  [ "$(cut -f2 "$O/status/$id.tsv" 2>/dev/null)" = 0 ] && { echo "skip (done): $id"; exit 0; }
  mkdir -p "$O/union/$tier"; u=$O/union/$tier/$prefix.rds; rm -f "$u"
  run=$O/runs/f$id; rm -rf "$run"; mkdir -p "$run"
  for e in "$wd"/*; do
    n=$(basename "$e")
    if [ -L "$e" ]; then cp -P "$e" "$run/"; elif [ -d "$e" ]; then mkdir -p "$run/$n"; cp -P "$e"/* "$run/$n/"; fi
  done
  rc=0; ns=0
  for f in "$run"/*_bcftools_variants.txt.gz "$run"/*_bcftools_snpblacklist.txt.gz; do
    [ -L "$f" ] || continue
    n=$(basename "$f")
    [ -s "$pile/$n" ] || { echo "no BAQ-off $n in $pile" > "$run/union.log"; rc=4; }
    rm -f "$f"; ln -s "$pile/$n" "$f"; ns=$(( ns + 1 ))
  done
  if [ $rc -eq 0 ]; then
    sweep_cmd "$wd" "--sweep_mode union --sweep_union $u" > "$run/union.sh"
    sweep_cmd "$wd" "--sweep_mode apply --sweep_union $u --sweep_settings_file $W/settings.tsv --sweep_settings centre --sweep_out $O/results/{setting}/$rel" > "$run/apply.sh"
    ( cd "$run" && echo "swapped $ns inputs" > union.log && bash union.sh >> union.log 2>&1 && bash apply.sh > apply.log 2>&1 ); rc=$?
  fi
  cat "$run"/*.log > "$O/status/$id.log"
  printf "%s\t%s\n" "$id" "$rc" > "$O/status/$id.tsv"
  [ $rc -eq 0 ] && rm -rf "$run"
  exit $rc
  ;;

score)
  cd "$B"
  if [ "$SGE_TASK_ID" -le 3 ]; then
    tier=${TIERS[$(( SGE_TASK_ID - 1 ))]}
    GT_UNIT_KEEP=$V8/unit_cov/celltype_keep_pwcov.txt GT_UNIT_RULE=gate_cov GT_FP_RULE=ag \
    GT_V7_TRUTH=$V9/truth/gt_v9_truth_sites.tsv.gz GT_UNIT_COV=$V9/unit_cov/ko_truth_unit_cov.tsv.gz \
      $RS human/scripts/gt_score_v7.R "$tier" "$O/results/centre/${PUB[$tier]}" "$O/scoring_pwcov/$tier" '_step5_nomm\.tsv\.gz$'
  else
    r=${O#$PM/}/results/centre
    GT_REDI_EXCLUDE_BED=$PM/pgp1_snps/PGP1_hu43860C_hg38_snvs.bed GT_CALL_RE='_step5_nomm\.tsv\.gz$' \
    GT_RES_DIRS="arm=$r/${PUB[arm]},sample=$r/${PUB[sample]},celltype=$r/${PUB[celltype]}" \
    GT_REDI_OUT=$O/rediportal $RS human/scripts/gt_score_rediportal_v8.R
  fi
  ;;

compare)
  cd "$B"; $RS human/scripts/gt_rank_v7.R "$O/scoring_pwcov" > "$O/scoring_pwcov/rank.log" 2>&1
  $RS human/scripts/gt_v9_baq_compare.R
  ;;

submit)
  cd "$B"; mkdir -p "$O/sub"; jf=""
  for m in $(awk -F'\t' 'NR>1{print $8}' "$O/manifest.tsv" | sort -nu); do
    f=$O/sub/filter_${m}G.tsv; { head -1 "$O/manifest.tsv"; awk -F'\t' -v m="$m" 'NR>1 && $8==m' "$O/manifest.tsv"; } > "$f"
    k=$(( $(wc -l < "$f") - 1 ))
    j=$(qsub -terse -N baqoff9_f$m -t 1-$k -tc 200 -l h_vmem=${m}G -v STAGE=filter,M=$f "$SELF" | cut -d. -f1)
    jf=$jf${jf:+,}$j
  done
  js=$(qsub -terse -N baqoff9_score -t 1-3 -l h_vmem=48G -l h_rt=12:00:00 -hold_jid "$jf" -v STAGE=score "$SELF" | cut -d. -f1)
  jm=$(qsub -terse -N baqoff9_cmp -l h_vmem=8G -l h_rt=1:00:00 -hold_jid "$js" -v STAGE=compare "$SELF")
  echo "filter $jf | score $js | compare $jm"
  ;;

*) echo "unknown STAGE '$STAGE'" >&2; exit 1 ;;
esac
