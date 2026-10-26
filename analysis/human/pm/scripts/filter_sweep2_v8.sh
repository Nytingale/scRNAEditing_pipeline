#!/bin/bash
# filter_sweep2_v8.sh   STAGE=manifest|submit|union|apply|check|score|rank   (2026-10-21)
# Downstream-filter sweep, second design (user 2026-10-21): coverage 5 / 10 / 15, variant reads 2 / 3 / 4, editing
# level 1 / 2 / 3%, one at a time from the main values (10 reads, 3 variant reads, 1%); dbSNP, multi-allelic and
# JACUSA2-score removal and the nearby-mismatch window filter stay on (not swept). Every threshold is applied to
# every caller (REDItools3 included, r3floor_rerun_v8.sh). Same NTC tasks as filter_sweep_v8.sh (its manifest.tsv).
#   union   BCFtools / REDItools3: the filter_sweep_v8.sh union tables (parse min_cov 1; no caller floor).
#           GATK / JACUSA2 / cellSNP-lite / scAllele: the task's input is swapped for the replayed caller output
#           (recall_cov5_v8.sh: caller min_cov 5, scAllele --AC 2), then parsed as filter_sweep_v8.sh.
#   apply   human/scripts/downstream_filter_sweep.R per setting (W/settings.tsv) -> W/results/<setting>/<rel_dir>/
#   check   centre step 4 / step 5 vs the published tables -> W/check_centre.tsv
#   score   array 1..3 x (settings + 1): gt_score_v7.R -> W/<SCN>/<setting>/mismatch/<tier>/ (last = window_off: centre
#           step 4, no window filter; all_off is scored at step 4 too); SCN / GT_UNIT_RULE /
#           GT_UNIT_KEEP from -v (headline: scoring_pwcov, gate_cov; power only: scoring_pw, gate)
#   rank    gt_rank_v7.R per setting, then gt_v8_filter_sweep_compare.R (W, SCN)
#   extra settings outside the Fig 6 sweep (W/settings_extra.tsv, e.g. alt1 = Fig 8 Relaxed): apply with
#           -v SETF=W/settings_extra.tsv, score with -v ONE=<setting> (tasks 1-3) into a separate SCN
#           (submit_alt1.sh)
#   STAGE=manifest bash filter_sweep2_v8.sh; HOLD=<recall gather + r3floor swap ids> STAGE=submit bash filter_sweep2_v8.sh
#$ -S /bin/bash
#$ -N fsweep2
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=48:00:00
set -uo pipefail
export LC_ALL=C
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V8=$PM/gt_eval_v8
W1=$V8/filter_sweep
W=$V8/filter_sweep2
REC=$W/recall
SELF=$PM/scripts/filter_sweep2_v8.sh
SWEEP_R=$B/human/scripts/downstream_filter_sweep.R
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
STAGE=${STAGE:-${1:-}}
CHUNK=${CHUNK:-1}
TIERS=(arm sample celltype)
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
mapfile -t SETTINGS < <(awk -F'\t' 'NR>1{print $1}' "$W/settings.tsv")
SCN=${SCN:-scoring_pwcov}
RECALLED=" gatk jacusa2 cellsnp scallele "

rows() {  # manifest rows for this array element
  local m=$1 first=$(( (SGE_TASK_ID - 1) * CHUNK + 2 ))
  sed -n "${first},$(( first + CHUNK - 1 ))p" "$m"
}
sweep_cmd() {  # <work_dir> <extra args>: the task's own command with the sweep script
  awk -v r="$SWEEP_R" -v x="$2" '
    /^Rscript .*downstream_filter(_scallele)?\.R/ { sub(/[^ ]*downstream_filter(_scallele)?\.R/, r) }
    { lines[++n] = $0; if ($0 ~ /[^[:space:]]/) last = n }
    END { sub(/[ \\]*$/, " \\", lines[last]); for (i = 1; i <= last; i++) print lines[i]; print "    " x }
  ' "$1/.command.sh"
}

case "$STAGE" in
manifest)
  set -e; mkdir -p "$W/sub"
  # re-called callers parse a larger (min_cov 5) input: twice the filter_sweep_v8.sh memory, at most 128G
  awk -F'\t' -v OFS='\t' -v r="$RECALLED" 'NR>1 && index(r, " "$5" ") {$8 = ($8 * 2 > 128 ? 128 : $8 * 2)} 1' \
    "$W1/manifest.tsv" > "$W/manifest.tsv"
  # apply rows: arm = one setting per row, other tiers = every setting of the task in one row
  printf "id\tsettings\tmem\n" > "$W/apply_manifest.tsv"
  all=$(IFS=,; echo "${SETTINGS[*]}")
  awk -F'\t' -v all="$all" 'NR>1{ if ($2 == "arm") { n = split(all, s, ","); for (i = 1; i <= n; i++) print $1"\t"s[i]"\t"$8 }
                                  else print $1"\t"all"\t"$8 }' "$W/manifest.tsv" >> "$W/apply_manifest.tsv"
  awk -F'\t' 'NR>1{n[$2" "$5]++} END{for (k in n) print k, n[k]}' "$W/manifest.tsv" | sort
  echo "apply rows: $(( $(wc -l < "$W/apply_manifest.tsv") - 1 ))"
  ;;

submit)
  # HOLD: the recall_cov5_v8.sh arrays (union of re-called callers) and r3floor swap (check compares to published)
  cd "$B"; mkdir -p "$W/sub"; hu=""; ha=""
  for m in $(awk -F'\t' 'NR>1{print $8}' "$W/manifest.tsv" | sort -nu); do
    f=$W/sub/union_${m}G.tsv; { head -1 "$W/manifest.tsv"; awk -F'\t' -v m="$m" 'NR>1 && $8==m' "$W/manifest.tsv"; } > "$f"
    n=$(( $(wc -l < "$f") - 1 ))
    j=$(qsub -terse -N fsw2_u$m -t 1-$n -tc 300 -l h_vmem=${m}G ${HOLD:+-hold_jid $HOLD} -v STAGE=union,M=$f "$SELF" | cut -d. -f1)
    hu=$hu${hu:+,}$j
  done
  for m in $(awk -F'\t' 'NR>1{print $3}' "$W/apply_manifest.tsv" | sort -nu); do
    f=$W/sub/apply_${m}G.tsv; { head -1 "$W/apply_manifest.tsv"; awk -F'\t' -v m="$m" 'NR>1 && $3==m' "$W/apply_manifest.tsv"; } > "$f"
    n=$(( $(wc -l < "$f") - 1 ))
    j=$(qsub -terse -N fsw2_a$m -t 1-$n -tc 300 -l h_vmem=${m}G -hold_jid "$hu" -v STAGE=apply,M=$f "$SELF" | cut -d. -f1)
    ha=$ha${ha:+,}$j
  done
  jc=$(qsub -terse -N fsw2_check -l h_vmem=4G -hold_jid "$ha" -v STAGE=check "$SELF")
  ns=$(( 3 * (${#SETTINGS[@]} + 1) )); K=$V8/unit_cov; out="union $hu | apply $ha | check $jc"
  for v in "scoring_pwcov gate_cov celltype_keep_pwcov.txt" "scoring_pw gate celltype_keep_min20_pw.txt"; do
    read -r scn rule keep <<< "$v"
    js=$(qsub -terse -N fsw2_s_${scn#scoring_} -t 1-$ns -l h_vmem=48G -l h_rt=12:00:00 -hold_jid "$ha" \
           -v STAGE=score,SCN=$scn,GT_UNIT_RULE=$rule,GT_UNIT_KEEP=$K/$keep "$SELF" | cut -d. -f1)
    jr=$(qsub -terse -N fsw2_r_${scn#scoring_} -l h_vmem=16G -hold_jid "$js" -v STAGE=rank,SCN=$scn "$SELF")
    out="$out | $scn score $js rank $jr"
  done
  echo "$out"
  ;;

union)
  mkdir -p "$W/status/union" "$W/runs"
  rows "${M:-$W/manifest.tsv}" | while IFS=$'\t' read -r id tier wd pub caller prefix rel mem; do
    st=$W/status/union/$id.tsv
    if [ -s "$st" ] && [ "$(cut -f2 "$st")" = 0 ]; then echo "skip (done): $id"; continue; fi
    mkdir -p "$W/union/$tier"; u=$W/union/$tier/$prefix.rds; rm -f "$u"
    if [[ "$RECALLED" != *" $caller "* ]]; then
      s1=$W1/status/union/$id.tsv
      [ "$(cut -f2 "$s1" 2>/dev/null)" = 0 ] || { printf "%s\t9\tNA\n" "$id" > "$st"; continue; }
      [ -s "$W1/union/$tier/$prefix.rds" ] && ln -s "$W1/union/$tier/$prefix.rds" "$u"
      printf "%s\t0\t%s\n" "$id" "$(cut -f3 "$s1")" > "$st"; continue
    fi
    run=$W/runs/u$id; rm -rf "$run"; mkdir -p "$run"
    for e in "$wd"/*; do
      n=$(basename "$e")
      if [ -L "$e" ]; then cp -P "$e" "$run/"
      elif [ -d "$e" ]; then mkdir -p "$run/$n"; cp -P "$e"/* "$run/$n/"
      fi
    done
    i=$(grep -oP -- '--input \K\S+' "$wd/.command.sh"); t=$(readlink "$wd/$i")
    new=$REC/calls/$(basename "$(dirname "$(dirname "$t")")")/$(basename "$(dirname "$t")")/$(basename "$t")
    rc=0
    if [ -s "$new" ]; then
      ln -sfn "$new" "$run/$i"
      sweep_cmd "$wd" "--sweep_mode union --sweep_union $u" > "$run/run.sh"
      ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=$?
    else echo "no re-called input $new" > "$run/run.log"; rc=4; fi
    n=NA; [ -s "$u" ] && n=$(sed -n 's/.*\[sweep\] union table: \([0-9]*\) sites.*/\1/p' "$run/run.log")
    [ $rc -eq 0 ] && [ ! -s "$u" ] && n=0
    cp "$run/run.log" "$W/status/union/$id.log"
    printf "%s\t%s\t%s\n" "$id" "$rc" "$n" > "$st"
    [ $rc -eq 0 ] && rm -rf "$run"
  done
  ;;

apply)
  mkdir -p "$W/status/apply" "$W/runs"
  rows "${M:-$W/apply_manifest.tsv}" | while IFS=$'\t' read -r id settings mem; do
    key=$id.${settings//,/_}; [ ${#key} -gt 80 ] && key=$id.all
    st=$W/status/apply/$key.tsv
    if [ -s "$st" ] && [ "$(cut -f3 "$st")" = 0 ]; then echo "skip (done): $key"; continue; fi
    IFS=$'\t' read -r _ tier wd pub caller prefix rel mem < <(awk -F'\t' -v i="$id" '$1==i' "$W/manifest.tsv")
    u=$W/union/$tier/$prefix.rds
    us=$(cut -f2,3 "$W/status/union/$id.tsv" 2>/dev/null)
    run=$W/runs/a$key; rm -rf "$run"; mkdir -p "$run"
    if [ "$us" = $'0\t0' ]; then
      for s in ${settings//,/ }; do d=$W/results/$s/$rel; mkdir -p "$d"; printf "%s\n0\n" "$s" > "$d/$prefix.done"; done
      rc=0
    elif [ ! -s "$u" ]; then
      echo "no union table for $id ($us)" > "$run/run.log"; rc=3
    else
      sweep_cmd "$wd" "--sweep_mode apply --sweep_union $u --sweep_settings_file ${SETF:-$W/settings.tsv} --sweep_settings $settings --sweep_out $W/results/{setting}/$rel" > "$run/run.sh"
      ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=$?
    fi
    cp "$run/run.log" "$W/status/apply/$key.log"
    printf "%s\t%s\t%s\n" "$id" "$settings" "$rc" > "$st"
    [ $rc -eq 0 ] && rm -rf "$run"
  done
  ;;

check)
  printf "id\ttier\tcaller\tprefix\tstep\tpublished\tsweep\tidentical\n" > "$W/check_centre.tsv"
  awk -F'\t' 'NR>1' "$W/manifest.tsv" | while IFS=$'\t' read -r id tier wd pub caller prefix rel mem; do
    for s in step4_cov_edit_filtered step5_nomm; do
      p=$PM/$rel/${prefix}_$s.tsv.gz; q=$W/results/centre/$rel/${prefix}_$s.tsv.gz
      mp=absent; mq=absent
      [ -e "$p" ] && mp=$(zcat "$p" | md5sum | cut -c1-32)
      [ -e "$q" ] && mq=$(zcat "$q" | md5sum | cut -c1-32)
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$id" "$tier" "$caller" "$prefix" "$s" "$mp" "$mq" "$([ "$mp" = "$mq" ] && echo TRUE || echo FALSE)"
    done
  done >> "$W/check_centre.tsv"
  awk -F'\t' 'NR>1{n[$3" "$5" "$8]++} END{for (k in n) print k, n[k]}' "$W/check_centre.tsv" | sort
  ;;

score)
  i=$(( SGE_TASK_ID - 1 )); tier=${TIERS[$(( i % 3 ))]}; k=$(( i / 3 ))
  if [ -n "${ONE:-}" ]; then s=$ONE; src=$ONE; re='_step5_nomm\.tsv\.gz$'   # extra setting (SETF), tasks 1-3
  elif [ $k -lt ${#SETTINGS[@]} ]; then s=${SETTINGS[$k]}; src=$s; re='_step5_nomm\.tsv\.gz$'
  else s=window_off; src=centre; re='_step4_cov_edit_filtered\.tsv\.gz$'; fi
  # all_off: no coverage / variant-read / editing / JACUSA2-score threshold and no window filter (step 4)
  [ "$s" = all_off ] && re='_step4_cov_edit_filtered\.tsv\.gz$'
  echo "setting=$s tier=$tier SCN=$SCN GT_UNIT_RULE=${GT_UNIT_RULE:-cov} GT_UNIT_KEEP=${GT_UNIT_KEEP:-}"
  cd "$B"
  GT_UNIT_KEEP=${GT_UNIT_KEEP:?} GT_V7_TRUTH=$V8/truth/gt_v7_truth_sites.tsv.gz \
  GT_UNIT_COV=$V8/unit_cov/ko_truth_unit_cov.tsv.gz GT_UNIT_RULE=${GT_UNIT_RULE:?} \
    $RS human/scripts/gt_score_v7.R "$tier" "$W/results/$src/${PUB[$tier]}" "$W/$SCN/$s/mismatch/$tier" "$re"
  ;;

rank)
  cd "$B"
  for d in "$W"/$SCN/*/mismatch; do $RS human/scripts/gt_rank_v7.R "$d" > "$d/rank.log" 2>&1 || echo "rank failed: $d"; done
  W=$W SCN=$SCN $RS human/scripts/gt_v8_filter_sweep_compare.R
  ;;

*) echo "unknown STAGE '$STAGE'" >&2; exit 1 ;;
esac
