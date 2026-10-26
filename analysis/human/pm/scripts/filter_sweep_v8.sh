#!/bin/bash
# filter_sweep_v8.sh   STAGE=manifest|union|apply|link|check|score|rank   (2026-10-19)
# Downstream-filter sweep (user 2026-10-19): one downstream filter at a time moved from the main values
# (min_cov 10, min_alt 3, min_edit_level 1%, multi-mismatch / dbSNP / JACUSA2-score removal on; settings.tsv) or
# turned off, on every NTC DOWNSTREAM_R task of the v8 arm, sample and 11-cell-type runs; each setting scored and
# ranked like the main v8 analysis (call state mismatch = step 5; centre step 4 = no window filter).
# Each task is replayed from its own .command.sh (gt_eval_v8/nomm/manifest.tsv), so caller, BAM type, strand mode
# and SNP blacklist are the task's own; human/scripts/downstream_filter_sweep.R does the parse + strand assignment
# once (union) and the filters + window search per setting (apply).
#   manifest  W/manifest.tsv: id, tier, work_dir, published_pre100, caller, prefix, rel_dir, mem (G, from the trace)
#             W/apply_manifest.tsv: id, settings (arm: one setting per row; sample / cell type: all per row)
#   union     array over manifest rows (M=<subset> CHUNK=n)  -> W/union/<tier>/<prefix>.rds, W/status/union/
#   apply     array over apply_manifest rows                 -> W/results/<setting>/<rel_dir>/, W/status/apply/
#   link      settings that cannot change a caller (settings.tsv callers) -> symlinks to the centre files
#   check     centre step 4 / step 5 vs the published tables (md5 of the uncompressed table) -> W/check_centre.tsv
#   score     array 1..3 x (settings + 1): gt_score_v7.R -> W/scoring/<setting>/mismatch/<tier>/ (window_off: centre
#             step 4); then rank: gt_rank_v7.R per setting, human/scripts/gt_v8_filter_sweep_compare.R
#   bash filter_sweep_v8.sh manifest; bash filter_sweep_v8.sh submit
#$ -S /bin/bash
#$ -N fsweep
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
W=$V8/filter_sweep
SELF=$PM/scripts/filter_sweep_v8.sh
SWEEP_R=$B/human/scripts/downstream_filter_sweep.R
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
STAGE=${STAGE:-${1:-}}
CHUNK=${CHUNK:-1}
TIERS=(arm sample celltype)
declare -A PUB=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
mapfile -t SETTINGS < <(awk -F'\t' 'NR>1{print $1}' "$W/settings.tsv")
SCN=${SCN:-scoring}

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
submit)
  # one union and one apply array per memory class; apply waits for every union array (SGE holds ignore exit codes:
  # check W/status before trusting link / score)
  cd "$B"; mkdir -p "$W/sub"
  # HOLD_UNION=<job ids>: union arrays already submitted; submit apply onwards only
  hu=${HOLD_UNION:-}; ha=""
  [ -z "$hu" ] && for m in $(awk -F'\t' 'NR>1{print $8}' "$W/manifest.tsv" | sort -nu); do
    f=$W/sub/union_${m}G.tsv; { head -1 "$W/manifest.tsv"; awk -F'\t' -v m="$m" 'NR>1 && $8==m' "$W/manifest.tsv"; } > "$f"
    n=$(( $(wc -l < "$f") - 1 ))
    j=$(qsub -terse -N fsw_u$m -t 1-$n -tc 300 -l h_vmem=${m}G -v STAGE=union,M=$f "$SELF" | cut -d. -f1); hu=$hu${hu:+,}$j
  done
  for m in $(awk -F'\t' 'NR>1{print $3}' "$W/apply_manifest.tsv" | sort -nu); do
    f=$W/sub/apply_${m}G.tsv; { head -1 "$W/apply_manifest.tsv"; awk -F'\t' -v m="$m" 'NR>1 && $3==m' "$W/apply_manifest.tsv"; } > "$f"
    n=$(( $(wc -l < "$f") - 1 ))
    j=$(qsub -terse -N fsw_a$m -t 1-$n -tc 300 -l h_vmem=${m}G -hold_jid "$hu" -v STAGE=apply,M=$f "$SELF" | cut -d. -f1); ha=$ha${ha:+,}$j
  done
  jl=$(qsub -terse -N fsw_link -l h_vmem=4G -hold_jid "$ha" -v STAGE=link "$SELF")
  jc=$(qsub -terse -N fsw_check -l h_vmem=4G -hold_jid "$jl" -v STAGE=check "$SELF")
  js=$(qsub -terse -N fsw_score -t 1-$(( 3 * (${#SETTINGS[@]} + 1) )) -l h_vmem=48G -l h_rt=12:00:00 -hold_jid "$jl" \
         -v STAGE=score "$SELF" | cut -d. -f1)
  jr=$(qsub -terse -N fsw_rank -l h_vmem=16G -hold_jid "$js" -v STAGE=rank "$SELF")
  echo "union $hu | apply $ha | link $jl | check $jc | score $js | rank $jr"
  ;;

manifest)
  set -e
  keep=$(grep -v '^$' "$V8/unit_cov/celltype_keep_min20.txt" | paste -sd'|')
  for t in "${TIERS[@]}"; do
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["process"]~/^DOWNSTREAM_R/ {print $h["hash"]"\t"$h["peak_rss"]}' \
      "$PM/${PUB[$t]}/pipeline_trace.txt"
  done | awk -F'\t' '{v=$2; g=0; if (v~/GB/) g=v+0; else if (v~/MB/) g=v/1024; else if (v~/KB/) g=v/1048576; print $1"\t"g}' \
    | sort -u > /tmp/fsweep_rss.tsv
  printf "id\ttier\twork_dir\tpublished_pre100\tcaller\tprefix\trel_dir\tmem\n" > "$W/manifest.tsv"
  awk -F'\t' -v keep="^(${keep})/" -v pm="$PM/" 'NR==FNR{r[$1]=$2; next} FNR==1{next}
    {rel = substr($3, length(pm) + 1); unit = rel; sub(/^[^/]+\//, "", unit)
     if ($1 == "celltype" && unit !~ keep) next
     split($2, p, "/"); h = p[length(p)-1] "/" substr(p[length(p)], 1, 6)
     if (!(h in r)) { print "no trace row for " $2 > "/dev/stderr"; exit 1 }
     split(rel, q, "/04_filtered/"); split(q[2], c, "/")
     pre = $3; sub(/^.*\//, "", pre); sub(/_pre100\.tsv\.gz$/, "", pre)
     d = rel; sub(/\/[^/]+$/, "", d)
     g = r[h] * 1.5 + 4; m = 8; while (m < g) m *= 2
     printf "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n", ++n, $1, $2, $3, c[1], pre, d, m}' \
    /tmp/fsweep_rss.tsv "$V8/nomm/manifest.tsv" >> "$W/manifest.tsv"
  awk -F'\t' 'NR>1{n[$2]++; m[$8]++} END{for (k in n) print "tier", k, n[k]; for (k in m) print "mem", k"G", m[k]}' "$W/manifest.tsv"
  # apply rows: arm = one setting per row, other tiers = every setting of the task in one row
  printf "id\tsettings\tmem\n" > "$W/apply_manifest.tsv"
  awk -F'\t' 'NR==FNR{if (FNR > 1) {s[++ns] = $1; c[$1] = $9}; next} FNR==1{next}
    {ok = ""
     for (i = 1; i <= ns; i++) if (c[s[i]] == "all" || index("," c[s[i]] ",", "," $5 ",")) {
       if ($2 == "arm") print $1"\t"s[i]"\t"$8; else ok = ok (ok == "" ? "" : ",") s[i] }
     if ($2 != "arm") print $1"\t"ok"\t"$8}' "$W/settings.tsv" "$W/manifest.tsv" >> "$W/apply_manifest.tsv"
  echo "apply rows: $(( $(wc -l < "$W/apply_manifest.tsv") - 1 ))"
  ;;

union)
  mkdir -p "$W/status/union" "$W/runs"
  rows "${M:-$W/manifest.tsv}" | while IFS=$'\t' read -r id tier wd pub caller prefix rel mem; do
    st=$W/status/union/$id.tsv
    if [ -s "$st" ] && [ "$(cut -f2 "$st")" = 0 ]; then echo "skip (done): $id"; continue; fi
    run=$W/runs/u$id; rm -rf "$run"; mkdir -p "$run" "$W/union/$tier"
    for e in "$wd"/*; do
      n=$(basename "$e")
      if [ -L "$e" ]; then cp -P "$e" "$run/"
      elif [ -d "$e" ]; then mkdir -p "$run/$n"; cp -P "$e"/* "$run/$n/"
      fi
    done
    u=$W/union/$tier/$prefix.rds; rm -f "$u"
    sweep_cmd "$wd" "--sweep_mode union --sweep_union $u" > "$run/run.sh"
    ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=$?
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
      sweep_cmd "$wd" "--sweep_mode apply --sweep_union $u --sweep_settings_file $W/settings.tsv --sweep_settings $settings --sweep_out $W/results/{setting}/$rel" > "$run/run.sh"
      ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=$?
    fi
    cp "$run/run.log" "$W/status/apply/$key.log"
    printf "%s\t%s\t%s\n" "$id" "$settings" "$rc" > "$st"
    [ $rc -eq 0 ] && rm -rf "$run"
  done
  ;;

link)
  # settings whose filter does not exist for a caller: the centre call set
  awk -F'\t' 'NR==FNR{if (FNR > 1) c[$1] = $9; next} FNR==1{next}
    {for (s in c) if (c[s] != "all" && !index("," c[s] ",", "," $5 ",")) print s"\t"$7"\t"$6}' \
    "$W/settings.tsv" "$W/manifest.tsv" | while IFS=$'\t' read -r s rel prefix; do
    d=$W/results/$s/$rel; mkdir -p "$d"
    for f in "$W/results/centre/$rel/$prefix"{_step4_cov_edit_filtered.tsv.gz,_step5_nomm.tsv.gz,.done}; do
      [ -e "$f" ] && ln -sfn "$f" "$d/"
    done
  done
  echo "linked"
  ;;

check)
  printf "id\ttier\tprefix\tstep\tpublished\tsweep\tidentical\n" > "$W/check_centre.tsv"
  awk -F'\t' 'NR>1' "$W/manifest.tsv" | while IFS=$'\t' read -r id tier wd pub caller prefix rel mem; do
    for s in step4_cov_edit_filtered step5_nomm; do
      p=$PM/$rel/${prefix}_$s.tsv.gz; q=$W/results/centre/$rel/${prefix}_$s.tsv.gz
      mp=absent; mq=absent
      [ -e "$p" ] && mp=$(zcat "$p" | md5sum | cut -c1-32)
      [ -e "$q" ] && mq=$(zcat "$q" | md5sum | cut -c1-32)
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$id" "$tier" "$prefix" "$s" "$mp" "$mq" "$([ "$mp" = "$mq" ] && echo TRUE || echo FALSE)"
    done
  done >> "$W/check_centre.tsv"
  awk -F'\t' 'NR>1{n[$4" "$7]++} END{for (k in n) print k, n[k]}' "$W/check_centre.tsv"
  ;;

score)
  i=$(( SGE_TASK_ID - 1 )); tier=${TIERS[$(( i % 3 ))]}; k=$(( i / 3 ))
  if [ $k -lt ${#SETTINGS[@]} ]; then s=${SETTINGS[$k]}; src=$s; re='_step5_nomm\.tsv\.gz$'
  else s=window_off; src=centre; re='_step4_cov_edit_filtered\.tsv\.gz$'; fi
  # -v SCN=scoring_pw,GT_UNIT_RULE=gate,GT_UNIT_KEEP=$V8/unit_cov/celltype_keep_pw.txt: gate-matched recall denominator
  echo "setting=$s tier=$tier SCN=$SCN GT_UNIT_RULE=${GT_UNIT_RULE:-cov}"
  cd "$B"
  GT_UNIT_KEEP=${GT_UNIT_KEEP:-$V8/unit_cov/celltype_keep_min20.txt} GT_V7_TRUTH=$V8/truth/gt_v7_truth_sites.tsv.gz \
  GT_UNIT_COV=$V8/unit_cov/ko_truth_unit_cov.tsv.gz GT_UNIT_RULE=${GT_UNIT_RULE:-cov} \
    $RS human/scripts/gt_score_v7.R "$tier" "$W/results/$src/${PUB[$tier]}" "$W/$SCN/$s/mismatch/$tier" "$re"
  ;;

rank)
  cd "$B"
  for d in "$W"/$SCN/*/mismatch; do $RS human/scripts/gt_rank_v7.R "$d" > "$d/rank.log" 2>&1 || echo "rank failed: $d"; done
  SCN=$SCN $RS human/scripts/gt_v8_filter_sweep_compare.R
  ;;

*) echo "unknown STAGE '$STAGE'" >&2; exit 1 ;;
esac
