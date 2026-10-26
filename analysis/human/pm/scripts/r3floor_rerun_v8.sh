#!/bin/bash
# r3floor_rerun_v8.sh   STAGE=manifest | run (SGE array, M=<manifest subset>) | check | swap   (2026-10-21)
# REDItools3 minimum coverage (user 2026-10-21: the downstream filter applies >= 10 reads, >= 3 variant reads and
# >= 1% editing to every caller). The REDItools3 parse of downstream_filter{,_v8,_mmonly}.R had no min_cov cut, so
# REDItools3 calls with < 10 reads reached every table; the cut is now in the parse (backups:
# nextflow/_superseded/2026-10-21_reditools3_min_cov/, human/pm/scripts/_superseded_2026-10-21_reditools3_min_cov/).
# Every REDItools3 NTC task of gt_eval_v8/mmonly/manifest.tsv (arm, sample, cell type, split B) is replayed from its
# own .command.sh in each call state the scorers read:
#   mm         downstream_filter_mmonly.R            -> <prefix>_step4_cov_edit_filtered.tsv.gz, <prefix>_step5_nomm.tsv.gz
#   both       nextflow/bin/downstream_filter.R       -> <prefix>_{pre100.tsv.gz,filtered.txt.gz,summary.txt,stepwise.tsv,
#                                                        metrics.tsv}
#   clustered  as both, --mismatch_window_filter 0   -> the same files in results_gt{,_persample,_celltype}_v8_nomm/
#              (arm / sample / cell type only, as nomm_rerun_v8.sh)
#   manifest  P/manifest.tsv  tier, work_dir, published_pre100, command, mem (G: 2 x trace peak_rss, 16-128)
#             P/sub/run_<mem>G.tsv: one array per memory class
#   run       -> P/results/<published path relative to human/pm>/..., P/status/<tier>__<prefix>.tsv  (exit per state)
#   check     the new step-4 table must equal the published one restricted to cov >= 10 -> P/check_step4.tsv
#   swap      moves each replaced published file to gt_eval_v8/_superseded/reditools3_nofloor_2026-10-21/<same path>
#             and copies the new one in (exit 100 if any row is not done or check failed)
#   STAGE=manifest bash r3floor_rerun_v8.sh; STAGE=submit bash r3floor_rerun_v8.sh
#$ -S /bin/bash
#$ -N r3floor_v8
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
set -uo pipefail
export LC_ALL=C
PM=${SCRNA_DIR:-.}/human/pm
BIN=${SCRNA_DIR:-.}/nextflow/bin
SELF=$PM/scripts/r3floor_rerun_v8.sh
P=$PM/gt_eval_v8/r3floor
SUP=$PM/gt_eval_v8/_superseded/reditools3_nofloor_2026-10-21
M=${M:-$P/manifest.tsv}
STAGE=${STAGE:-run}
declare -A TRACE=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8
                  [split]=results_gt_split_B_v8)

# published file -> staged file (same path under P/results)
stage_of() { echo "$P/results/${1#"$PM"/}"; }
# published pre100 -> the clustered-state (results_*_v8_nomm) path
nomm_of() { local rel=${1#"$PM"/}; echo "$PM/${rel%%/*}_nomm/${rel#*/}"; }
# every published file this row replaces
row_files() {  # <tier> <pub>
  local d=$(dirname "$2") b=$(basename "$2" _pre100.tsv.gz) x
  printf "%s\n" "$d/${b}_step4_cov_edit_filtered.tsv.gz" "$d/${b}_step5_nomm.tsv.gz"
  for x in pre100.tsv.gz filtered.txt.gz summary.txt stepwise.tsv metrics.tsv; do
    echo "$d/${b}_$x"; [ "$1" != split ] && echo "$(dirname "$(nomm_of "$2")")/${b}_$x"
  done
}

case "$STAGE" in
manifest)
  set -e; mkdir -p "$P/sub"
  printf "tier\twork_dir\tpublished_pre100\tcommand\tmem\n" > "$M"
  awk -F'\t' 'NR>1 && $3 ~ /\/reditools3\//' "$PM/gt_eval_v8/mmonly/manifest.tsv" |
  while IFS=$'\t' read -r tier wd pub cmd; do
    h=$(basename "$(dirname "$wd")")/$(basename "$wd" | cut -c1-6)
    rss=$(awk -F'\t' -v h="$h" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["hash"]==h{print $c["peak_rss"]}' \
          "$PM/${TRACE[$tier]}/pipeline_trace.txt" | tail -1)
    mem=$(echo "$rss" | awk '{g=$1; if($2=="MB")g/=1024; if($2=="KB")g/=1048576; m=16; while(m<2*g && m<128)m*=2; print m}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$tier" "$wd" "$pub" "$cmd" "$mem" >> "$M"
  done
  for m in $(awk -F'\t' 'NR>1{print $5}' "$M" | sort -nu); do
    { head -1 "$M"; awk -F'\t' -v m="$m" 'NR>1 && $5==m' "$M"; } > "$P/sub/run_${m}G.tsv"
  done
  cut -f1,5 "$M" | tail -n +2 | sort | uniq -c
  echo "manifest: $M ($(($(wc -l < "$M") - 1)) tasks)"
  ;;
submit)
  cd "$PM/../.."; ids=""
  for f in "$P"/sub/run_*G.tsv; do
    m=$(basename "$f" .tsv); m=${m#run_}; m=${m%G}; n=$(( $(wc -l < "$f") - 1 ))
    j=$(qsub -terse -N r3f_$m -t 1-$n -tc 200 -l h_vmem=${m}G -v STAGE=run,M=$f "$SELF" | cut -d. -f1); ids=$ids${ids:+,}$j
  done
  c=$(qsub -terse -N r3f_check -l h_vmem=8G -l h_rt=4:00:00 -hold_jid "$ids" -v STAGE=check "$SELF")
  echo "run $ids | check $c"
  ;;
run)
  export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
  mkdir -p "$P/status" "$P/runs"
  sed -n "$(( SGE_TASK_ID + 1 ))p" "$M" | while IFS=$'\t' read -r tier wd pub cmd mem; do
    base=$(basename "$pub" _pre100.tsv.gz); d=$(dirname "$pub")
    st=$P/status/${tier}__${base}.tsv
    states=(mm both); [ "$tier" != split ] && states+=(clustered)
    declare -A rc=()
    grep -q -- '--mismatch_window_filter 1' "$cmd" || { echo "no --mismatch_window_filter 1 in $cmd" >&2; exit 2; }
    for s in "${states[@]}"; do
      run=$P/runs/$tier/$base/$s; rm -rf "$run"; mkdir -p "$run"
      for e in "$wd"/*; do
        n=$(basename "$e")
        if [ -L "$e" ]; then cp -P "$e" "$run/"
        elif [ -d "$e" ]; then mkdir -p "$run/$n"; cp -P "$e"/* "$run/$n/"
        elif [[ "$n" == *_reditools3.txt.gz ]]; then ln -s "$e" "$run/$n"
        fi
      done
      case $s in
        mm)        sed "s#${BIN}/downstream_filter\.R#${PM}/scripts/downstream_filter_mmonly.R#" "$cmd" > "$run/run.sh" ;;
        both)      cp "$cmd" "$run/run.sh" ;;
        clustered) sed 's/--mismatch_window_filter 1/--mismatch_window_filter 0/' "$cmd" > "$run/run.sh" ;;
      esac
      [ "$s" = mm ] && ! grep -q "_mmonly\.R" "$run/run.sh" && { rc[$s]=3; continue; }
      ( cd "$run" && bash run.sh > run.log 2>&1 ); rc[$s]=$?
      [ "${rc[$s]}" -eq 0 ] || continue
      case $s in
        mm)
          dest=$(dirname "$(stage_of "$pub")"); mkdir -p "$dest"
          for x in step4_cov_edit_filtered step5_variable_regions_removed; do
            f=$run/${base}_filtered_${x}.tsv.gz; o=${x/step5_variable_regions_removed/step5_nomm}
            if [ -s "$f" ]; then cp "$f" "$dest/${base}_${o}.tsv.gz"; else : > "$dest/${base}_${o}.tsv.gz"; fi
          done
          cp "$run/run.log" "$dest/${base}_mm_run.log" ;;
        both|clustered)
          if [ "$s" = both ]; then dest=$(dirname "$(stage_of "$pub")"); else dest=$(dirname "$(stage_of "$(nomm_of "$pub")")"); fi
          mkdir -p "$dest"
          for x in pre100.tsv.gz filtered.txt.gz summary.txt stepwise.tsv metrics.tsv; do
            [ -e "$run/${base}_$x" ] && cp "$run/${base}_$x" "$dest/"
          done
          cp "$run/run.log" "$dest/${base}_${s}_run.log" ;;
      esac
      rm -rf "$run"
    done
    printf "%s\t%s\t%s\t%s\t%s\n" "$tier" "$base" "${rc[mm]:-NA}" "${rc[both]:-NA}" "${rc[clustered]:-NA}" > "$st"
  done
  ;;
check)
  # step 4 is upstream of every window pass: new = published restricted to cov >= 10 (site, varType, strand)
  export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
  printf "tier\tprefix\tn_old\tn_old_cov10\tn_new\tidentical\n" > "$P/check_step4.tsv"
  tail -n +2 "$M" | while IFS=$'\t' read -r tier wd pub cmd mem; do
    b=$(basename "$pub" _pre100.tsv.gz); old=$(dirname "$pub")/${b}_step4_cov_edit_filtered.tsv.gz
    new=$(dirname "$(stage_of "$pub")")/${b}_step4_cov_edit_filtered.tsv.gz
    [ -e "$new" ] || { printf "%s\t%s\tNA\tNA\tNA\tmissing\n" "$tier" "$b" >> "$P/check_step4.tsv"; continue; }
    key='NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} {print $h["site"]"\t"$h["varType"]"\t"$h["strand"]"\t"$h["cov"]}'
    zcat -f "$old" | awk -F'\t' "$key" | sort > /tmp/r3f_old_$$
    zcat -f "$new" | awk -F'\t' "$key" | sort > /tmp/r3f_new_$$
    awk -F'\t' '$4 >= 10' /tmp/r3f_old_$$ > /tmp/r3f_old10_$$
    same=$(cmp -s /tmp/r3f_old10_$$ /tmp/r3f_new_$$ && echo yes || echo no)
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tier" "$b" "$(wc -l < /tmp/r3f_old_$$)" "$(wc -l < /tmp/r3f_old10_$$)" \
      "$(wc -l < /tmp/r3f_new_$$)" "$same" >> "$P/check_step4.tsv"
  done
  rm -f /tmp/r3f_*_$$
  awk -F'\t' 'NR>1{c[$6]++} END{for(k in c) print k, c[k]}' "$P/check_step4.tsv"
  ;;
swap)
  set -e; bad=0; n=0
  while IFS=$'\t' read -r tier wd pub cmd mem; do
    b=$(basename "$pub" _pre100.tsv.gz); st=$P/status/${tier}__${b}.tsv
    want=$([ "$tier" = split ] && echo "0 0 NA" || echo "0 0 0")
    [ -s "$st" ] && [ "$(cut -f3-5 "$st" | tr '\t' ' ')" = "$want" ] || { echo "NOT DONE: $tier $b" >&2; bad=$((bad + 1)); }
  done < <(tail -n +2 "$M")
  nc=$(awk -F'\t' 'NR>1 && $6!="yes"' "$P/check_step4.tsv" | wc -l)
  [ "$bad" -eq 0 ] && [ "$nc" -eq 0 ] || { echo "ERROR: $bad rows not done, $nc step-4 checks failed; nothing swapped" >&2; exit 100; }
  while IFS=$'\t' read -r tier wd pub cmd mem; do
    while read -r f; do
      new=$(stage_of "$f"); [ -e "$new" ] || continue
      rel=${f#"$PM"/}
      if [ -e "$f" ] && [ ! -e "$SUP/$rel" ]; then mkdir -p "$SUP/$(dirname "$rel")"; mv "$f" "$SUP/$rel"; fi
      cp "$new" "$f"; n=$((n + 1))
      # mmonly_swap_v8.sh copies gt_eval_v8/mmonly/results over the published step-5 tables
      mm=$PM/gt_eval_v8/mmonly/results/$rel
      if [[ "$f" == *_step5_nomm.tsv.gz ]] && [ -e "$mm" ]; then
        s=$SUP/mmonly_results/$rel; [ -e "$s" ] || { mkdir -p "$(dirname "$s")"; mv "$mm" "$s"; }
        cp "$new" "$mm"
      fi
    done < <(row_files "$tier" "$pub")
  done < <(tail -n +2 "$M")
  echo "swapped $n files (old copies: $SUP)"
  ;;
esac
