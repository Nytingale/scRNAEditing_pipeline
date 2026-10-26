#!/bin/bash
# mmonly_rerun_v8.sh   STAGE=manifest | (SGE array, CHUNK manifest rows per element)   (2026-10-15)
# Main call state with the nearby-mismatch window chosen on pass 1 alone (user 2026-10-15: the joint w1 x w2 search,
# scored after both passes, let the clustered-sites pass steer w1). Every NTC DOWNSTREAM_R / DOWNSTREAM_R_SCALLELE
# task of the v8 arm, sample, cell-type and split-B runs, and the 8 JACUSA2 replicate-mode downstream runs, is re-run
# from its own .command.sh with nextflow/bin/downstream_filter{,_scallele}.R replaced by the copies
# human/pm/scripts/downstream_filter{,_scallele}_mmonly.R (--clustered_filter 0 by default: w2 = 0, w1 searched).
# nextflow/bin/ and the published trees are not touched here (swap: mmonly_swap_v8.sh).
#   manifest  gt_eval_v8/mmonly/manifest.tsv  tier, work_dir, published_pre100, command
#   array     -> gt_eval_v8/mmonly/results/<published path relative to human/pm>/<prefix>_{step5_nomm.tsv.gz,summary.txt,
#                stepwise.tsv,metrics.tsv,run.log}
#   status    gt_eval_v8/mmonly/status/<tier>__<prefix>.tsv  tier, dest, exit, best_w1, best_w2, n_AG_step5, old_w1, old_w2
#
#   STAGE=manifest bash mmonly_rerun_v8.sh
#   n=$(( ($(wc -l < gt_eval_v8/mmonly/manifest.tsv) - 1 + 1) / 2 )); qsub -t 1-$n -tc 400 -v CHUNK=2 mmonly_rerun_v8.sh
#$ -S /bin/bash
#$ -N mmonly_v8
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 4
#$ -l h_vmem=6G
#$ -l h_rt=24:00:00
set -uo pipefail
export LC_ALL=C
PM=${SCRNA_DIR:-.}/human/pm
BIN=${SCRNA_DIR:-.}/nextflow/bin
P=$PM/gt_eval_v8/mmonly
M=${M:-$P/manifest.tsv}
STAGE=${STAGE:-run}

if [ "$STAGE" = "manifest" ]; then
  set -e
  mkdir -p "$P"
  printf "tier\twork_dir\tpublished_pre100\tcommand\n" > "$M"
  # arm / sample / celltype: the same 2,212 NTC tasks as the clustered-state replay
  tail -n +2 "$PM/gt_eval_v8/nomm/manifest.tsv" | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $2"/.command.sh"}' >> "$M"
  # split B
  r=$PM/results_gt_split_B_v8
  find "$r" -path '*NTC_ctrl*/04_filtered/*' -name '*_pre100.tsv.gz' | awk -F/ '{print $NF"\t"$0}' | sort > /tmp/mmonly8_pub_split.tsv
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
    $h["process"]~/^DOWNSTREAM_R/ && $h["tag"]~/(^|_)NTC_ctrl(__[^:]*)?:/ && $h["status"]~/COMPLETED|CACHED/ {print $h["hash"]}' \
    "$r/pipeline_trace.txt" | sort -u |
  while read -r h; do
    d=$(ls -d "$PM/work_gt_split_B_v8/$h"*)
    p=$(grep -oP -- '--output_pre100 \K\S+' "$d/.command.sh")
    printf "%s\t%s\n" "$p" "$d"
  done | sort > /tmp/mmonly8_task_split.tsv
  join -t $'\t' /tmp/mmonly8_task_split.tsv /tmp/mmonly8_pub_split.tsv | awk -F'\t' -v OFS='\t' '{print "split", $2, $3, $2"/.command.sh"}' >> "$M"
  n_t=$(wc -l < /tmp/mmonly8_task_split.tsv); n_p=$(wc -l < /tmp/mmonly8_pub_split.tsv); n_j=$(awk -F'\t' '$1=="split"' "$M" | wc -l)
  echo "split: tasks=$n_t published=$n_p joined=$n_j"
  [ "$n_t" -eq "$n_p" ] && [ "$n_j" -eq "$n_t" ] || { echo "manifest mismatch for split" >&2; exit 1; }
  # JACUSA2 replicate mode (replicates_jacusa2_v8.sh): run dir with the summed replicate input; command = arm task
  for B in recal recalq20 sampleRmDup sampleUmiDup; do for s in annotation caller; do
    tag="NTC_ctrl:jacusa2:${B}:${s}:snp_off"
    h=$(awk -F'\t' -v t="$tag" '$4 == "DOWNSTREAM_R" && $5 == t && $7 == "COMPLETED" {print $2}' "$PM/results_gt_v8/pipeline_trace.txt" | tail -1)
    dwd=$(ls -d "$PM/work_gt_v8/$h"*)
    d=$PM/gt_eval_v8/replicates/runs/jacusa2_rep/$B/downstream_$s
    [ -d "$d" ] || { echo "missing $d" >&2; exit 1; }
    pub=$PM/gt_eval_v8/replicates/results_jacusa2rep_arm/NTC_ctrl/04_filtered/jacusa2/$B/$s/snp_off/NTC_ctrl_jacusa2_${B}_${s}_snpoff_pre100.tsv.gz
    printf "jacrep\t%s\t%s\t%s\n" "$d" "$pub" "$dwd/.command.sh" >> "$M"
  done; done
  cut -f1 "$M" | tail -n +2 | sort | uniq -c
  echo "manifest: $M ($(($(wc -l < "$M") - 1)) tasks)"
  exit 0
fi

export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
CHUNK=${CHUNK:-2}
first=$(( (SGE_TASK_ID - 1) * CHUNK + 2 ))
last=$(( first + CHUNK - 1 ))
mkdir -p "$P/status" "$P/runs"
wins() { sed -n 's/.*Best window combination: w1=\([0-9]*\), w2=\([0-9]*\).*/\1 \2/p' "$1" 2>/dev/null | head -1; }
sed -n "${first},${last}p" "$M" | while IFS=$'\t' read -r tier wd pub cmd; do
  base=$(basename "$pub" _pre100.tsv.gz)
  st=$P/status/${tier}__${base}.tsv
  if [ -s "$st" ] && [ "$(cut -f3 "$st")" = 0 ]; then echo "skip (done): $tier $base"; continue; fi
  run=$P/runs/$tier/$base
  rm -rf "$run"; mkdir -p "$run"
  for e in "$wd"/*; do
    n=$(basename "$e")
    if [ -L "$e" ]; then cp -P "$e" "$run/"
    elif [ -d "$e" ]; then mkdir -p "$run/$n"; cp -P "$e"/* "$run/$n/"
    fi
  done
  rc=0
  grep -q -- '--mismatch_window_filter 1' "$cmd" || { echo "no --mismatch_window_filter 1 in $cmd" >&2; rc=2; }
  sed -e "s#${BIN}/downstream_filter_scallele\.R#${PM}/scripts/downstream_filter_scallele_mmonly.R#" \
      -e "s#${BIN}/downstream_filter\.R#${PM}/scripts/downstream_filter_mmonly.R#" "$cmd" > "$run/run.sh"
  grep -q "_mmonly\.R" "$run/run.sh" || { echo "filter script not substituted in $cmd" >&2; rc=3; }
  if [ "$rc" -eq 0 ]; then ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=$?; fi
  rel=${pub#"$PM"/}
  dest=$P/results/$(dirname "$rel")
  w1=NA; w2=NA; nag=NA
  read -r ow1 ow2 < <(wins "$(dirname "$pub")/${base}_summary.txt")
  if [ "$rc" -eq 0 ]; then
    mkdir -p "$dest"
    for f in "$run/${base}"_{summary.txt,stepwise.tsv,metrics.tsv}; do [ -e "$f" ] && cp "$f" "$dest/"; done
    cp "$run/run.log" "$dest/${base}_run.log"
    read -r w1 w2 < <(wins "$run/run.log")
    s5=$run/${base}_filtered_step5_variable_regions_removed.tsv.gz
    if [ -s "$s5" ]; then
      cp "$s5" "$dest/${base}_step5_nomm.tsv.gz"
      nag=$(zcat "$s5" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="varType")v=i;if($i=="site")c=i};next}$v=="AG"{s[$c]=1}END{print length(s)+0}')
    else : > "$dest/${base}_step5_nomm.tsv.gz"; nag=0; fi
    rm -rf "$run"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tier" "$dest" "$rc" "${w1:-NA}" "${w2:-NA}" "$nag" "${ow1:-NA}" "${ow2:-NA}" > "$st"
done
