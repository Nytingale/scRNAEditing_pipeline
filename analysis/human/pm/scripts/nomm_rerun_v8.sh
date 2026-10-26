#!/bin/bash
# nomm_rerun_v8.sh   STAGE=manifest | (SGE array, CHUNK manifest rows per element)   (2026-10-15)
# "clustered" call state for the v8 caller runs (as nomm_rerun_v6.sh for v6): every NTC DOWNSTREAM_R /
# DOWNSTREAM_R_SCALLELE task of the arm, sample and cell-type v8 runs is re-run from its own .command.sh with
# --mismatch_window_filter 0 (pass 1 off, w1 = 0; only the clustered-sites pass w2 is searched). The published
# results_gt_*_v8/ trees are not touched.
#   manifest  gt_eval_v8/nomm/manifest.tsv  tier, work_dir, published_pre100 (trace hash -> work dir; must join 1:1)
#   array     -> human/pm/results_gt{,_persample,_celltype}_v8_nomm/<unit>/04_filtered/<caller>/<bam>/<strand>/<snp>/
#                <prefix>_{pre100.tsv.gz,filtered.txt.gz,summary.txt,stepwise.tsv,metrics.tsv,run.log}
#   status    gt_eval_v8/nomm/status/<tier>__<prefix>.tsv  tier, dest, exit, best_w1, best_w2, n_AG_pre100
#
#   STAGE=manifest bash nomm_rerun_v8.sh
#   n=$(( ($(wc -l < gt_eval_v8/nomm/manifest.tsv) - 1 + 9) / 10 )); qsub -t 1-$n -tc 80 -v CHUNK=10 nomm_rerun_v8.sh
#$ -S /bin/bash
#$ -N nomm_v8
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
set -uo pipefail
export LC_ALL=C
PM=${SCRNA_DIR:-.}/human/pm
P=$PM/gt_eval_v8/nomm
# M override: a subset manifest (same header) for retries, e.g. M=$P/manifest_retry.tsv CHUNK=1 -l h_vmem=64G
# (16G was too little for some tasks: without the nearby-mismatch pass, strand assignment sees 2-12 M sites)
M=${M:-$P/manifest.tsv}
STAGE=${STAGE:-run}

if [ "$STAGE" = "manifest" ]; then
  set -e
  mkdir -p "$P"
  printf "tier\twork_dir\tpublished_pre100\n" > "$M"
  declare -A RES=([arm]=results_gt_v8 [sample]=results_gt_persample_v8 [celltype]=results_gt_celltype_v8)
  declare -A WORK=([arm]=work_gt_v8 [sample]=work_gt_persample_v8 [celltype]=work_gt_celltype_v8)
  for t in arm sample celltype; do
    r=$PM/${RES[$t]}
    find "$r" -path '*NTC_ctrl*/04_filtered/*' -name '*_pre100.tsv.gz' | awk -F/ '{print $NF"\t"$0}' | sort > /tmp/nomm8_pub_$t.tsv
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
      $h["process"]~/^DOWNSTREAM_R/ && $h["tag"]~/(^|_)NTC_ctrl(__[^:]*)?:/ && $h["status"]~/COMPLETED|CACHED/ {print $h["hash"]}' \
      "$r/pipeline_trace.txt" | sort -u |
    while read -r h; do
      d=$(ls -d "$PM/${WORK[$t]}/$h"*)
      p=$(grep -oP -- '--output_pre100 \K\S+' "$d/.command.sh")
      printf "%s\t%s\n" "$p" "$d"
    done | sort > /tmp/nomm8_task_$t.tsv
    n_t=$(wc -l < /tmp/nomm8_task_$t.tsv); n_p=$(wc -l < /tmp/nomm8_pub_$t.tsv)
    join -t $'\t' /tmp/nomm8_task_$t.tsv /tmp/nomm8_pub_$t.tsv | awk -F'\t' -v t=$t '{print t"\t"$2"\t"$3}' >> "$M"
    n_j=$(awk -F'\t' -v t=$t '$1==t' "$M" | wc -l)
    echo "$t: tasks=$n_t published=$n_p joined=$n_j"
    [ "$n_t" -eq "$n_p" ] && [ "$n_j" -eq "$n_t" ] || { echo "manifest mismatch for $t" >&2; exit 1; }
  done
  echo "manifest: $M ($(($(wc -l < "$M") - 1)) tasks)"
  exit 0
fi

export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
CHUNK=${CHUNK:-10}
first=$(( (SGE_TASK_ID - 1) * CHUNK + 2 ))
last=$(( first + CHUNK - 1 ))
mkdir -p "$P/status" "$P/runs"
sed -n "${first},${last}p" "$M" | while IFS=$'\t' read -r tier wd pub; do
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
  grep -q -- '--mismatch_window_filter 1' "$wd/.command.sh" || { echo "no --mismatch_window_filter 1 in $wd" >&2; rc=2; }
  sed 's/--mismatch_window_filter 1/--mismatch_window_filter 0/' "$wd/.command.sh" > "$run/run.sh"
  ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=${rc:-$?}
  rel=${pub#"$PM"/}
  res=${rel%%/*}
  dest=$PM/${res}_nomm/$(dirname "${rel#*/}")
  w1=NA; w2=NA; nag=NA
  if [ "$rc" -eq 0 ]; then
    mkdir -p "$dest"
    for f in "$run/${base}"_{pre100.tsv.gz,filtered.txt.gz,summary.txt,stepwise.tsv,metrics.tsv}; do
      [ -e "$f" ] && cp "$f" "$dest/"
    done
    cp "$run/run.log" "$dest/${base}_run.log"
    read -r w1 w2 < <(sed -n 's/.*Best window combination: w1=\([0-9]*\), w2=\([0-9]*\).*/\1 \2/p' "$run/run.log" | head -1)
    if [ -s "$dest/${base}_pre100.tsv.gz" ]; then
      nag=$(zcat "$dest/${base}_pre100.tsv.gz" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="varType")v=i;if($i=="site")c=i};next}$v=="AG"{s[$c]=1}END{print length(s)+0}')
    else : > "$dest/${base}_pre100.tsv.gz"; nag=0; fi
    rm -rf "$run"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tier" "$dest" "$rc" "${w1:-NA}" "${w2:-NA}" "$nag" > "$P/status/${tier}__${base}.tsv"
  unset rc
done
