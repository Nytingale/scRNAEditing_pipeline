#!/bin/bash
# nomm_rerun_v6.sh   (SGE array; one element = CHUNK manifest rows; 2026-10-01)
# Separate result set without the nearby-mismatch window filter (pass 1, w1 fixed at 0; only the
# clustered-sites pass w2 is searched). Re-runs each published NTC downstream task's own
# .command.sh (manifest = gt_eval_v6/prewindow/manifest.tsv: tier, work_dir, published_pre100)
# in a scratch dir with --mismatch_window_filter 0, then copies the outputs to the same relative
# path under results_*_v6_nomm/ (the published results_*_v6/ trees are not touched).
# Output: human/pm/results_gt{,_persample,_celltype}_v6_nomm/<unit>/04_filtered/<caller>/<bam>/<strand>/<snp>/
#         <prefix>_{pre100.tsv.gz,filtered.txt.gz,summary.txt,stepwise.tsv,metrics.tsv,run.log}
# Status: gt_eval_v6/nomm/status/<tier>__<prefix>.tsv  tier, dest, exit, best_w1, best_w2, n_AG_pre100
#
#   n=$(( ($(wc -l < manifest.tsv) - 1 + CHUNK - 1) / CHUNK ))
#   qsub -t 1-$n -tc 60 -v CHUNK=10 nomm_rerun_v6.sh
#$ -S /bin/bash
#$ -N nomm_v6
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
set -uo pipefail
PM=${SCRNA_DIR:-.}/human/pm
P=$PM/gt_eval_v6/nomm
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
CHUNK=${CHUNK:-10}
MANIFEST=${MANIFEST:-$PM/gt_eval_v6/prewindow/manifest.tsv}
first=$(( (SGE_TASK_ID - 1) * CHUNK + 2 ))
last=$(( first + CHUNK - 1 ))
mkdir -p "$P/status" "$P/runs"
sed -n "${first},${last}p" "$MANIFEST" | while IFS=$'\t' read -r tier wd pub; do
  base=$(basename "$pub" _pre100.tsv.gz)
  run=$P/runs/$tier/$base
  rm -rf "$run"; mkdir -p "$run"
  for e in "$wd"/*; do
    n=$(basename "$e")
    if [ -L "$e" ]; then cp -P "$e" "$run/"
    elif [ -d "$e" ]; then mkdir -p "$run/$n"; cp -P "$e"/* "$run/$n/"
    fi
  done
  awk '{print} /^Rscript .*downstream_filter/{print "    --mismatch_window_filter 0 \\"}' "$wd/.command.sh" > "$run/run.sh"
  ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=$?
  rel=${pub#"$PM"/}
  res=${rel%%/*}
  dest=$PM/${res}_nomm/$(dirname "${rel#*/}")
  w1=NA; w2=NA; nag=NA
  if [ $rc -eq 0 ]; then
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
done
