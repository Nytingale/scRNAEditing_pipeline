#!/bin/bash
# prewindow_rerun_v6.sh   (SGE array; one element = CHUNK manifest rows)
# Re-runs each NTC downstream task's own .command.sh (from its Nextflow work dir) in a fresh
# directory under gt_eval_v6/prewindow/runs/, adding --save_intermediates 4 so the site table
# after the coverage/editing filter and before the mismatch-window filters (step 4) is saved.
# Only the staged input symlinks are copied, so nothing is written into the work dirs.
# The regenerated pre100 must equal the published one (same code and inputs), else the row fails.
# Output: gt_eval_v6/prewindow/results_<tier>/<unit>/04_filtered/<caller>/<bam>/<strand>/<snp>/
#         <prefix>_step4_cov_edit_filtered.tsv.gz   (empty file when step 4 has no sites)
# Status: gt_eval_v6/prewindow/status/<task>.tsv
#   tier, run_dir, exit, pre100_identical, n_step4, n_step4_AG, n_step4_AG_in_stepwise
#
#   n=$(( ($(wc -l < manifest.tsv) - 1 + CHUNK - 1) / CHUNK ))
#   qsub -t 1-$n -tc 60 -v CHUNK=10 prewindow_rerun_v6.sh
#$ -S /bin/bash
#$ -N prewindow_v6
#$ -q slice
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00
set -uo pipefail
PM=${SCRNA_DIR:-.}/human/pm
P=$PM/gt_eval_v6/prewindow
export PATH=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin:$PATH
CHUNK=${CHUNK:-10}
MANIFEST=${MANIFEST:-$P/manifest.tsv}   # same columns; e.g. a retry subset
first=$(( (SGE_TASK_ID - 1) * CHUNK + 2 ))
last=$(( first + CHUNK - 1 ))
mkdir -p "$P/status"
rc_all=0
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
  awk '{print} /^Rscript .*downstream_filter/{print "    --save_intermediates 4 \\"}' "$wd/.command.sh" > "$run/run.sh"
  ( cd "$run" && bash run.sh > run.log 2>&1 ); rc=$?
  same=NA; n4=NA; ag4=NA; ag4_stepwise=NA
  if [ $rc -eq 0 ]; then
    [ "$(zcat "$run/${base}_pre100.tsv.gz" | md5sum)" = "$(zcat "$pub" | md5sum)" ] && same=TRUE || same=FALSE
    rel=${pub#"$PM"/results_*/}
    dest=$P/results_$tier/$(dirname "$rel")
    mkdir -p "$dest"
    s4=$run/${base}_filtered_step4_cov_edit_filtered.tsv.gz
    out4=$dest/${base}_step4_cov_edit_filtered.tsv.gz
    ag4_stepwise=$(awk -F'\t' '$1=="AG" && $5=="Cov_and_edit_filtered"{print $2}' "$run/${base}_stepwise.tsv")
    if [ -s "$s4" ]; then
      cp "$s4" "$out4"
      read -r n4 ag4 < <(zcat "$s4" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="varType")v=i;next}{n++; if($v=="AG")a++}END{print n+0, a+0}')
    else : > "$out4"; n4=0; ag4=0; fi
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tier" "$run" "$rc" "$same" "$n4" "$ag4" "$ag4_stepwise" > "$P/status/${tier}__${base}.tsv"
done
