#!/bin/bash
# ding_3prime_bias_tso_chain.sh   (2026-10-31) 3' bias figure (response Fig 16) on the TSO-clipped Ding BAMs.
#   STAGE=submit bash human/scripts/ding_3prime_bias_tso_chain.sh   (login node)
#     ding_ss2_pbmc1_merge.sh -> ding_dedup_read_counts.sh -> this script STAGE=sub (K from the new 10x umiDup count)
#   STAGE=umimd bash ...   (login node; 2026-11-14, 10x UMI-aware MarkDuplicates): ding_dedup_read_counts.sh -> STAGE=sub
#   STAGE=sub (held job): ding_ss2_dedup_subsample.sh -v K -> ding_3prime_bias_pileup_sub.sh (92) ->
#     ding_3prime_bias_sub_count.sh (3) -> DING_3prime_bias_v5_source_data.R + DING_3prime_bias_v5_plot.R
#   K = the 10x UMI-aware MarkDuplicates BAM (02_bams_umimd/PBMC1.bam) primary mapped count made even (was umi_tools)
#$ -S /bin/bash
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -j y
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -N ding_tso_chain
#$ -l h_vmem=2G
#$ -l h_rt=1:00:00
set -euo pipefail
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
B=${SCRNA_DIR:-.}
H=$B/human
S=$H/scripts
SELF=$S/ding_3prime_bias_tso_chain.sh
MODE=${STAGE:?STAGE=submit|umimd|sub}; unset STAGE
case "$MODE" in
submit)
  mkdir -p $H/jzlxv3/three_prime_bias/pileup
  j1=$(qsub -terse $S/ding_ss2_pbmc1_merge.sh)
  j2=$(qsub -terse -hold_jid $j1 $S/ding_dedup_read_counts.sh)
  j3=$(qsub -terse -hold_jid $j2 -v STAGE=sub $SELF)
  echo "pbmc1 merge $j1 | dedup counts $j2 | chain $j3"
  ;;
umimd)
  j2=$(qsub -terse $S/ding_dedup_read_counts.sh)
  j3=$(qsub -terse -hold_jid $j2 -v STAGE=sub $SELF)
  echo "dedup counts $j2 | chain $j3"
  ;;
sub)
  n=$(awk -F'\t' -v b=$H/jzlxv3/02_bams_umimd/PBMC1.bam '$1 == b {print $2}' $H/jzlxv3/three_prime_bias/read_counts.tsv)
  [ -n "$n" ] || { echo "ERROR: no 10x UMI-aware MarkDuplicates count" >&2; exit 1; }
  K=$(( n / 2 * 2 )); echo "10x UMI-aware MarkDuplicates primary mapped $n -> K = $K"
  js=$(qsub -terse -t 1-4 -v K=$K $S/ding_ss2_dedup_subsample.sh | cut -d. -f1)
  jp=$(qsub -terse -t 1-92 -hold_jid $js $S/ding_3prime_bias_pileup_sub.sh | cut -d. -f1)
  jc=$(qsub -terse -t 1-3 -hold_jid $jp $S/ding_3prime_bias_sub_count.sh | cut -d. -f1)
  P=$B/publication/three_prime_bias/scripts
  jf=$(qsub -terse -q slice -N ding_v5_fig -hold_jid $jc -o $B/logs/ -j y -l h_vmem=16G -l h_rt=2:00:00 -b y \
    "cd $B; ${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript $P/DING_3prime_bias_v5_source_data.R && ${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript $P/DING_3prime_bias_v5_plot.R")
  echo "subsample $js | pileup $jp | count $jc | figure $jf"
  ;;
*) echo "ERROR: STAGE=submit|umimd|sub" >&2; exit 1 ;;
esac
