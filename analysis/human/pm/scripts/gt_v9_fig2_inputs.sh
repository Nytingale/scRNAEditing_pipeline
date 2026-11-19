#!/bin/bash
# gt_v9_fig2_inputs.sh STAGE=spec|specgather|kopile|ko|redi|keep|submit   (2026-11-16)
# Response figure 2 (GTDEF_v8_universe_celltype_robustness) inputs on the seed-8787 umi_tools BAMs and the v9 truth.
#   spec        array over gt_eval_v6/scan/regions.tsv = run_gt_truth_rebuild_v6.sh spectrumtask (merged NTC / ADAR1 KO,
#               every aligned base, BQ/MQ 30 and 20) -> gt_eval_v6/spectrum/chunks/
#   specgather  run_gt_truth_rebuild_v6.sh spectrumgather -> gt_eval_v6/spectrum/mismatch_spectrum.tsv
#   kopile      merged ADAR1 KO pileup at the candidates, BQ/MQ 30 (run_gt_truth_rebuild_v6.sh truthpile for that unit)
#               -> gt_eval_v6/pileups/bq30/pileup_merged_ADAR1_KO.tsv (panel d)
#   ko / redi   headline scorings (gate_cov, FP = A-to-G) of all 23 NTC cell types (GT_UNIT_KEEP empty), for their
#               covered positives only -> gt_eval_v9/celltype_all23/{ko,redi}/
#   keep        gt_v9_celltype_covered_positives.R -> gt_eval_v9/unit_cov/celltype_{covered_positives,keep}_pwcov.*
#   submit      login node: qsub only
#$ -S /bin/bash
#$ -N gt_v9_fig2
# -q <queue>: add as a "#$" line if your cluster needs one
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/human/pm/gt_eval_v9/logs/
#$ -j y
#$ -l h_vmem=8G
#$ -l h_rt=96:00:00
set -euo pipefail
B=${SCRNA_DIR:-.}
PM=$B/human/pm
V6=$PM/gt_eval_v6
V8=$PM/gt_eval_v8
V9=$PM/gt_eval_v9
SELF=$PM/scripts/gt_v9_fig2_inputs.sh
REB=$PM/scripts/run_gt_truth_rebuild_v6.sh
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
C=$V9/celltype_all23
STAGE=${STAGE:-${1:-}}
cd "$B"

case $STAGE in
spec)       bash "$REB" spectrumtask ;;
specgather) bash "$REB" spectrumgather ;;
kopile)
  export PATH="${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:${GT_PILEUP_BIN:-/opt/conda/envs/gt_pileup/bin}:/usr/local/bin:/usr/bin:/bin:${PATH}"
  bash "$B/human/scripts/gt_site_pileup.sh" --bam "$PM/02_gt_bams_umi_v6/merged_ADAR1_KO.bam" --bed "$V6/candidates/candidates_all.bed" \
    --fasta "$B/human/ref_annots/hg38.analysisSet.fa" --out "$V6/pileups/bq30/pileup_merged_ADAR1_KO.tsv" \
    --libtype forward --min-bq 30 --min-mq 30
  ;;
ko)
  GT_UNIT_RULE=gate_cov GT_FP_RULE=ag GT_V9_DIR=$V9 GT_UNIT_KEEP= GT_V7_TRUTH=$V9/truth/gt_v9_truth_sites.tsv.gz \
    GT_UNIT_COV=$V9/unit_cov/ko_truth_unit_cov.tsv.gz \
    $RS human/scripts/gt_score_v7.R celltype "$PM/results_gt_celltype_v8" "$C/ko/celltype" '_step5_nomm\.tsv\.gz$'
  ;;
redi)
  R=$V9/rediportal_v9
  GT_V9_DIR=$V9 GT_REDI_UDIR=$R/universe GT_REDI_REVCOV=$R/unit_cov_rev GT_UNIT_KEEP= \
    GT_REDI_EXCLUDE_BED=$PM/pgp1_snps/PGP1_hu43860C_hg38_snvs.bed GT_CALL_RE='_step5_nomm\.tsv\.gz$' \
    GT_REDI_PART=score GT_RES_DIRS="celltype=results_gt_celltype_v8" GT_REDI_OUT=$C/redi \
    $RS human/scripts/gt_score_rediportal_v9.R
  ;;
keep) $RS human/scripts/gt_v9_celltype_covered_positives.R 30 ;;
submit)
  mkdir -p "$V9/logs" "$V6/spectrum/chunks"
  nreg=$(wc -l < "$V6/scan/regions.tsv")
  j1=$(qsub -terse -N v9_spec -t 1-"$nreg" -v STAGE=spec "$SELF" | cut -d. -f1)
  j2=$(qsub -terse -N v9_specg -hold_jid "$j1" -v STAGE=specgather "$SELF")
  j3=$(qsub -terse -N v9_kopile -pe threaded 4 -l h_vmem=8G -v STAGE=kopile "$SELF")
  j4=$(qsub -terse -N v9_ct_ko -l h_vmem=48G -v STAGE=ko "$SELF")
  j5=$(qsub -terse -N v9_ct_redi -l h_vmem=48G -v STAGE=redi "$SELF")
  j6=$(qsub -terse -N v9_ct_keep -hold_jid "$j4,$j5" -v STAGE=keep "$SELF")
  echo "spec $j1 (1-$nreg) gather $j2 kopile $j3 ko $j4 redi $j5 keep $j6"
  ;;
*) echo "STAGE?" >&2; exit 1 ;;
esac
