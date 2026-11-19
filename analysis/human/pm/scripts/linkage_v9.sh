#!/bin/bash
# Linked-SNV analysis (R3.9, Response Figure 22) on the v9 KO positives: union of the six v9 KO ground truths
# (gt_eval_v9/truth/gt_v9_positive_sites.tsv.gz, any t) + the merged NTC arm GIREMI candidate A-to-G sites
# (test_giremi_ntc_arm/cands.txt, as in v8). Counters read the unprocessed per-teratoma BAMs (02_gt_bams; no
# deduplication), so only the site list differs from percell_gt_linkage_v8.
#   STAGE=sites    build percell_gt_linkage_v9/linkage_sites.tsv (login node: awk only)
#   STAGE=submit   span / cell / allele counters on the four NTC arms (tasks 5,10,15,20), then STAGE=summary
#   STAGE=summary  percell_gt_linkage_cell_summary.R, TERATOMA_v8_linkage_ased_thresholds.R,
#                  TERATOMA_v8_linkage_snp_{source_data,plot}.R with LINK_DIR=percell_gt_linkage_v9
#$ -N link9
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -l h_vmem=48G
set -eo pipefail
BASE=${SCRNA_DIR:-.}
PM=$BASE/human/pm
export LINK_DIR=$PM/percell_gt_linkage_v9
RS=${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript
cd "$BASE"

case "${STAGE:?}" in
sites)
  mkdir -p "$LINK_DIR"
  zcat "$PM/gt_eval_v9/truth/gt_v9_positive_sites.tsv.gz" |
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
                $h["class_t25"]=="positive"||$h["class_t50"]=="positive"||$h["class_t75"]=="positive"{print $h["site"]}' |
    sort -u > "$LINK_DIR/ko_positive_sites.txt"
  awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{ko[$1]=1;next} {k=$1":"$2; if(($3"_"$4)=="A_G"||($3"_"$4)=="T_C") c[k]=1}
       END{print "site","ko_positive","candidate"; for(k in ko) print k,1,(k in c)?1:0; for(k in c) if(!(k in ko)) print k,0,1}' \
    "$LINK_DIR/ko_positive_sites.txt" "$PM/test_giremi_ntc_arm/cands.txt" > "$LINK_DIR/linkage_sites.tsv"
  awk 'NR>1{a+=$2;b+=$3} END{print NR-1" sites, "a" KO positives, "b" candidates"}' "$LINK_DIR/linkage_sites.tsv"
  ;;
submit)
  [ -s "$LINK_DIR/linkage_sites.tsv" ] || { echo "run STAGE=sites first" >&2; exit 1; }
  j1=$(qsub -terse -t 5-20:5 -v LINK_DIR="$LINK_DIR" human/pm/scripts/percell_gt_linkage_v8.sh | cut -d. -f1)
  j2=$(qsub -terse -t 5-20:5 -v LINK_DIR="$LINK_DIR" human/pm/scripts/percell_gt_linkage_cell_v8.sh | cut -d. -f1)
  j3=$(qsub -terse -t 5-20:5 -v LINK_DIR="$LINK_DIR" human/pm/scripts/percell_gt_linkage_allele_v8.sh | cut -d. -f1)
  j4=$(qsub -terse -hold_jid "$j1,$j2,$j3" -v STAGE=summary human/pm/scripts/linkage_v9.sh)
  echo "span $j1 cell $j2 allele $j3 summary $j4"
  ;;
summary)
  # readlen tables (all 20 GT samples) do not depend on the sites: the source data script reads them from LINK_DIR
  for f in "$PM"/percell_gt_linkage_v8/*.readlen.tsv; do
    [ -s "$LINK_DIR/$(basename "$f")" ] || cp -p "$f" "$LINK_DIR/"
  done
  S=publication/teratoma_gt/scripts
  $RS human/pm/scripts/percell_gt_linkage_cell_summary.R
  $RS $S/TERATOMA_v8_linkage_ased_thresholds.R
  $RS $S/TERATOMA_v8_linkage_snp_source_data.R
  $RS $S/TERATOMA_v8_linkage_snp_plot.R
  ;;
esac
