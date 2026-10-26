#!/bin/bash
# Per-cell reads linking a site to a heterozygous PGP1 SNP (R3.9 linkage figure), one array task per per-teratoma
# GT sample (T1-T4 x 5 arms = 20; unprocessed reads, 02_gt_bams). Sites: v8 KO positives (union, 5,358) and the
# merged NTC arm GIREMI candidate A-to-G sites (percell_gt_linkage_v8/linkage_sites.tsv).
# Usage:  qsub -t 1-20 human/pm/scripts/percell_gt_linkage_v8.sh
#$ -N pc_link8
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 1
#$ -l h_vmem=32G
#$ -l h_rt=48:00:00
# Do NOT pass -V (SGE JSV / truncated PATH).

set -eo pipefail
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

BASE="${SCRNA_DIR:-.}"
PM="${BASE}/human/pm"
OUT="${PM}/percell_gt_linkage_v8"
SITES="${OUT}/linkage_sites.tsv"
[ -s "$SITES" ] || { echo "missing $SITES" >&2; exit 1; }

SAMPLES=($(for t in T1 T2 T3 T4; do for a in AAVS1_ctrl ADAR1_KO ADARB1_KO ADARB2_KO NTC_ctrl; do echo "${t}_${a}"; done; done))
S="${SAMPLES[$((SGE_TASK_ID - 1))]}"
BAM="${PM}/02_gt_bams/${S}.bam"
echo "=== $(date) === ${S} (${BAM})"

python "${PM}/scripts/percell_gt_linkage.py" \
    --bam "$BAM" \
    --sample "$S" \
    --cells "${PM}/gt_barcodes/${S}.tsv" \
    --sites "$SITES" \
    --vcf "${PM}/pgp1_snps/PGP1_hu43860C_hg38_snvs.vcf.gz" \
    --fasta "${BASE}/human/ref_annots/hg38.analysisSet.fa" \
    --out_prefix "${OUT}/${S}"

echo "=== $(date) === done ${S}"
