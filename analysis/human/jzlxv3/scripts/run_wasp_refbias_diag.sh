#!/bin/bash
# ============================================================================
# WASP reference-bias DIAGNOSTIC (Reviewer 2, Major 4) — human PBMC (jzlxv3)
# ============================================================================
# Quantifies the magnitude of STAR reference-mapping bias on A-to-I editing
# WITHOUT re-running the full editing pipeline. Uses the STAR WASP outputs that
# already exist:
#   ALL reads (biased estimate)   : 01_starsolobams_wasp/           (vW-tagged)
#   WASP-passed (bias-corrected)  : 01_starsolobams_wasp_filtered/  (vW==1|no-vW)
# STAR's WASP flips the allele at each REDIportal variant, re-maps, and flags
# reads whose mapping is allele-dependent (vW>1). Those are precisely the
# reference-biased reads. Comparing editing levels at REDIportal A-to-I sites
# between the two sets measures how much reference bias distorts editing.
#
# This is a DIAGNOSTIC only (human/mouse, REDIportal-guided); the principled
# cross-species correction is scAllele (local reassembly), already a caller.
#
# Pipeline:
#   1) REDIportal A-to-I sites (main chroms) -> BED
#   2) merge the 4 SRR BAMs in each set -> merged_all.bam / merged_pass.bam
#   3) strand-aware pileup at those sites (project engine gt_site_pileup.sh)
#   4) read-loss headline counts
#   5) R/ggplot2 analysis -> figures + summary TSV
#
# Usage:  qsub human/jzlxv3/scripts/run_wasp_refbias_diag.sh
#$ -N wasp_refbias
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 6
#$ -l h_vmem=10G
#$ -l h_rt=24:00:00
# Do NOT pass -V (SGE JSV / truncated PATH).

set -euo pipefail
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

THREADS=6
BASE="${SCRNA_DIR:-.}"
ENG="${BASE}/human/scripts"
REF="${BASE}/human/ref_annots/hg38.analysisSet.fa"
CHRLIST="${BASE}/human/ref_annots/chrlist_main_chrs.txt"
REDI="${BASE}/human/ref_annots/REDIportal_hg38_v2.vcf.gz"
DATA="${BASE}/human/jzlxv3"
ALLDIR="${DATA}/01_starsolobams_wasp"
PASSDIR="${DATA}/01_starsolobams_wasp_filtered"
OUT="${DATA}/wasp_refbias"
mkdir -p "$OUT"

echo "=== $(date) === WASP reference-bias diagnostic (jzlxv3 PBMC) ==="

# ── Step 1: REDIportal A-to-I sites on main chromosomes -> BED ──────────────
BED="${OUT}/reditportal_AI_mainchr.bed"
if [ ! -s "$BED" ]; then
    echo "--- $(date) --- building REDIportal A-to-I BED (main chroms) ---"
    zcat "$REDI" \
      | awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2-1, $2}' \
      | awk 'NR==FNR{keep[$1]=1; next} ($1 in keep)' "$CHRLIST" - \
      | sort -k1,1 -k2,2n > "$BED"
fi
echo "  sites: $(wc -l < "$BED")"

# ── Step 2: merge the 4 SRR BAMs per set (coordinate-sorted; cached) ────────
for SET in all pass; do
    if [ "$SET" = "all" ]; then SRC="$ALLDIR"; else SRC="$PASSDIR"; fi
    MB="${OUT}/merged_${SET}.bam"
    if [ ! -s "$MB" ]; then
        echo "--- $(date) --- merging ${SET} set from ${SRC} ---"
        # shellcheck disable=SC2086
        samtools merge -f -@ "$THREADS" "$MB" ${SRC}/*Aligned.sortedByCoord.out.bam
        samtools index -@ "$THREADS" "$MB"
    fi
done

# ── Step 3: strand-aware editing-level pileup at REDIportal sites ───────────
# 10x 3' libraries are forward (second-strand). -B (no BAQ), MQ>=30, BQ>=30.
for SET in all pass; do
    PU="${OUT}/pileup_${SET}.tsv"
    if [ ! -s "$PU" ]; then
        echo "--- $(date) --- pileup (${SET}) ---"
        bash "${ENG}/gt_site_pileup.sh" \
            --bam "${OUT}/merged_${SET}.bam" \
            --bed "$BED" \
            --fasta "$REF" \
            --out "$PU" \
            --libtype forward --min-bq 30 --min-mq 30 --max-depth 1000000
    fi
done

# ── Step 4: read-loss headline (global) ────────────────────────────────────
echo "--- $(date) --- read-loss counts ---"
A=$(samtools view -c -@ "$THREADS" "${OUT}/merged_all.bam")
P=$(samtools view -c -@ "$THREADS" "${OUT}/merged_pass.bam")
{
    echo -e "metric\tvalue"
    echo -e "reads_all\t${A}"
    echo -e "reads_wasp_pass\t${P}"
    echo -e "reads_removed\t$((A - P))"
    echo -e "frac_removed\t$(awk -v a="$A" -v p="$P" 'BEGIN{ if(a>0) printf "%.6f", (a-p)/a; else print "NA" }')"
} > "${OUT}/readloss.tsv"
cat "${OUT}/readloss.tsv"

# ── Step 5: R/ggplot2 analysis ─────────────────────────────────────────────
echo "--- $(date) --- R analysis ---"
Rscript "${DATA}/scripts/wasp_refbias_analysis.R" "$OUT"

echo "=== $(date) === DONE. Outputs in ${OUT} ==="
