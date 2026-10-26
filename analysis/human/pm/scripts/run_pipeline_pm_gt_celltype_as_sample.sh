#!/bin/bash
# Per-CELLTYPE teratoma GT discovery, run as "celltype-as-sample".
# Each cell type (pooled across the samples of one `group` -- for teratoma, the
# four teratomas of one CRISPR arm) is a pre-split BAM treated as a sample, so it
# flows through the SAMPLE-level pipeline and receives the IDENTICAL bam_type
# matrix, including BQSR (recal/recalq20) per celltype BAM. Pre-split-then-dedup
# == split-then-dedup, because dedup is cell-aware and cell type never splits a
# dedup group. Mirrors run_pipeline_pm_gt_persample.sh exactly except for the
# samplesheet / outdir / work-dir.
#
# The cell-type minimum is NOT hardcoded here: it is applied by
# validate_samplesheet.R --min-cells when the split plan is emitted (default 0 =
# keep all) and recorded in publication/CELLTYPE_SCOPE.tsv.
#
# Depends on split_gt_bams_celltype.sh; submit both via
# submit_celltype_as_sample.sh, which holds this driver on the split array.
# Combines the per-unit fragments into the samplesheet, then launches nextflow --
# but only after the completeness guard below passes.
#$ -N nf_pm_ct
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/
#$ -pe threaded 1
#$ -l h_vmem=64G
#$ -l h_rt=336:00:00

set -eo pipefail
export NXF_OPTS="-Xms1g -Xmx10g"

BASEDIR="${SCRNA_DIR:-.}"
NFDIR="${BASEDIR}/nextflow"
REFDIR="${BASEDIR}/human/ref_annots"
DATADIR="${BASEDIR}/human/pm"
SAMPLESHEET="${DATADIR}/samplesheet_pm_gt_celltype_as_sample.csv"
FRAGDIR="${FRAGDIR:-${DATADIR}/_celltype_frags}"
PLANDIR="${PLANDIR:-$FRAGDIR}"
RUNDIR="${DATADIR}/nf_run_gt_celltype_as_sample"
mkdir -p "$RUNDIR"

# ── COMPLETENESS GUARD ───────────────────────────────────────────────────────
# >>> Do NOT relax this to a row count. <<<
# -hold_jid releases on job COMPLETION, not success. The old check was
# `[ "$N" -ge 1 ]` -- a non-empty samplesheet -- which passes when SOME units
# succeeded. If the NTC arm OOMs while the other four finish, that check builds a
# perfectly valid-looking 26-row samplesheet and launches the full GT run WITH NO
# CONTROL ARM: NTC is 17 of teratoma's 43 cell types and 9540 of ~14.5k cells.
# The failure is silent and the output is a benchmark with nothing to compare
# against, which is far worse than a crash.
#
# So: every unit in _expected.txt must have both a _frag_<unit>.csv and the
# _ok_<unit> sentinel that split_gt_bams_celltype.sh touches only after running
# to completion under `set -e`.
EXPECTED="${PLANDIR}/_expected.txt"
[ -s "$EXPECTED" ] || { echo "ERROR: missing or empty ${EXPECTED}; was the plan emitted?" >&2; exit 1; }

MISSING=0
while read -r UNIT; do
    [ -n "$UNIT" ] || continue
    if [ ! -f "${FRAGDIR}/_ok_${UNIT}" ]; then
        echo "ERROR: unit '${UNIT}' has no completion sentinel (${FRAGDIR}/_ok_${UNIT})" >&2
        MISSING=$((MISSING + 1))
    elif [ ! -s "${FRAGDIR}/_frag_${UNIT}.csv" ]; then
        echo "ERROR: unit '${UNIT}' is marked complete but its fragment is missing or empty" >&2
        MISSING=$((MISSING + 1))
    fi
done < "$EXPECTED"

if [ "$MISSING" -gt 0 ]; then
    echo "ERROR: ${MISSING} of $(wc -l < "$EXPECTED") unit(s) did not complete the split." >&2
    echo "       REFUSING to launch -- a partial unit set silently produces a benchmark" >&2
    echo "       with a missing arm. Fix the failed split task(s) and resubmit." >&2
    exit 1
fi
echo "=== $(date) === completeness guard passed: all $(wc -l < "$EXPECTED") unit(s) present"

# ── Combine per-unit fragments into the samplesheet ─────────────────────────
# Assembled unit by unit in _expected.txt order rather than by a _frag_*.csv
# glob, so the sheet can only ever contain units the guard above verified.
# `group` carries the unit through, keeping this sheet spec-conformant and
# re-validatable in its own right.
echo "sample_id,bam,celltype_annot,called_barcodes,platform,fastq_r1,fastq_r2,group" > "$SAMPLESHEET"
while read -r UNIT; do
    [ -n "$UNIT" ] || continue
    cat "${FRAGDIR}/_frag_${UNIT}.csv" >> "$SAMPLESHEET"
done < "$EXPECTED"
N=$(($(wc -l < "$SAMPLESHEET") - 1))
[ "$N" -ge 1 ] || { echo "ERROR: samplesheet empty (split not done?)" >&2; exit 1; }
echo "=== $(date) === samplesheet: ${N} celltype-samples across $(wc -l < "$EXPECTED") unit(s)"

# The generated sheet must itself satisfy the spec before a multi-day run starts.
Rscript "${NFDIR}/bin/validate_samplesheet.R" "$SAMPLESHEET" \
    --dataset "${DATASET:-pm_gt_celltype}" --tsv-out "${BASEDIR}/publication/CELLTYPE_SCOPE.tsv" \
    || { echo "ERROR: generated samplesheet failed validation" >&2; exit 1; }

cd "$RUNDIR"
export PS1="${PS1:-}"
export PATH="/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:${PATH}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

echo "=== $(date) === Launching teratoma GT pipeline (CELLTYPE-as-SAMPLE) ==="

# >>> PER-CELL IS OFF. params.percell_analyses DEFAULTS TO TRUE. <<<
# SINTO_SPLIT_PERCELL + PERCELL_MPILEUP are memory-heavy, sit on NO
# ground-truth critical path, and keep the Nextflow job in state `r` for many
# hours after every benchmark table is written. On the 09-16 merged run per-cell
# ran 14.7 h past the point where all 5 arms were at 64/64, which made the job's
# own completion signal useless and forced the phase-2 scorer to poll data
# readiness instead. These are caller-comparison runs, so it is disabled.

nextflow run "${NFDIR}/main.nf" \
    -profile sge,conda \
    --input "${SAMPLESHEET}" \
    --outdir "${DATADIR}/results_gt_celltype_as_sample" \
    --genome_fasta "${REFDIR}/hg38.analysisSet.fa" \
    --genome_gtf "${REFDIR}/hg38.ncbiRefSeq.gtf" \
    --snp_file "${BASEDIR}/annotations/SNP_files/dbSNP156_common_GRanges.RData" \
    --strand_annots "${BASEDIR}/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --alu_annots "${BASEDIR}/annotations/annotation_files/extended_annots_GRCh38_analysis_set.RData" \
    --rediportal "${REFDIR}/REDIportal_hg38_v2.vcf.gz" \
    --known_sites_vcf "${REFDIR}/dbSNP156_GRCh38.vcf.gz" \
    --refseq_bed "${BASEDIR}/review/ref_annots/hg38_RefSeq.bed" \
    --main_chroms "${REFDIR}/chrlist_main_chrs.txt" \
    --callers "bcftools,reditools3,gatk,jacusa2,scallele" \
    --scallele_conda "${CONDA_ENVS:-/opt/conda/envs}/scAllele" \
    --dedup "none,coordinate,umi" \
    --strand_modes "annotation,hybrid,caller" \
    --bam_types "sample,sampleRmDup,sampleUmiDup,recal,recalq20" \
    --qual 30 \
    --min_cov 10 \
    --requant_editing false \
    --percell_analyses false \
    -work-dir "${DATADIR}/work_gt_celltype_as_sample" \
    -resume

echo "=== $(date) === Celltype-as-sample pipeline complete ==="
