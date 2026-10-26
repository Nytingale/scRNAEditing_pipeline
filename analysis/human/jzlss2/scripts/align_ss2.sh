#!/bin/bash
# ============================================================================
# STAR alignment for Smart-seq2 dataset (JZLSS2) — PER CELL
# ============================================================================
# Smart-seq2: each SRR is a single cell. No cell barcodes, no UMIs.
# Paired-end: two cDNA reads (R1 + R2). Index reads (_3, _4) are ignored.
# Read length: ~50bp → sjdbOverhang=49
#
# Verification: confirms exactly 2 non-index reads exist per SRR, and that
# neither contains polyT runs (which would indicate barcode reads).
#
# Usage: bash align_ss2.sh <SRR_accession>
#
# WASP MODE (reference-bias diagnostic, Reviewer 2 Major 4). Setting WASP_VCF
# adds --waspOutputMode SAMtag and the vA/vG/vW tags, so reads whose alignment
# is allele-dependent come out flagged vW:i:2-7 (vW:i:1 = passed, no tag =
# overlaps no variant). One script, not a copy, so the paired-end cDNA read
# detection below cannot drift between the two alignments being compared:
#
#   WASP_VCF=<vcf> OUTDIR=<dir> bash align_ss2.sh <SRR>
#
# Both variables are optional and default to the plain alignment.
# ============================================================================

set -euo pipefail

REFDIR="${SCRNA_DIR:-.}/human/ref_annots"
FQDIR="${SCRNA_DIR:-.}/human/jzlss2/fastq_files"
OUTDIR="${OUTDIR:-${SCRNA_DIR:-.}/human/jzlss2/01_starsolobams}"
STARDB="${REFDIR}/stardb_hg38_49"
NTHREADS=${NTHREADS:-4}
WASP_VCF="${WASP_VCF:-}"

SAMPLE=${1:?Usage: bash align_ss2.sh <SRR_accession>}

# WASP tags ride alongside the standard attributes; without a VCF, STAR must
# not be given vA/vG/vW at all or it errors.
SAM_ATTRS="NH HI nM AS MD"
WASP_OPTS=()
if [ -n "$WASP_VCF" ]; then
    if [ ! -f "$WASP_VCF" ]; then
        echo "ERROR: WASP VCF not found: $WASP_VCF" >&2
        exit 1
    fi
    SAM_ATTRS="NH HI nM AS MD vA vG vW"
    WASP_OPTS=(--waspOutputMode SAMtag --varVCFfile "$WASP_VCF")
fi

# ---- Verify paired-end cDNA reads ----
detect_pe_reads() {
  local sample="$1"
  local fqdir="$2"
  local long_files=() long_lens=()

  for f in "$fqdir"/${sample}_*.fastq.gz; do
    [[ -f "$f" ]] || continue
    local len=$(zcat "$f" | head -2 | tail -1 | wc -c)
    len=$((len - 1))
    if [[ $len -le 12 ]]; then
      continue  # skip index reads
    fi
    long_files+=("$f")
    long_lens+=("$len")
  done

  if [[ ${#long_files[@]} -lt 2 ]]; then
    echo "ERROR: expected >=2 cDNA reads for SS2 $sample, found ${#long_files[@]}" >&2
    return 1
  fi

  # Sanity check: neither read should have polyT (would indicate barcode contamination)
  for i in 0 1; do
    local pt=$(zcat "${long_files[$i]}" | head -40 | awk '
      NR%4==2 {
        total++; seq=$0; max_t=0; cur_t=0;
        for(j=1;j<=length(seq);j++) {
          if(substr(seq,j,1)=="T") cur_t++;
          else { if(cur_t>max_t) max_t=cur_t; cur_t=0 }
        }
        if(cur_t>max_t) max_t=cur_t;
        if(max_t>=15) polyt++
      }
      END { print polyt+0 }')
    if [[ $pt -gt 5 ]]; then
      echo "WARNING: $(basename ${long_files[$i]}) has polyT signatures — may be barcode read, not cDNA" >&2
    fi
  done

  echo "${long_files[0]} ${long_files[1]}"
}

READ_PAIR=$(detect_pe_reads "$SAMPLE" "$FQDIR") || exit 1
R1=$(echo "$READ_PAIR" | awk '{print $1}')
R2=$(echo "$READ_PAIR" | awk '{print $2}')

echo "=== $(date) === Read detection for $SAMPLE ==="
echo "  R1: $(basename $R1)"
echo "  R2: $(basename $R2)"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Idempotent: an array of 768 tasks will be partially re-run after any task
# hits the wall clock or the memory cap, and re-aligning a finished cell wastes
# a 31 GB index load. quickcheck, not -s, so a BAM truncated by a killed task
# is redone rather than trusted.
FINAL_BAM="${OUTDIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
if [ -s "$FINAL_BAM" ] && samtools quickcheck "$FINAL_BAM" 2>/dev/null; then
    echo "=== $(date) === SKIP $SAMPLE: $FINAL_BAM already complete ==="
    exit 0
fi

echo "=== $(date) === Aligning $SAMPLE (Smart-seq2, single cell) ==="
[ -n "$WASP_VCF" ] && echo "    WASP VCF: $WASP_VCF"

STAR \
    --genomeDir "$STARDB" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMattributes $SAM_ATTRS \
    --outFilterType BySJout \
    --outBAMsortingBinsN 50 \
    --limitBAMsortRAM 10000000000 \
    --genomeLoad NoSharedMemory \
    ${WASP_OPTS[@]+"${WASP_OPTS[@]}"} \
    --outFileNamePrefix "${SAMPLE}_" \
    --runThreadN "$NTHREADS"

echo "=== $(date) === Indexing BAM ==="
samtools index "${SAMPLE}_Aligned.sortedByCoord.out.bam"

if [ -n "$WASP_VCF" ]; then
    echo "=== $(date) === WASP tag summary ==="
    samtools view "${SAMPLE}_Aligned.sortedByCoord.out.bam" \
      | grep -oP 'vW:i:\d+' | sort | uniq -c | sort -rn || true
fi

echo "=== $(date) === Done: $SAMPLE ==="
