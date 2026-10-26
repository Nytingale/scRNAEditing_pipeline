#!/bin/bash
# One-shot submitter for the celltype-as-sample GT discovery.
#
#   0) validate the samplesheet and emit the split plan   (runs HERE, on the
#      login node, so a bad sheet fails in seconds instead of inside an array
#      task -- and so the array size is known before anything is queued)
#   1) split_gt_bams_celltype.sh          array -t 1-<n units>, sized from the plan
#   2) run_pipeline_pm_gt_celltype_as_sample.sh   -hold_jid on the split; checks
#      the completion sentinels, builds the samplesheet, runs the pipeline
#
# NOTHING here is teratoma-specific except the DEFAULTS below, all overridable:
#
#   SHEET      samplesheet to derive units from
#   DATASET    name recorded in publication/CELLTYPE_SCOPE.tsv
#   MIN_CELLS  cell-type minimum. DEFAULT 0 = keep every cell type. A non-zero
#              value is a filter that is NOT part of the original analysis, so it
#              is opt-in and always recorded in the scope TSV. It cannot be a
#              global constant: a teratoma-derived floor of 50 would erase jzlcel
#              outright (0 of its 48 cell types reach it).
#
# Usage:
#   ./submit_celltype_as_sample.sh
#   SHEET=.../samplesheet_jzl10x.csv DATASET=jzl10x ./submit_celltype_as_sample.sh
set -eo pipefail

QSUB=/cm/shared/apps/sge/2011.11p1/bin/linux-x64/qsub
BASEDIR="${SCRNA_DIR:-.}"
SCRIPTS="${BASEDIR}/human/pm/scripts"

SHEET="${SHEET:-${BASEDIR}/human/pm/samplesheet_pm_gt_persample.csv}"
DATASET="${DATASET:-pm_gt}"
MIN_CELLS="${MIN_CELLS:-0}"

DATADIR="$(dirname "$SHEET")"
PLANDIR="${PLANDIR:-${DATADIR}/_celltype_frags}"
OUTBAM="${OUTBAM:-${DATADIR}/02_gt_bams_celltype}"
OUTBC="${OUTBC:-${DATADIR}/gt_barcodes_celltype}"
SCOPE_TSV="${SCOPE_TSV:-${BASEDIR}/publication/CELLTYPE_SCOPE.tsv}"

mkdir -p "${BASEDIR}/logs" "$PLANDIR"

export PS1="${PS1:-}"
source ~/.bashrc 2>/dev/null || true
mamba activate scRNAEditing

# ── 0. Validate + emit the split plan ───────────────────────────────────────
# --emit-plan writes nothing unless the sheet validates clean, so a failure here
# cannot leave a partial plan for the array to pick up. Stale _plan_/_cells_/
# _keep_ files from a previous run with different groups are cleared by the
# validator; the sentinels are cleared below.
echo "=== $(date) === validating ${SHEET} (dataset=${DATASET}, min_cells=${MIN_CELLS})"
Rscript "${BASEDIR}/nextflow/bin/validate_samplesheet.R" "$SHEET" \
    --dataset "$DATASET" --min-cells "$MIN_CELLS" \
    --tsv-out "$SCOPE_TSV" --emit-plan "$PLANDIR"

EXPECTED="${PLANDIR}/_expected.txt"
[ -s "$EXPECTED" ] || { echo "ERROR: validator wrote no ${EXPECTED}" >&2; exit 1; }
N=$(wc -l < "$EXPECTED")
[ "$N" -ge 1 ] || { echo "ERROR: no units in ${EXPECTED}" >&2; exit 1; }

# Sentinels and fragments from any previous run MUST go, or the driver's
# completeness check would pass on stale evidence.
rm -f "${PLANDIR}"/_ok_* "${PLANDIR}"/_frag_*.csv

echo "=== units (${N}) ==="
nl -ba "$EXPECTED"

# ── 1. Split array, sized from the plan ─────────────────────────────────────
OUT=$($QSUB -t "1-${N}" \
      -v PLANDIR="$PLANDIR",OUTBAM="$OUTBAM",OUTBC="$OUTBC",FRAGDIR="$PLANDIR" \
      "${SCRIPTS}/split_gt_bams_celltype.sh")
echo "$OUT"
JID=$(echo "$OUT" | grep -oE '[0-9]+' | head -1)
[ -n "$JID" ] || { echo "ERROR: could not parse split job id" >&2; exit 1; }
echo "split job id = ${JID}  (array 1-${N})"

# ── 2. Driver, held on the whole split array ────────────────────────────────
# -hold_jid releases on COMPLETION, not success, which is exactly why the driver
# verifies a sentinel per unit rather than trusting the hold.
$QSUB -hold_jid "$JID" \
      -v PLANDIR="$PLANDIR",OUTBAM="$OUTBAM",OUTBC="$OUTBC",FRAGDIR="$PLANDIR" \
      "${SCRIPTS}/run_pipeline_pm_gt_celltype_as_sample.sh"
echo "Submitted launch driver (held on split ${JID})."
