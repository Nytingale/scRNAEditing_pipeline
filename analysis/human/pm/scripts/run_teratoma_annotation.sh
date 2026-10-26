#!/bin/bash
#$ -N ter_annot
#$ -pe threaded 8
#$ -l h_vmem=12G
#$ -l h_rt=48:00:00
#$ -j y
#$ -o ${SCRNA_DIR:-.}/human/pm/logs/teratoma_annotation.log
#$ -cwd
#$ -o ${SCRNA_DIR:-.}/logs/
#$ -e ${SCRNA_DIR:-.}/logs/

set -euo pipefail
export PS1="${PS1:-}"

RSCRIPT="${CONDA_ENVS:-/opt/conda/envs}/scRNAEditing/bin/Rscript"
CRDIR="${SCRNA_DIR:-.}/human/pm/cellranger_out"

echo "=== Teratoma cell type annotation ==="
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "Cores: ${NSLOTS:-unset}"
echo "Rscript: $RSCRIPT"

if [ ! -x "$RSCRIPT" ]; then
    echo "ERROR: conda Rscript not found: $RSCRIPT" >&2
    exit 1
fi

echo "=== Preflight CellRanger outs ==="
for s in T1 T2 T3 T4; do
    mtx="${CRDIR}/${s}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
    if [ ! -s "$mtx" ]; then
        mtx="${CRDIR}/${s}/outs/filtered_feature_bc_matrix/matrix.mtx"
    fi
    if [ ! -s "$mtx" ]; then
        echo "ERROR: $s filtered matrix missing (CellRanger not finished?): $mtx" >&2
        exit 1
    fi
    echo "  $s matrix OK ($(ls -lh "$mtx" | awk '{print $5}'))"
done

export OMP_NUM_THREADS="${NSLOTS:-4}"
export MKL_NUM_THREADS="${NSLOTS:-4}"
export OPENBLAS_NUM_THREADS=1

"$RSCRIPT" ${SCRNA_DIR:-.}/human/pm/scripts/teratoma_annotation.R

echo "=== Done: $(date) ==="
