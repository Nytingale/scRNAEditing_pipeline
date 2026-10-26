#!/bin/bash
# Submit STARsolo alignment jobs for teratoma (PM) dataset
# 8 SRR accessions → 8 parallel alignment jobs
# After all finish, merge job combines per-teratoma replicates

SCRIPTDIR="${SCRNA_DIR:-.}/human/pm/scripts"
LOGDIR="${SCRNA_DIR:-.}/human/pm/logs"
mkdir -p "$LOGDIR"

SRRS="SRR26992449 SRR26992450 SRR26992452 SRR26992453 SRR26992455 SRR26992456 SRR26992458 SRR26992459"

ALIGN_JOBIDS=""
for SRR in $SRRS; do
    JOB=$(qsub -N "star_${SRR}" \
        -pe threaded 4 \
        -l h_rt=24:00:00 \
        -l h_vmem=10G \
        -j y \
        -o "${LOGDIR}/star_${SRR}.log" \
        -cwd \
        -V \
        -b y \
        "source ~/.bashrc && mamba activate scRNAEditing && bash ${SCRIPTDIR}/align_teratoma.sh ${SRR}")

    JOBID=$(echo "$JOB" | awk '{print $3}')
    ALIGN_JOBIDS="${ALIGN_JOBIDS:+${ALIGN_JOBIDS},}${JOBID}"
    echo "Submitted: $JOB"
done

# Merge job depends on all alignment jobs
qsub -N "merge_pm" \
    -hold_jid "$ALIGN_JOBIDS" \
    -pe threaded 4 \
    -l h_rt=8:00:00 \
    -l h_vmem=4G \
    -j y \
    -o "${LOGDIR}/merge_teratoma.log" \
    -cwd \
    -V \
    -b y \
    "source ~/.bashrc && mamba activate scRNAEditing && bash ${SCRIPTDIR}/merge_teratoma.sh"

echo ""
echo "Alignment jobs: $ALIGN_JOBIDS"
echo "Merge job will run after all alignments complete."
echo "Monitor with: qstat -u \$USER"
echo "Logs in: $LOGDIR"
