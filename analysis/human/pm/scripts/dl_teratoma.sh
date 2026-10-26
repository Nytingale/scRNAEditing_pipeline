#!/bin/bash
#$ -l h_rt=48:00:00
#$ -pe threaded 4
#$ -j y

source ~/.bashrc
mamba activate scRNAEditing
export PATH=${SOFTWARE_DIR:-/opt/software}/sratoolkit.3.0.7-centos_linux64/bin:$PATH

cd "$WKDIR"
mkdir -p fastq_files
cd fastq_files

echo "=== $(date) === Prefetching ${SRRID} ==="
prefetch --max-size 80G "$SRRID"

echo "=== $(date) === Dumping FASTQ ==="
parallel-fastq-dump --sra-id "$SRRID" --threads 4 --outdir ./ --split-files --gzip

echo "=== $(date) === Cleaning up SRA cache ==="
rm -rf "$SRRID"

echo "=== $(date) === Done ==="
ls -lh ${SRRID}*.fastq.gz 2>/dev/null
