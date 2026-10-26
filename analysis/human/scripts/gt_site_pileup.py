#!/usr/bin/env python3
"""gt_site_pileup.py  (dataset-agnostic GT evaluator, step 2/3)

Independent re-quantification of editing at discovered sites. For one arm BAM,
count the A/C/G/T reads at each union site with base- and mapping-quality
filters, so downstream editing levels are computed directly from the reads and
do NOT depend on any variant caller's cov_ref/cov_alt or strand handling.

Input : --sites  gt_eval_union_sites.tsv  (chr, pos[1-based], strand, site; header)
Output: TSV  site, chr, pos, A, C, G, T, cov   (cov = A+C+G+T after QC)

Quality filters mirror the pipeline defaults (base & mapping quality >= 30) and
count only primary, mapped, non-duplicate, non-QCfail alignments.

Usage:
  python gt_site_pileup.py --bam merged_ARM.bam --sites union_sites.tsv \
         --out pileup_ARM.tsv [--min-bq 30] [--min-mq 30]
"""
import argparse
import csv
import sys

import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sites", required=True, help="TSV: chr,pos(1-based),strand,site with header")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-bq", type=int, default=30, help="min base quality")
    ap.add_argument("--min-mq", type=int, default=30, help="min mapping quality")
    args = ap.parse_args()

    min_mq = args.min_mq

    def keep(read):
        return (
            (not read.is_unmapped)
            and (not read.is_secondary)
            and (not read.is_supplementary)
            and (not read.is_qcfail)
            and (not read.is_duplicate)
            and (read.mapping_quality >= min_mq)
        )

    # Load sites
    sites = []
    with open(args.sites) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            sites.append((row["site"], row["chr"], int(row["pos"])))
    if not sites:
        sys.exit(f"No sites in {args.sites}")

    bam = pysam.AlignmentFile(args.bam, "rb")
    refs = set(bam.references)

    n = 0
    n_missing_contig = 0
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["site", "chr", "pos", "A", "C", "G", "T", "cov"])
        for site, chrom, pos in sites:
            if chrom not in refs:
                n_missing_contig += 1
                w.writerow([site, chrom, pos, 0, 0, 0, 0, 0])
                continue
            # count_coverage returns 4 arrays (A,C,G,T) over [start,stop)
            cov = bam.count_coverage(
                chrom, pos - 1, pos,
                quality_threshold=args.min_bq,
                read_callback=keep,
            )
            a, c, g, t = (int(cov[0][0]), int(cov[1][0]), int(cov[2][0]), int(cov[3][0]))
            w.writerow([site, chrom, pos, a, c, g, t, a + c + g + t])
            n += 1
            if n % 20000 == 0:
                sys.stderr.write(f"  ...{n} sites piled up\n")
    bam.close()
    sys.stderr.write(
        f"Done: {n} sites written to {args.out}"
        + (f" ({n_missing_contig} sites on contigs absent from BAM)\n" if n_missing_contig else "\n")
    )


if __name__ == "__main__":
    main()
