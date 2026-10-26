#!/usr/bin/env python3
"""Per-cell ref / edited read counts at teratoma ground-truth sites.

One pass per per-teratoma GT BAM (02_gt_bams/<T>_<arm>.bam), restricted to the
v5 TRUE + FP sites, grouping reads by CB. Read filters follow the pipeline's
Unprocessed ('sample') tier: uniquely mapped (NH == 1), not secondary /
supplementary / QC-fail, base and mapping quality >= 30, NO dedup (the most
favourable case for per-cell depth). The edited base is taken from the genome
reference: A>G (ref A) or T>C (ref T, minus-strand A>G). pysam is used because
per-barcode pileup counting has no efficient samtools/R equivalent.

Output (gz TSV, only cell x site pairs with >= 1 counted read):
  sample  barcode  site  class  ref_reads  alt_reads
"""
import argparse
import gzip
import sys
from collections import defaultdict

import pysam

EDIT_ALT = {"A": "G", "T": "C"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cells", required=True, help="GT barcode TSV (barcode in col 1)")
    ap.add_argument("--sites", required=True, help="TSV: site(chr:pos) class")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--min_q", type=int, default=30)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with open(a.cells) as fh:
        cells = {ln.split("\t")[0].strip() for ln in fh if ln.strip()}
    sites = []
    with open(a.sites) as fh:
        next(fh)
        for ln in fh:
            site, cls = ln.rstrip("\n").split("\t")[:2]
            chrom, pos = site.rsplit(":", 1)
            sites.append((chrom, int(pos), site, cls))

    fa = pysam.FastaFile(a.fasta)
    bam = pysam.AlignmentFile(a.bam)
    contigs = set(bam.references)
    n_rows = 0
    with gzip.open(a.out, "wt") as out:
        out.write("sample\tbarcode\tsite\tclass\tref_reads\talt_reads\n")
        for chrom, pos, site, cls in sites:
            if chrom not in contigs:
                continue
            ref = fa.fetch(chrom, pos - 1, pos).upper()
            alt = EDIT_ALT.get(ref)
            if alt is None:
                continue
            counts = defaultdict(lambda: [0, 0])
            for col in bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="nofilter",
                                  min_base_quality=0, max_depth=100_000_000,
                                  ignore_orphans=False, ignore_overlaps=False):
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip or pr.query_position is None:
                        continue
                    r = pr.alignment
                    if r.is_secondary or r.is_supplementary or r.is_qcfail:
                        continue
                    if r.mapping_quality < a.min_q:
                        continue
                    if r.has_tag("NH") and r.get_tag("NH") != 1:
                        continue
                    if not r.has_tag("CB"):
                        continue
                    cb = r.get_tag("CB")
                    if cb not in cells:
                        continue
                    if r.query_qualities[pr.query_position] < a.min_q:
                        continue
                    b = r.query_sequence[pr.query_position].upper()
                    if b == ref:
                        counts[cb][0] += 1
                    elif b == alt:
                        counts[cb][1] += 1
            for cb, (nr, na) in counts.items():
                out.write(f"{a.sample}\t{cb}\t{site}\t{cls}\t{nr}\t{na}\n")
                n_rows += 1
    print(f"{a.sample}: {len(sites)} sites, {len(cells)} cells, {n_rows} cell-site rows", file=sys.stderr)


if __name__ == "__main__":
    main()
