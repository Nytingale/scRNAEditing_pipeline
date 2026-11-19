#!/usr/bin/env python3
"""Per-site read and molecule counts of each base (Reviewer 2, major 1; user 2026-10-22).

  molecule_counts_v8.py --bam B.bam --fasta G.fa --sites S.tsv.gz --min_bq Q --min_mapq M --out O.tsv.gz [--threads N]

sites: tab-separated chr, pos (1-based), no header, one row per position.
A molecule is one (CB cell barcode, UB UMI) pair; a read without CB, or with UB missing or '-', counts as a molecule
of its own. Reads: MAPQ >= M, not unmapped / secondary / supplementary / QC-fail / duplicate, base quality >= Q, an
aligned base at the site (no deletion or splice gap).
Output: chr, pos, ref (genome), depth, n_A n_C n_G n_T (reads), m_A m_C m_G m_T (molecules), capped (pileup depth
cap reached; counts are then from the first MAX_DEPTH reads).
"""
import argparse
import gzip
import multiprocessing as mp
from collections import defaultdict

import pysam

MAX_DEPTH = 200000
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
BASES = "ACGT"


def count_chrom(job):
    bam_path, fasta_path, chrom, positions, min_bq, min_mapq = job
    bam = pysam.AlignmentFile(bam_path)
    fa = pysam.FastaFile(fasta_path)
    out = []
    for pos in positions:
        p0 = pos - 1
        reads = dict.fromkeys(BASES, 0)
        mols = defaultdict(set)
        depth = 0
        for col in bam.pileup(chrom, p0, p0 + 1, truncate=True, stepper="nofilter", min_base_quality=0,
                              ignore_overlaps=False, ignore_orphans=False, max_depth=MAX_DEPTH):
            if col.reference_pos != p0:
                continue
            depth = col.nsegments
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                aln = pr.alignment
                if aln.flag & SKIP_FLAGS or aln.mapping_quality < min_mapq:
                    continue
                q = pr.query_position
                if aln.query_qualities[q] < min_bq:
                    continue
                b = aln.query_sequence[q]
                if b not in reads:
                    continue
                reads[b] += 1
                cb = aln.get_tag("CB") if aln.has_tag("CB") else None
                ub = aln.get_tag("UB") if aln.has_tag("UB") else None
                key = (cb, ub) if cb and ub and ub != "-" else ("read", aln.query_name)
                mols[b].add(key)
        ref = fa.fetch(chrom, p0, p0 + 1).upper()
        out.append("\t".join([chrom, str(pos), ref, str(depth)] + [str(reads[b]) for b in BASES] +
                             [str(len(mols[b])) for b in BASES] + [str(int(depth >= MAX_DEPTH))]))
    return chrom, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--min_bq", type=int, required=True)
    ap.add_argument("--min_mapq", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    by_chrom = defaultdict(set)
    with gzip.open(a.sites, "rt") as f:
        for line in f:
            c, p = line.split("\t")[:2]
            by_chrom[c].add(int(p))
    contigs = set(pysam.AlignmentFile(a.bam).references)
    missing = sorted(set(by_chrom) - contigs)
    if missing:
        raise SystemExit(f"contigs not in BAM: {missing[:5]}")
    jobs = sorted(((a.bam, a.fasta, c, sorted(ps), a.min_bq, a.min_mapq) for c, ps in by_chrom.items()),
                  key=lambda j: -len(j[3]))
    with mp.Pool(a.threads) as pool:
        res = dict(pool.imap_unordered(count_chrom, jobs))
    with gzip.open(a.out, "wt") as o:
        o.write("chr\tpos\tref\tdepth\t" + "\t".join(f"n_{b}" for b in BASES) + "\t" +
                "\t".join(f"m_{b}" for b in BASES) + "\tcapped\n")
        for c in sorted(res):
            for line in res[c]:
                o.write(line + "\n")
    print(f"{sum(len(v) for v in res.values())} sites, {len(res)} contigs -> {a.out}")


if __name__ == "__main__":
    main()
