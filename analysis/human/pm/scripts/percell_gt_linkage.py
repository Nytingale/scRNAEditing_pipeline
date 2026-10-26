#!/usr/bin/env python3
"""Per-cell reads linking a candidate editing site to a heterozygous donor SNP (R3.9).

One pass per per-teratoma GT BAM (02_gt_bams/<T>_<arm>.bam), chromosome by chromosome. A read links a site and a
SNP when it has an aligned base at both, each with base quality >= 30, the site base is the reference or the
edited base (A>G for ref A, T>C for ref T) and the SNP base is the SNP reference or alternate allele. Read filters
follow percell_gt_site_counts.py (Unprocessed tier): uniquely mapped (NH == 1), not secondary / supplementary /
QC-fail, MAPQ >= 30, CB in the GT barcode list, no dedup; distinct UB per cell are also counted (UMI-level).
SNPs: PGP1 GT 0/1, FILTER PASS, single-base SNVs; sites at any PGP1 variant position are dropped. pysam is used
because per-read checks of pairs of positions have no efficient samtools / R equivalent.

Outputs (<out_prefix>.*):
  span.tsv.gz     sample barcode site n_reads n_umis       (cell x site pairs with >= 1 linking read)
  cov.tsv.gz      sample barcode site n_reads n_umis       (KO-positive sites only; reads at the site, any SNP)
  readlen.tsv     sample read_length n                     (aligned query length of counted reads)
"""
import argparse
import bisect
import gzip
import sys
from collections import Counter, defaultdict

import pysam

EDIT_ALT = {"A": "G", "T": "C"}


def load_snps(vcf):
    het, allpos = defaultdict(dict), defaultdict(set)
    with gzip.open(vcf, "rt") as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            f = ln.split("\t", 10)
            chrom, pos, ref, alt = f[0], int(f[1]), f[3].upper(), f[4].upper()
            allpos[chrom].add(pos)
            if f[6] == "PASS" and f[9].split(":")[0].strip() == "0/1" and len(ref) == 1 and len(alt) == 1:
                het[chrom][pos - 1] = (ref, alt)
    return het, allpos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--sites", required=True, help="TSV: site(chr:pos) ko_positive candidate")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--min_q", type=int, default=30)
    ap.add_argument("--out_prefix", required=True)
    a = ap.parse_args()

    with open(a.cells) as fh:
        cells = {ln.split("\t")[0].strip() for ln in fh if ln.strip()}
    het, allpos = load_snps(a.vcf)
    fa = pysam.FastaFile(a.fasta)

    sites = defaultdict(list)  # chrom -> [(pos0, site, ref, alt, ko)]
    n_drop = 0
    with open(a.sites) as fh:
        next(fh)
        for ln in fh:
            site, ko, _ = ln.rstrip("\n").split("\t")
            chrom, pos = site.rsplit(":", 1)
            pos = int(pos)
            if pos in allpos.get(chrom, ()):
                n_drop += 1
                continue
            ref = fa.fetch(chrom, pos - 1, pos).upper()
            if ref not in EDIT_ALT:
                continue
            sites[chrom].append((pos - 1, site, ref, EDIT_ALT[ref], ko == "1"))
    for c in sites:
        sites[c].sort()
    print(f"{a.sample}: {sum(map(len, sites.values()))} sites ({n_drop} at PGP1 variant positions dropped), "
          f"{sum(map(len, het.values()))} het SNVs, {len(cells)} cells", file=sys.stderr)

    bam = pysam.AlignmentFile(a.bam)
    contigs = set(bam.references)
    span = defaultdict(lambda: [0, set()])
    cov = defaultdict(lambda: [0, set()])
    rlen = Counter()
    for chrom in sorted(sites):
        if chrom not in contigs:
            continue
        sl = sites[chrom]
        spos = [s[0] for s in sl]
        snp = het.get(chrom, {})
        npos = sorted(snp)
        for r in bam.fetch(chrom):
            if r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_unmapped:
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
            blocks = r.get_blocks()
            hits_s, hits_n = [], []
            for b0, b1 in blocks:
                i = bisect.bisect_left(spos, b0)
                while i < len(spos) and spos[i] < b1:
                    hits_s.append(i)
                    i += 1
                j = bisect.bisect_left(npos, b0)
                while j < len(npos) and npos[j] < b1:
                    hits_n.append(npos[j])
                    j += 1
            if not hits_s:
                continue
            rlen[r.query_length] += 1
            q = r.query_qualities
            seq = r.query_sequence
            qmap = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
            ub = r.get_tag("UB") if r.has_tag("UB") else None
            snp_ok = False
            for p in hits_n:
                qp = qmap.get(p)
                if qp is not None and q[qp] >= a.min_q and seq[qp].upper() in snp[p]:
                    snp_ok = True
                    break
            for i in hits_s:
                p0, site, ref, alt, ko = sl[i]
                qp = qmap.get(p0)
                if qp is None or q[qp] < a.min_q or seq[qp].upper() not in (ref, alt):
                    continue
                if ko:
                    c = cov[(cb, site)]
                    c[0] += 1
                    if ub:
                        c[1].add(ub)
                if snp_ok:
                    c = span[(cb, site)]
                    c[0] += 1
                    if ub:
                        c[1].add(ub)
        print(f"{a.sample}: {chrom} done, {len(span)} linked cell x site pairs", file=sys.stderr, flush=True)

    for name, d in (("span", span), ("cov", cov)):
        with gzip.open(f"{a.out_prefix}.{name}.tsv.gz", "wt") as out:
            out.write("sample\tbarcode\tsite\tn_reads\tn_umis\n")
            for (cb, site), (n, u) in d.items():
                out.write(f"{a.sample}\t{cb}\t{site}\t{n}\t{len(u)}\n")
    with open(f"{a.out_prefix}.readlen.tsv", "w") as out:
        out.write("sample\tread_length\tn\n")
        for k in sorted(rlen):
            out.write(f"{a.sample}\t{k}\t{rlen[k]}\n")
    print(f"{a.sample}: done; {len(span)} linked pairs, {len(cov)} KO-positive coverage pairs", file=sys.stderr)


if __name__ == "__main__":
    main()
