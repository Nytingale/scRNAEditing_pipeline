#!/usr/bin/env python3
"""Allele-resolved linkage of candidate editing sites to heterozygous donor SNVs (R3.9, bulk ASED criteria).

Companion of percell_gt_linkage.py / percell_gt_linkage_cell.py (same read filters, sites and SNVs). One pass per
per-teratoma GT BAM, chromosome by chromosome, at molecule level (CB + UB + STARsolo GX gene; the reads of one
molecule are combined, as both mates of a fragment are in bulk paired-end data). A molecule is counted for a site
and a SNV when its reads have one unambiguous base at each (site: reference or edited base; SNV: reference or
alternate allele; BQ >= 30). Counts are kept per cell; pooling over cells (pseudobulks, as in bulk allele-specific
editing analyses: Park et al. 2017, Genome Biol; Zhou et al. 2018, Nucleic Acids Res) is done downstream.
pysam is used because per-read checks of pairs of positions have no efficient samtools / R equivalent.

Output <out_prefix>.allele_cell.tsv.gz:
  sample barcode site snv snv_allele unedited edited      (molecules; snv_allele = ref | alt)
"""
import argparse
import bisect
import gzip
import sys
from collections import defaultdict

import pysam

from percell_gt_linkage import EDIT_ALT, load_snps

SNV_WINDOW = 5000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--min_q", type=int, default=30)
    ap.add_argument("--out_prefix", required=True)
    a = ap.parse_args()

    with open(a.cells) as fh:
        cells = {ln.split("\t")[0].strip() for ln in fh if ln.strip()}
    het, allpos = load_snps(a.vcf)
    fa = pysam.FastaFile(a.fasta)
    sites = defaultdict(list)
    with open(a.sites) as fh:
        next(fh)
        for ln in fh:
            site = ln.split("\t", 1)[0]
            chrom, pos = site.rsplit(":", 1)
            pos = int(pos)
            if pos in allpos.get(chrom, ()):
                continue
            ref = fa.fetch(chrom, pos - 1, pos).upper()
            if ref in EDIT_ALT:
                sites[chrom].append((pos - 1, site, ref, EDIT_ALT[ref]))
    for c in sites:
        sites[c].sort()

    bam = pysam.AlignmentFile(a.bam)
    contigs = set(bam.references)
    out = gzip.open(f"{a.out_prefix}.allele_cell.tsv.gz", "wt")
    out.write("sample\tbarcode\tsite\tsnv\tsnv_allele\tunedited\tedited\n")
    n_rows = 0
    for chrom in sorted(sites):
        if chrom not in contigs:
            continue
        sl = sites[chrom]
        spos = [s[0] for s in sl]
        snp = het.get(chrom, {})
        def near_site(p):
            i = bisect.bisect_left(spos, p)
            return (i < len(spos) and spos[i] - p <= SNV_WINDOW) or (i > 0 and p - spos[i - 1] <= SNV_WINDOW)

        npos = [p for p in sorted(snp) if near_site(p)]
        mol_site = defaultdict(lambda: defaultdict(set))  # (cb, ub, gene) -> site_idx -> {base}
        mol_snv = defaultdict(lambda: defaultdict(set))   # (cb, ub, gene) -> snv pos0 -> {base}
        for r in bam.fetch(chrom):
            if r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_unmapped:
                continue
            if r.mapping_quality < a.min_q or (r.has_tag("NH") and r.get_tag("NH") != 1):
                continue
            if not (r.has_tag("CB") and r.has_tag("UB") and r.has_tag("GX")):
                continue
            cb = r.get_tag("CB")
            if cb not in cells:
                continue
            gene = r.get_tag("GX")
            if gene == "-":
                continue
            hits_s, hits_n = [], []
            for b0, b1 in r.get_blocks():
                i = bisect.bisect_left(spos, b0)
                while i < len(spos) and spos[i] < b1:
                    hits_s.append(i)
                    i += 1
                j = bisect.bisect_left(npos, b0)
                while j < len(npos) and npos[j] < b1:
                    hits_n.append(npos[j])
                    j += 1
            if not hits_s and not hits_n:
                continue
            key = (cb, r.get_tag("UB"), gene)
            q = r.query_qualities
            seq = r.query_sequence
            qmap = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
            for p in hits_n:
                qp = qmap.get(p)
                if qp is not None and q[qp] >= a.min_q and seq[qp].upper() in snp[p]:
                    mol_snv[key][p].add(seq[qp].upper())
            for i in hits_s:
                p0, _, ref, alt = sl[i]
                qp = qmap.get(p0)
                if qp is not None and q[qp] >= a.min_q and seq[qp].upper() in (ref, alt):
                    mol_site[key][i].add(seq[qp].upper())
        counts = defaultdict(lambda: [0, 0])  # (cb, site_idx, snv, allele) -> [unedited, edited]
        for key, sd in mol_site.items():
            nd = mol_snv.get(key)
            if not nd:
                continue
            for i, sb in sd.items():
                if len(sb) != 1:
                    continue
                edited = int(next(iter(sb)) == sl[i][3])
                for p, nb in nd.items():
                    if len(nb) != 1:
                        continue
                    allele = "ref" if next(iter(nb)) == snp[p][0] else "alt"
                    counts[(key[0], i, p, allele)][edited] += 1
        for (cb, i, p, allele), (u, e) in counts.items():
            out.write(f"{a.sample}\t{cb}\t{sl[i][1]}\t{chrom}:{p + 1}\t{allele}\t{u}\t{e}\n")
            n_rows += 1
        print(f"{a.sample}: {chrom} done, {len(mol_site)} site molecules, {n_rows} rows", file=sys.stderr, flush=True)
    out.close()
    print(f"{a.sample}: done; {n_rows} site x SNV allele rows", file=sys.stderr)


if __name__ == "__main__":
    main()
