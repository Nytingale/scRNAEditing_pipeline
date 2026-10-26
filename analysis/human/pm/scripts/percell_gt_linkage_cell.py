#!/usr/bin/env python3
"""Single-cell linkage of candidate editing sites to heterozygous donor SNVs beyond the same read (R3.9).

Companion of percell_gt_linkage.py (same read filters, sites and SNVs). One pass per per-teratoma GT BAM,
chromosome by chromosome, at UMI level (CB + UB = one molecule; gene = STARsolo GX tag, GeneFull):
  site molecules   UMIs with >= 1 read with the reference or edited base at the site (BQ >= 30)
  molecule-linked  site molecules of which any read (same CB + UB + gene) has a SNV allele (BQ >= 30)
  cell-gene SNVs   per cell and gene: UMIs with a SNV allele, per SNV and allele
pysam is used because per-read checks of pairs of positions have no efficient samtools / R equivalent.

Output <out_prefix>.cell.tsv.gz, one row per cell x site with >= 1 site molecule:
  sample barcode site gene site_umis mol_linked_umis gene_snv_umis gene_snvs gene_snvs_ge2 gene_snvs_mono
gene_snvs_ge2: SNVs of the gene with >= 2 UMIs in the cell; gene_snvs_mono: of those, SNVs with one allele only.
"""
import argparse
import bisect
import gzip
import sys
from collections import Counter, defaultdict

import pysam

from percell_gt_linkage import EDIT_ALT, load_snps


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
    out = gzip.open(f"{a.out_prefix}.cell.tsv.gz", "wt")
    out.write("sample\tbarcode\tsite\tgene\tsite_umis\tmol_linked_umis\tgene_snv_umis\tgene_snvs\tgene_snvs_ge2\t"
              "gene_snvs_mono\n")
    n_rows = 0
    for chrom in sorted(sites):
        if chrom not in contigs:
            continue
        sl = sites[chrom]
        spos = [s[0] for s in sl]
        snp = het.get(chrom, {})
        npos = sorted(snp)
        site_mol = defaultdict(set)                      # (cb, site_idx) -> {(ub, gene)}
        site_gene = defaultdict(Counter)                 # site_idx -> gene counts
        mol_snp = set()                                  # (cb, ub, gene) with a SNV allele
        cg_snp = defaultdict(lambda: defaultdict(dict))  # (cb, gene) -> snp -> {ub: allele}
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
            ub = r.get_tag("UB")
            q = r.query_qualities
            seq = r.query_sequence
            qmap = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
            for p in hits_n:
                qp = qmap.get(p)
                if qp is None or q[qp] < a.min_q:
                    continue
                b = seq[qp].upper()
                if b in snp[p]:
                    mol_snp.add((cb, ub, gene))
                    cg_snp[(cb, gene)][p].setdefault(ub, b)
            for i in hits_s:
                p0, _, ref, alt = sl[i]
                qp = qmap.get(p0)
                if qp is None or q[qp] < a.min_q or seq[qp].upper() not in (ref, alt):
                    continue
                site_mol[(cb, i)].add((ub, gene))
                site_gene[i][gene] += 1
        for (cb, i), mols in site_mol.items():
            gene = site_gene[i].most_common(1)[0][0]
            mols = {ub for ub, g in mols if g == gene}
            if not mols:
                continue
            linked = sum((cb, ub, gene) in mol_snp for ub in mols)
            snvs = cg_snp.get((cb, gene), {})
            umis = sum(len(u) for u in snvs.values())
            ge2 = [u for u in snvs.values() if len(u) >= 2]
            mono = sum(len(set(u.values())) == 1 for u in ge2)
            out.write(f"{a.sample}\t{cb}\t{sl[i][1]}\t{gene}\t{len(mols)}\t{linked}\t{umis}\t{len(snvs)}\t"
                      f"{len(ge2)}\t{mono}\n")
            n_rows += 1
        print(f"{a.sample}: {chrom} done, {n_rows} rows", file=sys.stderr, flush=True)
    out.close()
    print(f"{a.sample}: done; {n_rows} cell x site rows", file=sys.stderr)


if __name__ == "__main__":
    main()
