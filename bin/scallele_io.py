#!/usr/bin/env python
"""Helpers for the Nextflow SCALLELE process (scAllele env: pysam, no bcftools).

Modes:
  faidx FASTA
  bam_chroms BAM OUT
  intersect_chroms BAM REQ OUT
  mapped_reads BAM CHR   -> prints mapped read count on CHR (0 if absent), from the index
  vcf VCF MIN_COV [BLACKLIST_OUT]
        -> TSV on stdout: chr pos ref alt cov_ref cov_alt cov (PASS+LowQual)
        -> if BLACKLIST_OUT given, also append PASS-only SNV sites
           (chr\tpos) there = scAllele's sample-specific genetic-variant set.
"""
from __future__ import print_function

import sys


def cmd_faidx(fasta):
    import pysam
    pysam.faidx(fasta)


def cmd_bam_chroms(bam_path, out_path):
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    with open(out_path, "w") as out:
        out.write("\n".join(bam.references) + "\n")


def cmd_mapped_reads(bam_path, chrom):
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    n = 0
    for s in bam.get_index_statistics():
        if s.contig == chrom:
            n = s.mapped
            break
    print(n)


def cmd_intersect_chroms(bam_path, req_path, out_path):
    import pysam
    req = [ln.strip() for ln in open(req_path) if ln.strip()]
    bam = pysam.AlignmentFile(bam_path, "rb")
    present = set(bam.references)
    keep = [c for c in req if c in present]
    if not keep:
        sys.exit("ERROR: no requested chromosomes present in BAM")
    with open(out_path, "w") as out:
        out.write("\n".join(keep) + "\n")
    print("scAllele chroms: " + ",".join(keep), flush=True)


# scAllele VCF FILTER is BINARY: write_vcf_file() collapses every internal
# category (LowAC/LowDP/LowAB/LowBaseQual/LowReadPos/Overploidy/LowQual) to
# "LowQual", and emits "PASS" only if the GLM passed. So the VCF contains just
# PASS or LowQual.
# CRITICAL (editing discovery): scAllele's shipped GLM was trained with
# RNA-editing sites REMOVED (gqv_glm.process_dataframe drops variant_type
# == 'RNA_EDITING?'), so it scores A-to-I (A>G on +, T>C on -) editing as
# "not a genetic variant" -> QUAL < 3.01 -> FILTER = LowQual, NEVER PASS.
# A PASS-only filter would delete essentially all A-to-I editing. For discovery
# we therefore keep PASS + LowQual (= every candidate SNV scAllele emitted) and
# let downstream_filter.R apply the SAME uniform editing filters used for all
# callers (coverage, alt-count, edit level, SNP removal, strand, A>G
# enrichment). scAllele's contribution is local DBG reassembly (reduced
# reference bias in the ALT allele counts), not its genetic-variant classifier.
SCALLELE_KEEP_FILTERS = frozenset(("PASS", "LowQual", ".", ""))


def cmd_vcf(vcf_path, min_cov, blacklist_out=None,
            keep_filters=SCALLELE_KEEP_FILTERS):
    """SNVs only. FORMAT GT:GQ:DP:AC:RC:AB:RCL (first sample column).

    scAllele deletes its own VCF header (rm -f prefix.* after write_vcf_header),
    so the file is header-less data rows. We skip any '#'/blank lines anyway.
    Columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1[...]

    The main SNV TSV (stdout) keeps PASS + LowQual = every candidate SNV, so
    scAllele is benchmarked fairly like the other callers. If blacklist_out is
    given, PASS-only SNV sites (chr\tpos) are ALSO written there: because the
    GLM was trained WITHOUT editing, a PASS call is scAllele's confident
    genetic-variant (SNP) call, so this is a SAMPLE-SPECIFIC SNP blacklist that
    supplements static dbSNP. downstream_filter.R --snp_blacklist subtracts it
    (editing-safe: A-to-I lands in LowQual, never PASS).
    """
    min_cov = int(min_cov)
    bases = set("ACGT")
    bl_fh = open(blacklist_out, "a") if blacklist_out else None
    try:
        with open(vcf_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                chrom, pos, _vid, ref, alt, _qual, filt = f[:7]
                if len(ref) != 1 or len(alt) != 1:
                    continue
                if ref not in bases or alt not in bases:
                    continue
                if filt not in keep_filters:
                    continue
                keys = f[8].split(":")
                vals = f[9].split(":")
                d = dict(zip(keys, vals))

                def _int(key, default=0):
                    raw = d.get(key, ".")
                    if raw in (".", "", None):
                        return default
                    try:
                        return int(float(raw))
                    except ValueError:
                        return default

                dp = _int("DP")
                ac = _int("AC")
                rc = _int("RC")
                if dp < min_cov:
                    continue
                # PASS = scAllele's genetic-variant call -> blacklist.
                if bl_fh is not None and filt == "PASS":
                    bl_fh.write("{}\t{}\n".format(chrom, pos))
                sys.stdout.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        chrom, pos, ref, alt, rc, ac, dp
                    )
                )
    finally:
        if bl_fh is not None:
            bl_fh.close()


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    mode = argv[1]
    if mode == "faidx":
        cmd_faidx(argv[2])
    elif mode == "bam_chroms":
        cmd_bam_chroms(argv[2], argv[3])
    elif mode == "intersect_chroms":
        cmd_intersect_chroms(argv[2], argv[3], argv[4])
    elif mode == "mapped_reads":
        cmd_mapped_reads(argv[2], argv[3])
    elif mode == "vcf":
        blacklist_out = argv[4] if len(argv) > 4 else None
        cmd_vcf(argv[2], argv[3], blacklist_out)
    else:
        sys.exit("unknown mode: " + mode)


if __name__ == "__main__":
    main(sys.argv)
