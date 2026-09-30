#!/usr/bin/env python
"""scAllele VCF -> SNV table for SCALLELE_CHROM (bin/scallele_chrom.sh).

Usage: scallele_vcf.py VCF MIN_COV [READ_CLUSTERS]
  -> TSV on stdout (PASS + LowQual, SNVs only):
     chr pos ref alt cov_ref cov_alt cov filter varType_tag strand
  READ_CLUSTERS = scAllele's <prefix>.read_cluster_info.tab (bam_file index chrom
  strand start end max_coverage). With it, `strand` is the read-cluster strand of
  EVERY SNV; without it, only RNA_EDITING? rows carry a strand (legacy behaviour).

Read by bin/downstream_filter_scallele.R. Kept separate from scallele_io.py
because Nextflow hashes the content of bin/ scripts named in a task script, and
scallele_io.py is also called by REDITOOLS3 and JACUSA2_CALL1.
"""
from __future__ import print_function

import sys


# scAllele VCF FILTER is BINARY: write_vcf_file() collapses every internal
# category (LowAC/LowDP/LowAB/LowBaseQual/LowReadPos/Overploidy/LowQual) to
# "LowQual", and emits "PASS" only if the hard filters pass and the GLM gives
# QUAL >= 3.01 (gqv_glm.update_vars).
# PASS is NOT used to select editing calls: the GLM was trained on genetic
# variants with RNA_EDITING? candidates removed from training
# (gqv_glm.process_dataframe), so filtering on PASS would score editing with a
# classifier that never saw editing. PASS is not "genetic variant" either
# (well-supported edits do PASS), so it is not used as a SNP blacklist.
# Every SNV (PASS + LowQual) is kept; FILTER is recorded, never used.
#
# varType_tag = scAllele's INFO varType. RNA_EDITING? is a deterministic label
# (gqv_utils): A>G in a '+' read cluster or T>C in a '-' read cluster, the
# cluster strand coming from the reads via --strandedness. The VCF itself only
# implies that strand for tagged rows ('+' for A>G, '-' for T>C).
# For untagged rows the strand comes from the read-cluster table: FORMAT RCL is
# the read-cluster index, and the SNV lies inside that cluster's span. The table
# is 0-based and the VCF 1-based, so the match is (index, start <= pos-1 <= end);
# indices repeat within a table, the span makes the match unique.
# Validated on Replogle chr7 / chr10 (2026-09-29): every SNV matched exactly one
# cluster, and the cluster strand reproduced all 94,747 RNA_EDITING? tags with no
# extra ones. Tagged rows keep the tag-implied strand (identical to before); a
# tagged row whose cluster strand contradicts its tag, or that has no cluster,
# is a hard error (exit 3), so no editing call can be silently re-stranded or lost.
SCALLELE_KEEP_FILTERS = frozenset(("PASS", "LowQual", ".", ""))
EDIT_TAG = "RNA_EDITING?"


def read_clusters(path):
    """index -> list of (start, end, strand) from scAllele's read_cluster_info.tab."""
    clusters = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(header)}
        for k in ("index", "strand", "start", "end"):
            if k not in col:
                sys.exit("scallele_vcf.py: read-cluster table lacks column '{}'".format(k))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            clusters.setdefault(int(f[col["index"]]), []).append(
                (int(f[col["start"]]), int(f[col["end"]]), f[col["strand"]]))
    return clusters


def cluster_strand(clusters, rcl, pos0):
    hits = set(s for (a, b, s) in clusters.get(rcl, ()) if a <= pos0 <= b)
    if len(hits) == 1:
        return hits.pop()
    return None if not hits else "."


def cmd_vcf(vcf_path, min_cov, clusters_path=None, keep_filters=SCALLELE_KEEP_FILTERS):
    """SNVs only. FORMAT GT:GQ:DP:AC:RC:AB:RCL (first sample column).

    scAllele deletes its own VCF header (rm -f prefix.* after write_vcf_header),
    so the file is header-less data rows. We skip any '#'/blank lines anyway.
    Columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1[...]
    """
    min_cov = int(min_cov)
    bases = set("ACGT")
    clusters = read_clusters(clusters_path) if clusters_path else None
    n_out = n_unmatched = 0
    with open(vcf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, _vid, ref, alt, _qual, filt, info = f[:8]
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
            tag = "."
            for kv in info.split(";"):
                if kv.startswith("varType="):
                    tag = kv[len("varType="):] or "."
                    break
            if tag == EDIT_TAG and (ref, alt) == ("A", "G"):
                strand = "+"
            elif tag == EDIT_TAG and (ref, alt) == ("T", "C"):
                strand = "-"
            else:
                strand = "."
            if clusters is not None:
                cs = cluster_strand(clusters, _int("RCL", -1), int(pos) - 1)
                if tag == EDIT_TAG:
                    if cs is None or cs not in (strand, "."):
                        sys.stderr.write(
                            "scallele_vcf.py: RNA_EDITING? {}:{} {}>{} has read-cluster strand {} "
                            "(RCL {}); expected {}\n".format(chrom, pos, ref, alt, cs, d.get("RCL"), strand))
                        sys.exit(3)
                elif cs is None:
                    n_unmatched += 1
                else:
                    strand = cs
            n_out += 1
            sys.stdout.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    chrom, pos, ref, alt, rc, ac, dp, filt or ".", tag, strand
                )
            )
    if clusters is not None:
        sys.stderr.write("scallele_vcf.py: {} SNVs written; {} untagged SNVs without a read cluster "
                         "(strand '.')\n".format(n_out, n_unmatched))


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    cmd_vcf(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else None)
