#!/bin/bash
# obim_RefSeq.bed (--refseq_bed): BED12, one row per transcript of obim.ncbiRefSeq.gtf (from 01), sorted by
# chromosome and start; thickStart / thickEnd = CDS span (= txEnd when non-coding), as in the human and mouse BEDs.
# The original study's obim_RefSeq.bed was gtf2bed output (one row per GTF record), not BED12.
#
# Usage: REFDIR=<dir> bash 03_refseq_bed.sh
set -euo pipefail
REFDIR="${REFDIR:-$PWD}"
cd "$REFDIR"
[ -f obim.ncbiRefSeq.gtf ] || { echo "ERROR: run 01_download_genome_and_gtf.sh first" >&2; exit 1; }

awk 'BEGIN {FS = "\t"; OFS = "\t"}
    $3 != "exon" && $3 != "CDS" {next}
    {
        if (!match($9, /transcript_id "[^"]+"/)) next
        t = substr($9, RSTART + 15, RLENGTH - 16)
        s = $4 - 1; e = $5
        if ($3 == "CDS") {
            if (!(t in cs) || s < cs[t]) cs[t] = s
            if (!(t in ce) || e > ce[t]) ce[t] = e
            next
        }
        chr[t] = $1; str[t] = $7; n[t]++; xs[t, n[t]] = s; xe[t, n[t]] = e
    }
    END {
        for (t in n) {
            k = n[t]
            for (i = 1; i <= k; i++) { a[i] = xs[t, i]; b[i] = xe[t, i] }
            for (i = 2; i <= k; i++) {
                x = a[i]; y = b[i]
                for (j = i - 1; j >= 1 && a[j] > x; j--) { a[j + 1] = a[j]; b[j + 1] = b[j] }
                a[j + 1] = x; b[j + 1] = y
            }
            ts = a[1]; te = b[k]; sz = ""; st = ""
            for (i = 1; i <= k; i++) { sz = sz (b[i] - a[i]) ","; st = st (a[i] - ts) ","; if (b[i] > te) te = b[i] }
            th = (t in cs) ? cs[t] : te; tt = (t in ce) ? ce[t] : te
            print chr[t], ts, te, t, 0, str[t], th, tt, 0, k, sz, st
        }
    }' obim.ncbiRefSeq.gtf | LC_ALL=C sort -k1,1 -k2,2n > obim_RefSeq.bed
echo "$(wc -l < obim_RefSeq.bed) transcripts -> $REFDIR/obim_RefSeq.bed"
