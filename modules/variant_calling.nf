/*
 * Variant calling: bcftools, REDItools3, GATK, JACUSA2, Cellsnp-lite, scAllele
 */

// `--flag false` arrives as the string "false", which is truthy in Groovy.
def toBool(v) {
    if (v instanceof Boolean) return v
    return (v?.toString()?.trim()?.toLowerCase() in ['true', '1', 'yes', 'y', 't', 'on'])
}

process BCFTOOLS_MPILEUP {
    tag { "${meta.id}:${bam_label}" }
    label 'process_calling'
    publishDir { "${params.outdir}/${meta.id}/03_variants/bcftools/${bam_label}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), val(bam_label)
    path fasta
    path main_chroms

    output:
    tuple val(meta), val('bcftools'), val(bam_label),
          path("${meta.id}_${bam_label}_bcftools_variants.txt.gz"), emit: variants
    tuple val(meta), val('bcftools'), val(bam_label),
          path("${meta.id}_${bam_label}_bcftools_snpblacklist.txt.gz"), emit: snp_blacklist

    script:
    def qual = (meta.qual_override != null) ? meta.qual_override : params.qual
    def min_cov = params.min_cov
    def baq = toBool(params.bcftools_baq) ? 'baq' : 'nobaq'
    """
    #!/bin/bash
    set -euo pipefail
    # Work-dir-local FASTA .fai BEFORE the xargs workers start. Otherwise each
    # worker's mpileup auto-builds ${fasta}.fai beside the staged symlink, and
    # workers that read a half-written index get REF=N for whole chromosomes.
    if [ ! -f ${fasta}.fai ]; then
        REALFA=\$(readlink -f ${fasta})
        if [ -n "\$REALFA" ] && [ -f "\${REALFA}.fai" ]; then
            cp -f "\${REALFA}.fai" ${fasta}.fai
        else
            samtools faidx ${fasta}
        fi
    fi

    # Determine which chromosomes to process
    if [ "${main_chroms.name}" != "NO_FILE" ] && [ -f "${main_chroms}" ]; then
        cp ${main_chroms} _chroms.txt
    else
        # Fall back to all chromosomes in the BAM header
        samtools view -H ${bam} | grep '@SQ' | cut -f2 | sed 's/SN://' > _chroms.txt
    fi

    # Per-chromosome pileup SCATTERED across ${task.cpus} cores. bcftools mpileup's
    # pileup engine is single-threaded (--threads only helps BCF compression, which
    # here was piped uncompressed -> no-op), so the old serial per-chrom loop ran on
    # ~1 core. We now run one worker per chromosome in parallel; each worker does a
    # SINGLE mpileup that feeds BOTH the variant table and the SNP-blacklist
    # genotyping (was two full pileups). ~cpus x faster + halves the pileup work,
    # keeping bcftools' benchmark runtime honest vs the windowed REDItools3.
    grep -vE '^[[:space:]]*\$' _chroms.txt | xargs -P ${task.cpus} -I{} \\
        bcftools_mpileup_chrom.sh {} ${bam} ${fasta} ${qual} ${min_cov} ${baq}

    # Concatenate per-chromosome outputs in chromosome (file) order.
    > _all_variants.txt
    > _snp_blacklist.txt
    while read CHR; do
        [ -z "\$CHR" ] && continue
        [ -f "_var_\${CHR}.txt" ] && cat "_var_\${CHR}.txt" >> _all_variants.txt
        [ -f "_bl_\${CHR}.txt" ]  && cat "_bl_\${CHR}.txt"  >> _snp_blacklist.txt
    done < _chroms.txt

    gzip -c _all_variants.txt > ${meta.id}_${bam_label}_bcftools_variants.txt.gz
    gzip -c _snp_blacklist.txt > ${meta.id}_${bam_label}_bcftools_snpblacklist.txt.gz
    rm -f _all_variants.txt _snp_blacklist.txt _chroms.txt _var_*.txt _bl_*.txt
    """
}


process REDITOOLS3 {
    tag { "${meta.id}:${bam_label}" }
    label 'process_calling'
    publishDir { "${params.outdir}/${meta.id}/03_variants/reditools3/${bam_label}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), val(bam_label), val(strandedness)
    path fasta
    path main_chroms

    output:
    tuple val(meta), val('reditools3'), val(bam_label),
          path("${meta.id}_${bam_label}_reditools3.txt.gz"), emit: variants

    script:
    def strand_flag = strandedness == 'forward' ? '-s 1' :
                      strandedness == 'reverse' ? '-s 2' : '-s 0'
    def qual = (meta.qual_override != null) ? meta.qual_override : params.qual
    """
    #!/bin/bash
    set -euo pipefail
    # Work-dir-local FASTA .fai BEFORE workers start. REDItools3 3.7's 12 threads
    # otherwise race to auto-create the index at startup → spurious
    # "Reference name chrN not found in FASTA file" KeyError (jzl10x resume10).
    # Nextflow often stages the FASTA as a symlink; pysam looks for ${fasta}.fai
    # beside that name, not beside the real file, so copy the existing .fai in.
    if [ ! -f ${fasta}.fai ]; then
        REALFA=\$(readlink -f ${fasta})
        if [ -n "\$REALFA" ] && [ -f "\${REALFA}.fai" ]; then
            cp -f "\${REALFA}.fai" ${fasta}.fai
        else
            samtools faidx ${fasta}
        fi
    fi

    # Restrict to chromosomes present in BOTH the BAM header AND the reference
    # FASTA, and (when supplied) intersected with --main_chroms — for the SAME
    # main-chromosome scope as BCFTOOLS/SCALLELE (fair cross-caller comparison;
    # downstream_filter.R STEP 2 clips every caller to main_chroms regardless).
    # REDItools3 iterates every BAM @SQ contig (region_args.py reads
    # bam.references, ignoring -B) and looks up its reference sequence, so a BAM
    # contig absent from the FASTA — e.g. Cell Ranger's Ensembl-style
    # 'KI270728.1' vs the analysis-set 'chr16_KI270728v1_random' — is a fatal
    # KeyError. Looping per-chrom with -g over the intersection avoids this.
    if [ "${main_chroms.name}" != "NO_FILE" ] && [ -f "${main_chroms}" ]; then
        cp ${main_chroms} _chroms_req.txt
    else
        scallele_io.py bam_chroms ${bam} _chroms_req.txt
    fi
    scallele_io.py intersect_chroms ${bam} _chroms_req.txt _chroms_bam.txt
    awk 'NR==FNR{f[\$1]=1; next} (\$1 in f)' <(cut -f1 ${fasta}.fai) _chroms_bam.txt > _chroms.txt

    # Match method_scripts/run_reditools3.pbs (REDItools3 3.7):
    #   python -m reditools analyze -r REF -g CHR -w WIN -o OUT -s STRAND -V -q Q -bq Q -t N BAM
    # BAM is positional; there is no -i. One invocation per chromosome (-g),
    # concatenated keeping a single header row (fread-parsed downstream).
    # -w splits each chromosome into WINDOW-sized ranges so the ${task.cpus}
    # threads are actually used (REDItools parallelises over ranges, not within
    # one). Verified byte-identical to no-window output (md5-identical on chr1)
    # while ~5x faster — important for the runtime figures in the benchmark.
    OUT=${meta.id}_${bam_label}_reditools3.txt
    rm -f "\$OUT"
    first=1
    while read -r CHR; do
        [ -z "\$CHR" ] && continue
        python -m reditools analyze \\
            -r ${fasta} \\
            -g "\$CHR" \\
            -w 10000000 \\
            -o _rd_"\$CHR".txt \\
            ${strand_flag} \\
            -V \\
            -q ${qual} -bq ${qual} \\
            -t ${task.cpus} \\
            ${bam}
        if [ -s _rd_"\$CHR".txt ]; then
            if [ "\$first" -eq 1 ]; then
                cat _rd_"\$CHR".txt > "\$OUT"; first=0
            else
                tail -n +2 _rd_"\$CHR".txt >> "\$OUT"
            fi
        fi
        rm -f _rd_"\$CHR".txt
    done < _chroms.txt
    # Guarantee a header even if no chromosome produced sites.
    if [ ! -f "\$OUT" ]; then
        printf 'Region\\tPosition\\tReference\\tStrand\\tCoverage\\tMeanQ\\tBaseCount[A,C,G,T]\\tAllSubs\\tFrequency\\tgCoverage\\tgMeanQ\\tgBaseCount[A,C,G,T]\\tgAllSubs\\tgFrequency\\n' > "\$OUT"
    fi

    gzip "\$OUT"
    rm -f _chroms.txt _chroms_req.txt _chroms_bam.txt
    """
}


process GATK_HAPLOTYPECALLER {
    tag { "${meta.id}:${bam_label}" }
    label 'process_calling'
    errorStrategy { (task.exitStatus == null || task.exitStatus in [137, 140, 143, 104, 71]) && task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries    2
    publishDir { "${params.outdir}/${meta.id}/03_variants/gatk/${bam_label}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), val(bam_label)
    path fasta
    path main_chroms

    output:
    tuple val(meta), val('gatk'), val(bam_label),
          path("${meta.id}_${bam_label}_gatk_variants.txt.gz"), emit: variants
    tuple val(meta), val('gatk'), val(bam_label),
          path("${meta.id}_${bam_label}_gatk_snpblacklist.txt.gz"), emit: snp_blacklist

    script:
    // Soft clips off and mapping / base quality = QUAL, plus 8 read filters disabled
    // on non-recalibrated BAMs. recaldefault: GATK defaults on the recal BAM.
    def is_default = (bam_label == 'recaldefault')
    def is_recal   = bam_label.startsWith('recal') && !is_default
    def qual = (meta.qual_override != null) ? meta.qual_override : params.qual
    def gatk_filters
    if (is_default) {
        gatk_filters = "--dont-use-soft-clipped-bases true"
    } else if (is_recal) {
        gatk_filters = """--dont-use-soft-clipped-bases true \\
        --minimum-mapping-quality ${qual} \\
        --min-base-quality-score ${qual}"""
    } else {
        gatk_filters = """--dont-use-soft-clipped-bases true \\
        --disable-read-filter GoodCigarReadFilter \\
        --disable-read-filter MappedReadFilter \\
        --disable-read-filter MappingQualityAvailableReadFilter \\
        --disable-read-filter NonZeroReferenceLengthAlignmentReadFilter \\
        --disable-read-filter NotDuplicateReadFilter \\
        --disable-read-filter NotSecondaryAlignmentReadFilter \\
        --disable-read-filter PassesVendorQualityCheckReadFilter \\
        --disable-read-filter WellformedReadFilter \\
        --minimum-mapping-quality ${qual} \\
        --min-base-quality-score ${qual}"""
    }
    def min_cov = params.min_cov
    // HaplotypeCaller is single-threaded: one JVM per main chromosome across all cpus.
    // Heaps come from the enforced cgroup (~4 GB per cpu), 75% split across workers;
    // a retry uses fewer, bigger JVMs at the same allocation.
    def scatter_n  = Math.max(1, (task.cpus / task.attempt).intValue())
    def hc_heap_gb = Math.max(1, ((task.cpus * 4 * 0.75) / scatter_n).intValue())
    // One line, so it can be echoed to the workers' argv file (a heredoc would be
    // broken by Nextflow's script indent stripping).
    def filter_tokens = gatk_filters.replaceAll(/\\\s*\n\s*/, ' ').trim()
    """
    #!/bin/bash
    set -euo pipefail

    # Ensure FASTA dict and fai exist
    if [ ! -f ${fasta}.fai ]; then samtools faidx ${fasta}; fi
    DICT=\$(echo ${fasta} | sed 's/.fa\$/.dict/' | sed 's/.fasta\$/.dict/')
    if [ ! -f \$DICT ]; then gatk CreateSequenceDictionary -R ${fasta}; fi

    # Restrict traversal to --main_chroms (same main-chromosome scope as the
    # other callers; downstream_filter.R STEP 2 clips to main_chroms anyway).
    # One contig name per line = one scatter unit. Falls back to the full
    # reference when --main_chroms is absent.
    if [ "${main_chroms.name}" != "NO_FILE" ] && [ -f "${main_chroms}" ]; then
        grep -vE '^[[:space:]]*\$' ${main_chroms} > _chroms.txt
    else
        cut -f1 ${fasta}.fai > _chroms.txt
    fi

    # Original jzl10x: rewrite as one sample RG before HaplotypeCaller. This is
    # load-bearing, not cosmetic: the merged pseudobulk BAMs carry one RG per
    # source cell/plate (267 RGs with 267 distinct SMs for the Replogle NTC arm),
    # which would otherwise make HaplotypeCaller emit a 267-sample VCF and break
    # the single-sample [%AD]/[%DP] query below.
    gatk --java-options "-Xmx${(task.cpus * 4 * 0.75).intValue()}G" AddOrReplaceReadGroups \\
        -I ${bam} \\
        -O _rg.bam \\
        -RGID ${meta.id} \\
        -RGLB lib1 \\
        -RGPL ILLUMINA \\
        -RGPU unit1 \\
        -RGSM ${meta.id}
    samtools index _rg.bam

    # Hand the caller's filter flags to the workers as a file; each worker splits it
    # back into argv (they cannot be passed through xargs as one token).
    echo '${filter_tokens}' > _gatk_filters.txt

    # Scatter the call across chromosomes; the worker fails the task if any single
    # chromosome fails, so a dropped contig can never silently shrink the results.
    xargs -P ${scatter_n} -I{} \\
        gatk_haplotypecaller_chrom.sh {} _rg.bam ${fasta} ${hc_heap_gb} _gatk_filters.txt < _chroms.txt

    # Concatenate in reference (.fai) order so the merged VCF is ordered exactly as
    # a single unscattered traversal would have emitted it.
    awk 'NR==FNR {want[\$1]=1; next} want[\$1] {print "_raw_"\$1".vcf.gz"}' \\
        _chroms.txt ${fasta}.fai > _vcf_list.txt
    bcftools concat -Oz -o _raw.vcf.gz -f _vcf_list.txt

    # Parse VCF to TSV: chr, pos, ref, alt, cov_ref, cov_alt, cov, sb (strand bias counts)
    # FORMAT is GT:AD:DP:GQ:PL:SB — SB = ref_fwd,ref_rev,alt_fwd,alt_rev
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\t[%DP]\\t[%SB]\\n' _raw.vcf.gz | \\
    awk -F'\\t' 'BEGIN {OFS=FS} {
        split(\$5, ad, ","); ref_ct=ad[1]; alt_ct=ad[2]; cov=\$6;
        if (cov >= ${min_cov}) print \$1, \$2, \$3, \$4, ref_ct, alt_ct, cov, \$7
    }' | \\
    gzip > ${meta.id}_${bam_label}_gatk_variants.txt.gz

    # Per-sample SNP blacklist (R2.M2/M4): HaplotypeCaller already genotypes,
    # so keep only homozygous-ALT (1/1) SNVs at DP>=min_cov. Editing-safe --
    # A-to-I never reaches clean homozygous ALT in pseudobulk; het (0/1) is NOT
    # included (indistinguishable from ~50% editing).
    bcftools view -i 'GT="AA"' _raw.vcf.gz 2>/dev/null | \\
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%DP]\\n' 2>/dev/null | \\
    awk -F'\\t' -v mc=${min_cov} 'length(\$3)==1 && length(\$4)==1 && (\$5+0)>=mc {print \$1"\\t"\$2}' | \\
    gzip > ${meta.id}_${bam_label}_gatk_snpblacklist.txt.gz

    rm -f _rg.bam _rg.bam.bai _raw.vcf.gz _raw.vcf.gz.tbi \\
          _raw_*.vcf.gz _raw_*.vcf.gz.tbi _chroms.txt _vcf_list.txt _gatk_filters.txt
    """
}


process JACUSA2_CALL1 {
    tag { "${meta.id}:${bam_label}" }
    label 'process_calling'
    // -Xmx is sized from the enforced cgroup (~8 GB per 2-thread slot), not from the
    // unenforced `memory` request: a heap above the cap is killed before it collects.
    // cpus / memory / time come from withLabel process_calling (benchmark parity).
    errorStrategy { (task.exitStatus == null || task.exitStatus in [137, 140, 143, 104, 71]) ? 'retry' : 'finish' }
    maxRetries    3
    publishDir { "${params.outdir}/${meta.id}/03_variants/jacusa2/${bam_label}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), val(bam_label), val(strandedness)
    path fasta
    path main_chroms

    output:
    tuple val(meta), val('jacusa2'), val(bam_label),
          path("${meta.id}_${bam_label}_jacusa2.out.gz"), emit: variants

    script:
    // -P library type (strandedness). Each library uses the strand RSeQC
    // inferred: forward=FR-SECONDSTRAND (10x 3'), reverse=FR-FIRSTSTRAND.
    def strand_flag = strandedness == 'forward' ? 'FR-SECONDSTRAND' :
                      strandedness == 'reverse' ? 'FR-FIRSTSTRAND' : 'UNSTRANDED'
    def qual = (meta.qual_override != null) ? meta.qual_override : params.qual
    def min_cov = params.min_cov
    """
    #!/bin/bash
    set -euo pipefail

    # Resolve the JACUSA2 jar from the active conda env (set -u safe); the
    # \$JACUSA2_JAR env var is not exported by the conda package.
    JACUSA2_JAR=\$(ls \${CONDA_PREFIX:-}/share/jacusa2-*/JACUSA_v*.jar 2>/dev/null | head -n1)
    if [ -z "\$JACUSA2_JAR" ]; then
        echo "ERROR: no JACUSA_v*.jar under \${CONDA_PREFIX:-<no conda env>}/share/jacusa2-*/ (install jacusa2 into the env)" >&2
        exit 1
    fi

    # JACUSA2 needs an .fai-indexed reference (staged FASTA is a symlink).
    [ -f ${fasta}.fai ] || samtools faidx ${fasta}

    # Restrict the scan to chromosomes present in the reference FASTA (and, when
    # supplied, --main_chroms) via a BED -- the SAME main-chromosome scope as
    # BCFTOOLS/SCALLELE (fair cross-caller comparison; downstream_filter.R
    # STEP 2 clips every caller to main_chroms regardless). JACUSA2 initialises
    # every BAM @SQ contig at startup and fetches its reference sequence, so a
    # BAM contig absent from the FASTA — e.g. Cell Ranger's Ensembl-style
    # 'KI270728.1' vs the analysis-set 'chr16_KI270728v1_random' — is a fatal
    # 'Unable to find entry for contig' SAMException. The -b BED confines the
    # scan to reference contigs and avoids this.
    if [ "${main_chroms.name}" != "NO_FILE" ] && [ -f "${main_chroms}" ]; then
        awk 'NR==FNR{L[\$1]=\$2; next} (\$1 in L){print \$1"\\t0\\t"L[\$1]}' \\
            ${fasta}.fai ${main_chroms} > _main.bed
    else
        # No --main_chroms: fall back to every reference contig also present in
        # the BAM header (still avoids FASTA-missing contigs).
        scallele_io.py bam_chroms ${bam} _bamchroms.txt
        awk 'NR==FNR{L[\$1]=\$2; next} (\$1 in L){print \$1"\\t0\\t"L[\$1]}' \\
            ${fasta}.fai _bamchroms.txt > _main.bed
    fi

    # Match method_scripts/run_jacusa2_call1.pbs (JACUSA2 call-1):
    #   call-1 -R REF -c 10 -m Q -q Q -p N -P STRAND -r OUT BAM
    # -P sets the library type; -m min mapping quality, -q min base quality,
    # -c min coverage; -b confines the scan to the main-chromosome BED. BAM is
    # positional.
    # -Xmx is sized from the ENFORCED cgroup (~4 GB per allocated cpu = 48 GB at the
    # standard 12), NOT from the unenforced `memory` request: the old 0.8*64 = 51 GB
    # put the heap ceiling above the cap and let the deepest BAM grow into the OOM
    # kill. 75% leaves room for non-heap JVM overhead (metaspace, stacks, GC).
    java -Xmx${(task.cpus * 4 * 0.75).intValue()}G -jar "\$JACUSA2_JAR" call-1 \\
        -R ${fasta} \\
        -c ${min_cov} \\
        -m ${qual} \\
        -q ${qual} \\
        -p ${task.cpus} \\
        -P ${strand_flag} \\
        -b _main.bed \\
        -r ${meta.id}_${bam_label}_jacusa2.out \\
        ${bam}

    gzip ${meta.id}_${bam_label}_jacusa2.out
    rm -f _main.bed _bamchroms.txt
    """
}


/*
 * JACUSA2 call-1 per chromosome + ordered gather (--jacusa2_scatter). A worker-thread
 * OOM hangs the whole-genome JVM; here -XX:+ExitOnOutOfMemoryError (exit 3) makes it a
 * retry. call-1 scores sites independently, so the gathered rows equal the
 * whole-genome rows. Runtime = sum over JACUSA2_CHROM tasks. JACUSA2_CALL1 is kept
 * unchanged for the teratoma runs (task cache).
 */
process JACUSA2_CHROM {
    tag { "${meta.id}:${bam_label}:${chrom}" }

    input:
    tuple val(meta), path(bam), path(bai), val(bam_label), val(strandedness), val(chrom)
    path fasta

    output:
    tuple val(meta), val(bam_label), val(chrom),
          path("${meta.id}__${bam_label}__${chrom}.jacusa2.out.gz"), emit: chrom_out

    script:
    def strand_flag = strandedness == 'forward' ? 'FR-SECONDSTRAND' :
                      strandedness == 'reverse' ? 'FR-FIRSTSTRAND' : 'UNSTRANDED'
    def qual = (meta.qual_override != null) ? meta.qual_override : params.qual
    def min_cov = params.min_cov
    def pref = "${meta.id}__${bam_label}__${chrom}"
    // memory and -l mem are set equal in nextflow.config, so task.memory is the enforced cap;
    // 75% leaves room for non-heap JVM overhead.
    def heap_gb = Math.max(1, (task.memory.toGiga() * 0.75).intValue())
    """
    #!/bin/bash
    set -euo pipefail

    JACUSA2_JAR=\$(ls \${CONDA_PREFIX:-}/share/jacusa2-*/JACUSA_v*.jar 2>/dev/null | head -n1)
    if [ -z "\$JACUSA2_JAR" ]; then
        echo "ERROR: no JACUSA_v*.jar under \${CONDA_PREFIX:-<no conda env>}/share/jacusa2-*/ (install jacusa2 into the env)" >&2
        exit 1
    fi
    [ -f ${fasta}.fai ] || samtools faidx ${fasta}

    # One-contig BED. A contig missing from the reference or the BAM header has
    # nothing to call: header-only output, as JACUSA2 writes for an empty contig.
    awk -v c=${chrom} '\$1 == c { print \$1"\\t0\\t"\$2 }' ${fasta}.fai > _chrom.bed
    if [ ! -s _chrom.bed ] || ! samtools view -H ${bam} | awk -v c=${chrom} '\$1 == "@SQ" && \$2 == "SN:"c { f = 1 } END { exit !f }'; then
        printf '#contig\\tstart\\tend\\tname\\tscore\\tstrand\\tbases11\\tinfo\\tfilter\\tref\\n' | gzip -c > ${pref}.jacusa2.out.gz
        echo "JACUSA2 ${pref}: ${chrom} not in reference / BAM header"
        exit 0
    fi

    java -Xmx${heap_gb}G -XX:+ExitOnOutOfMemoryError -jar "\$JACUSA2_JAR" call-1 \\
        -R ${fasta} \\
        -c ${min_cov} \\
        -m ${qual} \\
        -q ${qual} \\
        -p ${task.cpus} \\
        -P ${strand_flag} \\
        -b _chrom.bed \\
        -r _jacusa2.out \\
        ${bam}

    gzip -c _jacusa2.out > ${pref}.jacusa2.out.gz
    echo "JACUSA2 ${pref}: \$(grep -vc '^#' _jacusa2.out || true) rows"
    rm -f _jacusa2.out _chrom.bed
    """
}

process JACUSA2_GATHER {
    tag { "${meta.id}:${bam_label}" }
    label 'process_low'
    publishDir { "${params.outdir}/${meta.id}/03_variants/jacusa2/${bam_label}" }, mode: 'copy'

    input:
    tuple val(meta), val(bam_label), val(chroms), path(outs)
    val chrom_order

    output:
    tuple val(meta), val('jacusa2'), val(bam_label),
          path("${meta.id}_${bam_label}_jacusa2.out.gz"), emit: variants

    script:
    def order = chrom_order.join(' ')
    def pref  = "${meta.id}__${bam_label}__"
    """
    #!/bin/bash
    set -euo pipefail
    : > _banner.txt
    : > _rows.txt
    for CHR in ${order}; do
        f="${pref}\${CHR}.jacusa2.out.gz"
        [ -f "\$f" ] || { echo "ERROR: missing \$f" >&2; exit 1; }
        [ -s _banner.txt ] || { gzip -cd "\$f" | grep '^##' > _banner.txt || true; }
        gzip -cd "\$f" | { grep -v '^#' || true; } >> _rows.txt
    done
    { cat _banner.txt
      printf '#contig\\tstart\\tend\\tname\\tscore\\tstrand\\tbases11\\tinfo\\tfilter\\tref\\n'
      cat _rows.txt; } | gzip -c > ${meta.id}_${bam_label}_jacusa2.out.gz
    echo "JACUSA2_GATHER ${meta.id}:${bam_label}: \$(wc -l < _rows.txt) rows"
    rm -f _banner.txt _rows.txt
    """
}


/*
 * Cellsnp-lite mode 2b: bulk pileup of the pseudobulk BAM (--cellTAG None, --UMItag
 * None). No strand output (annotation strand mode only) and no base-quality filter
 * (--minMAPQ = QUAL). Output: cov_ref = DP - AD, cov_alt = AD, cov = DP + OTH; rows
 * with AD >= 1. Runs params.cellsnp_bin from its own env so the pipeline's samtools
 * stays on PATH.
 */
process CELLSNP_LITE {
    tag { "${meta.id}:${bam_label}" }
    label 'process_calling'
    publishDir { "${params.outdir}/${meta.id}/03_variants/cellsnp/${bam_label}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), val(bam_label)
    path fasta
    path main_chroms

    output:
    tuple val(meta), val('cellsnp'), val(bam_label),
          path("${meta.id}_${bam_label}_cellsnp_variants.txt.gz"), emit: variants

    script:
    def qual = (meta.qual_override != null) ? meta.qual_override : params.qual
    def min_cov = params.min_cov
    """
    #!/bin/bash
    set -euo pipefail
    if [ ! -f ${fasta}.fai ]; then
        REALFA=\$(readlink -f ${fasta})
        if [ -n "\$REALFA" ] && [ -f "\${REALFA}.fai" ]; then
            cp -f "\${REALFA}.fai" ${fasta}.fai
        else
            samtools faidx ${fasta}
        fi
    fi

    # Same main-chromosome scope as the other callers: --main_chroms (or every
    # BAM contig) intersected with the BAM header and the reference.
    if [ "${main_chroms.name}" != "NO_FILE" ] && [ -f "${main_chroms}" ]; then
        grep -vE '^[[:space:]]*\$' ${main_chroms} > _chroms_req.txt
    else
        samtools view -H ${bam} | awk -F'\\t' '\$1=="@SQ"{sub(/^SN:/,"",\$2); print \$2}' > _chroms_req.txt
    fi
    samtools view -H ${bam} | awk -F'\\t' '\$1=="@SQ"{sub(/^SN:/,"",\$2); print \$2}' > _chroms_bam.txt
    awk 'NR==FNR{f[\$1]=1; next} (\$1 in f)' _chroms_bam.txt _chroms_req.txt \\
        | awk 'NR==FNR{f[\$1]=1; next} (\$1 in f)' <(cut -f1 ${fasta}.fai) - > _chroms.txt
    CHROMS=\$(paste -sd, _chroms.txt)

    ${params.cellsnp_bin} \\
        -s ${bam} \\
        -I ${meta.id} \\
        -O _csnp \\
        --chrom "\$CHROMS" \\
        --cellTAG None \\
        --UMItag None \\
        -f ${fasta} \\
        --minCOUNT ${min_cov} \\
        --minMAF 0 \\
        --minMAPQ ${qual} \\
        -p ${task.cpus} \\
        --gzip

    gzip -cd _csnp/cellSNP.base.vcf.gz | awk -F'\\t' 'BEGIN{OFS="\\t"} !/^#/ {
        ad = dp = oth = 0; n = split(\$8, kv, ";")
        for (i = 1; i <= n; i++) { split(kv[i], p, "=")
            if (p[1] == "AD") ad = p[2] + 0; else if (p[1] == "DP") dp = p[2] + 0; else if (p[1] == "OTH") oth = p[2] + 0 }
        r = toupper(\$4); a = toupper(\$5)
        if (ad >= 1 && r ~ /^[ACGT]\$/ && a ~ /^[ACGT]\$/ && r != a) print \$1, \$2, r, a, dp - ad, ad, dp + oth
    }' | gzip > ${meta.id}_${bam_label}_cellsnp_variants.txt.gz

    rm -rf _csnp _chroms.txt _chroms_req.txt _chroms_bam.txt
    """
}


/*
 * scAllele per chromosome + ordered gather. Runtime = sum over SCALLELE_CHROM tasks
 * (human/scripts/nf_trace_scallele_runtime.R), not measured under process_calling.
 * All SNVs (PASS + LowQual) are kept: scAllele's GLM was trained on genetic variants
 * with editing candidates excluded, so PASS is not used. strand_mode 'caller' uses
 * the read-cluster strand (bin/scallele_vcf.py). No per-sample SNP blacklist. The raw
 * VCF and read-cluster table are published so parser changes need no rerun.
 */
process SCALLELE_CHROM {
    tag { "${meta.id}:${bam_label}:${chrom}" }
    label 'process_scallele_chrom'
    // Dedicated env: scAllele 0.0.9.7 needs scipy 1.15 / numpy 1.26;
    // scRNAEditing has scipy 1.17 / numpy 2.x and must not mix.
    conda params.scallele_conda

    input:
    tuple val(meta), path(bam), path(bai), val(bam_label), val(strandedness), val(chrom)
    path fasta

    output:
    tuple val(meta), val(bam_label), val(chrom),
          path("${meta.id}__${bam_label}__${chrom}.snvs.txt"),
          path("${meta.id}__${bam_label}__${chrom}.vcf.gz"),
          path("${meta.id}__${bam_label}__${chrom}.clusters.tsv.gz"), emit: chrom_out

    script:
    // Library type matches JACUSA2_CALL1 / RSeQC:
    //   10x 3' forward = second-strand; reverse = first-strand.
    def strand_flag = strandedness == 'forward' ? 'fr-secondstrand' :
                      strandedness == 'reverse' ? 'fr-firststrand' : 'fr-unstrand'
    def qual = (meta.qual_override != null) ? meta.qual_override : params.qual
    def pref = "${meta.id}__${bam_label}__${chrom}"
    """
    #!/bin/bash
    set -euo pipefail

    # Work-dir-local FASTA .fai (staged FASTA is often a symlink).
    if [ ! -f ${fasta}.fai ]; then
        REALFA=\$(readlink -f ${fasta})
        if [ -n "\$REALFA" ] && [ -f "\${REALFA}.fai" ]; then
            cp -f "\${REALFA}.fai" ${fasta}.fai
        else
            scallele_io.py faidx ${fasta}
        fi
    fi

    # A contig absent from / unmapped in this BAM has nothing to call.
    N=\$(scallele_io.py mapped_reads ${bam} ${chrom})
    if [ "\$N" -eq 0 ]; then
        : > ${pref}.snvs.txt; : | gzip -c > ${pref}.vcf.gz; : | gzip -c > ${pref}.clusters.tsv.gz
        echo "scAllele ${pref}: no reads on ${chrom}"
        exit 0
    fi

    # SCALLELE_STRICT=1: a scAllele failure fails the task (-> retry with more
    # memory) instead of silently yielding an empty contig.
    SCALLELE_STRICT=1 scallele_chrom.sh ${chrom} ${bam} ${fasta} ${params.min_cov} ${params.min_alt} \\
        ${params.min_edit_level} ${qual} ${strand_flag} ${params.scallele_run_mode}
    mv _snvs_${chrom}.txt ${pref}.snvs.txt
    mv _vcf_${chrom}.vcf.gz ${pref}.vcf.gz
    mv _rc_${chrom}.tsv.gz ${pref}.clusters.tsv.gz
    echo "scAllele ${pref}: \$(wc -l < ${pref}.snvs.txt) SNVs, \$(awk -F'\\t' '\$9 == "RNA_EDITING?"' ${pref}.snvs.txt | wc -l) RNA_EDITING?"
    """
}

process SCALLELE_GATHER {
    tag { "${meta.id}:${bam_label}" }
    label 'process_low'
    publishDir { "${params.outdir}/${meta.id}/03_variants/scallele/${bam_label}" }, mode: 'copy'

    input:
    tuple val(meta), val(bam_label), val(chroms), path(snvs), path(vcfs), path(clusters)
    val chrom_order

    output:
    tuple val(meta), val('scallele'), val(bam_label),
          path("${meta.id}_${bam_label}_scallele_variants.txt.gz"), emit: variants
    tuple val(meta), val(bam_label),
          path("${meta.id}_${bam_label}_scallele.vcf.gz"), emit: vcf
    tuple val(meta), val(bam_label),
          path("${meta.id}_${bam_label}_scallele_read_clusters.tsv.gz"), emit: clusters

    script:
    def order = chrom_order.join(' ')
    def pref  = "${meta.id}__${bam_label}__"
    """
    #!/bin/bash
    set -euo pipefail
    : > _all_snvs.txt
    : > _all.vcf
    : > _all_rc.tsv
    for CHR in ${order}; do
        s="${pref}\${CHR}.snvs.txt"; v="${pref}\${CHR}.vcf.gz"; r="${pref}\${CHR}.clusters.tsv.gz"
        [ -f "\$s" ] || { echo "ERROR: missing \$s" >&2; exit 1; }
        [ -f "\$v" ] || { echo "ERROR: missing \$v" >&2; exit 1; }
        [ -f "\$r" ] || { echo "ERROR: missing \$r" >&2; exit 1; }
        cat "\$s" >> _all_snvs.txt
        gzip -cd "\$v" >> _all.vcf
        # one header (bam_file index chrom strand start end max_coverage), then rows
        if [ -s _all_rc.tsv ]; then gzip -cd "\$r" | tail -n +2 >> _all_rc.tsv; else gzip -cd "\$r" >> _all_rc.tsv; fi
    done
    gzip -c _all_snvs.txt > ${meta.id}_${bam_label}_scallele_variants.txt.gz
    gzip -c _all.vcf > ${meta.id}_${bam_label}_scallele.vcf.gz
    gzip -c _all_rc.tsv > ${meta.id}_${bam_label}_scallele_read_clusters.tsv.gz
    rm -f _all_snvs.txt _all.vcf _all_rc.tsv
    """
}
