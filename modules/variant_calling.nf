/*
 * Variant calling modules — bcftools, REDItools3, GATK, JACUSA2, scAllele
 *
 * Each caller is run per-chromosome for parallelism where applicable.
 */

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
    """
    #!/bin/bash
    set -euo pipefail

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
        bcftools_mpileup_chrom.sh {} ${bam} ${fasta} ${qual} ${min_cov}

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
    // The call is scattered across chromosomes (see script), so the memory the task
    // needs scales with the scatter width. The escalating-cpus retry that lifts the
    // cgroup cap lives in nextflow.config (`withName: JACUSA2_CALL1|
    // GATK_HAPLOTYPECALLER`) — a config selector would override a `cpus` set here.
    errorStrategy { task.exitStatus in [137, 140, 143, 104, 71] && task.attempt <= 2 ? 'retry' : 'finish' }
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
    // Match the original GATK design (base/mapping quality = QUAL throughout;
    // q20 lanes set QUAL=20 via qual_override):
    //   non-recal (run_gatk_call_RGbams_nonrecal.pbs): disable 8 read filters
    //     so GATK can run on NON-recalibrated BAMs, + soft-clip off + mq/bq QUAL.
    //   recal / recalq20 (recalibrated BAM): soft-clip off + mq/bq QUAL, no
    //     filter disabling (recalibration makes the default filters usable).
    //   recaldefault (GATK-only supplementary run, run_gatk_call_RGbams_
    //     noRenaming_default.pbs): recal BAM with GATK default settings.
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
    // HaplotypeCaller is single-threaded per invocation (GATK4 dropped -nct), so
    // one whole-genome traversal of a deep pseudobulk BAM ran for days and hit the
    // queue wall limit. Scatter one JVM per main chromosome instead, across ALL
    // allocated cpus — same as BCFTOOLS_MPILEUP, so the per-caller runtime stays a
    // like-for-like comparison at the shared process_calling allocation.
    // Heaps are sized from the ENFORCED cgroup (~4 GB per allocated cpu), not the
    // unenforced `memory` request, with 75% split across the workers. 3 GB/JVM at
    // the standard 12 cpus is ample: HaplotypeCaller holds one active region at a
    // time and downsamples deep pileups (measured peak here was under 1 GB).
    // A retry de-parallelises instead of asking for more slots, so an OOM is met
    // with fewer, bigger JVMs while the allocation stays fixed at 12 cpus.
    def scatter_n  = Math.max(1, (task.cpus / task.attempt).intValue())
    def hc_heap_gb = Math.max(1, ((task.cpus * 4 * 0.75) / scatter_n).intValue())
    // Collapse the filter block's backslash-newline continuations to one line so it
    // can be written to the workers' argv file with a plain echo. (A heredoc would
    // be at the mercy of the indent stripping Nextflow applies to script blocks.)
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
    // JACUSA2's largest merged BAM (Replogle NTC:sample) OOM-killed (exit 137) with
    // epilog MEM_USED_G 48.01 / MEM_LIMIT_G 48.00.
    // CLUSTER MEMORY MODEL (measured): the enforced cap is a cgroup of ~8 GB per
    // SLOT, and the `threaded` PE allocates 2 threads/slot, so 12 threads = 6
    // slots = 48 GB. The Nextflow `memory` directive (emitted as -l h_rss,mem_free)
    // is IGNORED for enforcement (a 64 GB request still gave a 48 GB cgroup; an
    // 8-thread/4-slot task gave 32 GB).
    // ROOT CAUSE: -Xmx was derived from that unenforced request (0.8 * 64 = 51 GB),
    // so the JVM's heap ceiling sat ABOVE the 48 GB cgroup — the JVM had no reason
    // to collect before the kernel killed it, and NTC (deepest BAM) grew straight
    // into the cap. The other arms peaked at only 11.5-15.5 GB.
    // FIX: size -Xmx from the ENFORCED cgroup (see script), keeping the allocation
    // at 12 threads. Adding slots would also have worked, but every caller is
    // benchmarked at the same allocation, so raising it for one caller would
    // invalidate the per-caller runtime comparison and would run the deepest arm
    // under different resources than the KD arms.
    // NB: directives set HERE are overridden whenever nextflow.config carries a
    // matching withLabel/withName selector (verified) — this process takes its
    // cpus/memory/time from `withLabel: process_calling`.
    errorStrategy { task.exitStatus in [137, 140, 143, 104, 71] ? 'retry' : 'finish' }
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
        JACUSA2_JAR=/home/yng7/miniforge3/envs/scRNAEditing/share/jacusa2-2.1.17-0/JACUSA_v2.1.17.jar
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


process SCALLELE {
    tag { "${meta.id}:${bam_label}" }
    label 'process_calling'
    // Dedicated env: scAllele 0.0.9.7 needs scipy 1.15 / numpy 1.26;
    // scRNAEditing has scipy 1.17 / numpy 2.x and must not mix.
    conda params.scallele_conda
    // EDITING DISCOVERY: emit ALL scAllele SNVs (PASS+LowQual). scAllele's GLM
    // was trained WITHOUT editing, so A-to-I always lands in LowQual — a
    // PASS-only filter would delete the editing we want. downstream_filter.R
    // then applies the uniform editing filters (see bin/scallele_io.py).
    // SNP BLACKLIST: the PASS sites (scAllele's confident genetic variants) are
    // ALSO emitted as a per-sample SNP blacklist (snp_blacklist output) and
    // subtracted downstream when --scallele_snp_subtract (default true).
    // Chr-scatter can peak on the deepest contig; a retry de-parallelises (see
    // scatter_n) so an OOM is met with lower peak memory at the fixed allocation.
    errorStrategy { task.exitStatus in [137, 140, 143, 104, 71] ? 'retry' : 'finish' }
    maxRetries    3
    publishDir { "${params.outdir}/${meta.id}/03_variants/scallele/${bam_label}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), val(bam_label), val(strandedness)
    path fasta
    path main_chroms

    output:
    tuple val(meta), val('scallele'), val(bam_label),
          path("${meta.id}_${bam_label}_scallele_variants.txt.gz"), emit: variants
    tuple val(meta), val('scallele'), val(bam_label),
          path("${meta.id}_${bam_label}_scallele_snpblacklist.txt.gz"), emit: snp_blacklist

    script:
    // Library type matches JACUSA2_CALL1 / RSeQC:
    //   10x 3' forward = second-strand; reverse = first-strand.
    def strand_flag = strandedness == 'forward' ? 'fr-secondstrand' :
                      strandedness == 'reverse' ? 'fr-firststrand' : 'fr-unstrand'
    def qual = (meta.qual_override != null) ? meta.qual_override : params.qual
    def min_cov = params.min_cov
    def min_alt = params.min_alt
    def min_edit = params.min_edit_level
    def run_mode = params.scallele_run_mode
    // Concurrent contigs = cpus on attempt 1, halving-ish each retry to cut peak
    // memory (each scAllele worker is single-threaded, so this = worker count).
    def scatter_n = Math.max(1, (task.cpus / task.attempt).intValue())
    """
    #!/bin/bash
    set -euo pipefail

    # Work-dir-local FASTA .fai. Staged FASTA is often a symlink; pysam/HTSeq
    # look beside the staged name (same pattern as REDITOOLS3).
    if [ ! -f ${fasta}.fai ]; then
        REALFA=\$(readlink -f ${fasta})
        if [ -n "\$REALFA" ] && [ -f "\${REALFA}.fai" ]; then
            cp -f "\${REALFA}.fai" ${fasta}.fai
        else
            scallele_io.py faidx ${fasta}
        fi
    fi

    if [ "${main_chroms.name}" != "NO_FILE" ] && [ -f "${main_chroms}" ]; then
        cp ${main_chroms} _chroms_req.txt
    else
        scallele_io.py bam_chroms ${bam} _chroms_req.txt
    fi
    scallele_io.py intersect_chroms ${bam} _chroms_req.txt _chroms.txt

    # Chr-SCATTER (matches BCFTOOLS_MPILEUP / GATK): scAllele's internal threading
    # (-n) scales poorly, so instead of the old serial per-chromosome loop we run
    # one single-threaded (-n 1) scAllele per contig in parallel across the
    # allocated cpus. Each worker writes per-chr outputs (no shared-file races),
    # concatenated below in _chroms.txt order. A retry de-parallelises (fewer
    # concurrent contigs) rather than asking for more slots, so an OOM is met with
    # lower peak memory at the fixed allocation.
    grep -vE '^[[:space:]]*\$' _chroms.txt | xargs -P ${scatter_n} -I{} \\
        scallele_chrom.sh {} ${bam} ${fasta} ${min_cov} ${min_alt} ${min_edit} ${qual} ${strand_flag} ${run_mode}

    # Concatenate per-chromosome outputs in chromosome (file) order.
    : > _all_snvs.txt
    : > _snp_blacklist.txt   # PASS-only genetic-variant sites (chr,pos)
    while read CHR; do
        [ -z "\$CHR" ] && continue
        [ -f "_snvs_\${CHR}.txt" ] && cat "_snvs_\${CHR}.txt" >> _all_snvs.txt
        [ -f "_bl_\${CHR}.txt" ]   && cat "_bl_\${CHR}.txt"   >> _snp_blacklist.txt
    done < _chroms.txt

    gzip -c _all_snvs.txt > ${meta.id}_${bam_label}_scallele_variants.txt.gz
    gzip -c _snp_blacklist.txt > ${meta.id}_${bam_label}_scallele_snpblacklist.txt.gz
    rm -f _all_snvs.txt _snp_blacklist.txt _snvs_*.txt _bl_*.txt _chroms.txt _chroms_req.txt
    """
}
