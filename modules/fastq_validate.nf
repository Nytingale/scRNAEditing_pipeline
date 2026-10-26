/*
 * Checks that fastq_r1 / fastq_r2 are the barcode / cDNA reads (SRA numbering is
 * unreliable), from read length and polyT.
 */

def getBarcodeCore(platform) {
    def lengths = [
        '10xv2':   26,  // 16bp CB + 10bp UMI
        '10xv3':   28,  // 16bp CB + 12bp UMI
        'dropseq': 20,  // 12bp CB + 8bp UMI
        'celseq2': 12,  // 6bp UMI + 6bp CB
        'seqwell': 20,  // 12bp CB + 8bp UMI
        'dronc':   20,  // 12bp CB + 8bp UMI
        'sciseq':  18,  // 10bp CB + 8bp UMI
        'indrop':  14,  // 8bp CB + 6bp UMI
    ]
    lengths.getOrDefault(platform, 0)
}


process VALIDATE_FASTQ_READS {
    tag { "${meta.lane_id ?: meta.id}" }
    label 'process_low'
    errorStrategy 'terminate'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("validated_barcode.fastq.gz"), path("validated_cdna.fastq.gz"), emit: validated

    script:
    def platform = meta.platform
    def bc_core = getBarcodeCore(platform)
    def check_n = 20

    if (platform == 'smartseq2') {
        // Smart-seq2: no barcode read distinction, just pass through
        """
        ln -s ${r1} validated_barcode.fastq.gz
        ln -s ${r2} validated_cdna.fastq.gz
        """
    } else {
        """
        #!/bin/bash
        set -euo pipefail

        BC_CORE=${bc_core}
        CHECK_N=${check_n}
        R1="${r1}"
        R2="${r2}"

        # Get read lengths (nf_head avoids SIGPIPE 141 under pipefail)
        LEN_R1=\$(nf_head.sh 2 zcat "\$R1" | tail -1 | tr -d '\\n' | wc -c)
        LEN_R2=\$(nf_head.sh 2 zcat "\$R2" | tail -1 | tr -d '\\n' | wc -c)

        echo "[VALIDATE] ${meta.lane_id ?: meta.id}: R1=\$(basename \$R1) len=\${LEN_R1}bp, R2=\$(basename \$R2) len=\${LEN_R2}bp"

        BARCODE_FILE=""
        CDNA_FILE=""

        if [ "\$LEN_R1" -ne "\$LEN_R2" ]; then
            # Different lengths: shorter is barcode
            if [ "\$LEN_R1" -lt "\$LEN_R2" ]; then
                BARCODE_FILE="\$R1"; CDNA_FILE="\$R2"
                echo "[VALIDATE] Length-based: R1 is barcode (\${LEN_R1}bp), R2 is cDNA (\${LEN_R2}bp)"
            else
                BARCODE_FILE="\$R2"; CDNA_FILE="\$R1"
                echo "[VALIDATE] Length-based: R2 is barcode (\${LEN_R2}bp), R1 is cDNA (\${LEN_R1}bp)"
            fi
        else
            # Same length: detect polyT run after barcode region
            PT_R1=\$(nf_head.sh \$((CHECK_N*4)) zcat "\$R1" | awk -v start="\$BC_CORE" '
                NR%4==2 {
                    total++; seq=\$0; len=length(seq);
                    if(len<=start) next;
                    tail_region=substr(seq,start+1); max_t=0; cur_t=0;
                    for(i=1;i<=length(tail_region);i++){
                        if(substr(tail_region,i,1)=="T") cur_t++;
                        else { if(cur_t>max_t) max_t=cur_t; cur_t=0 }
                    }
                    if(cur_t>max_t) max_t=cur_t;
                    if(max_t>=6) polyt++
                }
                END { print polyt+0 }')

            PT_R2=\$(nf_head.sh \$((CHECK_N*4)) zcat "\$R2" | awk -v start="\$BC_CORE" '
                NR%4==2 {
                    total++; seq=\$0; len=length(seq);
                    if(len<=start) next;
                    tail_region=substr(seq,start+1); max_t=0; cur_t=0;
                    for(i=1;i<=length(tail_region);i++){
                        if(substr(tail_region,i,1)=="T") cur_t++;
                        else { if(cur_t>max_t) max_t=cur_t; cur_t=0 }
                    }
                    if(cur_t>max_t) max_t=cur_t;
                    if(max_t>=6) polyt++
                }
                END { print polyt+0 }')

            echo "[VALIDATE] PolyT detection: R1=\${PT_R1}/\${CHECK_N}, R2=\${PT_R2}/\${CHECK_N}"

            if [ "\$PT_R1" -gt "\$PT_R2" ] && [ "\$PT_R1" -gt 0 ]; then
                BARCODE_FILE="\$R1"; CDNA_FILE="\$R2"
                echo "[VALIDATE] PolyT-based: R1 is barcode, R2 is cDNA"
            elif [ "\$PT_R2" -gt "\$PT_R1" ] && [ "\$PT_R2" -gt 0 ]; then
                BARCODE_FILE="\$R2"; CDNA_FILE="\$R1"
                echo "[VALIDATE] PolyT-based: R2 is barcode, R1 is cDNA"
            else
                echo "ERROR: Cannot distinguish barcode from cDNA for ${meta.lane_id ?: meta.id}" >&2
                echo "  Both reads are \${LEN_R1}bp with similar polyT: R1=\${PT_R1}, R2=\${PT_R2}" >&2
                echo "  Platform ${platform} expects ${bc_core}bp barcode core" >&2
                exit 1
            fi
        fi

        # Validate barcode length >= core length
        BC_LEN=\$(nf_head.sh 2 zcat "\$BARCODE_FILE" | tail -1 | tr -d '\\n' | wc -c)
        if [ "\$BC_LEN" -lt "\$BC_CORE" ]; then
            echo "ERROR: Barcode read length (\${BC_LEN}bp) < expected core (\${BC_CORE}bp) for ${platform}" >&2
            exit 1
        fi

        ln -s "\$BARCODE_FILE" validated_barcode.fastq.gz
        ln -s "\$CDNA_FILE" validated_cdna.fastq.gz
        echo "[VALIDATE] PASSED: barcode=\$(basename \$BARCODE_FILE), cDNA=\$(basename \$CDNA_FILE)"
        """
    }
}
