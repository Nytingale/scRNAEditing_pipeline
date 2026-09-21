Barcode Whitelist Files for STARsolo
=====================================

Place barcode whitelist files here for STARsolo remapping (--remap true).
These are only needed if you are remapping from FASTQs — not needed for BAM input.

Files included:

  10x Chromium v2:   737K-august-2016.txt   (737,280 barcodes)
  10x Chromium v3:   3M-february-2018.txt   (6,794,880 barcodes)
  InDrop v2 CB1:     inDropCB1.txt          (384 gel barcodes, part 1)
  InDrop v2 CB2:     inDropCB2.txt          (384 gel barcodes, part 2)

Sources:
  - CellRanger whitelists: cellranger-x.x.x/lib/python/cellranger/barcodes/
  - STAR source: https://github.com/alexdobin/STAR/tree/master/extras

Other platforms (Drop-seq, CEL-Seq2, Seq-Well, DroNc-seq, sci-seq, InDrop v1)
do not require whitelists — STARsolo uses --soloCBwhitelist None.

InDrop v2 uses inDropCB1.txt and inDropCB2.txt by default. To override,
specify custom paths via the samplesheet 'whitelist_files' column.
