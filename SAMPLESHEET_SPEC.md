# Samplesheet / metadata specification

The samplesheet is the single metadata file a user must supply. Everything the
pipeline and the cell-type split need is derived from it — no dataset-specific
constants live in the scripts.

Validate with:

```bash
Rscript nextflow/bin/validate_samplesheet.R <samplesheet.csv> \
    --dataset <name> --tsv-out publication/CELLTYPE_SCOPE.tsv [--min-cells N]
```

---

## Columns

| column | required | meaning |
|---|---|---|
| `sample_id` | yes | The **dedup unit**. Rows sharing a `sample_id` are lanes/runs of one sample and are merged before dedup. Cell keys are unique *within* a `sample_id`, not across. |
| `bam` | yes | Aligned, coordinate-sorted BAM for that lane. A `.bai` must sit beside it. |
| `celltype_annot` | yes | 2-column TSV, no header: `cell_key <TAB> celltype`. |
| `called_barcodes` | yes | Cell keys passing QC. 1 or 2 columns; only the first is read. |
| `platform` | yes | One of `10xv2`, `10xv3`, `celseq2`, `smartseq2`. Determines topology and whether UMIs exist. |
| `group` | no | Pooling unit for cell-type analysis — e.g. a CRISPR arm or treatment condition. Cell types are pooled across `sample_id`s **within** a group. Absent or blank means one group named `all`. |
| `cell_id` | conditional | **Required for plate-based platforms**, where one BAM is one cell. Must match a key in `celltype_annot`. Ignored for barcoded platforms. |
| `fastq_r1`, `fastq_r2` | no | Optional provenance; may be empty. |

Unknown extra columns are permitted and ignored.

---

## The two topologies

`platform` selects between them. This is the only branch in the cell-type split.

**Barcoded** (`10xv2`, `10xv3`, `celseq2`) — one BAM holds many cells, keyed by
the `CB`/`CR` tag. `celltype_annot` keys are barcodes.

- Split each sample's BAM by cell type using its own annotation.
- Rewrite `CB`/`CR` to `<sample_id>_<barcode>` during the split.
- Merge same-celltype fragments across samples within a group.

The prefix is load-bearing, not cosmetic. 10x barcodes are reused across
samples, merged BAMs carry no `@RG`, and the pipeline dedups *after* it receives
the BAM. Since dedup is cell-aware (`umi_tools --per-cell`, MarkDuplicates
`--BARCODE_TAG`), two cells that collide on a barcode would be treated as one
and their independent molecules collapsed. Prefixing makes keys globally unique,
which makes dedup-after-merge produce exactly what per-sample dedup would.

**Plate-based** (`smartseq2`) — one BAM *is* one cell, named by `cell_id`. There
is no `CB` tag.

- No splitting: route each whole BAM to the cell type its `cell_id` maps to.
- Merge by cell type within the group.
- No prefixing — plate `cell_id`s are already globally unique, and there is no
  tag to rewrite.

## UMIs

`smartseq2` has no `UB` tag, so `--dedup umi` is invalid for it; `dedup.nf`
hard-errors rather than degrading. Only `none` and `coordinate` apply. Every
other platform above carries UMIs.

## Why `sample_id` is the dedup unit

A cell's reads are spread across its lanes, so deduping per lane would miss
cross-lane PCR duplicates and over-count. Deduping across samples would let
reused barcodes merge unrelated cells. The sample is the only granularity that
is correct on both counts.

---

## Rules the validator enforces

1. All required columns present; no blank values in them.
2. Every `bam` exists, is non-empty, and has a `.bai` beside it.
3. Every `celltype_annot` and `called_barcodes` exists and is non-empty;
   `celltype_annot` has at least 2 tab-separated columns.
4. `platform` is recognised.
5. A `sample_id` maps to exactly one `platform` and exactly one `group`.
6. Plate-based rows carry a non-blank `cell_id`; `cell_id` is unique across the
   whole sheet and is present in that row's `celltype_annot`.
7. Barcoded rows: cell keys are unique within each `sample_id`.

It also **reports** (does not fail on) the things that drive the split:

- groups, samples per group, and lanes per sample;
- cell types per group clearing `--min-cells`, and the cells retained;
- cross-sample cell-key collisions per group — the count that the barcode
  prefix exists to neutralise. Non-zero is expected and fine for barcoded data;
  non-zero for plate-based data is a genuine error, since `cell_id`s must be
  unique.

---

## The cell-type minimum is opt-in, per dataset, and recorded

`--min-cells` **defaults to 0 — keep every cell type.** A minimum cell count is
not part of the original analysis, so it is never applied silently.

It cannot be a global constant. Measured largest cell type per dataset:

| dataset | cell types | cells | largest type |
|---|---|---|---|
| pm_gt (NTC_ctrl) | 23 | 9540 | 2286 |
| jzl10x | 27 | 9806 | 1174 |
| jzlxv3 | 8 | 3222 | 962 |
| jzlss2 | 7 | 526 | 193 |
| jzlcel | 48 | 526 | **33** |

A teratoma-derived floor of 50 would erase `jzlcel` outright — 0 of its 48 cell
types reach it. Pick the value per dataset and pass it explicitly.

Whenever it is set, `--tsv-out` must be given so the choice lands in
`publication/CELLTYPE_SCOPE.tsv`, the supplementary table recording, per dataset
and group, the threshold applied and the cell types and cells it kept or
dropped. Rows are keyed on `dataset` and replaced on re-run, so the file is
idempotent and always reflects the current scope.
