# rsemverse

<!-- badges: start -->
[![R-CMD-check](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210459.svg)](https://doi.org/10.5281/zenodo.24210459)
<!-- badges: end -->

**RSEM results in R, one step at a time, with the pipe.** rsemverse reads
`*.genes.results` and `*.isoforms.results` files, checks them, and turns them into the
tables that plots, scores and downstream tools need. Every step takes a table and returns a
table, so steps chain with `|>` the way tidyverse functions do, and each step prints one line
about what it did.

```r
rsem_read("rsem", species = "human", expr = "count") |>   # files -> gene symbols x samples
  rsem_filter() |>                                        # drop lowly expressed genes
  rsem_normalize(method = "tmm") |>                       # between-sample normalization
  rsem_transform("log2") |>                               # log2(value + 1)
  rsem_sample_qc()                                        # PCA, correlation, per-sample numbers
```

<a id="install"></a>

## Installation

```r
install.packages("BiocManager")
BiocManager::install("DaikiKumakura/rsemverse")
BiocManager::install("org.Hs.eg.db")   # human; org.Mm.eg.db for mouse (only for gene symbols)
```

## Quick start

Three small RSEM files come with the package, so this runs as it is:

```r
library(rsemverse)
ex   <- system.file("extdata", "example", package = "rsemverse")
meta <- file.path(ex, "metadata.csv")

logcpm <- file.path(ex, "data") |>
  rsem_read(meta, species = "human", expr = "count") |>
  rsem_filter() |>
  rsem_subset("human", keep = "protein_coding", remove = c("mito", "ribo")) |>
  rsem_normalize(method = "tmm") |>
  rsem_transform("log2")
```

```text
rsem_read: 3 samples, 44 genes -> 43 rows (species = human, expr = count, org.Hs.eg.db 3.22.0)
rsem_filter (edger): kept 42 of 43 genes (CPM >= 1260 in >= 3 of 3 samples and total count >= 15); removed 1
rsem_subset: kept 40 of 42 genes; keep protein_coding; remove mito, ribo (...)
rsem_normalize (tmm, count): 40 genes x 3 samples; effective library sizes: S1 6587, S2 7145, S3 7827
rsem_transform: 40 genes x 3 samples; log2(value + 1)
```

The table at the end knows how it was made:

```r
rsem_steps(logcpm)
#> "read: count of genes"  "filter: edger, kept 42 of 43 genes"  "subset: kept 40 of 42 genes"
#> "normalize: tmm"        "transform: log2(value + 1)"
```

From here, a long table for ggplot2 (`rsem_long(logcpm, meta)`), a file for another tool
(`rsem_write(logcpm, "out.tsv")`) or the checks below.

## Documentation

[Package website](https://daikikumakura.github.io/rsemverse/) · [Function reference](https://daikikumakura.github.io/rsemverse/reference/) · [Changelog](https://daikikumakura.github.io/rsemverse/news/) · [Releases](https://github.com/DaikiKumakura/rsemverse/releases) · [R-universe](https://daikikumakura.r-universe.dev/rsemverse)

## How it fits together

- **One kind of table.** Genes are the rows (gene symbols), samples the columns. `rsem_read()`
  makes it; `rsem_filter()`, `rsem_subset()`, `rsem_normalize()`, `rsem_transform()` and the
  reference z-scores take it and return it.
- **It remembers.** The table carries what it holds (`count`, `tpm`, ...), the original gene IDs
  behind each row, the normalization factors and the steps. So `rsem_filter()` and
  `rsem_normalize()` need no `expr = ` after `rsem_read()`, `rsem_sample_qc()` knows when the
  values are already log2, and `rsem_steps()` lists what was done.
- **It does not guess.** Samples are named after the files, groups come from your sample table,
  rows with the same symbol are summed (never the maximum), IDs without a symbol keep their ID.
  Files with different genes, a different order or a different set of transcripts per gene stop the
  run instead of being merged by position.
- **Checks first.** `rsem_check()` lists every problem of all files at once before you read them.

## The functions

| Step | Function | What it does |
|---|---|---|
| **Read** | `rsem_read()` | files -> a table of one value (`count`, `tpm`, `fpkm`, `length`, `isopct`, posterior values) |
| | `rsem_import()` | several values at once (`count`, `tpm`, `length` by default) |
| **Prepare** | `rsem_filter()` | remove lowly expressed genes (edgeR rule, thresholds) or keep the most variable |
| | `rsem_subset()` | keep or remove gene groups: protein-coding, mitochondrial, ribosomal, rRNA, ... |
| | `rsem_normalize()` | CPM, TMM, TMMwsp, upper quartile, DESeq2 median of ratios, quantile |
| | `rsem_transform()` | log2, z-score, centering, ranks |
| | `rsem_fit_reference()` `rsem_apply_reference()` | z-scores against a fixed set of control samples |
| **Check the files** | `rsem_check()` | every problem in the files, all at once |
| | `rsem_reconcile()` | do gene files and isoform files agree (up to RSEM's rounding) |
| | `rsem_id_report()` | which gene IDs were merged into one row |
| **QC** | `rsem_sample_qc()` | PCA, correlation, relative log expression per sample |
| | `rsem_qc_summary()` | one table with the quality numbers of RSEM's `.cnt`, `.model` and `.theta` files |
| | `rsem_mapping_qc()` `rsem_alignment_profile()` `rsem_model_qc()` `rsem_noise_qc()` | the parts of that table, in full |
| | `rsem_biotype_profile()` `rsem_length_profile()` `rsem_detection_curve()` | share of each gene type, effective lengths, genes kept per threshold |
| **Transcripts** | `rsem_select_isoforms()` `rsem_isoform_usage()` | keep isoforms by expression; how a gene's TPM is split |
| **Hand over** | `rsem_long()` | one row per gene and sample, with your sample table, for ggplot2 |
| | `rsem_write()` | a table (or GCT for GSEA) with gene IDs and types |
| | `rsem_as_tximport()` | the list `DESeqDataSetFromTximport()` takes |
| | `rsem_as_summarizedexperiment()` | a `SummarizedExperiment` |
| | `rsem_steps()` | what was done to a table |

## Recipes

```r
# Heatmap of the most variable genes
rsem_read("rsem", meta, species = "human", expr = "count") |>
  rsem_filter() |> rsem_normalize(method = "tmm") |>
  rsem_filter(method = "variable", n_top = 500) |>      # chooses on the normalized values
  rsem_transform(c("log2", "zscore"))

# Treated samples against controls, gene by gene
tpm <- rsem_read("rsem", meta, species = "human", expr = "tpm") |> rsem_filter()
ref <- rsem_fit_reference(tpm, reference = "control", metadata = meta, group = "group")
tpm |> rsem_apply_reference(ref)

# Look at the files before reading them
rsem_check("rsem")                                    # problems, all at once
rsem_reconcile("rsem")                                # gene files against isoform files
rsem_qc_summary("rsem")                               # .cnt, .model, .theta, TPM totals

# Hand over
rsem_import("rsem", meta, assays = c("count", "tpm", "length")) |> rsem_as_tximport()   # DESeq2
rsem_read("rsem", meta, species = "human", expr = "tpm") |> rsem_write("tpm.gct", format = "gct")  # GSEA
```

Cell-type deconvolution takes the TPM table straight from `rsem_read(expr = "tpm")`
([bulkdeconv](https://daikikumakura.github.io/bulkdeconv/)); do not filter, normalize or
log-transform first.

## Your files

Put the RSEM files in a folder. Without a sample table, `rsem_read("rsem", metadata = NULL, ...)`
finds every `*.genes.results` file under it (`.gz` too) and names each sample after its file
(`S1.genes.results` -> `S1`), or after its folder with `sample_from = "parent_dir"`. It prints which
file became which sample. With a sample table, you choose the names and keep the groups:

```text
sample_id,ref,group,age
S1,S1.genes.results,control,45
S2,S2.genes.results,treated,52
```

`ref` is the file name in the folder; any other column (group, batch, age, ...) goes to
`rsem_long()` and the QC functions. Gene symbols come from `org.Hs.eg.db` or `org.Mm.eg.db`; for
other organisms use `id = "original"` (the IDs of the files) or `id = "annotation"` with a GTF.

## Differential expression

Test on the RSEM files, not on rounded counts: a gene's effective length differs between samples,
and DESeq2 and edgeR take that as an offset through `tximport`. `rsem_as_tximport()` gives the same
list as `tximport(type = "rsem", txIn = FALSE)`; genes with an effective length of 0 stop DESeq2, so
choose `zero_length = "drop"` or import transcript-level files. Use rsemverse's tables for filtering
decisions, plots and scores. For transcript-level tests that use how certain the read assignment is,
edgeR has `catchRSEM(files, ngibbs = )` for files made with `--calc-pme` or `--calc-ci`; give the
number of Gibbs samples that was used.

## Related tools

[bulkdeconv](https://daikikumakura.github.io/bulkdeconv/) estimates cell-type composition from the
TPM table; [rsemflow](https://daikikumakura.github.io/rsemflow/) runs DESeq2, GSEA and GSVA from the
command line on the same files.

<a id="checked-against"></a>

## Validation and limitations

- **edgeR, DESeq2, limma, tximport:** filtering and normalization (TMMwsp included) equal edgeR 4.8.2,
  DESeq2 1.50.2 and limma 3.66.0; `rsem_as_tximport()` equals `tximport(type = "rsem")`.
- **Real data:** six ENCODE K562 RSEM runs (59,526 genes each), with their transcript files and the GENCODE V29
  annotation (`validation/`).
- **RSEM itself:** RSEM 1.3.3 was run on a toy data set (`--calc-pme --calc-ci --append-names
  --estimate-rspd`); its output is in `tests/testthat/fixtures/rsem-real/` and is read by the tests.
- **Platforms:** R CMD check passes on Linux, macOS and Windows.

More: [function reference](https://daikikumakura.github.io/rsemverse/reference/) ·
[from hand-written scripts](https://daikikumakura.github.io/rsemverse/articles/from-scripts.html) ·
[workflow](https://daikikumakura.github.io/rsemverse/articles/workflow.html) · [NEWS](NEWS.md)

## Citation

Cite the archived version: <https://doi.org/10.5281/zenodo.24210459> (all versions). The metadata
is in [`CITATION.cff`](CITATION.cff).

## License

MIT © Daiki Kumakura
