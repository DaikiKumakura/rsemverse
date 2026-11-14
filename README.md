# rsemverse

<!-- badges: start -->
[![R-CMD-check](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Read RSEM gene-level results as a table of gene symbols by samples, then
filter, normalize, transform and subset it, or turn it into a long table for
plotting.

**Documentation:** [Function reference and examples](https://daikikumakura.github.io/rsemverse/).

## Usage

Put the RSEM files in a folder and list them in `metadata.csv`:

```text
metadata.csv
data/
  S1.genes.results
  S2.genes.results
```

```text
sample_id,ref
S1,S1.genes.results
S2,S2.genes.results
```

`ref` is the file name inside `data`. Other columns are ignored.

```r
library(rsemverse)

df <- rsem_read(
  data     = "data",
  metadata = "metadata.csv",
  species  = "human",   # "human" or "mouse"
  expr     = "tpm"      # "count" or "tpm"
)
```

`df` is a data frame: gene symbols in the rows, `sample_id` in the columns, in
the order of `metadata.csv`. Save it as you like, for example
`write.csv(df, "tpm.csv")`.

## What it does

* Removes the version of Ensembl IDs (`ENSG00000141510.18` → `ENSG00000141510`)
  and looks up the symbol in `org.Hs.eg.db` or `org.Mm.eg.db`.
* Adds up the rows that end up with the same symbol.
* Keeps a gene without a symbol under its original ID, so nothing is lost.

It stops with a message if the files do not match (missing file, different
genes, wrong species, and so on). Values are not rounded or normalized. See
`?rsem_read` for the details.

## Prepare the table

Every function takes the table and returns a table of the same layout (genes
in the rows, samples in the columns), so the steps can be chained. Each one
prints what it did, with numbers.

```r
counts <- rsem_read("data", "metadata.csv", species = "human", expr = "count")

expressed <- rsem_filter(counts, expr = "count")          # edgeR filterByExpr rule
pc        <- rsem_subset(expressed, "human",
                         keep = "protein_coding", remove = c("mito", "ribo"))
tmm       <- rsem_normalize(pc, expr = "count", method = "tmm")
log_tmm   <- rsem_transform(tmm, "log2")                  # log2(TMM-CPM + 1)
long      <- rsem_long(log_tmm, metadata = "metadata.csv")
```

| Function | What it does | Choices |
|---|---|---|
| `rsem_filter()` | Removes genes | `"edger"` (edgeR `filterByExpr()`, all samples as one group; default for counts), `"detected"` (value above a threshold in n or a share of samples; default for TPM: TPM >= 1 in 20% of samples), `"nonzero"`, `"variable"` (top 2000 by variance of log2) |
| `rsem_normalize()` | Normalizes between samples | counts: `"cpm"`, `"tmm"`, `"upper_quartile"`, `"deseq2"` (median of ratios), `"quantile"`; TPM: `"quantile"` |
| `rsem_transform()` | Changes the scale | `"log2"`, `"zscore"`, `"center"`, `"rank"`, applied in the order given |
| `rsem_subset()` | Keeps or removes gene groups | `"protein_coding"`, `"mito"`, `"ribo"`, `"hb"`, `"sex"`, `"ercc"` |
| `rsem_long()` | One row per gene and sample, with every column of `metadata.csv` | none: numbers stay numbers, text stays text |

The normalization methods give the same values as edgeR, DESeq2 and limma
(checked in the tests) without depending on them. The gene groups come from
the same `org.Hs.eg.db` / `org.Mm.eg.db` as the symbols, so no gene list is
needed. Differential expression tools (DESeq2, edgeR) need the raw counts, not
normalized values. See each function's help page for the details.

## Use the table with other tools

The table goes straight into [bulkdeconv](https://daikikumakura.github.io/bulkdeconv/),
which estimates cell-type composition or abundance from tumor bulk RNA-seq
(use `expr = "tpm"` for most methods):

```r
tpm <- rsem_read("data", "metadata.csv", species = "human", expr = "tpm")
res <- bulkdeconv::bulkdeconv(tpm, expr = "tpm", species = "human", method = "mcp_counter")
res$score   # samples x cell types
```

A method refuses a table with too few of its signature genes, so a small test
table may not work with every method.

[rsemflow](https://daikikumakura.github.io/rsemflow/) reads the RSEM files
itself, together with a sample metadata table, to run DESeq2, GSEA, GSVA and
PCA, so it does not need `rsem_read()` first.

## Installation

```r
install.packages("BiocManager")
BiocManager::install("DaikiKumakura/rsemverse")
BiocManager::install("org.Hs.eg.db")   # human; org.Mm.eg.db for mouse
```

## Related

- [rsemverse on the author's site](https://daikikumakura.github.io/software/rsemverse.html): summary, installation and status.
- [bulkdeconv](https://daikikumakura.github.io/software/bulkdeconv.html) turns the table into cell-type scores; [rsemflow](https://daikikumakura.github.io/software/rsemflow.html) produces the RSEM files.

## License

MIT © Daiki Kumakura
