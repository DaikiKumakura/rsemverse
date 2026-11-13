# rsemverse

<!-- badges: start -->
[![R-CMD-check](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Read RSEM gene-level results as a table of gene symbols by samples.

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

## License

MIT © Daiki Kumakura
