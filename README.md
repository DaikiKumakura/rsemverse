# rsemverse

<!-- badges: start -->
[![R-CMD-check](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210460.svg)](https://doi.org/10.5281/zenodo.24210460)
<!-- badges: end -->

**From RSEM gene results to analysis-ready gene tables in R, for human and
mouse bulk RNA-seq.**

rsemverse reads RSEM `*.genes.results` files into one table of gene symbols by
samples, then filters, selects, normalizes and transforms that table, and
finally turns it into a long table for plotting. Every function takes the
same table and returns the same table, so the steps chain in any order.

- **No gene lists to prepare.** Gene symbols and gene groups (protein-coding,
  mitochondrial, ribosomal, …) come from `org.Hs.eg.db` / `org.Mm.eg.db`.
- **No sample groups needed.** Every step works on the table alone; the sample
  table is used only to label samples.
- **Every step reports what it did,** with numbers: how many genes were kept,
  which threshold was used, the normalization factors.
- **Checked against the reference tools.** Filtering and normalization give
  the same results as edgeR, DESeq2 and limma, without depending on them.

**Documentation:** <https://daikikumakura.github.io/rsemverse/>

## Workflow

```text
  RSEM *.genes.results + metadata.csv
                 │
          rsem_read()         read: gene symbol × sample table (count or TPM)
                 │
          rsem_filter()       remove lowly expressed genes
          rsem_subset()       keep or remove gene groups
          rsem_normalize()    normalize between samples
          rsem_transform()    log2, z-score, ranks
                 │
          rsem_long()         one row per gene and sample, with the sample table
                 │
         ggplot2, dplyr, bulkdeconv, …
```

| Step | Function | Choices (default first) |
|---|---|---|
| Read | `rsem_read()` | `expr = "count"` or `"tpm"`; `species = "human"` or `"mouse"` |
| Filter | `rsem_filter()` | `"edger"` (counts) / `"detected"` (TPM), `"nonzero"`, `"variable"` |
| Select genes | `rsem_subset()` | `"protein_coding"`, `"mito"`, `"ribo"`, `"hb"`, `"sex"`, `"ercc"` |
| Normalize | `rsem_normalize()` | `"cpm"`, `"tmm"`, `"upper_quartile"`, `"deseq2"`, `"quantile"` |
| Transform | `rsem_transform()` | `"log2"`, `"zscore"`, `"center"`, `"rank"` |
| Long table | `rsem_long()` | none: every column of the sample table is added as it is |

## Installation

```r
install.packages("BiocManager")
BiocManager::install("DaikiKumakura/rsemverse")
BiocManager::install("org.Hs.eg.db")   # human; org.Mm.eg.db for mouse
```

## Quick start

The package includes three small RSEM files, so this runs as it is:

```r
library(rsemverse)
ex <- system.file("extdata", "example", package = "rsemverse")

counts <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
                    species = "human", expr = "count")
kept   <- rsem_filter(counts, expr = "count")
genes  <- rsem_subset(kept, "human", keep = "protein_coding", remove = c("mito", "ribo"))
logcpm <- rsem_transform(rsem_normalize(genes, expr = "count", method = "tmm"), "log2")
long   <- rsem_long(logcpm, metadata = file.path(ex, "metadata.csv"))
```

Each call prints one line:

```text
rsem_read: 3 samples, 44 genes -> 43 rows (species = human, expr = count, org.Hs.eg.db 3.22.0)
rsem_filter (edger): kept 42 of 43 genes (CPM >= 1260 in >= 3 of 3 samples and total count >= 15); removed 1
rsem_subset: kept 40 of 42 genes; keep protein_coding; remove mito, ribo (genes per group: ...)
rsem_normalize (tmm, count): 40 genes x 3 samples; effective library sizes: S1 6587, S2 7145, S3 7827
rsem_transform: 40 genes x 3 samples; log2(value + 1)
rsem_long: 120 rows (40 genes x 3 samples) with 1 sample column(s): ref (text)
```

## Your data

Put the RSEM files in a folder and list them in `metadata.csv`. `ref` is the
file name inside the folder. Add any other columns you like (group, batch,
age, …); `rsem_read()` ignores them and `rsem_long()` adds them to the long
table.

```text
metadata.csv          data/
                        S1.genes.results
                        S2.genes.results
```

```text
sample_id,ref,group,age
S1,S1.genes.results,control,45
S2,S2.genes.results,treated,52
```

## The table

Every function except `rsem_long()` returns the same kind of table:

```text
        S1     S2     S3      <- sample_id, in the order of metadata.csv
TP53    12.4   10.9   11.7
GAPDH   980.1  1012.5 955.0
...
^ gene symbols as row names
```

A data frame with the symbols in the first column, or a numeric matrix with
row names, is accepted too, so tables from elsewhere work as well.

## The steps

### Read: `rsem_read()`

Reads the expected count or the TPM of every file, removes the version of
Ensembl IDs (`ENSG00000141510.18` → `ENSG00000141510`), looks up the symbols,
and adds up rows that get the same symbol. A gene without a symbol keeps its
original ID, so nothing is lost: the column sums equal the input. It stops
with a message when files are missing, genes differ between files, or the
species does not match.

### Filter: `rsem_filter()`

| Method | Keeps a gene when | Typical use |
|---|---|---|
| `"edger"` (default for counts) | CPM reaches 10 counts in a median-sized sample, in all samples (or 10 + 70% of the rest above 10), and the total count is at least 15; edgeR `filterByExpr()` with all samples as one group | before normalization or differential expression |
| `"detected"` (default for TPM) | the value reaches `min_value` (TPM 1, count 10) in `min_samples` samples (20% of them) | TPM tables; "expressed in at least n samples" |
| `"nonzero"` | it is above zero in any sample | removing empty rows only |
| `"variable"` | it is among the `n_top` (2000) with the largest variance of log2 | before PCA, clustering, heatmaps |

### Select genes: `rsem_subset()`

`keep` takes a gene when it is in any listed group; `remove` then drops genes
in any listed group.

| Group | Genes | Why |
|---|---|---|
| `"protein_coding"` | gene type protein-coding in NCBI Gene | coding-only analyses, comparisons with proteomics |
| `"mito"` | the 37 genes of the mitochondrial genome | a high share points to damaged cells or degraded RNA |
| `"ribo"` | cytosolic ribosomal protein genes (`RPL*`, `RPS*`) | highly expressed; often removed before clustering |
| `"hb"` | hemoglobin genes (`HBA1`, `HBB`, …) | blood contamination |
| `"sex"` | chromosome Y genes and `XIST` | separate samples by sex |
| `"ercc"` | ERCC spike-ins | technical controls |

```r
rsem_subset(counts, "human", keep = "protein_coding", remove = c("mito", "ribo", "sex"))
```

### Normalize: `rsem_normalize()`

| Method | Corrects for | Values |
|---|---|---|
| `"cpm"` | sequencing depth | counts per million |
| `"tmm"` | depth and RNA composition (edgeR) | CPM with TMM factors |
| `"upper_quartile"` | depth, using the 75th percentile | CPM with upper-quartile factors |
| `"deseq2"` | depth and RNA composition (DESeq2) | counts divided by median-of-ratios size factors |
| `"quantile"` | the whole distribution (strong assumption) | every sample gets the same distribution |

TPM already corrects for depth and gene length, so for TPM only `"quantile"`
is offered. The factors are kept in `attr(x, "factors")`. Filter first: genes
with very low counts make the factors less stable.

### Transform: `rsem_transform()`

| Method | Result | Typical use |
|---|---|---|
| `"log2"` | `log2(value + 1)` | plots, distances, PCA |
| `"zscore"` | each gene: mean 0, standard deviation 1 | heatmaps |
| `"center"` | each gene: mean 0 | PCA on centered data |
| `"rank"` | rank within each sample, from 0 to 1 | rank-based scores, comparing differently scaled data |

Methods run in the order given: `rsem_transform(x, c("log2", "zscore"))`.

### Long table: `rsem_long()`

Returns `gene`, `sample_id`, `value` and every column of the sample table,
with no setting per column: numbers stay numbers, text stays text, and
`sample_id` keeps its leading zeros. Samples removed earlier are left out of
the sample table with a message.

```r
long <- rsem_long(logcpm, metadata = "metadata.csv")
ggplot2::ggplot(subset(long, gene == "TP53"), ggplot2::aes(group, value)) +
  ggplot2::geom_boxplot()
```

## Recipes

| Goal | Steps |
|---|---|
| Heatmap of variable genes | counts → `rsem_filter()` → `rsem_normalize("tmm")` → `rsem_filter("variable")` → `rsem_transform(c("log2", "zscore"))` |
| PCA | counts → `rsem_filter()` → `rsem_subset(remove = c("mito", "ribo", "sex"))` → `rsem_normalize("tmm")` → `rsem_transform("log2")` |
| Plot genes by group | any of the above → `rsem_long(metadata = "metadata.csv")` → ggplot2 |
| Cell-type deconvolution | `rsem_read(expr = "tpm")` straight into [bulkdeconv](https://daikikumakura.github.io/bulkdeconv/); do not filter, normalize or log-transform first |
| Differential expression | `rsem_read(expr = "count")` → `rsem_filter()` → DESeq2 (after `round()`) or edgeR; both need the raw counts, not normalized values |

For the heatmap, choose the genes on the normalized values and scale them
last (`"variable"` takes the log itself):

```r
counts <- rsem_read("data", "metadata.csv", species = "human", expr = "count")
tmm    <- rsem_normalize(rsem_filter(counts, "count"), "count", "tmm")
top    <- rsem_filter(tmm, "count", method = "variable", n_top = 500)
heat   <- rsem_transform(top, c("log2", "zscore"))
```

## Validation

- **Tests:** the filtering and normalization results equal those of edgeR
  4.8.2, DESeq2 1.50.2 and limma 3.66.0 on simulated counts.
- **Real data:** on six ENCODE K562 RSEM files (59,526 genes each),
  `rsem_read()` gave symbols to 36,385 genes, kept 23,141 under their ID and
  preserved the column sums. The edgeR rule kept the same 12,808 genes as
  `edgeR::filterByExpr()`, every normalization was within 3e-16 (relative) of
  the reference tools, and all 37 mitochondrial genes were found. See
  [`validation/`](validation/).
- **Platforms:** R CMD check passes on Linux, macOS and Windows in CI.

## Related tools

- [bulkdeconv](https://daikikumakura.github.io/bulkdeconv/) estimates cell-type
  composition from the TPM table returned by `rsem_read()`.
- [rsemflow](https://daikikumakura.github.io/rsemflow/) reads the same RSEM
  files to run DESeq2, GSEA, GSVA and PCA from the command line.

## Citation

Cite the archived version: <https://doi.org/10.5281/zenodo.24210460>. The
citation metadata is in [`CITATION.cff`](CITATION.cff).

## License

MIT © Daiki Kumakura
