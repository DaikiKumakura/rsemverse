# rsemverse

<!-- badges: start -->
[![R-CMD-check](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Documentation:** <https://daikikumakura.github.io/rsemverse/> ·
[Get started](https://daikikumakura.github.io/rsemverse/articles/rsemverse.html) ·
[Installation](https://daikikumakura.github.io/rsemverse/articles/installation.html) ·
[Tutorial](https://daikikumakura.github.io/rsemverse/articles/tutorial.html) ·
[How it works](https://daikikumakura.github.io/rsemverse/articles/how-it-works.html) ·
[Troubleshooting](https://daikikumakura.github.io/rsemverse/articles/troubleshooting.html) ·
[Real-data check](https://daikikumakura.github.io/rsemverse/articles/real-data.html)

Read RSEM gene-level results as a table of gene symbols by samples.

```r
library(rsemverse)

df <- rsem_read(
  data     = "data",
  metadata = "metadata.csv",
  species  = "human",   # "human" or "mouse"
  expr     = "tpm"      # "count" or "tpm"
)
```

`df` is a data frame with one row per gene symbol and one column per sample.
Save it as you like, for example `write.csv(df, "tpm.csv")`.

## Input

```text
metadata.csv
data/
  S1.genes.results
  S2.genes.results
  ...
```

`metadata.csv` has two columns, `sample_id` and `ref`. `ref` is the name of the
sample's RSEM file inside `data`:

```text
sample_id,ref
S1,S1.genes.results
S2,S2.genes.results
```

Other columns are ignored, and the output columns follow the row order of this
file. Each RSEM file needs `gene_id` and `expected_count` (for `expr = "count"`)
or `TPM` (for `expr = "tpm"`). All files must have the same genes in the same
order.

## What `rsem_read()` does

1. Checks the metadata and the files, and stops with a message on any problem:
   a missing file, a missing column, duplicated IDs, non-numeric or missing
   values, or different gene tables across files.
2. Checks that the gene IDs belong to the species you gave.
3. Removes the version of Ensembl IDs (`ENSG00000141510.18` becomes
   `ENSG00000141510`) and looks up the gene symbol in `org.Hs.eg.db` or
   `org.Mm.eg.db`. The OrgDb keys have no version, so this step is required.
4. Sums the rows that share a symbol, for counts and for TPM alike.

Details that are good to know:

* An Ensembl ID with several symbols gets the symbol of its lowest Entrez ID.
* A gene without a symbol keeps the `gene_id` from the file (spike-ins,
  numeric IDs, IDs ending in `_PAR_Y`, and Ensembl IDs the OrgDb does not
  know), so the column sums equal those of the input. If such an ID equals a
  gene symbol, the function stops instead of merging them.
* Rows appear in the order of first occurrence in the files.
* Values are not rounded, and TPM is not renormalized.
* A message reports the numbers of samples, genes, rows, genes with a symbol,
  merged symbols, and genes kept under their original ID, together with the
  version of the OrgDb package. Use `suppressMessages()` to silence it.

The symbols depend on the installed version of the OrgDb package, so record it
(it is in the message) when you report results.

## Checked on real data

`validation/encode-k562.R` downloads six ENCODE RSEM files (human K562 cells,
GRCh38, GENCODE V29; control and two shRNA knockdowns, two replicates each),
checks their MD5 sums, and runs `rsem_read()`. Result of November 10, 2026, with
org.Hs.eg.db 3.22.0:

| | |
| --- | --- |
| Genes per file | 59,526: 58,735 versioned Ensembl IDs, 97 ERCC spike-ins, 649 numeric IDs, 45 ending in `_PAR_Y` |
| Rows returned | 59,427 |
| With a symbol | 36,385 genes (89 symbols merged from several IDs) |
| Kept under the original ID | 23,141 genes: the 22,350 Ensembl IDs that have no symbol in the OrgDb, plus the 791 non-Ensembl IDs |
| Column sums | equal to the input, for counts and for TPM |
| Share of the total in the kept rows | 3.76–4.90% of the expected count and 5.89–7.99% of the TPM, per sample |
| Run time | 3.6 s for 6 samples |
| Knockdown targets (mean count, knockdown / control) | PTBP1 0.35 in the PTBP1 knockdown; SRSF1 0.19 in the SRSF1 knockdown |

About 39% of the genes in a GENCODE-based RSEM file have no symbol in the
OrgDb, but they hold only a few percent of the counts. They remain in the
table under their Ensembl ID; select the rows whose name starts with `ENSG` if
you want to drop them. These unannotated genes were not characterized further.

## Installation

```r
install.packages("BiocManager")
BiocManager::install("DaikiKumakura/rsemverse")
BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db"))  # the species you need
```

R 4.1 or newer. The package needs `data.table` and `AnnotationDbi`; install
only the OrgDb package of the species you use.

## Example

The package contains a synthetic example (3 samples, 44 genes with real human
Ensembl IDs, one symbol shared by two IDs, and two IDs without a symbol):

```r
ex <- system.file("extdata", "example", package = "rsemverse")
tpm <- rsem_read(file.path(ex, "data"), file.path(ex, "metadata.csv"),
                 species = "human", expr = "tpm")
```

The expression values are synthetic.

## License

MIT © Daiki Kumakura
