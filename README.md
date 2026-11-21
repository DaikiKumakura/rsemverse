# rsemverse

<!-- badges: start -->
[![R-CMD-check](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DaikiKumakura/rsemverse/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24210459.svg)](https://doi.org/10.5281/zenodo.24210459)
<!-- badges: end -->

**From RSEM results to analysis-ready gene tables in R, for human and
mouse bulk RNA-seq.**

rsemverse reads RSEM `*.genes.results` (and `*.isoforms.results`) files into
one table of gene symbols by samples, then filters, selects, normalizes and
transforms that table, and finally turns it into a long table for plotting or
a file for other tools. Every function takes the
same table and returns the same table, so the steps chain in any order.

- **No gene lists to prepare.** Gene symbols and gene groups (protein-coding,
  mitochondrial, ribosomal, …) come from `org.Hs.eg.db` / `org.Mm.eg.db`.
- **No sample groups needed.** Every step works on the table alone; the sample
  table is used only to label samples, and to name the reference group when
  you ask for z-scores against controls.
- **Every step reports what it did,** with numbers: how many genes were kept,
  which threshold was used, the normalization factors.
- **Checked against the reference tools.** Filtering and normalization give
  the same results as edgeR, DESeq2 and limma, without depending on them.

**Documentation:** <https://daikikumakura.github.io/rsemverse/>

## Workflow

```text
  RSEM *.genes.results (+ metadata.csv, or found in the folder)
                 │
          rsem_read()              read: gene symbol × sample table (count, TPM, effective length)
                 │
          rsem_detection_curve()   how many genes each threshold keeps
          rsem_filter()            remove lowly expressed genes
          rsem_subset()            keep or remove gene groups (OrgDb or GTF)
          rsem_biotype_profile()   share of rRNA, lncRNA, … in each sample
          rsem_check()             problems in the RSEM files, all at once, before reading
          rsem_mapping_qc()        alignable, filtered, unique-gene reads from the .cnt files
          rsem_sample_qc()         PCA, sample correlation, detected genes per sample
          rsem_normalize()         normalize between samples
          rsem_transform()         log2, z-score, ranks
          rsem_fit_reference()     z-scores against a fixed control group
          rsem_apply_reference()
                 │
          rsem_long()              one row per gene and sample, with the sample table
          rsem_write()             TSV with gene IDs and gene types
                 │
         ggplot2, dplyr, bulkdeconv, …
```

Transcripts: `rsem_read(level = "isoform")` and `rsem_select_isoforms()`.
Effective lengths: `rsem_read(expr = "length")` and `rsem_length_profile()`.

| Step | Function | Choices (default first) |
|---|---|---|
| Read | `rsem_read()` | `expr = "count"`, `"tpm"`, `"length"` (`"isopct"` for transcripts); `species = "human"` or `"mouse"`; `metadata = NULL` finds the files (`.gz` too); `id = "original"` for any organism; `id = "annotation"` takes names from a GTF; `--append-names` IDs are understood; `"pme_count"`, `"pme_tpm"`, `"count_sd"`, `"tpm_ci_lower"`, `"tpm_ci_upper"` for `--calc-pme` and `--calc-ci` output |
| Thresholds | `rsem_detection_curve()` | thresholds × number or share of samples |
| Filter | `rsem_filter()` | `"edger"` (counts) / `"detected"` (TPM), `"nonzero"`, `"variable"` |
| Select genes | `rsem_subset()` | `"protein_coding"`, `"mito"`, `"ribo"`, `"hb"`, `"sex"`, `"ercc"`; with a GTF also `"rrna"`, `"mito_rrna"`, `"mito_trna"`, `"snorna"`, `"snrna"`, `"mirna"`, `"lncrna"` |
| Gene types | `rsem_biotype_profile()` | share of each gene type of a GTF per sample; `total =` for the unfiltered denominator |
| Sample QC | `rsem_sample_qc()` | PCA on the most variable genes, sample correlation (`"spearman"`, `"pearson"`, `"kendall"`) on detected, all or variable genes, per-sample summary |
| Check files | `rsem_check()` | errors, warnings and notes about the files: columns, values, IDs, zero lengths, different references |
| Read statistics | `rsem_mapping_qc()` | alignable, unalignable, filtered, unique-gene and multi-gene fractions from `.cnt` files |
| Normalize | `rsem_normalize()` | `"cpm"`, `"tmm"`, `"tmmwsp"`, `"upper_quartile"`, `"deseq2"`, `"quantile"` |
| Transform | `rsem_transform()` | `"log2"`, `"zscore"`, `"center"`, `"rank"` |
| Reference z-scores | `rsem_fit_reference()`, `rsem_apply_reference()` | reference samples or group; `zero_sd = "na"` or `"drop"` |
| Long table | `rsem_long()` | none: every column of the sample table is added as it is |
| Write | `rsem_write()` | original gene IDs and GTF gene types next to the values; `.gz` compresses |

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

Without a sample table, `rsem_read("data", metadata = NULL, ...)` finds every
`*.genes.results` file under the folder and names each sample by its file
name (`S1.genes.results` → `S1`) or, with `sample_from = "parent_dir"`, by its
folder (`S1/rsem.genes.results` → `S1`). It prints the files it found and
stops if two files give the same sample ID. Groups are never guessed from
file names.

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
original ID, so nothing is lost: the column sums equal the input. Rows that
share a symbol are summed, never reduced to their maximum, and on the
original scale, before any log. The IDs behind each row stay in
`attr(x, "gene_ids")`. It stops with a message when files are missing, genes
differ between files, or the species does not match.

For RefSeq-based references, whose IDs (`gene1000`) say little,
`id = "annotation", annotation = "reference.gtf"` names the rows by the
`gene_name` of the annotation, also for IDs that carry the name that
`--append-names` appends (`<gene_id>_<gene_name>`); the part before an
underscore must be a `gene_id` of the annotation and the part after it that
gene's `gene_name`, so `gene1_2` is not taken for `gene1` and underscores in
IDs or names are safe. Posterior estimates written by
`--calc-pme` or `--calc-ci` are read with `expr = "pme_count"`, `"pme_tpm"`,
`"count_sd"`, `"tpm_ci_lower"` and `"tpm_ci_upper"`; the posterior means add
up like counts, while standard deviations and credibility bounds need
`id = "original"` because they cannot be summed.

`expr = "length"` reads the effective lengths and reports genes whose length
is 0: RSEM gives 0 to transcripts shorter than the mean fragment length, and
`tximport`/`DESeqDataSetFromTximport()` stop on such genes.
`rsem_length_profile()` shows how much each gene's length varies between
samples. `level = "isoform"` reads the transcript files, with the
transcript-to-gene table in `attr(x, "tx2gene")`.

### Filter: `rsem_filter()`

| Method | Keeps a gene when | Typical use |
|---|---|---|
| `"edger"` (default for counts) | CPM reaches 10 counts in a median-sized sample, in all samples (or 10 + 70% of the rest above 10), and the total count is at least 15; edgeR `filterByExpr()`. All samples are one group unless `group` (or `design`) is given: then the smallest group sets the number of samples | before normalization or differential expression |
| `"detected"` (default for TPM) | the value reaches `min_value` (TPM 1, count 10) in `min_samples` samples (20% of them) | TPM tables; "expressed in at least n samples" |
| `"nonzero"` | it is above zero in any sample | removing empty rows only |
| `"variable"` | it is among the `n_top` (2000) with the largest variance of log2 | before PCA, clustering, heatmaps |

Before choosing a threshold, `rsem_detection_curve()` counts the genes each
threshold would keep, for several sample requirements at once; each row is
exactly the number `rsem_filter(method = "detected")` keeps with those
settings.

```r
rsem_detection_curve(tpm, "tpm", thresholds = c(0.5, 1, 2, 5), min_samples = c(1, 0.2, 0.5))
```

With 3 controls and 10 treated samples, the all-in-one-group rule asks for
expression in all 13 samples; `rsem_filter(counts, "count", group = condition)`
asks for the size of the smallest group, as `edgeR::filterByExpr(group = )`
does. Groups only set how many samples a gene must be expressed in; no gene
is chosen by comparing groups.

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

`"ribo"` is ribosomal **protein** genes, not ribosomal RNA. The OrgDb does
not separate RNA gene types, so rRNA, mitochondrial rRNA and tRNA, snoRNA,
snRNA, miRNA and lncRNA come from a GTF file, ideally the one the RSEM
reference was built from:

```r
rsem_subset(tpm, remove = c("rrna", "mito_rrna"), annotation = "gencode.v29.annotation.gtf.gz")
rsem_biotype_profile(tpm, "gencode.v29.annotation.gtf.gz")   # share of each gene type per sample
```

If the GTF has no gene of a requested type, `rsem_subset()` stops and
`rsem_biotype_profile()` reports `NA` instead of 0: an annotation without rRNA
genes cannot tell how much rRNA there was.

### Normalize: `rsem_normalize()`

| Method | Corrects for | Values |
|---|---|---|
| `"cpm"` | sequencing depth | counts per million |
| `"tmm"` | depth and RNA composition (edgeR) | CPM with TMM factors |
| `"tmmwsp"` | as TMM, for data with many zeros (edgeR `"TMMwsp"`) | CPM with TMMwsp factors |
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

### Z-scores against controls: `rsem_fit_reference()`, `rsem_apply_reference()`

`"zscore"` centres each gene on the mean of all samples in the table, so the
values change whenever samples are added or removed. To measure every sample
against a fixed group (controls, vehicle), fit the mean and standard
deviation on that group only and apply them to any table on the same scale,
including new samples:

```r
ref <- rsem_fit_reference(tpm, "control", metadata = "metadata.csv", group = "group")  # log2(TPM + 1)
z   <- rsem_apply_reference(tpm, ref)
z_new <- rsem_apply_reference(tpm_validation, ref)   # same reference, nothing refitted
```

Genes with no variation among the controls have no z-score: `NA` by default,
or removed with `zero_sd = "drop"`. At least two reference samples are needed.

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

### Check the files: `rsem_check()`

Before reading a new data set, `rsem_check("rsem")` lists every problem in the
files at once instead of stopping at the first: missing columns, non-numeric,
missing, infinite or negative values, duplicated IDs, genes with effective
length 0, a TPM total that is not 1e6, and, compared with the first file, IDs
that are missing, in another order, or listed with a different set of
transcripts (files quantified against different references). Each finding has
a severity (`ERROR`, `WARNING`, `INFO`), the sample, the number of rows and a
description. `rsem_read()` itself stops on the same errors, and its
`validate_reference = TRUE` (the default) checks the `transcript_id(s)` column
of gene-level files as sets; `allow_reorder = TRUE` accepts files with the same
IDs in another order.

### Sample QC: `rsem_sample_qc()`

Returns the tables behind a first look at the samples: principal components
(`pca_scores`, `pca_variance`, `pca_loadings`), the sample correlation and
distance matrices, and per sample the total, the number of detected genes and
how well it correlates with the others. The PCA uses the `n_top` most variable
genes of `log2(value + 1)`, centred and not scaled; the correlation uses its
own gene set (`correlation_genes`: genes detected in any sample, all genes, or
the PCA genes), because restricting it to variable genes changes the values a
great deal. No sample is flagged or removed.

```r
qc <- rsem_sample_qc(filtered, metadata = "metadata.csv", n_top = 2000)
qc$pca_variance          # share of variance per component
ggplot2::ggplot(qc$pca_scores, ggplot2::aes(PC1, PC2, colour = group)) + ggplot2::geom_point()
qc$sample_summary        # total, detected genes, median correlation with the others
```

### Read statistics: `rsem_mapping_qc()`

Reads the `.cnt` file RSEM writes for each sample (`<name>.stat/<name>.cnt`)
and returns the alignable, unalignable and filtered fractions, the share of
alignable reads that map to one gene or to several, and the alignments per
read. A file whose totals do not add up stops the function. The statistics
describe the alignments given to RSEM, not the abundance estimates.

```r
rsem_mapping_qc("rsem", metadata = NULL)
```

### Write: `rsem_write()`

Writes a TSV (gzip-compressed for `.gz`) with the gene, the original gene IDs
behind it (several joined by `;` where they were summed into one symbol), the
GTF gene type if an annotation is given, and the samples:

```r
tpm <- rsem_read("data", "metadata.csv", species = "human", expr = "tpm")
rsem_write(tpm, "tpm.tsv.gz", annotation = "gencode.v29.annotation.gtf.gz")
```

### Transcripts: `rsem_select_isoforms()`

Keeps transcripts by TPM, by their share of the gene's most expressed
transcript, or only the most expressed one, combining samples by `"any"`,
`"all"`, `"prevalence"` or `"mean"` (the rules of Trinity's
`filter_low_expr_transcripts.pl`). The kept values are the original estimates,
not re-quantified; removing transcripts by expression can remove important
ones, so use it only for a stated purpose.

## Recipes

| Goal | Steps |
|---|---|
| Heatmap of variable genes | counts → `rsem_filter()` → `rsem_normalize("tmm")` → `rsem_filter("variable")` → `rsem_transform(c("log2", "zscore"))` |
| PCA | counts → `rsem_filter()` → `rsem_subset(remove = c("mito", "ribo", "sex"))` → `rsem_normalize("tmm")` → `rsem_sample_qc()` (it takes the log itself) |
| Plot genes by group | any of the above → `rsem_long(metadata = "metadata.csv")` → ggplot2 |
| Cell-type deconvolution | `rsem_read(expr = "tpm")` straight into [bulkdeconv](https://daikikumakura.github.io/bulkdeconv/); do not filter, normalize or log-transform first |
| Differential expression | filter with `rsem_filter(group = ...)`, but test on the RSEM files through `tximport`, not on rounded counts (see below) |
| Treated versus controls, sample by sample | `rsem_read(expr = "tpm")` → `rsem_filter()` → `rsem_fit_reference(group = ...)` → `rsem_apply_reference()` |
| rRNA and mitochondrial share per sample | `rsem_read()` → `rsem_biotype_profile(annotation = gtf)` |
| Share a table with IDs and gene types | `rsem_read()` → `rsem_write(annotation = gtf)` |

For the heatmap, choose the genes on the normalized values and scale them
last (`"variable"` takes the log itself):

```r
counts <- rsem_read("data", "metadata.csv", species = "human", expr = "count")
tmm    <- rsem_normalize(rsem_filter(counts, "count"), "count", "tmm")
top    <- rsem_filter(tmm, "count", method = "variable", n_top = 500)
heat   <- rsem_transform(top, c("log2", "zscore"))
```

## Differential expression

RSEM's expected counts are not integer counts of a fixed gene length: a
gene's effective length differs between samples when the isoform usage or the
fragment length distribution differs. DESeq2 and edgeR can model this with a
length offset, which `tximport` supplies. Rounding `rsem_read(expr = "count")`
and passing it to DESeq2 ignores the offset, and the tximport documentation
does not recommend using the estimated counts without it; the exceptions it
names are protocols whose counts do not depend on gene length, such as 3' tag
RNA-seq, so the choice rests on the library preparation and the model, not on
how similar the isoform usage looks. The route described in the tximport and
DESeq2 documentation reads the RSEM files directly:

```r
library(tximport)
library(DESeq2)
files <- setNames(file.path("data", samples$ref), samples$sample_id)
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
```

`DESeqDataSetFromTximport()` stops when a gene has effective length 0.
`rsem_read(expr = "length")` reports those genes and `rsem_length_profile()`
lists them; whether to remove them or to import the transcript-level files
(`txIn = TRUE`) is a decision for the analysis, and replacing 0 by 1 changes
the offset of those genes. Use rsemverse's tables for filtering decisions,
plots and scores, and the files themselves for the test.

## Validation

- **Tests:** the filtering and normalization results (TMMwsp included) equal
  those of edgeR 4.8.2, DESeq2 1.50.2 and limma 3.66.0 on simulated counts.
- **Real data:** on six ENCODE K562 RSEM files (59,526 genes each),
  `rsem_read()` gave symbols to 36,385 genes, kept 23,141 under their ID and
  preserved the column sums. The edgeR rule kept the same 12,808 genes as
  `edgeR::filterByExpr()`, every normalization was within 3e-16 (relative) of
  the reference tools, and all 37 mitochondrial genes were found. The
  functions added in 0.3.0 were checked on the same experiments, with their
  transcript files and the GENCODE V29 annotation. See
  [`validation/`](https://github.com/DaikiKumakura/rsemverse/tree/main/validation).
- **Platforms:** R CMD check passes on Linux, macOS and Windows in CI.

## Related tools

- [bulkdeconv](https://daikikumakura.github.io/bulkdeconv/) estimates cell-type
  composition from the TPM table returned by `rsem_read()`.
- [rsemflow](https://daikikumakura.github.io/rsemflow/) reads the same RSEM
  files to run DESeq2, GSEA and GSVA from the command line.

## Citation

Cite the archived version: <https://doi.org/10.5281/zenodo.24210459> (all versions; version 0.3.0 is <https://doi.org/10.5281/zenodo.24382507>). The
citation metadata is in [`CITATION.cff`](https://github.com/DaikiKumakura/rsemverse/blob/main/CITATION.cff).

## License

MIT © Daiki Kumakura
