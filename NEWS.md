# rsemverse 0.3.0

Steps that scripts starting from RSEM output often write by hand, done with
the checks those scripts tend to skip. Existing calls give the same results
as in 0.2.0.

## Reading

* `rsem_read(metadata = NULL)` finds the `*.genes.results` files under a
  folder and names samples by file name or, with `sample_from =
  "parent_dir"`, by folder. It prints the files found and stops on duplicated
  sample IDs.
* `rsem_read(expr = "length")` reads effective lengths and reports genes with
  length 0, the values that stop `tximport`/`DESeqDataSetFromTximport()`.
  Rows merged into one symbol get the TPM-weighted mean length.
* `rsem_read(level = "isoform")` reads `*.isoforms.results` (TPM, counts,
  lengths or `IsoPct`), with the transcript-to-gene table in
  `attr(x, "tx2gene")`.
* `rsem_read(id = "original")` keeps the IDs of the files and needs no OrgDb.
* The original IDs behind each row are kept in `attr(x, "gene_ids")`.

## New functions

* `rsem_detection_curve()`: genes kept for each detection threshold and
  sample requirement, equal to what `rsem_filter(method = "detected")` keeps.
* `rsem_fit_reference()`, `rsem_apply_reference()`, `rsem_reference_scale()`:
  per-gene z-scores against a fixed reference group, fitted once and applied
  to any samples. Genes with no variation in the reference give `NA` or are
  removed, never smoothed.
* `rsem_biotype_profile()`: share of each GTF gene type per sample, with `NA`
  for types the annotation does not contain.
* `rsem_write()`: TSV (or `.tsv.gz`) with the original gene IDs and GTF gene
  types next to the values.
* `rsem_select_isoforms()`: the transcript selection rules of Trinity's
  `filter_low_expr_transcripts.pl`, with an explicit rule for combining
  samples.
* `rsem_length_profile()`: variation of each gene's effective length across
  samples and the share of samples with length 0.

## Changed functions

* `rsem_subset()` takes a GTF file (`annotation`) for the groups `"rrna"`,
  `"mito_rrna"`, `"mito_trna"`, `"snorna"`, `"snrna"`, `"mirna"` and
  `"lncrna"`, and stops when the annotation has no gene of a requested type.
  Rows are matched to the GTF through the original gene IDs that
  `rsem_read()` keeps, so genes whose OrgDb symbol differs from the GTF gene
  name are not lost. The documentation now says that `"ribo"` means
  ribosomal protein genes.
* `rsem_normalize(method = "tmmwsp")`: TMMwsp factors, equal to
  `edgeR::calcNormFactors(method = "TMMwsp")` (edgeR 4.8.2) in the tests.

## Documentation

* Two articles: a workflow from RSEM files to z-scores against controls, and
  a table of hand-written RSEM steps with their rsemverse equivalents.
* `validation/encode-k562-0.3.R` checks the new functions on the six ENCODE
  K562 experiments, their transcript files and the GENCODE V29 annotation.

# rsemverse 0.2.0

New functions that take the table from `rsem_read()` and return a table of
the same layout, so they can be chained. None of them needs a gene list or a
sample grouping from the user.

* `rsem_filter()` removes genes: the edgeR `filterByExpr()` rule with all
  samples as one group (default for counts), a value threshold in a number or
  share of samples (default for TPM: TPM >= 1 in 20% of samples), genes that
  are zero everywhere, or all but the most variable genes.
* `rsem_normalize()` normalizes between samples: CPM, TMM, upper quartile and
  DESeq2 median of ratios for counts, and quantile normalization. The tests
  check that the results equal those of edgeR 4.8.2, DESeq2 1.50.2 and limma
  3.66.0; the package does not depend on them.
* `rsem_transform()` applies log2 (with a pseudocount), per-gene z-scores or
  centering, and within-sample ranks, in the order given.
* `rsem_subset()` keeps or removes built-in gene groups: protein-coding genes
  (gene type from the OrgDb), mitochondrial-genome genes, cytosolic ribosomal
  protein genes, hemoglobin genes, chromosome Y genes with XIST, and ERCC
  spike-ins. Mitochondrial rRNA and tRNA genes, which have no Ensembl link in
  the OrgDb and keep their Ensembl ID in `rsem_read()`, are found through a
  built-in list of the 37 mitochondrial Ensembl IDs.
* `rsem_long()` returns one row per gene and sample with every column of the
  sample table added as it is (numbers stay numeric, other columns stay text).

On the six ENCODE K562 files in `validation/`, every function ran on all
59,427 rows, and the normalized values were within 3e-16 (relative) of edgeR,
DESeq2 and limma. See `validation/README.md`.

# rsemverse 0.1.0

First version.

* `rsem_read()` reads RSEM `*.genes.results` files listed in a metadata CSV
  (`sample_id`, `ref`) and returns expected counts or TPM per gene symbol and
  sample, using `org.Hs.eg.db` (human) or `org.Mm.eg.db` (mouse). Ensembl
  versions are removed before the lookup, and rows that share a symbol are
  summed.
