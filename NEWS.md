# rsemverse 0.5.1

* `rsem_read(id = "annotation")` is stricter: an ID with an appended name is
  accepted only as `<gene_id>_<gene_name>`, with the part before an
  underscore a `gene_id` of the annotation and the part after it exactly that
  gene's `gene_name`. Before, an ID such as `gene1_2` was read as the gene
  `gene1` whenever the annotation had `gene1`. IDs that do not fit keep their
  ID.
* More tests for `rsem_check()` and `rsem_read(id = "annotation")`: the order
  of transcripts inside a gene, the TPM total tolerance, zero effective
  length with a positive count (reported, never changed), duplicated or
  unnamed genes in the annotation.

# rsemverse 0.5.0

Protection against comparing files that should not be compared, more of
RSEM's own output, and a check that reports problems before reading.

## Changes

* `rsem_read()` checks that every gene-level file lists the same transcripts
  for each gene (`transcript_id(s)`, compared as sets): files quantified
  against different references now stop the function
  (`validate_reference = TRUE`, the default; skipped when the column is
  missing). Isoform files were already checked for the same gene per
  transcript.
* `rsem_read(allow_reorder = TRUE)` accepts files with the same IDs in a
  different order and puts them in the order of the first file; the default
  still stops, and the message says when only the order differs.
* `rsem_read(id = "annotation", annotation = )` names rows by the `gene_name`
  of a GTF file or data frame, also for IDs with the name that
  `--append-names` appends (as in RefSeq-based references: `gene1000_TP53`).
  An ID is shortened only to an ID that the annotation has, so underscores in
  IDs or names are never cut at the wrong place. Rows without a name keep
  their ID.
* `rsem_read()` reads RSEM's posterior estimates: `expr = "pme_count"`,
  `"pme_tpm"` (posterior means, which add up over rows like counts),
  `"count_sd"`, `"tpm_ci_lower"` and `"tpm_ci_upper"` (these cannot be summed
  and need `id = "original"` or `level = "isoform"`). A file without the
  column says whether to use `--calc-pme` or `--calc-ci`.
* `rsem_sample_qc()` stops when the table is already log2-transformed (the
  history that `rsem_transform()` records now also survives `rsem_filter()`,
  `rsem_subset()`, `rsem_normalize()` and the reference z-scores) instead of
  transforming twice; use `transform = "none"`.
* README: the differential expression section no longer allows the
  length-offset shortcut because the isoform usage looks similar; the
  tximport documentation names protocols without a gene-length bias, such as
  3' tag RNA-seq, as the exception.

## New function

* `rsem_check()` reports problems in the RSEM files (columns, values, IDs,
  zero effective lengths, the TPM total, IDs and transcripts that differ from
  the first file) for all files at once, with a severity for each finding.

# rsemverse 0.4.0

Fixes found when checking the 0.3.0 functions against each other, and three
additions.

## Fixes

* The gene IDs behind each row (`attr(x, "gene_ids")`) and the
  transcript-to-gene table (`attr(x, "tx2gene")`) are now carried through
  `rsem_filter()`, `rsem_subset()`, `rsem_normalize()`, `rsem_transform()`
  and `rsem_apply_reference()`, for the rows that remain. Before, they were
  dropped, so `rsem_write()` after a filter had no original IDs and GTF gene
  types were matched by symbol only.
* `rsem_read()` understands the IDs of `rsem-calculate-expression
  --append-names` (`ENSG00000141510.18_TP53`): the version and the appended
  name are removed before the symbol lookup. `_PAR_Y` IDs keep their ID as
  before.
* `rsem_read(id = "original")` no longer needs `species` and no longer checks
  that the IDs are Ensembl IDs, so results for any organism or a de novo
  transcriptome can be read.
* Gzip-compressed result files (`*.genes.results.gz`) are found and read,
  with a sample table or without, and without extra packages.
* Transcript files whose gene differs for the same transcript between samples
  (a different reference) now stop `rsem_read(level = "isoform")`.
* `rsem_biotype_profile()` gives `NA` shares for a sample whose total is 0,
  and takes `total` for the denominator, so rRNA or mitochondrial shares can
  be reported against the unfiltered total; the documentation says that the
  default denominator is the total of the genes in `x`.
* README: the differential expression recipe no longer suggests rounding
  expected counts as the route to DESeq2; a new section describes the
  `tximport` route, the effective-length offset and the zero-length genes.
* The tests that depended on the symbols of one org.Hs.eg.db release compute
  their expectations from the installed release.

## Additions

* `rsem_filter(method = "edger", group = , design = )` sets the number of
  samples a gene must be expressed in from the groups or the design, as
  `edgeR::filterByExpr()` does (checked against edgeR 4.8.2 in the tests).
* `rsem_sample_qc()` returns the PCA (scores, variance, loadings), the sample
  correlation and distance matrices, and a per-sample summary (total,
  detected genes, median correlation with the other samples). The PCA uses
  the most variable genes and the correlation its own gene set; nothing is
  flagged or removed. Equal to `stats::prcomp()` and `stats::cor()` in the
  tests.
* `rsem_mapping_qc()` reads the `.cnt` statistics files of RSEM: alignable,
  unalignable and filtered fractions, unique-gene and multi-gene fractions of
  the alignable reads, and alignments per read. It stops on files whose
  totals do not add up.

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
