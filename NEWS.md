# rsemverse 0.8.0

Made for the pipe: steps chain without repeating what the table holds, and the table says how it was made.
Also the gaps that two external reviews of 0.7.2 and our own tests showed.

## The pipe

* `rsem_read()` records what the table holds (`count`, `tpm`, ...; the attribute `"expr"`) and the first
  step. `rsem_filter()`, `rsem_detection_curve()` and `rsem_normalize()` take it as the default of `expr`, so
  `rsem_read(..., expr = "count") |> rsem_filter() |> rsem_normalize(method = "tmm")` needs no second
  `expr = `. A normalized or log-transformed table has no `expr`; `rsem_filter(method = "variable")`
  and `"nonzero"` do not need one.
* `rsem_filter()`, `rsem_subset()`, `rsem_normalize()`, `rsem_transform()` and the reference z-scores
  add a line to the attribute `"steps"`, and the table keeps the normalization `"factors"` through them
  (before, a transformation lost them). `rsem_steps()` returns the steps.
* `rsem_sample_qc()` takes `transform = NULL` as the default: log2, or none when the table is already
  log2-transformed. Asking for `"log2"` on such a table still stops.

## New functions

* `rsem_qc_summary()` joins the quality numbers of the `.cnt`, `.model` and `.theta` files and of the gene
  results by sample, and warns when the read type of a `.cnt` file differs from the model type of the `.model`.
* `rsem_isoform_usage()` gives, per gene and sample, the share of the dominant transcript, the number of
  expressed transcripts, the entropy and the effective number of isoforms; no TPM gives `NA`.
* `rsem_as_summarizedexperiment()` wraps `rsem_import()` values, the sample table and the gene IDs.

## Fixes

* GTF names and types: the row type decides (gene, then transcript, then the others), also when the
  rows are in other chunks of the file.
* `rsem_model_qc()` stops on an optional read length flag other than 0 or 1 and on bounds that are not
  whole numbers.
* `rsem_alignment_profile()` stops on negative, fractional or repeated numbers of alignments and on a
  histogram that does not match `N0` and `N2`, and warns when the alignments in it are not `nHits`
  (they are in RSEM's own `.cnt`).
* `rsem_reconcile()` also compares `posterior_mean_count`, `pme_TPM` and `pme_FPKM` (the values that add
  up), and `max_diff_gene` is the gene of `max_abs_diff`.
* `rsem_check()` reports an empty or missing `transcript_id(s)` (gene files) or `gene_id` (transcript files).

## Documentation

* The README is rewritten: one page, the pipe first, a table of the functions, recipes that were run on real files.

# rsemverse 0.7.2

Gaps in 0.7.1 found by testing against RSEM's own output.

## New function

* `rsem_noise_qc()` reads the `.theta` files: the noise share (first value of the third line) and the
  uncorrected value of the second line. In a run of RSEM with 600 added unalignable reads of 2,750, the
  second-line value is `600 / 2750` and the third-line value `600 / 2700` (the 50 reads of a
  transcript with effective length 0 left out), so the third line is the corrected one. The share is
  not the unalignable fraction of `rsem_mapping_qc()` in general.

## Changes

* A GTF with gene rows for some genes only is read completely: a gene without a gene row, or with a
  gene row without a name or type, gets them from its transcript rows and then its other rows (before,
  the exon rows were used only when there was no gene row at all). The rows of genes that are already
  complete are not parsed, so GENCODE V29 takes 24 seconds (16 of them reading the file) and gives the same
  58,721 genes.
* `rsem_model_qc()` stops on a strand probability outside 0 to 1, on any probability outside 0 to 1 in
  the length distribution or the read start positions, and on a number of bins that is not a positive
  integer; its documentation says that `forward_prob` is the strand probability the run was given, not
  one learned from the reads.

# rsemverse 0.7.1

Checked against the output of RSEM itself.

## New function

* `rsem_model_qc()` reads the `.model` files: the fragment (or read) length distribution with its mean,
  standard deviation and mode, the strand probability, and the read start position distribution when RSEM
  estimated it (`--estimate-rspd`).

## Changes

* `rsem_read(id = "annotation")` and the GTF functions take names and types from the transcript or exon rows
  when the GTF has no gene rows (RSEM needs only exon rows).

## Checked against RSEM's own output

RSEM 1.3.3 was built and run (`--calc-pme --calc-ci --append-names`, and `--estimate-rspd`) on a toy data
set; its files are in `tests/testthat/fixtures/rsem-real/`. All new columns, the `.cnt` file, the `.model`
file and the names after `--append-names` are read from them in the tests. Two things the
documentation did not say: the coefficient of quartile variation of a gene without reads is written as 0
(not nan), and the number of RSPD bins and their probabilities follow the RSPD flag in the `.model` file
only when the RSPD was estimated.

Not done: the `.theta` file (its noise fraction) and the sequencing error sections of the `.model`.

# rsemverse 0.7.0

Found by comparing the package with the RSEM documentation and source and with public
workflows (100 searches, November 23, 2026).

## New functions

* `rsem_import()` reads several values of the same files at once (expected count, TPM and
  effective length by default) into tables with the same genes and samples, and
  `rsem_as_tximport()` turns them into the list that `DESeq2::DESeqDataSetFromTximport()`
  takes, without changing a value. On the six ENCODE K562 runs it equals
  `tximport::tximport(type = "rsem", txIn = FALSE)` (largest difference 0 for counts and
  lengths). Genes with an effective length of 0 stop it unless `zero_length = "drop"` or
  `"keep"` is chosen; rows that combine several gene IDs are refused.
* `rsem_alignment_profile()` returns the histogram of alignments per read that follows the
  first three lines of a `.cnt` file (the number of reads with 0, 1, 2, ... alignments).

## Changes

* `rsem_read()` also reads `FPKM`, `pme_FPKM`, `FPKM_ci_lower_bound`, `FPKM_ci_upper_bound`,
  `TPM_coefficient_of_quartile_variation`, `FPKM_coefficient_of_quartile_variation` and, for
  transcripts, `IsoPct_from_pme_TPM` (`expr = "fpkm"`, `"pme_fpkm"`, `"fpkm_ci_lower"`,
  `"fpkm_ci_upper"`, `"tpm_cqv"`, `"fpkm_cqv"`, `"isopct_pme"`). The interval, quartile
  variation and percentage cannot be summed and need `id = "original"`. A quartile variation
  that RSEM could not compute (written as nan) is `NA`. FPKM equals `1e3 / l_bar * TPM` of the
  documentation on the ENCODE files (median ratio 1.0000).
* `rsem_check()` reports posterior columns with missing, infinite or negative values and an
  interval whose lower bound is above its upper bound, and, for transcript files, a
  transcript whose sequence length differs from the first file (the gene-level length is an
  abundance-weighted mean and may differ). On the ENCODE transcript files no length differs.
* `rsem_reconcile()` has the column `transcript_mapping_checked`, and says when the gene file
  has no `transcript_id(s)` column, so that the transcripts could not be compared.
* `rsem_mapping_qc()` stops on an infinite value, a read type other than 0 to 3 and more
  uncertain reads than alignable reads, and no longer prints `Inf` when no sample has
  alignable reads.
* `rsem_sample_qc()` adds `rle_median` and `rle_iqr` to `sample_summary`: relative log
  expression (the deviation of a sample's log2 values from the median of each gene), `NA`
  unless the values are log2.
* `rsem_write(format = "gct")` writes the GCT 1.2 format of GSEA and ssGSEA.


# rsemverse 0.6.2

## Changes

* `rsem_long()` and `rsem_write()` keep missing values. The z-scores of
  `rsem_apply_reference()` are `NA` for a gene with no spread in the reference
  group; the two functions refused such a table, so the step from reference z-scores to
  a long table or a file failed. The other functions still stop on missing values.
* `rsem_reconcile()` stops on negative values (RSEM writes none: a pair of files with the
  same negative value in both would have counted as in agreement), on empty or missing IDs,
  and on a file without the `gene_id` or `transcript_id` column.

# rsemverse 0.6.1

Gaps in the checks of 0.6.0 and two older ones, found by reading the code.

## Changes

* `rsem_reconcile()` also compares the transcripts: the `transcript_id(s)` that the
  gene file lists for a gene must be the transcripts of that gene in the isoform file
  (new column `transcript_mapping_mismatch`). A transcript that is renamed leaves the
  sums unchanged, so the values alone could not show a different reference. It now
  stops on missing, infinite or non-numeric values, on duplicated IDs, and on an
  infinite `tolerance_extra`, instead of returning `NA`.
* `rsem_read()` stops on infinite or negative values in the column it reads (RSEM
  writes none; they mean a damaged file). Before, only missing values stopped it.
* `rsem_fit_reference()` and `rsem_apply_reference()` stop when the table is already
  log2-transformed and the reference applies `log2` itself, as `rsem_sample_qc()`
  does, instead of taking the logarithm twice.
* Not changed: the `rsem_check()` warning for the TPM total still uses 1,000. The
  largest difference that rounding can cause (0.005 per row) is a worst case; the
  typical error is about 0.003 times the square root of the number of rows (3 for a
  million rows), and a total off by 1,500 is not explained by rounding in practice.

# rsemverse 0.6.0

Two checks that use the structure of RSEM's own output, and a report of merged IDs.

## New functions

* `rsem_reconcile()` compares `*.genes.results` and `*.isoforms.results` of the
  same samples: the gene value of `expected_count`, `TPM` and `FPKM` must equal the
  sum over its transcripts, within the rounding that two-decimal output allows
  (`0.005 * (n + 1)` for `n` transcripts). On six ENCODE K562 runs no gene is beyond
  rounding (largest difference 0.05); with the isoform file of another sample, 20,405
  of 59,526 genes are.
* `rsem_id_report()` lists the rows of a table from `rsem_read()` that combine several
  gene IDs, with the IDs. On the ENCODE files, 89 rows.

## Changes

* The `rsem_check()` warning for a TPM total away from 1e6 now states the difference and
  the largest difference that two-decimal rounding of that many rows can cause.
  The threshold of 1,000 is unchanged: on the ENCODE files the total is within 0.1 of 1e6.

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
