#' rsemverse: from RSEM results to analysis-ready gene tables
#'
#' rsemverse turns RSEM `*.genes.results` (and `*.isoforms.results`) files
#' into one table of gene symbols by samples and prepares that table for
#' analysis. Every function takes the same table and returns the same table,
#' so the steps chain in any order, and each prints what it did.
#'
#' @section Workflow:
#' 1. [rsem_read()]: read counts, TPM or effective lengths into a gene symbol
#'    by sample table; find the files without a sample table if wanted.
#' 2. [rsem_detection_curve()] and [rsem_filter()]: see how many genes each
#'    detection threshold keeps, then remove lowly expressed genes (edgeR
#'    `filterByExpr()` rule, detection thresholds) or keep the most variable.
#' 3. [rsem_subset()] and [rsem_biotype_profile()]: keep or remove gene groups
#'    taken from the OrgDb or from a GTF file (rRNA, snoRNA, lncRNA, ...),
#'    and see what share of each sample they take. [rsem_mapping_qc()] reads
#'    the read assignment statistics of the RSEM `.cnt` files, and
#'    [rsem_sample_qc()] gives PCA, sample correlation and a per-sample summary.
#' 4. [rsem_normalize()]: normalize between samples (CPM, TMM, TMMwsp, upper
#'    quartile, DESeq2 median of ratios, quantile).
#' 5. [rsem_transform()] and [rsem_fit_reference()] /
#'    [rsem_apply_reference()]: log2, z-score, centering, ranks, and z-scores
#'    against a fixed reference group.
#' 6. [rsem_long()] and [rsem_write()]: one row per gene and sample for
#'    ggplot2, or a table with gene IDs and annotation for other tools.
#'
#' Transcripts: [rsem_read()] with `level = "isoform"` and
#' [rsem_select_isoforms()]. Effective lengths: [rsem_length_profile()].
#'
#' No gene list and no sample grouping are needed except where a function
#' compares with a reference group. Filtering and normalization give the same
#' results as edgeR, DESeq2 and limma without depending on them.
#'
#' @keywords internal
"_PACKAGE"
