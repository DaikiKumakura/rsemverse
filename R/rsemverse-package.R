#' rsemverse: RSEM results in R, one step at a time, with the pipe
#'
#' rsemverse reads RSEM `*.genes.results` and `*.isoforms.results` files, checks them, and turns them
#' into tables for plots, scores and downstream tools. Every step takes a table and returns a table,
#' so the steps chain with `|>`, and each prints one line about what it did. The table remembers what
#' it holds (`expr`), the original gene IDs, the normalization factors and the steps, so a step needs
#' no repeated `expr = `; [rsem_steps()] lists what was done.
#'
#' ```
#' rsem_read("rsem", species = "human", expr = "count") |>
#'   rsem_filter() |> rsem_normalize(method = "tmm") |> rsem_transform("log2") |> rsem_sample_qc()
#' ```
#'
#' @section The functions:
#' * Read: [rsem_read()] (one value), [rsem_import()] (several values).
#' * Prepare: [rsem_filter()], [rsem_subset()], [rsem_normalize()], [rsem_transform()],
#'   [rsem_fit_reference()] with [rsem_apply_reference()].
#' * Check the files: [rsem_check()], [rsem_reconcile()], [rsem_id_report()].
#' * QC: [rsem_sample_qc()], [rsem_qc_summary()] (and its parts [rsem_mapping_qc()],
#'   [rsem_alignment_profile()], [rsem_model_qc()], [rsem_noise_qc()]), [rsem_biotype_profile()],
#'   [rsem_length_profile()], [rsem_detection_curve()].
#' * Transcripts: [rsem_select_isoforms()], [rsem_isoform_usage()].
#' * Hand over: [rsem_long()], [rsem_write()], [rsem_as_tximport()], [rsem_as_summarizedexperiment()],
#'   [rsem_steps()].
#'
#' No gene list and no sample grouping are needed except where a function compares with a reference
#' group. Filtering and normalization give the same results as edgeR, DESeq2 and limma without
#' depending on them.
#'
#' @keywords internal
"_PACKAGE"
