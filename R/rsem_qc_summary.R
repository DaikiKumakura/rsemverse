#' One table of the quality numbers of every sample
#'
#' Puts the numbers that RSEM wrote about each sample next to each other: how the reads were assigned
#' (`.cnt`, [rsem_mapping_qc()]), what model was learned (`.model`, [rsem_model_qc()]), the noise share
#' (`.theta`, [rsem_noise_qc()]) and a few numbers of the gene results. It repeats no calculation of
#' those functions; it joins their results by sample and says whether the files agree.
#'
#' The `.cnt`, `.model` and `.theta` files are looked for under `data` (RSEM writes them in
#' `<name>.stat/`); a kind that is not there is left out with a message, so the table also works for
#' results that were copied without their `.stat` folders. The samples are named after the files, as in
#' [rsem_read()], so `A.genes.results` and `A.stat/A.cnt` are the sample `A`.
#'
#' @inheritParams rsem_mapping_qc
#' @return A data frame with one row per sample (the samples of the gene results): `sample_id`;
#'   from the gene results `genes`, `tpm_total`, `genes_zero_length` (effective length 0) and
#'   `genes_detected` (TPM of at least 1); from the `.cnt` file `n_total`, `alignable_fraction`,
#'   `unique_gene_fraction`, `uncertain_fraction`, `alignments_per_read` and `read_type`; from the
#'   `.model` file `model_type`, `forward_prob`, `length_mean`, `length_sd` and `rspd_estimated`; from the
#'   `.theta` file `noise_fraction`; and `files_agree`, `FALSE` when the read type of the `.cnt` file is
#'   not the model type of the `.model` file (`NA` when one of them is missing).
#' @examples
#' \dontrun{
#' "rsem" |> rsem_qc_summary()
#' }
#' @export
rsem_qc_summary <- function(data = "data", metadata = NULL, sample_from = "filename") {
  sample_from <- .choice(sample_from, c("filename", "parent_dir"), "sample_from")
  samples <- suppressMessages(
    if (is.null(metadata)) .discover_samples(data, "gene", sample_from) else .read_samples(data, metadata)
  )
  rows <- lapply(seq_len(nrow(samples)), function(i) {
    x <- .fread_rsem(samples$path[[i]], select = c("TPM", "effective_length"))
    data.frame(
      sample_id = samples$sample_id[[i]], genes = nrow(x), tpm_total = sum(x$TPM), genes_zero_length = sum(x$effective_length == 0),
      genes_detected = sum(x$TPM >= 1), stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  part <- function(label, f) {
    r <- tryCatch(suppressMessages(f()), error = function(e) NULL)
    if (is.null(r)) message(sprintf("rsem_qc_summary: no %s files found under '%s', left out", label, data))
    r
  }
  join <- function(out, tab, keep) {
    if (is.null(tab)) return(out)
    tab <- tab[, c("sample_id", keep), drop = FALSE]
    merge(out, tab, by = "sample_id", all.x = TRUE, sort = FALSE)
  }
  cnt <- part(".cnt", function() rsem_mapping_qc(data, sample_from = sample_from))
  mod <- part(".model", function() rsem_model_qc(data, sample_from = sample_from)$sample_summary)
  th <- part(".theta", function() rsem_noise_qc(data, sample_from = sample_from))
  out <- join(out, cnt, c("n_total", "alignable_fraction", "unique_gene_fraction", "uncertain_fraction", "alignments_per_read", "read_type"))
  out <- join(out, mod, c("model_type", "forward_prob", "length_mean", "length_sd", "rspd_estimated"))
  out <- join(out, th, "noise_fraction")
  out$files_agree <- if (all(c("read_type", "model_type") %in% names(out))) {
    ifelse(is.na(out$read_type) | is.na(out$model_type), NA, out$read_type == out$model_type)
  } else {
    NA
  }
  out <- out[match(samples$sample_id, out$sample_id), , drop = FALSE]
  rownames(out) <- NULL
  message(sprintf(
    "rsem_qc_summary: %d samples; %s", nrow(out),
    paste(c("genes", if (!is.null(cnt)) ".cnt", if (!is.null(mod)) ".model", if (!is.null(th)) ".theta"), collapse = ", ")
  ))
  if (isFALSE(all(out$files_agree %in% TRUE | is.na(out$files_agree)))) {
    warning("The read type of the .cnt file differs from the model type of the .model file for: ",
            paste(out$sample_id[which(out$files_agree %in% FALSE)], collapse = ", "), call. = FALSE)
  }
  out
}
