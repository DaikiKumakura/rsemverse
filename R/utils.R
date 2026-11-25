# Shared helpers for the functions that work on a gene by sample table.

# Turn the table into a numeric matrix with genes as row names. Accepts the
# data frame rsem_read() returns (symbols as row names), a data frame whose
# first column holds the gene names, or a numeric matrix with row names.
.expr_matrix <- function(x, nonnegative = TRUE, arg = "x", allow_na = FALSE) {
  if (is.data.frame(x)) {
    if (ncol(x) >= 1L && !is.numeric(x[[1L]])) {
      genes <- as.character(x[[1L]])
      x <- x[-1L]
    } else {
      genes <- rownames(x)
    }
    numeric_col <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_col)) {
      .stopf(
        "`%s` must have numeric sample columns; not numeric: %s",
        arg, paste(names(x)[!numeric_col], collapse = ", ")
      )
    }
    m <- as.matrix(x)
    if (!length(m)) m <- matrix(numeric(0), nrow(x), ncol(x), dimnames = list(NULL, names(x)))
    rownames(m) <- genes
  } else if (is.matrix(x) && is.numeric(x)) {
    m <- x
  } else {
    .stopf("`%s` must be a data frame or a numeric matrix, such as the table rsem_read() returns.", arg)
  }
  storage.mode(m) <- "double"

  if (!ncol(m)) .stopf("`%s` has no sample columns.", arg)
  if (!nrow(m)) .stopf("`%s` has no genes.", arg)
  genes <- rownames(m)
  if (is.null(genes) || anyNA(genes) || any(!nzchar(genes))) {
    .stopf("`%s` needs gene names, as row names or as the first column.", arg)
  }
  if (anyDuplicated(genes)) {
    .stopf("`%s` has duplicated gene names: %s", arg, paste(utils::head(unique(genes[duplicated(genes)]), 5L), collapse = ", "))
  }
  samples <- colnames(m)
  if (is.null(samples) || anyNA(samples) || any(!nzchar(samples))) {
    .stopf("`%s` needs sample names as column names.", arg)
  }
  if (anyDuplicated(samples)) {
    .stopf("`%s` has duplicated sample names: %s", arg, paste(unique(samples[duplicated(samples)]), collapse = ", "))
  }
  if (!allow_na && anyNA(m)) .stopf("`%s` contains missing values.", arg)
  if (any(is.infinite(m))) .stopf("`%s` contains infinite values.", arg)
  if (nonnegative && any(m < 0, na.rm = TRUE)) {
    .stopf("`%s` contains negative values; it must be counts or TPM, not log-transformed or scaled values.", arg)
  }
  m
}

# The same layout as rsem_read(): a data frame with genes as row names.
# `from` is the table the values came from: the gene IDs behind each row
# (attr "gene_ids") and the transcript-to-gene table (attr "tx2gene") that
# rsem_read() keeps are carried over for the rows that remain.
.as_table <- function(m, from = NULL, step = NULL, keep_expr = TRUE) {
  res <- as.data.frame(m, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(from)) return(res)
  # What the table is (expr) and what was done to it (steps, factors) go along the pipe; the factors are
  # per sample, so they stay valid when rows are removed.
  if (keep_expr && !is.null(attr(from, "expr", exact = TRUE))) attr(res, "expr") <- attr(from, "expr", exact = TRUE)
  if (!is.null(attr(from, "factors", exact = TRUE))) attr(res, "factors") <- attr(from, "factors", exact = TRUE)
  st <- c(attr(from, "steps", exact = TRUE), step)
  if (length(st)) attr(res, "steps") <- st
  tf <- attr(from, "transform", exact = TRUE)
  if (!is.null(tf)) attr(res, "transform") <- tf
  g <- attr(from, "gene_ids", exact = TRUE)
  if (is.data.frame(g)) {
    g <- g[g$gene %in% rownames(res), , drop = FALSE]
    rownames(g) <- NULL
    attr(res, "gene_ids") <- g
  }
  tg <- attr(from, "tx2gene", exact = TRUE)
  if (is.data.frame(tg)) {
    tg <- tg[tg$transcript_id %in% rownames(res), , drop = FALSE]
    rownames(tg) <- NULL
    attr(res, "tx2gene") <- tg
  }
  res
}

.fmt <- function(n) {
  formatC(n, format = "d", big.mark = ",")
}

.number <- function(x, name, lower = -Inf, integer = FALSE, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < lower || x > upper ||
      (integer && x != round(x))) {
    .stopf(
      "`%s` must be a single %s of at least %s%s.", name, if (integer) "whole number" else "number",
      format(lower), if (is.finite(upper)) paste0(" and at most ", format(upper)) else ""
    )
  }
  x
}

# Number of samples a gene must pass: a value below 1 is a proportion of the
# samples (at least one sample), otherwise a count.
.needed_samples <- function(min_samples, n, name = "min_samples") {
  min_samples <- .number(min_samples, name, lower = 0)
  needed <- if (min_samples < 1) max(1, ceiling(min_samples * n - 1e-9)) else min_samples
  if (needed > n) .stopf("`%s` (%s) is larger than the number of samples (%d).", name, format(min_samples), n)
  needed
}

# The value that a table made by rsem_read() holds ("count", "tpm", ...), for the `expr` default of the functions
# that follow it in a pipe.
.need_expr <- function(expr) {
  if (is.null(expr)) .stopf("`expr` must be \"count\" or \"tpm\": the table does not record what it holds (it was changed by a transformation or did not come from rsem_read()); say which.")
  expr
}
