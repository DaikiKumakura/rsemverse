# Shared helpers for the functions that work on a gene by sample table.

# Turn the table into a numeric matrix with genes as row names. Accepts the
# data frame rsem_read() returns (symbols as row names), a data frame whose
# first column holds the gene names, or a numeric matrix with row names.
.expr_matrix <- function(x, nonnegative = TRUE, arg = "x") {
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
  if (anyNA(m)) .stopf("`%s` contains missing values.", arg)
  if (any(!is.finite(m))) .stopf("`%s` contains infinite values.", arg)
  if (nonnegative && any(m < 0)) {
    .stopf("`%s` contains negative values; it must be counts or TPM, not log-transformed or scaled values.", arg)
  }
  m
}

# The same layout as rsem_read(): a data frame with genes as row names.
.as_table <- function(m) {
  as.data.frame(m, check.names = FALSE, stringsAsFactors = FALSE)
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
