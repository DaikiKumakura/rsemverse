#' What was done to a table
#'
#' The functions of this package add a line to the table they return (the attribute `"steps"`),
#' so that a table at the end of a pipe tells how it was made.
#'
#' @param x A table made by [rsem_read()] or one of the functions that take and return a table.
#' @return A character vector with one element per step, oldest first; `character(0)` for a table
#'   that was not made by this package.
#' @examples
#' x <- data.frame(S1 = c(5, 0), S2 = c(7, 1), row.names = c("g1", "g2"))
#' rsem_steps(x)
#' @export
rsem_steps <- function(x) {
  st <- attr(x, "steps", exact = TRUE)
  if (is.null(st)) character(0) else st
}
