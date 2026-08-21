#' @keywords internal
#' @noRd
.strip_ucell_suffix <- function(x) {
  sub("_UCell$", "", x)
}

#' @keywords internal
#' @noRd
.msg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
}

#' @keywords internal
#' @noRd
.check_matrix_like <- function(x, arg_name = "input") {
  ok <- is.matrix(x) || methods::is(x, "Matrix") || methods::is(x, "dgCMatrix")
  if (!ok) {
    stop(sprintf(
      "`%s` must be a matrix, dgCMatrix, or Seurat object (genes x cells raw counts).",
      arg_name
    ), call. = FALSE)
  }
  invisible(TRUE)
}
