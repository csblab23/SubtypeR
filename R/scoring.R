#' Score gene signatures per cell with UCell
#'
#' Thin, validating wrapper around [UCell::ScoreSignatures_UCell()]. Signatures
#' are intersected with genes actually present in `norm_matrix`; any signature
#' left with fewer than `min_genes` genes is dropped with a warning rather than
#' silently scored on a near-empty gene set.
#'
#' @param norm_matrix A normalized expression matrix, genes in rows, cells in
#'   columns (e.g. `Seurat::GetAssayData(obj, assay = "SCT", slot = "data")`).
#' @param signatures Named list of character vectors (gene symbols per
#'   signature), as returned by [read_gmt()] or [load_builtin_signatures()].
#' @param min_genes Minimum number of signature genes that must be present in
#'   `norm_matrix` for a signature to be scored. Default `2`.
#' @param maxRank Passed to `UCell::ScoreSignatures_UCell`. Default `1500`.
#' @param ncores Number of cores for UCell. Default `1`.
#' @param w_neg Passed to `UCell::ScoreSignatures_UCell`. Default `1`.
#' @param verbose Print progress / filtering messages. Default `TRUE`.
#'
#' @return A data.frame of UCell scores, cells in rows, one column per scored
#'   signature (column names carry UCell's `"_UCell"` suffix).
#' @export
score_subtypes <- function(norm_matrix,
                            signatures,
                            min_genes = 2,
                            maxRank = 1500,
                            ncores = 1,
                            w_neg = 1,
                            verbose = TRUE) {

  .check_matrix_like(norm_matrix, "norm_matrix")

  if (!is.list(signatures) || is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    stop("`signatures` must be a named list of character vectors (gene symbols).", call. = FALSE)
  }

  present_genes <- rownames(norm_matrix)

  filtered <- lapply(signatures, function(genes) intersect(genes, present_genes))
  n_kept <- vapply(filtered, length, integer(1))

  dropped <- names(filtered)[n_kept < min_genes]
  if (length(dropped) > 0) {
    warning(sprintf(
      "Dropping signature(s) with < %d genes present in the expression matrix: %s",
      min_genes, paste(dropped, collapse = ", ")
    ), call. = FALSE)
  }

  filtered <- filtered[n_kept >= min_genes]

  if (length(filtered) == 0) {
    stop("No signature retained >= min_genes genes present in norm_matrix. Check gene symbol matching (case, aliases, species).", call. = FALSE)
  }

  .msg(verbose, sprintf(
    "Scoring %d signature(s) with UCell: %s",
    length(filtered), paste(names(filtered), collapse = ", ")
  ))

  scores <- UCell::ScoreSignatures_UCell(
    matrix = norm_matrix,
    features = filtered,
    maxRank = maxRank,
    ncores = ncores,
    w_neg = w_neg
  )

  as.data.frame(scores)
}
