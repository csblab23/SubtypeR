#' Assign a per-cell subtype label from UCell scores
#'
#' For each signature, a per-signature threshold is computed as the
#' `quantile_prob`-th quantile of that signature's score across all cells
#' (default: 75th percentile, i.e. Q3). A cell "passes" a signature if its
#' score for that signature exceeds the threshold. Cells are then labeled:
#' \itemize{
#'   \item a single signature name, if the cell passes exactly one signature;
#'   \item `"Mixed"`, if the cell passes more than one signature;
#'   \item `"Unknown"`, if the cell passes none.
#' }
#'
#' @param scores A data.frame or matrix of UCell scores, cells in rows,
#'   signatures in columns (as returned by [score_subtypes()]).
#' @param quantile_prob Quantile probability used as the pass/fail threshold
#'   per signature. Default `0.75` (Q3, matching the original workflow).
#'   Must be in `(0, 1)`.
#'
#' @return A data.frame with one row per cell, containing: the original score
#'   columns, `n_signatures_passed`, and `Subtype` (signature names with any
#'   `"_UCell"` suffix stripped, or `"Mixed"` / `"Unknown"`).
#' @export
assign_subtypes <- function(scores, quantile_prob = 0.75) {

  if (!(is.data.frame(scores) || is.matrix(scores))) {
    stop("`scores` must be a data.frame or matrix (cells x signatures).", call. = FALSE)
  }
  if (!is.numeric(quantile_prob) || length(quantile_prob) != 1 ||
      quantile_prob <= 0 || quantile_prob >= 1) {
    stop("`quantile_prob` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  scores_df <- as.data.frame(scores)
  score_cols <- colnames(scores_df)

  if (is.null(rownames(scores_df))) {
    rownames(scores_df) <- paste0("cell_", seq_len(nrow(scores_df)))
  }

  thresholds <- apply(scores_df, 2, stats::quantile, probs = quantile_prob, na.rm = TRUE)

  pass_matrix <- sweep(as.matrix(scores_df), 2, thresholds, FUN = ">")
  n_pass <- rowSums(pass_matrix)

  subtype <- rep("Unknown", nrow(scores_df))

  one_pass_idx <- which(n_pass == 1)
  if (length(one_pass_idx) > 0) {
    subtype[one_pass_idx] <- score_cols[max.col(pass_matrix[one_pass_idx, , drop = FALSE], ties.method = "first")]
  }

  subtype[n_pass > 1] <- "Mixed"

  subtype <- ifelse(subtype %in% c("Mixed", "Unknown"), subtype, .strip_ucell_suffix(subtype))

  out <- tibble::tibble(
    Cell = rownames(scores_df),
    scores_df,
    n_signatures_passed = n_pass,
    Subtype = subtype
  )

  attr(out, "thresholds") <- thresholds
  attr(out, "quantile_prob") <- quantile_prob

  as.data.frame(out)
}
