#' SubtypeR: UCell-Based Single-Cell Molecular Subtype Assignment
#'
#' @description
#' SubtypeR assigns a molecular subtype to every cell in a single-cell
#' RNA-seq dataset using UCell gene-signature scoring. See [RunSubtypeR()]
#' for the main entry point.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

# Silence R CMD check NOTES for non-standard evaluation column names used
# inside dplyr/ggplot2 pipelines throughout the package.
utils::globalVariables(c(
  "Signature", "Score", "med", "Cell", "Subtype", "Freq", ".data"
))
