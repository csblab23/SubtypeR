#' Read a .gmt gene signature file
#'
#' Parses a standard `.gmt` file (as used by MSigDB / GSEA) into a named list
#' of character vectors, where each element is one signature's gene set.
#' Each line of the file must be tab-separated in the form:
#' `signature_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...`
#'
#' @param file Path to a `.gmt` file.
#'
#' @return A named list of character vectors (gene symbols per signature).
#' @export
#'
#' @examples
#' \dontrun{
#' sigs <- read_gmt("my_signatures.gmt")
#' }
read_gmt <- function(file) {
  if (!file.exists(file)) {
    stop(sprintf("GMT file not found: %s", file), call. = FALSE)
  }

  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) == 0) {
    stop("GMT file is empty.", call. = FALSE)
  }

  split_lines <- strsplit(lines, "\t")

  bad <- which(vapply(split_lines, length, integer(1)) < 3)
  if (length(bad) > 0) {
    stop(sprintf(
      "Malformed GMT line(s) %s: each line needs name, description, and >=1 gene, tab-separated.",
      paste(bad, collapse = ", ")
    ), call. = FALSE)
  }

  sig_names <- vapply(split_lines, `[`, character(1), 1)
  if (any(duplicated(sig_names))) {
    dups <- unique(sig_names[duplicated(sig_names)])
    stop(sprintf(
      "Duplicate signature name(s) in GMT file: %s",
      paste(dups, collapse = ", ")
    ), call. = FALSE)
  }

  genes <- lapply(split_lines, function(x) {
    g <- x[-c(1, 2)]
    g <- unique(g[nzchar(trimws(g))])
    g
  })
  names(genes) <- sig_names

  empty <- names(genes)[vapply(genes, length, integer(1)) == 0]
  if (length(empty) > 0) {
    stop(sprintf(
      "Signature(s) with no genes after parsing: %s",
      paste(empty, collapse = ", ")
    ), call. = FALSE)
  }

  genes
}

#' List built-in signature panels available in "auto" mode
#'
#' @return A character vector of built-in cancer-type panel names that can be
#'   passed to `RunSubtypeR(mode = "auto", cancer_type = ...)`.
#' @export
#'
#' @examples
#' list_builtin_panels()
list_builtin_panels <- function() {
  c("HGSOC")
}

#' Load a built-in signature panel
#'
#' @param cancer_type Name of a built-in panel. See [list_builtin_panels()].
#'
#' @return A named list of character vectors (gene symbols per signature).
#' @export
#'
#' @examples
#' sigs <- load_builtin_signatures("HGSOC")
#' names(sigs)
load_builtin_signatures <- function(cancer_type = "HGSOC") {
  available <- list_builtin_panels()

  cancer_type_match <- available[match(toupper(cancer_type), toupper(available))]
  if (is.na(cancer_type_match)) {
    stop(sprintf(
      "No built-in panel for cancer_type = '%s'. Available panels: %s. Use mode = 'manual' with your own .gmt file for other cancer types.",
      cancer_type, paste(available, collapse = ", ")
    ), call. = FALSE)
  }

  gmt_path <- system.file(
    "extdata", paste0(cancer_type_match, "_signatures.gmt"),
    package = "SubtypeR"
  )

  if (!nzchar(gmt_path) || !file.exists(gmt_path)) {
    stop(sprintf("Could not locate bundled panel file for '%s'.", cancer_type_match), call. = FALSE)
  }

  read_gmt(gmt_path)
}
