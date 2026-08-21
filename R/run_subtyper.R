#' Assign single-cell molecular subtypes with UCell
#'
#' `RunSubtypeR()` is the main entry point of the package, analogous to how
#' `cytotrace2()` assigns a stemness score to every cell. Given raw counts (or
#' an existing Seurat object), it runs SCTransform (retaining all genes),
#' scores a panel of gene signatures per cell with UCell, and assigns each
#' cell a single subtype label based on a per-signature quantile threshold.
#'
#' Two modes are supported:
#' \itemize{
#'   \item `mode = "auto"`: uses a built-in signature panel bundled with the
#'     package, selected via `cancer_type` (currently `"HGSOC"`; see
#'     [list_builtin_panels()]).
#'   \item `mode = "manual"`: uses a user-supplied `.gmt` file (`gmt_file`),
#'     applicable to any cancer type or custom gene-signature set.
#' }
#'
#' @param input Raw counts matrix / dgCMatrix (genes x cells), or a Seurat
#'   object with raw counts in its default assay.
#' @param mode Either `"auto"` (use a built-in panel) or `"manual"` (use a
#'   user-supplied `.gmt` file). Default `"auto"`.
#' @param cancer_type Built-in panel name, used only when `mode = "auto"`.
#'   Default `"HGSOC"`. See [list_builtin_panels()].
#' @param gmt_file Path to a `.gmt` file, required when `mode = "manual"`.
#' @param assay Assay name to use if `input` is a raw matrix. Default `"RNA"`.
#' @param min_cells,min_features Passed to [Seurat::CreateSeuratObject()] when
#'   `input` is a raw matrix. Defaults `3` and `200`.
#' @param vars_to_regress Optional metadata columns to regress out in
#'   SCTransform. Default `NULL`.
#' @param min_genes_per_signature Minimum number of a signature's genes that
#'   must be present in the data for that signature to be scored. Default `2`.
#' @param maxRank,w_neg,ncores Passed through to
#'   `UCell::ScoreSignatures_UCell()`. Defaults `1500`, `1`, `1`.
#' @param assign_quantile Quantile probability used as the per-signature
#'   pass/fail threshold for subtype calling. Default `0.75` (Q3).
#' @param seed Random seed for reproducibility. Default `1448145`.
#' @param add_to_seurat If `TRUE` (default) and preprocessing produced a
#'   Seurat object, per-signature scores and the `Subtype` call are added to
#'   its metadata.
#' @param verbose Print progress messages. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{seurat_object}{The (SCTransform-normalized) Seurat object, with
#'     `Subtype` and per-signature UCell scores added to `meta.data` if
#'     `add_to_seurat = TRUE`.}
#'   \item{results}{A data.frame, one row per cell, with per-signature UCell
#'     scores, `n_signatures_passed`, and `Subtype`.}
#'   \item{signatures_used}{The named list of gene signatures actually scored
#'     (after intersecting with genes present in the data).}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Automatic mode, built-in HGSOC panel
#' res <- RunSubtypeR(raw_counts, mode = "auto", cancer_type = "HGSOC")
#'
#' # Manual mode, any cancer type via a user .gmt file
#' res <- RunSubtypeR(raw_counts, mode = "manual", gmt_file = "my_signatures.gmt")
#'
#' table(res$results$Subtype)
#' plot_subtype_scores(res$results)
#' }
RunSubtypeR <- function(input,
                         mode = c("auto", "manual"),
                         cancer_type = "HGSOC",
                         gmt_file = NULL,
                         assay = "RNA",
                         min_cells = 3,
                         min_features = 200,
                         vars_to_regress = NULL,
                         min_genes_per_signature = 2,
                         maxRank = 1500,
                         w_neg = 1,
                         ncores = 1,
                         assign_quantile = 0.75,
                         seed = 1448145,
                         add_to_seurat = TRUE,
                         verbose = TRUE) {

  mode <- match.arg(mode)
  set.seed(seed)

  # ---- 1. Resolve gene signatures -----------------------------------------
  if (mode == "auto") {
    .msg(verbose, sprintf("[Mode: auto] Loading built-in '%s' signature panel...", cancer_type))
    signatures <- load_builtin_signatures(cancer_type)
  } else {
    if (is.null(gmt_file)) {
      stop("`gmt_file` must be supplied when mode = 'manual'.", call. = FALSE)
    }
    .msg(verbose, sprintf("[Mode: manual] Reading signatures from %s ...", gmt_file))
    signatures <- read_gmt(gmt_file)
  }
  .msg(verbose, sprintf("Loaded %d signature(s): %s", length(signatures), paste(names(signatures), collapse = ", ")))

  # ---- 2. Preprocess: raw counts -> SCTransform-normalized Seurat object --
  seurat_obj <- preprocess_counts(
    input = input,
    assay = assay,
    min_cells = min_cells,
    min_features = min_features,
    vars_to_regress = vars_to_regress,
    seed = seed,
    verbose = verbose
  )

  norm_matrix <- Seurat::GetAssayData(seurat_obj, assay = "SCT", slot = "data")

  # ---- 3. Score signatures with UCell -------------------------------------
  scores <- score_subtypes(
    norm_matrix = norm_matrix,
    signatures = signatures,
    min_genes = min_genes_per_signature,
    maxRank = maxRank,
    ncores = ncores,
    w_neg = w_neg,
    verbose = verbose
  )

  # ---- 4. Assign per-cell subtype -----------------------------------------
  .msg(verbose, sprintf("Assigning subtypes at quantile threshold = %.2f ...", assign_quantile))
  results <- assign_subtypes(scores, quantile_prob = assign_quantile)

  .msg(verbose, "Subtype composition:")
  if (isTRUE(verbose)) print(table(results$Subtype))

  # ---- 5. Attach to Seurat object -----------------------------------------
  if (isTRUE(add_to_seurat)) {
    meta <- results[, !(colnames(results) %in% "Cell"), drop = FALSE]
    rownames(meta) <- results$Cell
    meta <- meta[colnames(seurat_obj), , drop = FALSE]
    seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = meta)
  }

  list(
    seurat_object = seurat_obj,
    results = results,
    signatures_used = signatures
  )
}
