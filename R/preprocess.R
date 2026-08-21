#' Preprocess raw counts into a normalized Seurat object
#'
#' Accepts either a raw count matrix (genes x cells) or an existing Seurat
#' object and runs SCTransform normalization. Critically, `return.only.var.genes`
#' is forced to `FALSE` so that every gene present in the input is retained in
#' the normalized (`data`) slot -- this is required so that UCell can score
#' signature genes that may not be highly variable.
#'
#' @param input A raw counts matrix / dgCMatrix (genes x cells), or a Seurat
#'   object already containing raw counts in its default assay.
#' @param assay Name to use for the raw-count assay if `input` is a matrix.
#'   Ignored if `input` is already a Seurat object. Default `"RNA"`.
#' @param min_cells Passed to [Seurat::CreateSeuratObject()] when `input` is a
#'   matrix. Default `3`.
#' @param min_features Passed to [Seurat::CreateSeuratObject()] when `input`
#'   is a matrix. Default `200`.
#' @param vars_to_regress Optional character vector of metadata columns to
#'   regress out during SCTransform (e.g. `"percent.mt"`). Default `NULL`.
#' @param seed Random seed passed through to SCTransform for reproducibility.
#' @param verbose Print progress messages. Default `TRUE`.
#'
#' @return A Seurat object with an `"SCT"` assay, normalized `data` slot
#'   containing all genes (not just variable features).
#' @export
preprocess_counts <- function(input,
                               assay = "RNA",
                               min_cells = 3,
                               min_features = 200,
                               vars_to_regress = NULL,
                               seed = 1448145,
                               verbose = TRUE) {

  if (methods::is(input, "Seurat")) {
    .msg(verbose, "Input is already a Seurat object; using its default assay counts.")
    seurat_obj <- input
  } else {
    .check_matrix_like(input, "input")
    .msg(verbose, "Creating Seurat object from raw count matrix...")
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = input,
      assay = assay,
      min.cells = min_cells,
      min.features = min_features
    )
  }

  .msg(verbose, "Running SCTransform (return.only.var.genes = FALSE, all genes retained)...")
  seurat_obj <- Seurat::SCTransform(
    seurat_obj,
    assay = Seurat::DefaultAssay(seurat_obj),
    new.assay.name = "SCT",
    return.only.var.genes = FALSE,
    vars.to.regress = vars_to_regress,
    seed.use = seed,
    verbose = verbose
  )

  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  seurat_obj
}
