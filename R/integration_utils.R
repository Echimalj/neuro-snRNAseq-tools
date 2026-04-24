#' Integration and clustering utilities for Seurat objects
#'
#' Helper functions for Harmony integration, integrated UMAP generation,
#' neighbor graph construction, clustering, and resolution testing.
#'
#' @keywords internal
NULL

#' Run Harmony integration and clustering
#'
#' @param seu Seurat object.
#' @param group_var Character. Metadata column to integrate over.
#' @param assay Character. Assay to use. Default is `"SCT"`.
#' @param reduction Character. Input reduction for Harmony.
#' @param harmony_reduction Character. Output Harmony reduction name.
#' @param dims Integer vector. Dimensions to use.
#' @param resolution Numeric. Clustering resolution.
#' @param verbose Logical.
#'
#' @return Seurat object with Harmony reduction, UMAP, neighbors, and clusters.
#' @export
run_harmony_clustering <- function(seu,
                                   group_var,
                                   assay = "SCT",
                                   reduction = "pca",
                                   harmony_reduction = "harmony",
                                   dims = 1:20,
                                   resolution = 0.5,
                                   verbose = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Package 'harmony' is required.", call. = FALSE)
  }

  if (!group_var %in% colnames(seu@meta.data)) {
    stop("group_var not found in metadata: ", group_var, call. = FALSE)
  }

  Seurat::DefaultAssay(seu) <- assay

  seu <- harmony::RunHarmony(
    object = seu,
    group.by.vars = group_var,
    reduction = reduction,
    assay.use = assay,
    reduction.save = harmony_reduction,
    verbose = verbose
  )

  seu <- Seurat::FindNeighbors(
    object = seu,
    reduction = harmony_reduction,
    dims = dims,
    verbose = verbose
  )

  seu <- Seurat::FindClusters(
    object = seu,
    resolution = resolution,
    verbose = verbose
  )

  seu <- Seurat::RunUMAP(
    object = seu,
    reduction = harmony_reduction,
    dims = dims,
    verbose = verbose
  )

  return(seu)
}

#' Test multiple clustering resolutions
#'
#' @param seu Seurat object.
#' @param reduction Character. Reduction to use.
#' @param dims Integer vector. Dimensions to use.
#' @param resolutions Numeric vector of clustering resolutions.
#' @param verbose Logical.
#'
#' @return Seurat object with multiple clustering columns added to metadata.
#' @export
test_clustering_resolutions <- function(seu,
                                        reduction = "harmony",
                                        dims = 1:20,
                                        resolutions = c(0.1, 0.3, 0.5, 0.7, 0.9),
                                        verbose = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  seu <- Seurat::FindNeighbors(
    object = seu,
    reduction = reduction,
    dims = dims,
    verbose = verbose
  )

  seu <- Seurat::FindClusters(
    object = seu,
    resolution = resolutions,
    verbose = verbose
  )

  return(seu)
}

#' Plot clustree for tested resolutions
#'
#' @param seu Seurat object.
#' @param prefix Character. Metadata prefix for clustering columns.
#'
#' @return clustree plot.
#' @export
plot_clustree <- function(seu,
                          prefix = "SCT_snn_res.") {
  if (!requireNamespace("clustree", quietly = TRUE)) {
    stop("Package 'clustree' is required.", call. = FALSE)
  }

  clustree::clustree(seu, prefix = prefix)
}
