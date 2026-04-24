#' Preprocessing utilities for Seurat objects
#'
#' Helper functions for SCTransform normalization, PCA, UMAP, neighbor graph
#' construction, and clustering.
#'
#' @keywords internal
NULL

#' Run SCTransform preprocessing and initial clustering
#'
#' @param seu Seurat object.
#' @param vst_flavor Character. SCTransform flavor. Default is `"v2"`.
#' @param npcs Integer. Number of PCs to compute.
#' @param dims Integer vector. PCs to use for UMAP/neighbors.
#' @param resolution Numeric. Clustering resolution.
#' @param assay Character. Assay to use as input.
#' @param new_assay_name Character. Name of the SCT assay.
#' @param vars_to_regress Optional character vector of metadata variables to regress.
#' @param verbose Logical.
#'
#' @return Processed Seurat object.
#' @export
run_sct_pipeline <- function(seu,
                             vst_flavor = "v2",
                             npcs = 22,
                             dims = 1:20,
                             resolution = 0.5,
                             assay = "RNA",
                             new_assay_name = "SCT",
                             vars_to_regress = NULL,
                             verbose = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("sctransform", quietly = TRUE)) {
    stop("Package 'sctransform' is required.", call. = FALSE)
  }

  Seurat::DefaultAssay(seu) <- assay

  seu <- Seurat::SCTransform(
    object = seu,
    assay = assay,
    new.assay.name = new_assay_name,
    vst.flavor = vst_flavor,
    vars.to.regress = vars_to_regress,
    verbose = verbose
  )

  Seurat::DefaultAssay(seu) <- new_assay_name

  seu <- Seurat::RunPCA(
    object = seu,
    assay = new_assay_name,
    npcs = npcs,
    verbose = verbose
  )

  seu <- Seurat::RunUMAP(
    object = seu,
    reduction = "pca",
    dims = dims,
    verbose = verbose
  )

  seu <- Seurat::FindNeighbors(
    object = seu,
    reduction = "pca",
    dims = dims,
    verbose = verbose
  )

  seu <- Seurat::FindClusters(
    object = seu,
    resolution = resolution,
    verbose = verbose
  )

  return(seu)
}

#' Plot processed Seurat clusters
#'
#' @param seu Seurat object.
#' @param reduction Character. Reduction to plot.
#' @param group_by Optional metadata column.
#' @param label Logical. Whether to label clusters.
#' @param repel Logical. Whether labels should repel.
#'
#' @return ggplot object.
#' @export
plot_clusters <- function(seu,
                          reduction = "umap",
                          group_by = NULL,
                          label = TRUE,
                          repel = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  Seurat::DimPlot(
    object = seu,
    reduction = reduction,
    group.by = group_by,
    label = label,
    repel = repel
  )
}
