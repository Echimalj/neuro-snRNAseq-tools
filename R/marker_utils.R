#' Marker detection utilities for Seurat objects
#'
#' Helper functions for preparing the RNA assay for marker discovery and running
#' cluster-level marker detection.
#'
#' @keywords internal
NULL

#' Prepare RNA assay for marker detection
#'
#' Normalizes, identifies variable features, scales the RNA assay, and joins
#' layers when needed for Seurat v5 workflows.
#'
#' @param seu Seurat object.
#' @param assay Character. Assay to prepare.
#' @param features Optional vector of features to scale. If NULL, all genes.
#' @param run_join_layers Logical. Whether to run `JoinLayers()`.
#' @param verbose Logical.
#'
#' @return Seurat object with prepared RNA assay.
#' @export
prepare_rna_for_markers <- function(seu,
                                    assay = "RNA",
                                    features = NULL,
                                    run_join_layers = TRUE,
                                    verbose = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  Seurat::DefaultAssay(seu) <- assay

  seu <- Seurat::NormalizeData(
    object = seu,
    assay = assay,
    verbose = verbose
  )

  seu <- Seurat::FindVariableFeatures(
    object = seu,
    assay = assay,
    verbose = verbose
  )

  if (is.null(features)) {
    features <- rownames(seu[[assay]])
  }

  seu <- Seurat::ScaleData(
    object = seu,
    assay = assay,
    features = features,
    verbose = verbose
  )

  if (run_join_layers && "JoinLayers" %in% getNamespaceExports("Seurat")) {
    seu <- Seurat::JoinLayers(seu)
  }

  return(seu)
}

#' Find cluster markers
#'
#' Wrapper around `Seurat::FindAllMarkers()`.
#'
#' @param seu Seurat object.
#' @param assay Character. Assay to use.
#' @param only_pos Logical. Whether to return only positive markers.
#' @param logfc_threshold Numeric. Minimum log fold-change threshold.
#' @param min_pct Numeric. Minimum detection percentage.
#'
#' @return Data frame of marker genes.
#' @export
find_cluster_markers <- function(seu,
                                 assay = "RNA",
                                 only_pos = TRUE,
                                 logfc_threshold = 0.25,
                                 min_pct = 0.1) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  Seurat::DefaultAssay(seu) <- assay

  markers <- Seurat::FindAllMarkers(
    object = seu,
    assay = assay,
    only.pos = only_pos,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct
  )

  return(markers)
}

#' Save marker table
#'
#' @param markers Data frame of marker results.
#' @param file Character. Output file path.
#' @param sep Character. Delimiter.
#'
#' @return Invisibly returns output file path.
#' @export
save_marker_table <- function(markers,
                              file,
                              sep = "\t") {
  write.table(
    markers,
    file = file,
    sep = sep,
    quote = FALSE,
    row.names = FALSE
  )

  invisible(file)
}
