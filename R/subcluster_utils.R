#' Subclustering utilities for Seurat objects
#'
#' Helper functions for targeted subclustering and correlation-based merging.
#'
#' @keywords internal
NULL

#' Subcluster selected identities
#'
#' @param seu Seurat object.
#' @param clusters_to_subcluster Character vector of active identity labels to subcluster.
#' @param output_col Metadata column to store subcluster labels.
#' @param original_col Metadata column to store original identities.
#' @param resolution Clustering resolution for subclustering.
#' @param dims Dimensions to use.
#' @param npcs Number of PCs to compute.
#' @param assay Input assay.
#' @param vst_flavor SCTransform flavor.
#'
#' @return Seurat object with subcluster labels added.
#' @export
subcluster_selected_idents <- function(seu,
                                       clusters_to_subcluster,
                                       output_col = "subcluster",
                                       original_col = "original_cluster",
                                       resolution = 0.1,
                                       dims = 1:20,
                                       npcs = 20,
                                       assay = "RNA",
                                       vst_flavor = "v2") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("sctransform", quietly = TRUE)) {
    stop("Package 'sctransform' is required.", call. = FALSE)
  }

  seu[[original_col]] <- as.character(Seurat::Idents(seu))
  seu[[output_col]] <- as.character(Seurat::Idents(seu))

  missing_clusters <- setdiff(clusters_to_subcluster, levels(Seurat::Idents(seu)))
  if (length(missing_clusters) > 0) {
    stop(
      "These clusters are not present in active identities: ",
      paste(missing_clusters, collapse = ", "),
      call. = FALSE
    )
  }

  for (cl in clusters_to_subcluster) {
    message("Subclustering: ", cl)

    data_sub <- subset(seu, idents = cl)

    Seurat::DefaultAssay(data_sub) <- assay

    data_sub <- Seurat::SCTransform(
      data_sub,
      vst.flavor = vst_flavor,
      verbose = FALSE
    )

    data_sub <- Seurat::RunPCA(
      data_sub,
      npcs = npcs,
      verbose = FALSE
    )

    data_sub <- Seurat::RunUMAP(
      data_sub,
      reduction = "pca",
      dims = dims,
      verbose = FALSE
    )

    data_sub <- Seurat::FindNeighbors(
      data_sub,
      reduction = "pca",
      dims = dims,
      verbose = FALSE
    )

    data_sub <- Seurat::FindClusters(
      data_sub,
      resolution = resolution,
      verbose = FALSE
    )

    new_labels <- paste0(cl, "_sub", data_sub$seurat_clusters)

    seu@meta.data[colnames(data_sub), output_col] <- new_labels

    rm(data_sub)
    gc()
  }

  return(seu)
}

#' Set identities from metadata column
#'
#' @param seu Seurat object.
#' @param metadata_col Metadata column to use as active identities.
#'
#' @return Seurat object.
#' @export
set_idents_from_metadata <- function(seu, metadata_col) {
  if (!metadata_col %in% colnames(seu@meta.data)) {
    stop("metadata_col not found: ", metadata_col, call. = FALSE)
  }

  Seurat::Idents(seu) <- seu[[metadata_col]][, 1]
  return(seu)
}

#' Compute average-expression Pearson correlation between identities
#'
#' @param seu Seurat object.
#' @param assay Assay to use.
#' @param slot Slot/layer to use.
#' @param method Correlation method.
#'
#' @return Pearson correlation matrix.
#' @export
compute_identity_correlation <- function(seu,
                                         assay = "RNA",
                                         slot = "data",
                                         method = "pearson") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  avg_exp <- Seurat::AggregateExpression(
    seu,
    assays = assay,
    return.seurat = TRUE
  )

  mat <- Seurat::GetAssayData(
    avg_exp,
    assay = assay,
    slot = slot
  )

  cor_mat <- stats::cor(mat, method = method)

  return(cor_mat)
}

#' Save identity correlation matrix
#'
#' @param cor_mat Correlation matrix.
#' @param file Output file.
#'
#' @return Invisibly returns file path.
#' @export
save_correlation_matrix <- function(cor_mat, file) {
  write.table(
    cor_mat,
    file = file,
    sep = "\t",
    quote = FALSE
  )

  invisible(file)
}

#' Plot identity correlation heatmap
#'
#' @param cor_mat Correlation matrix.
#' @param main Plot title.
#'
#' @return pheatmap object.
#' @export
plot_identity_correlation <- function(cor_mat,
                                      main = "Pearson Correlation Between Subclusters") {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.", call. = FALSE)
  }

  pheatmap::pheatmap(
    cor_mat,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "average",
    main = main
  )
}

#' Merge identity labels using a named mapping
#'
#' @param seu Seurat object.
#' @param merge_map Named character vector. Names are old labels, values are new labels.
#' @param metadata_col Optional metadata column to store merged labels.
#' @param set_idents Logical. Whether to set merged labels as active identities.
#'
#' @return Seurat object with merged labels.
#' @export
merge_identity_labels <- function(seu,
                                  merge_map,
                                  metadata_col = "merged_celltype",
                                  set_idents = TRUE) {
  current_labels <- as.character(Seurat::Idents(seu))

  merged_labels <- current_labels
  idx <- current_labels %in% names(merge_map)
  merged_labels[idx] <- unname(merge_map[current_labels[idx]])

  seu[[metadata_col]] <- merged_labels

  if (set_idents) {
    Seurat::Idents(seu) <- seu[[metadata_col]][, 1]
  }

  return(seu)
}
