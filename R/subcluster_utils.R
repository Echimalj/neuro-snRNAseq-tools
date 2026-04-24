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

#' Save a Seurat subset by identity
#'
#' @param seu Seurat object.
#' @param idents Character vector of identities to subset.
#' @param file Output RDS file.
#'
#' @return Subsetted Seurat object.
#' @export
save_subset_by_idents <- function(seu,
                                  idents = NULL,
                                  file = NULL,
                                  group_col = NULL,
                                  groups = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  # If idents is NULL → keep all identities
  seu_sub <- if (is.null(idents)) {
    seu
  } else {
    subset(seu, idents = idents)
  }

  # Optional filtering by group
  if (!is.null(group_col) && !is.null(groups)) {
    if (!group_col %in% colnames(seu_sub@meta.data)) {
      stop("group_col not found in metadata: ", group_col, call. = FALSE)
    }

    keep_cells <- rownames(seu_sub@meta.data)[
      seu_sub@meta.data[[group_col]] %in% groups
    ]

    seu_sub <- subset(seu_sub, cells = keep_cells)
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seu_sub, file = file)
    message("Saved subset: ", file)
  }

  return(seu_sub)
}

#' Run standalone subcluster workflow
#'
#' @param seu_sub Subsetted Seurat object.
#' @param assay Input assay.
#' @param resolution Clustering resolution.
#' @param dims Dimensions to use.
#' @param npcs Number of PCs.
#' @param vst_flavor SCTransform flavor.
#'
#' @return Processed Seurat subset.
#' @export
run_standalone_subcluster_workflow <- function(seu_sub,
                                               assay = "RNA",
                                               resolution = 0.6,
                                               dims = 1:20,
                                               npcs = 20,
                                               vst_flavor = "v2") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("sctransform", quietly = TRUE)) {
    stop("Package 'sctransform' is required.", call. = FALSE)
  }

  Seurat::DefaultAssay(seu_sub) <- assay

  seu_sub <- Seurat::SCTransform(
    seu_sub,
    vst.flavor = vst_flavor,
    verbose = FALSE
  )

  seu_sub <- Seurat::RunPCA(
    seu_sub,
    npcs = npcs,
    verbose = FALSE
  )

  seu_sub <- Seurat::RunUMAP(
    seu_sub,
    reduction = "pca",
    dims = dims,
    verbose = FALSE
  )

  seu_sub <- Seurat::FindNeighbors(
    seu_sub,
    reduction = "pca",
    dims = dims,
    verbose = FALSE
  )

  seu_sub <- Seurat::FindClusters(
    seu_sub,
    resolution = resolution,
    verbose = FALSE
  )

  return(seu_sub)
}

#' Save UMAP plots for a Seurat object
#'
#' @param seu Seurat object.
#' @param file Output EPS/PDF/PNG file.
#' @param split_by Optional metadata column to split by.
#' @param width Plot width.
#' @param height Plot height.
#' @param label Logical. Whether to label clusters.
#' @param pt_size Point size.
#'
#' @return Invisibly returns file path.
#' @export
save_umap_plot <- function(seu,
                           file,
                           split_by = NULL,
                           width = 6,
                           height = 6,
                           label = FALSE,
                           pt_size = 0.5,
                           reduction = "umap") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  grDevices::cairo_ps(filename = file, width = width, height = height)

  print(
    Seurat::DimPlot(
      seu,
      reduction = reduction,
      label = label,
      split.by = split_by,
      pt.size = pt_size
    )
  )

  grDevices::dev.off()

  invisible(file)
}

#' Prepare subset RNA assay and find markers
#'
#' @param seu_sub Seurat object.
#' @param assay Assay to use.
#' @param only_pos Logical. Passed to FindAllMarkers.
#' @param output_file Optional output file.
#'
#' @return Marker table.
#' @export
find_subset_cluster_markers <- function(seu_sub,
                                        assay = "RNA",
                                        only_pos = TRUE,
                                        output_file = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  Seurat::DefaultAssay(seu_sub) <- assay

  seu_sub <- Seurat::NormalizeData(seu_sub, verbose = FALSE)
  seu_sub <- Seurat::FindVariableFeatures(seu_sub, verbose = FALSE)
  seu_sub <- Seurat::ScaleData(seu_sub, features = rownames(seu_sub), verbose = FALSE)

  if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
    seu_sub <- Seurat::JoinLayers(seu_sub)
  }

  markers <- Seurat::FindAllMarkers(seu_sub, only.pos = only_pos)

  if (!is.null(output_file)) {
    write.table(
      markers,
      file = output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  return(markers)
}

#' Run propeller on a subclustered subset
#'
#' @param seu_sub Seurat object.
#' @param sample_col Sample metadata column.
#' @param group_col Group metadata column.
#' @param output_file Optional CSV output file.
#'
#' @return Propeller result.
#' @export
run_subset_propeller <- function(seu_sub,
                                 sample_col = "orig.ident",
                                 group_col = "Genotype",
                                 output_file = NULL) {
  if (!requireNamespace("speckle", quietly = TRUE)) {
    stop("Package 'speckle' is required.", call. = FALSE)
  }

  prop <- speckle::propeller(
    clusters = Seurat::Idents(seu_sub),
    sample = seu_sub[[sample_col]][, 1],
    group = seu_sub[[group_col]][, 1]
  )

  if (!is.null(output_file)) {
    write.csv(prop, file = output_file)
  }

  return(prop)
}
