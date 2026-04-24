#' Annotation utilities for Seurat clusters
#'
#' Helper functions for applying manually curated cluster labels.
#'
#' @keywords internal
NULL

#' Apply manual cluster labels
#'
#' @param seu Seurat object.
#' @param cluster_labels Named character vector. Names should match current identity levels.
#' @param new_metadata_col Optional metadata column to store renamed labels.
#' @param set_idents Logical. Whether to set renamed labels as active identities.
#'
#' @return Seurat object with updated labels.
#' @export
apply_cluster_labels <- function(seu,
                                 cluster_labels,
                                 new_metadata_col = "celltype",
                                 set_idents = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  current_levels <- levels(Seurat::Idents(seu))
  missing_labels <- setdiff(current_levels, names(cluster_labels))

  if (length(missing_labels) > 0) {
    stop(
      "Missing labels for these identity levels: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }

  renamed <- plyr::mapvalues(
    x = as.character(Seurat::Idents(seu)),
    from = names(cluster_labels),
    to = unname(cluster_labels)
  )

  seu[[new_metadata_col]] <- renamed

  if (set_idents) {
    Seurat::Idents(seu) <- seu[[new_metadata_col]][, 1]
  }

  return(seu)
}

#' Create named cluster label vector
#'
#' @param labels Character vector of new labels.
#' @param cluster_levels Character vector of cluster levels.
#'
#' @return Named character vector.
#' @export
make_cluster_label_vector <- function(labels, cluster_levels) {
  if (length(labels) != length(cluster_levels)) {
    stop("labels and cluster_levels must have the same length.", call. = FALSE)
  }

  names(labels) <- cluster_levels
  labels
}

#' Collapse detailed cell class labels into broad cell classes
#'
#' @param x Character vector of cell labels.
#' @param keep Character vector of labels to keep unchanged.
#'
#' @return Character vector of collapsed cell class labels.
#' @export
collapse_cellclass_labels <- function(x,
                                      keep = character()) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required.", call. = FALSE)
  }

  dplyr::case_when(
    x %in% keep ~ x,
    stringr::str_detect(x, "^ExNeuron") ~ "ExNeuron",
    stringr::str_detect(x, "^InhNeuron") ~ "InhNeuron",
    stringr::str_detect(x, "^Oligodendrocytes") ~ "Oligodendrocytes",
    stringr::str_detect(x, "^OPC") ~ "OPC",
    stringr::str_detect(x, "^Astrocytes") ~ "Astrocytes",
    stringr::str_detect(x, "^Microglia") ~ "Microglia",
    stringr::str_detect(x, "^VLMC") ~ "VLMC",
    x %in% c("SMC", "Pericytes", "Fibroblast", "Endothelial") ~ x,
    TRUE ~ "Other"
  )
}

#' Add collapsed cell class metadata to a Seurat object
#'
#' @param seu Seurat object.
#' @param input_col Metadata column containing detailed cell labels.
#' @param output_col Metadata column to store collapsed labels.
#' @param keep Character vector of labels to keep unchanged.
#' @param set_idents Logical. Whether to set collapsed labels as active identities.
#'
#' @return Seurat object with collapsed cell class metadata.
#' @export
add_collapsed_cellclass <- function(seu,
                                    input_col = "cellclass",
                                    output_col = "collapsed_cellclass",
                                    keep = character(),
                                    set_idents = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  if (!input_col %in% colnames(seu@meta.data)) {
    stop("input_col not found in metadata: ", input_col, call. = FALSE)
  }

  new_labels <- collapse_cellclass_labels(
    x = as.character(seu[[input_col]][, 1]),
    keep = keep
  )

  names(new_labels) <- colnames(seu)

  seu <- Seurat::AddMetaData(
    object = seu,
    metadata = new_labels,
    col.name = output_col
  )

  if (set_idents) {
    Seurat::Idents(seu) <- seu[[output_col]][, 1]
  }

  return(seu)
}
