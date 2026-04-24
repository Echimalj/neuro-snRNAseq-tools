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
