#' Quality control utilities for Seurat objects
#'
#' Helper functions for adding QC metrics, visualizing QC features, filtering
#' cells, and summarizing cell retention after filtering.
#'
#' Common mitochondrial gene patterns:
#' - Human: "^MT-"
#' - Mouse: "^mt-"
#'
#' @keywords internal
NULL

#' Add mitochondrial percentage to a Seurat object
#'
#' @param seu Seurat object.
#' @param mt_pattern Character. Regex pattern for mitochondrial genes.
#'   Use `"^MT-"` for human and `"^mt-"` for mouse.
#' @param col_name Character. Metadata column name for mitochondrial percentage.
#'
#' @return Seurat object with mitochondrial percentage metadata.
#' @export
add_mito_percent <- function(seu,
                             mt_pattern = "^MT-",
                             col_name = "percent.mt") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  seu[[col_name]] <- Seurat::PercentageFeatureSet(
    object = seu,
    pattern = mt_pattern
  )

  return(seu)
}

#' Add ribosomal percentage to a Seurat object
#'
#' @param seu Seurat object.
#' @param ribo_pattern Character. Regex pattern for ribosomal genes.
#' @param col_name Character. Metadata column name.
#'
#' @return Seurat object with ribosomal percentage metadata.
#' @export
add_ribo_percent <- function(seu,
                             ribo_pattern = "^RP[SL]",
                             col_name = "percent.ribo") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  seu[[col_name]] <- Seurat::PercentageFeatureSet(
    object = seu,
    pattern = ribo_pattern
  )

  return(seu)
}

#' Add hemoglobin percentage to a Seurat object
#'
#' @param seu Seurat object.
#' @param hb_pattern Character. Regex pattern for hemoglobin genes.
#' @param col_name Character. Metadata column name.
#'
#' @return Seurat object with hemoglobin percentage metadata.
#' @export
add_hb_percent <- function(seu,
                           hb_pattern = "^HB[^(P)]",
                           col_name = "percent.hb") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  seu[[col_name]] <- Seurat::PercentageFeatureSet(
    object = seu,
    pattern = hb_pattern
  )

  return(seu)
}

#' Add standard QC metrics
#'
#' Adds mitochondrial percentage and optionally ribosomal/hemoglobin percentage.
#'
#' @param seu Seurat object.
#' @param species Character. Either `"human"` or `"mouse"`.
#' @param add_ribo Logical. Whether to add ribosomal percentage.
#' @param add_hb Logical. Whether to add hemoglobin percentage.
#'
#' @return Seurat object with QC metadata.
#' @export
add_qc_metrics <- function(seu,
                           species = c("human", "mouse"),
                           add_ribo = FALSE,
                           add_hb = FALSE) {
  species <- match.arg(species)

  mt_pattern <- if (species == "human") "^MT-" else "^mt-"
  ribo_pattern <- if (species == "human") "^RP[SL]" else "^Rp[sl]"
  hb_pattern <- if (species == "human") "^HB[^(P)]" else "^Hb[^(p)]"

  seu <- add_mito_percent(seu, mt_pattern = mt_pattern)

  if (add_ribo) {
    seu <- add_ribo_percent(seu, ribo_pattern = ribo_pattern)
  }

  if (add_hb) {
    seu <- add_hb_percent(seu, hb_pattern = hb_pattern)
  }

  return(seu)
}

#' Plot basic QC violin plots
#'
#' @param seu Seurat object.
#' @param features Character vector of QC features to plot.
#' @param group_by Optional metadata column for grouping.
#' @param split_by Optional metadata column for splitting.
#' @param ncol Integer. Number of columns.
#' @param pt_size Numeric. Point size.
#'
#' @return A ggplot/patchwork object.
#' @export
plot_basic_qc <- function(seu,
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                          group_by = NULL,
                          split_by = NULL,
                          ncol = 3,
                          pt_size = 0.1) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  missing_features <- setdiff(features, colnames(seu@meta.data))

  if (length(missing_features) > 0) {
    stop(
      "These QC features are missing from metadata: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }

  Seurat::VlnPlot(
    object = seu,
    features = features,
    group.by = group_by,
    split.by = split_by,
    ncol = ncol,
    pt.size = pt_size
  )
}

#' Filter cells using basic QC thresholds
#'
#' @param seu Seurat object.
#' @param min_features Numeric. Minimum detected features.
#' @param max_features Optional numeric. Maximum detected features.
#' @param min_counts Optional numeric. Minimum RNA counts.
#' @param max_counts Optional numeric. Maximum RNA counts.
#' @param max_mt Numeric. Maximum mitochondrial percentage.
#' @param mt_col Character. Metadata column containing mitochondrial percentage.
#'
#' @return Filtered Seurat object.
#' @export
filter_cells_basic <- function(seu,
                               min_features = 200,
                               max_features = NULL,
                               min_counts = NULL,
                               max_counts = NULL,
                               max_mt = 1,
                               mt_col = "percent.mt") {
  required_cols <- c("nFeature_RNA", "nCount_RNA", mt_col)
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata columns for filtering: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  keep <- seu$nFeature_RNA > min_features &
    seu[[mt_col]][, 1] < max_mt

  if (!is.null(max_features)) {
    keep <- keep & seu$nFeature_RNA < max_features
  }

  if (!is.null(min_counts)) {
    keep <- keep & seu$nCount_RNA > min_counts
  }

  if (!is.null(max_counts)) {
    keep <- keep & seu$nCount_RNA < max_counts
  }

  seu_filtered <- subset(seu, cells = colnames(seu)[keep])

  return(seu_filtered)
}

#' Count cells before and after filtering
#'
#' @param before Seurat object before filtering.
#' @param after Seurat object after filtering.
#' @param group_col Character. Metadata column to summarize.
#'
#' @return Data frame with before/after cell counts and retained fraction.
#' @export
summarize_filtering <- function(before,
                                after,
                                group_col = "orig.ident") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  if (!group_col %in% colnames(before@meta.data)) {
    stop("group_col not found in before object metadata: ", group_col, call. = FALSE)
  }

  if (!group_col %in% colnames(after@meta.data)) {
    stop("group_col not found in after object metadata: ", group_col, call. = FALSE)
  }

  before_df <- as.data.frame(table(before@meta.data[[group_col]]))
  colnames(before_df) <- c(group_col, "n_before")

  after_df <- as.data.frame(table(after@meta.data[[group_col]]))
  colnames(after_df) <- c(group_col, "n_after")

  out <- dplyr::left_join(before_df, after_df, by = group_col)
  out$n_after[is.na(out$n_after)] <- 0
  out$fraction_retained <- out$n_after / out$n_before

  return(out)
}
