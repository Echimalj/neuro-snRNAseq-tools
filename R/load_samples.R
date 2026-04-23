#' Load a single 10x sample with optional SoupX correction
#'
#' Utilities for loading 10x Genomics outputs into Seurat objects, either
#' directly from the filtered matrix or after ambient RNA correction with SoupX.
#'
#' Recommended sample sheet columns:
#' - sample_id
#' - sample_path
#' - optional: use_soupx
#'
#' @keywords internal
NULL

#' Load one sample using SoupX
#'
#' This function loads a Cell Ranger output directory with SoupX, estimates
#' ambient RNA contamination, adjusts counts, and creates a Seurat object.
#'
#' @param sample_path Character. Path to a Cell Ranger output directory
#'   compatible with `SoupX::load10X()`.
#' @param sample_id Character. Sample identifier used as the Seurat project name.
#' @param min_cells Integer. Passed to `Seurat::CreateSeuratObject()`.
#' @param min_features Integer. Passed to `Seurat::CreateSeuratObject()`.
#'
#' @return A Seurat object.
#' @export
load_sample_with_soupx <- function(sample_path,
                                   sample_id,
                                   min_cells = 3,
                                   min_features = 200) {
  if (!requireNamespace("SoupX", quietly = TRUE)) {
    stop("Package 'SoupX' is required for load_sample_with_soupx().", call. = FALSE)
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for load_sample_with_soupx().", call. = FALSE)
  }

  if (!dir.exists(sample_path)) {
    stop("Sample path does not exist: ", sample_path, call. = FALSE)
  }

  soup <- SoupX::load10X(sample_path)
  soup <- SoupX::autoEstCont(soup)
  adjusted <- SoupX::adjustCounts(soup)

  seu <- Seurat::CreateSeuratObject(
    counts = adjusted,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_features
  )

  rm(soup, adjusted)
  gc()

  return(seu)
}

#' Load one sample without SoupX
#'
#' This function reads a standard 10x filtered matrix directory using
#' `Seurat::Read10X()` and creates a Seurat object.
#'
#' @param sample_path Character. Path to a directory readable by
#'   `Seurat::Read10X(data.dir = ...)`, usually a `filtered_feature_bc_matrix`
#'   directory.
#' @param sample_id Character. Sample identifier used as the Seurat project name.
#' @param min_cells Integer. Passed to `Seurat::CreateSeuratObject()`.
#' @param min_features Integer. Passed to `Seurat::CreateSeuratObject()`.
#'
#' @return A Seurat object.
#' @export
load_sample_basic <- function(sample_path,
                              sample_id,
                              min_cells = 3,
                              min_features = 200) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for load_sample_basic().", call. = FALSE)
  }

  if (!dir.exists(sample_path)) {
    stop("Sample path does not exist: ", sample_path, call. = FALSE)
  }

  counts <- Seurat::Read10X(data.dir = sample_path)

  seu <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_features
  )

  rm(counts)
  gc()

  return(seu)
}

#' Load one sample with or without SoupX
#'
#' Wrapper that dispatches to either `load_sample_with_soupx()` or
#' `load_sample_basic()`.
#'
#' @param sample_path Character. Path to sample directory.
#' @param sample_id Character. Sample identifier.
#' @param min_cells Integer. Minimum cells threshold for Seurat object creation.
#' @param min_features Integer. Minimum features threshold for Seurat object creation.
#' @param use_soupx Logical. If `TRUE`, use SoupX-based loading.
#'
#' @return A Seurat object.
#' @export
load_sample <- function(sample_path,
                        sample_id,
                        min_cells = 3,
                        min_features = 200,
                        use_soupx = TRUE) {
  if (!is.logical(use_soupx) || length(use_soupx) != 1 || is.na(use_soupx)) {
    stop("'use_soupx' must be a single TRUE/FALSE value.", call. = FALSE)
  }

  if (use_soupx) {
    return(load_sample_with_soupx(
      sample_path = sample_path,
      sample_id = sample_id,
      min_cells = min_cells,
      min_features = min_features
    ))
  }

  return(load_sample_basic(
    sample_path = sample_path,
    sample_id = sample_id,
    min_cells = min_cells,
    min_features = min_features
  ))
}

#' Load multiple samples from a sample sheet
#'
#' Iterates over a sample sheet and returns a named list of Seurat objects.
#'
#' Required columns in `sample_sheet`:
#' - `sample_id`
#' - `sample_path`
#'
#' Optional columns:
#' - `use_soupx`
#'
#' If `use_soupx` is supplied as a function argument, it overrides the
#' sample-sheet column for all samples.
#'
#' @param sample_sheet Data frame containing at least `sample_id` and `sample_path`.
#' @param min_cells Integer. Minimum cells threshold for Seurat object creation.
#' @param min_features Integer. Minimum features threshold for Seurat object creation.
#' @param use_soupx Optional logical scalar. If `NULL`, the function will look
#'   for a `use_soupx` column in `sample_sheet`. If not present, defaults to `TRUE`.
#'
#' @return A named list of Seurat objects.
#' @export
load_samples_from_sheet <- function(sample_sheet,
                                    min_cells = 3,
                                    min_features = 200,
                                    use_soupx = NULL) {
  required_cols <- c("sample_id", "sample_path")
  missing_cols <- setdiff(required_cols, colnames(sample_sheet))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required sample_sheet columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(use_soupx)) {
    if (!is.logical(use_soupx) || length(use_soupx) != 1 || is.na(use_soupx)) {
      stop("'use_soupx' must be NULL or a single TRUE/FALSE value.", call. = FALSE)
    }
  }

  seurat_list <- lapply(seq_len(nrow(sample_sheet)), function(i) {
    sample_use_soupx <- if (!is.null(use_soupx)) {
      use_soupx
    } else if ("use_soupx" %in% colnames(sample_sheet)) {
      as.logical(sample_sheet$use_soupx[i])
    } else {
      TRUE
    }

    load_sample(
      sample_path = sample_sheet$sample_path[i],
      sample_id = sample_sheet$sample_id[i],
      min_cells = min_cells,
      min_features = min_features,
      use_soupx = sample_use_soupx
    )
  })

  names(seurat_list) <- sample_sheet$sample_id
  return(seurat_list)
}
