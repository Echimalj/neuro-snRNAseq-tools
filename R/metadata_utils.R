#' Metadata utilities for Seurat sample lists
#'
#' Helper functions for annotating individual Seurat objects from a sample sheet
#' and merging them into a single object.
#'
#' Recommended sample sheet columns:
#' - sample_id
#' - condition
#' - batch
#' - species
#'
#' Optional columns:
#' - donor_id
#' - sex
#' - age
#' - genotype
#' - diagnosis
#'
#' @keywords internal
NULL

#' Check required sample sheet columns
#'
#' @param sample_sheet Data frame containing sample metadata.
#' @param required_cols Character vector of required column names.
#'
#' @return Invisibly returns TRUE if all required columns are present.
#' @export
check_sample_sheet_columns <- function(sample_sheet,
                                       required_cols = c("sample_id")) {
  missing_cols <- setdiff(required_cols, colnames(sample_sheet))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required sample_sheet columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Annotate a list of Seurat objects from a sample sheet
#'
#' Adds sample-level metadata from a sample sheet to each Seurat object in a
#' named list. The list names should match the `sample_id` column.
#'
#' @param seurat_list Named list of Seurat objects.
#' @param sample_sheet Data frame containing sample metadata.
#' @param sample_id_col Character. Column in `sample_sheet` containing sample IDs.
#' @param metadata_cols Optional character vector of columns to add. If NULL,
#'   all columns except `sample_path` are added.
#'
#' @return Named list of annotated Seurat objects.
#' @export
annotate_samples_from_sheet <- function(seurat_list,
                                        sample_sheet,
                                        sample_id_col = "sample_id",
                                        metadata_cols = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required.", call. = FALSE)
  }

  check_sample_sheet_columns(sample_sheet, required_cols = sample_id_col)

  if (is.null(names(seurat_list)) || any(names(seurat_list) == "")) {
    stop("seurat_list must be a named list. Names should match sample IDs.", call. = FALSE)
  }

  sample_ids <- sample_sheet[[sample_id_col]]

  missing_in_sheet <- setdiff(names(seurat_list), sample_ids)
  if (length(missing_in_sheet) > 0) {
    stop(
      "These seurat_list names are missing from the sample sheet: ",
      paste(missing_in_sheet, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(metadata_cols)) {
    metadata_cols <- setdiff(colnames(sample_sheet), "sample_path")
  }

  metadata_cols <- unique(c(sample_id_col, metadata_cols))
  missing_metadata_cols <- setdiff(metadata_cols, colnames(sample_sheet))

  if (length(missing_metadata_cols) > 0) {
    stop(
      "These metadata columns are missing from the sample sheet: ",
      paste(missing_metadata_cols, collapse = ", "),
      call. = FALSE
    )
  }

  for (sid in names(seurat_list)) {
    row_idx <- which(sample_sheet[[sample_id_col]] == sid)

    if (length(row_idx) != 1) {
      stop(
        "Expected exactly one row for sample_id '", sid,
        "', but found ", length(row_idx), ".",
        call. = FALSE
      )
    }

    for (col in metadata_cols) {
      seurat_list[[sid]][[col]] <- sample_sheet[[col]][row_idx]
    }
  }

  return(seurat_list)
}

#' Merge a named list of Seurat objects
#'
#' Merges a list of Seurat objects and uses the list names as cell ID prefixes.
#'
#' @param seurat_list Named list of Seurat objects.
#' @param project_name Character. Project name for the merged Seurat object.
#'
#' @return A merged Seurat object.
#' @export
merge_seurat_samples <- function(seurat_list,
                                 project_name = "merged_seurat") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  if (length(seurat_list) < 2) {
    stop("seurat_list must contain at least two Seurat objects.", call. = FALSE)
  }

  if (is.null(names(seurat_list)) || any(names(seurat_list) == "")) {
    stop("seurat_list must be named before merging.", call. = FALSE)
  }

  merged <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = project_name
  )

  return(merged)
}

#' Summarize cells per sample
#'
#' Returns a table of cell counts per sample or metadata column.
#'
#' @param seu Seurat object.
#' @param group_col Character. Metadata column to summarize.
#'
#' @return A data frame with cell counts.
#' @export
summarize_cells_by_group <- function(seu,
                                     group_col = "orig.ident") {
  if (!group_col %in% colnames(seu@meta.data)) {
    stop("Column not found in metadata: ", group_col, call. = FALSE)
  }

  out <- as.data.frame(table(seu@meta.data[[group_col]]))
  colnames(out) <- c(group_col, "n_cells")

  return(out)
}
