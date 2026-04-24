#' Speckle/propeller utilities for cell proportion testing
#'
#' Helper functions for running propeller on all cells or selected identity groups.
#'
#' @keywords internal
NULL

#' Run propeller on a Seurat object
#'
#' @param seu Seurat object.
#' @param cluster_col Optional metadata column for clusters. If NULL, active identities are used.
#' @param sample_col Metadata column containing sample IDs.
#' @param group_col Metadata column containing biological groups.
#'
#' @return propeller result.
#' @export
run_propeller <- function(seu,
                          cluster_col = NULL,
                          sample_col = "orig.ident",
                          group_col = "condition") {
  if (!requireNamespace("speckle", quietly = TRUE)) {
    stop("Package 'speckle' is required.", call. = FALSE)
  }

  required_cols <- c(sample_col, group_col)
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  clusters <- if (is.null(cluster_col)) {
    Seurat::Idents(seu)
  } else {
    if (!cluster_col %in% colnames(seu@meta.data)) {
      stop("cluster_col not found in metadata: ", cluster_col, call. = FALSE)
    }
    seu[[cluster_col]][, 1]
  }

  speckle::propeller(
    clusters = clusters,
    sample = seu[[sample_col]][, 1],
    group = seu[[group_col]][, 1]
  )
}

#' Run propeller on selected identities
#'
#' @param seu Seurat object.
#' @param idents Character vector of identities to keep.
#' @param sample_col Metadata column containing sample IDs.
#' @param group_col Metadata column containing biological groups.
#'
#' @return propeller result.
#' @export
run_propeller_for_idents <- function(seu,
                                     idents,
                                     sample_col = "orig.ident",
                                     group_col = "condition") {
  seu_sub <- subset(seu, idents = idents)

  run_propeller(
    seu = seu_sub,
    cluster_col = NULL,
    sample_col = sample_col,
    group_col = group_col
  )
}

#' Run propeller for multiple identity groups
#'
#' @param seu Seurat object.
#' @param identity_sets Named list of character vectors.
#' @param sample_col Metadata column containing sample IDs.
#' @param group_col Metadata column containing biological groups.
#'
#' @return Named list of propeller results.
#' @export
run_propeller_sets <- function(seu,
                               identity_sets,
                               sample_col = "orig.ident",
                               group_col = "condition") {
  if (is.null(names(identity_sets)) || any(names(identity_sets) == "")) {
    stop("identity_sets must be a named list.", call. = FALSE)
  }

  results <- lapply(identity_sets, function(x) {
    run_propeller_for_idents(
      seu = seu,
      idents = x,
      sample_col = sample_col,
      group_col = group_col
    )
  })

  return(results)
}

#' Save propeller results
#'
#' @param propeller_results Propeller result or named list of propeller results.
#' @param output_dir Output directory.
#' @param prefix File prefix.
#'
#' @return Invisibly returns output paths.
#' @export
save_propeller_results <- function(propeller_results,
                                   output_dir = "results",
                                   prefix = "CellTypeProportions") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.list(propeller_results) && !is.data.frame(propeller_results)) {
    files <- vapply(names(propeller_results), function(nm) {
      file <- file.path(output_dir, paste0(prefix, "_", nm, ".csv"))
      write.csv(propeller_results[[nm]], file = file, row.names = TRUE)
      file
    }, character(1))

    return(invisible(files))
  }

  file <- file.path(output_dir, paste0(prefix, ".csv"))
  write.csv(propeller_results, file = file, row.names = TRUE)

  invisible(file)
}

#' Run Speckle ANOVA-style posthoc pairwise contrasts
#'
#' Useful after a multi-group propeller/ANOVA result when you want pairwise
#' comparisons between groups.
#'
#' @param seu Seurat object.
#' @param sample_col Metadata column containing biological replicate/sample IDs.
#' @param group_col Metadata column containing groups.
#' @param cluster_col Optional metadata column for clusters. If NULL, active identities are used.
#' @param group_levels Character vector specifying group order.
#' @param contrasts_list Optional named list of contrasts. If NULL, all pairwise contrasts are generated.
#' @param robust Logical. Passed to speckle::propeller.ttest().
#' @param trend Logical. Passed to speckle::propeller.ttest().
#' @param global_fdr Logical. Whether to add global BH correction across all posthoc tests.
#'
#' @return Data frame of pairwise posthoc results.
#' @export
run_propeller_posthoc <- function(seu,
                                  sample_col = "orig.ident",
                                  group_col = "condition",
                                  cluster_col = NULL,
                                  group_levels = NULL,
                                  contrasts_list = NULL,
                                  robust = TRUE,
                                  trend = TRUE,
                                  global_fdr = TRUE) {
  if (!requireNamespace("speckle", quietly = TRUE)) {
    stop("Package 'speckle' is required.", call. = FALSE)
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  required_cols <- c(sample_col, group_col)
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  clusters <- if (is.null(cluster_col)) {
    Seurat::Idents(seu)
  } else {
    if (!cluster_col %in% colnames(seu@meta.data)) {
      stop("cluster_col not found in metadata: ", cluster_col, call. = FALSE)
    }
    seu[[cluster_col]][, 1]
  }

  sample <- seu[[sample_col]][, 1]
  group <- seu[[group_col]][, 1]

  if (is.null(group_levels)) {
    group_levels <- unique(as.character(group))
  }

  group <- factor(group, levels = group_levels)

  prop_list <- speckle::getTransformedProps(
    clusters = clusters,
    sample = sample
  )

  samps <- colnames(prop_list$TransformedProps)

  sample_meta <- data.frame(
    sample = sample,
    group = group
  ) |>
    dplyr::distinct(.data$sample, .keep_all = TRUE) |>
    dplyr::filter(.data$sample %in% samps) |>
    dplyr::mutate(sample = factor(.data$sample, levels = samps)) |>
    dplyr::arrange(.data$sample)

  sample_meta$group <- factor(sample_meta$group, levels = group_levels)

  design <- stats::model.matrix(~ 0 + group, data = sample_meta)
  colnames(design) <- levels(sample_meta$group)

  if (is.null(contrasts_list)) {
    pairwise <- utils::combn(group_levels, 2, simplify = FALSE)

    contrasts_list <- lapply(pairwise, function(x) {
      contrast_string <- paste0(x[1], " - ", x[2])
      contrast_string
    })

    names(contrasts_list) <- vapply(
      pairwise,
      function(x) paste0(x[1], "_vs_", x[2]),
      character(1)
    )
  }

  contrast_matrix <- limma::makeContrasts(
    contrasts = unlist(contrasts_list),
    levels = design
  )

  posthoc_list <- lapply(colnames(contrast_matrix), function(cn) {
    res <- speckle::propeller.ttest(
      prop.list = prop_list,
      design = design,
      contrasts = contrast_matrix[, cn],
      robust = robust,
      trend = trend,
      sort = FALSE
    )

    res$Comparison <- cn
    res
  })

  posthoc <- dplyr::bind_rows(posthoc_list, .id = "Contrast_ID")

  if (global_fdr && "P.Value" %in% colnames(posthoc)) {
    posthoc$FDR_global <- stats::p.adjust(posthoc$P.Value, method = "BH")
  }

  posthoc <- tibble::rownames_to_column(posthoc, "cluster_raw")

  posthoc$cluster <- sub("\\.\\.\\..*$", "", posthoc$cluster_raw)

  return(posthoc)
}
