#' Doublet detection utilities for Seurat objects
#'
#' Helper functions for estimating 10x multiplet rates, selecting PCs,
#' running DoubletFinder per sample, annotating doublets, and retaining singlets.
#'
#' @keywords internal
NULL

#' Estimate 10x multiplet rate from recovered cell number
#'
#' Uses approximate 10x Genomics multiplet-rate expectations based on recovered
#' cell numbers.
#'
#' @param n_cells Integer. Number of recovered cells.
#'
#' @return Numeric multiplet rate.
#' @export
estimate_10x_multiplet_rate <- function(n_cells) {
  if (!is.numeric(n_cells) || length(n_cells) != 1 || is.na(n_cells)) {
    stop("'n_cells' must be a single numeric value.", call. = FALSE)
  }

  multiplet_rates_10x <- data.frame(
    Multiplet_rate = c(
      0.004, 0.008, 0.016, 0.023, 0.031, 0.039,
      0.046, 0.054, 0.061, 0.069, 0.076
    ),
    Loaded_cells = c(
      800, 1600, 3200, 4800, 6400, 8000,
      9600, 11200, 12800, 14400, 16000
    ),
    Recovered_cells = c(
      500, 1000, 2000, 3000, 4000, 5000,
      6000, 7000, 8000, 9000, 10000
    )
  )

  if (n_cells <= min(multiplet_rates_10x$Recovered_cells)) {
    return(min(multiplet_rates_10x$Multiplet_rate))
  }

  if (n_cells >= max(multiplet_rates_10x$Recovered_cells)) {
    return(max(multiplet_rates_10x$Multiplet_rate))
  }

  rate <- multiplet_rates_10x |>
    dplyr::filter(Recovered_cells <= n_cells) |>
    dplyr::slice(which.max(Recovered_cells)) |>
    dplyr::pull(Multiplet_rate)

  return(as.numeric(rate))
}

#' Choose significant PCs from a PCA reduction
#'
#' Uses a heuristic based on cumulative variance and elbow-like drops in
#' percent standard deviation.
#'
#' @param seu Seurat object with a PCA reduction.
#' @param reduction Character. Reduction name, usually `"pca"`.
#' @param cumulative_cutoff Numeric. Cumulative variance cutoff.
#' @param percent_stdv_cutoff Numeric. Individual PC percent standard deviation cutoff.
#' @param delta_cutoff Numeric. Minimum drop between consecutive PCs.
#' @param default_pcs Integer. Fallback number of PCs if heuristic fails.
#'
#' @return Integer number of PCs.
#' @export
choose_significant_pcs <- function(seu,
                                   reduction = "pca",
                                   cumulative_cutoff = 90,
                                   percent_stdv_cutoff = 5,
                                   delta_cutoff = 0.1,
                                   default_pcs = 20) {
  if (!reduction %in% names(seu@reductions)) {
    stop("Reduction not found in Seurat object: ", reduction, call. = FALSE)
  }

  stdv <- seu[[reduction]]@stdev

  if (length(stdv) < 2 || all(is.na(stdv))) {
    warning("Could not estimate significant PCs. Returning default_pcs.")
    return(default_pcs)
  }

  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)

  co1 <- which(cumulative > cumulative_cutoff & percent_stdv < percent_stdv_cutoff)[1]

  pc_drops <- percent_stdv[1:(length(percent_stdv) - 1)] -
    percent_stdv[2:length(percent_stdv)]

  co2 <- sort(which(pc_drops > delta_cutoff), decreasing = TRUE)[1] + 1

  pcs <- suppressWarnings(min(co1, co2, na.rm = TRUE))

  if (!is.finite(pcs) || is.na(pcs)) {
    pcs <- default_pcs
  }

  pcs <- min(pcs, length(stdv))
  pcs <- max(pcs, 2)

  return(as.integer(pcs))
}

#' Preprocess a sample for DoubletFinder
#'
#' Runs a lightweight Seurat preprocessing workflow required by DoubletFinder.
#'
#' @param seu Seurat object.
#' @param assay Character. Assay to use for preprocessing.
#' @param resolution Numeric. Clustering resolution for homotypic doublet modeling.
#'
#' @return Preprocessed Seurat object.
#' @export
preprocess_for_doubletfinder <- function(seu,
                                         assay = "RNA",
                                         resolution = 0.1) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  Seurat::DefaultAssay(seu) <- assay

  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, nfeatures.print = 10, verbose = FALSE)

  n_pcs <- choose_significant_pcs(seu, reduction = "pca")

  seu <- Seurat::RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, dims = 1:n_pcs, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = resolution, verbose = FALSE)

  attr(seu, "doubletfinder_n_pcs") <- n_pcs

  return(seu)
}

#' Find optimal pK for DoubletFinder
#'
#' @param seu Preprocessed Seurat object.
#' @param pcs Integer vector of PCs to use.
#' @param sct Logical. Passed to DoubletFinder functions.
#'
#' @return Numeric pK value.
#' @export
find_optimal_pk <- function(seu,
                            pcs,
                            sct = FALSE) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    stop("Package 'DoubletFinder' is required.", call. = FALSE)
  }

  sweep_list <- DoubletFinder::paramSweep(seu, PCs = pcs, sct = sct)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep_list)
  bcmvn <- DoubletFinder::find.pK(sweep_stats)

  optimal_pk <- bcmvn |>
    dplyr::filter(BCmetric == max(BCmetric, na.rm = TRUE)) |>
    dplyr::pull(pK) |>
    as.character() |>
    as.numeric()

  if (length(optimal_pk) == 0 || is.na(optimal_pk[1])) {
    stop("Could not determine optimal pK.", call. = FALSE)
  }

  return(optimal_pk[1])
}

#' Run DoubletFinder on a single Seurat sample
#'
#' @param seu_sample Seurat object corresponding to one sample.
#' @param multiplet_rate Optional numeric multiplet rate. If NULL, estimated
#'   from recovered cells using `estimate_10x_multiplet_rate()`.
#' @param assay Character. Assay to use.
#' @param sct Logical. Whether DoubletFinder should use SCT mode.
#' @param resolution Numeric. Clustering resolution used for homotypic modeling.
#'
#' @return Data frame with `cell_id` and `doublet_finder` columns.
#' @export
run_doubletfinder_custom <- function(seu_sample,
                                     multiplet_rate = NULL,
                                     assay = "RNA",
                                     sct = FALSE,
                                     resolution = 0.1) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    stop("Package 'DoubletFinder' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  if (is.null(multiplet_rate)) {
    multiplet_rate <- estimate_10x_multiplet_rate(ncol(seu_sample))
  }

  sample_name <- unique(seu_sample$orig.ident)
  message("Running DoubletFinder for sample: ", paste(sample_name, collapse = ", "))
  message("Using multiplet rate: ", multiplet_rate)

  sample <- preprocess_for_doubletfinder(
    seu = seu_sample,
    assay = assay,
    resolution = resolution
  )

  n_pcs <- attr(sample, "doubletfinder_n_pcs")
  pcs <- 1:n_pcs

  optimal_pk <- find_optimal_pk(
    seu = sample,
    pcs = pcs,
    sct = sct
  )

  annotations <- sample@meta.data$seurat_clusters
  homotypic_prop <- DoubletFinder::modelHomotypic(annotations)

  n_exp_poi <- round(multiplet_rate * ncol(sample))
  n_exp_poi_adj <- round(n_exp_poi * (1 - homotypic_prop))

  sample <- DoubletFinder::doubletFinder(
    seu = sample,
    PCs = pcs,
    pK = optimal_pk,
    nExp = n_exp_poi_adj,
    sct = sct
  )

  df_cols <- grep("DF.classifications", colnames(sample@meta.data), value = TRUE)

  if (length(df_cols) == 0) {
    stop("DoubletFinder classification column was not found.", call. = FALSE)
  }

  final_df_col <- df_cols[length(df_cols)]
  sample@meta.data$doublet_finder <- sample@meta.data[[final_df_col]]

  res <- sample@meta.data["doublet_finder"]
  res <- tibble::rownames_to_column(res, "cell_id")

  return(res)
}

#' Annotate doublets by sample
#'
#' Splits a merged Seurat object by sample, runs DoubletFinder on each sample,
#' and adds a `doublet_finder` metadata column to the original object.
#'
#' @param seu Merged Seurat object.
#' @param split_by Character. Metadata column used to split samples.
#' @param multiplet_rate Optional numeric scalar or named numeric vector. If a
#'   named vector is supplied, names should match sample IDs.
#' @param assay Character. Assay to use.
#' @param sct Logical. Whether to run DoubletFinder in SCT mode.
#' @param resolution Numeric. Clustering resolution for homotypic modeling.
#'
#' @return Seurat object with `doublet_finder` metadata column.
#' @export
annotate_doublets_by_sample <- function(seu,
                                        split_by = "orig.ident",
                                        multiplet_rate = NULL,
                                        assay = "RNA",
                                        sct = FALSE,
                                        resolution = 0.1) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  if (!split_by %in% colnames(seu@meta.data)) {
    stop("split_by column not found in metadata: ", split_by, call. = FALSE)
  }

  sample_list <- Seurat::SplitObject(seu, split.by = split_by)

  doublet_metadata <- lapply(names(sample_list), function(sample_id) {
    sample_rate <- NULL

    if (!is.null(multiplet_rate)) {
      if (length(multiplet_rate) == 1 && is.null(names(multiplet_rate))) {
        sample_rate <- multiplet_rate
      } else if (!is.null(names(multiplet_rate)) && sample_id %in% names(multiplet_rate)) {
        sample_rate <- multiplet_rate[[sample_id]]
      } else {
        stop(
          "multiplet_rate must be NULL, a scalar, or a named vector matching sample IDs.",
          call. = FALSE
        )
      }
    }

    run_doubletfinder_custom(
      seu_sample = sample_list[[sample_id]],
      multiplet_rate = sample_rate,
      assay = assay,
      sct = sct,
      resolution = resolution
    )
  })

  doublet_metadata <- dplyr::bind_rows(doublet_metadata)
  rownames(doublet_metadata) <- doublet_metadata$cell_id
  doublet_metadata$cell_id <- NULL

  seu <- Seurat::AddMetaData(seu, metadata = doublet_metadata)

  return(seu)
}

#' Keep only singlets
#'
#' @param seu Seurat object with a `doublet_finder` metadata column.
#' @param doublet_col Character. Metadata column containing singlet/doublet labels.
#' @param singlet_label Character. Label used for singlets.
#'
#' @return Filtered Seurat object containing only singlets.
#' @export
keep_singlets <- function(seu,
                          doublet_col = "doublet_finder",
                          singlet_label = "Singlet") {
  if (!doublet_col %in% colnames(seu@meta.data)) {
    stop("doublet_col not found in metadata: ", doublet_col, call. = FALSE)
  }

  keep_cells <- rownames(seu@meta.data)[seu@meta.data[[doublet_col]] == singlet_label]
  subset(seu, cells = keep_cells)
}

#' Summarize DoubletFinder calls
#'
#' @param seu Seurat object with doublet classifications.
#' @param group_col Character. Metadata column for grouping.
#' @param doublet_col Character. Metadata column with doublet labels.
#'
#' @return Data frame with singlet/doublet counts per group.
#' @export
summarize_doublets <- function(seu,
                               group_col = "orig.ident",
                               doublet_col = "doublet_finder") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  required_cols <- c(group_col, doublet_col)
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  out <- seu@meta.data |>
    dplyr::count(.data[[group_col]], .data[[doublet_col]], name = "n_cells") |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::mutate(fraction = n_cells / sum(n_cells)) |>
    dplyr::ungroup()

  return(out)
}
