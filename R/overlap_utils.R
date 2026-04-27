#' Gene set overlap utilities
#'
#' Helper functions for Fisher exact overlap testing between DEG sets and
#' external gene sets.
#'
#' @keywords internal
NULL

#' Fisher exact test for DEG overlap
#'
#' @param deg_set1 First gene set.
#' @param deg_set2 Second gene set.
#' @param universe Background gene universe.
#'
#' @return Tibble with overlap statistics.
#' @export
fisher_deg_overlap <- function(deg_set1,
                               deg_set2,
                               universe) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  deg_set1 <- unique(stats::na.omit(as.character(deg_set1)))
  deg_set2 <- unique(stats::na.omit(as.character(deg_set2)))
  universe <- unique(stats::na.omit(as.character(universe)))

  deg_set1 <- intersect(deg_set1, universe)
  deg_set2 <- intersect(deg_set2, universe)

  a <- length(intersect(deg_set1, deg_set2))
  b <- length(setdiff(deg_set1, deg_set2))
  c <- length(setdiff(deg_set2, deg_set1))
  d <- length(setdiff(universe, union(deg_set1, deg_set2)))

  cont <- matrix(c(a, b, c, d), nrow = 2)

  test <- stats::fisher.test(cont)

  expected <- (length(deg_set1) * length(deg_set2)) / length(universe)
  direction <- if (a >= expected) "enriched" else "depleted"

  tibble::tibble(
    overlap = a,
    expected = expected,
    fisher_p = test$p.value,
    odds_ratio = unname(test$estimate),
    direction = direction,
    set1_size = length(deg_set1),
    set2_size = length(deg_set2),
    universe_size = length(universe)
  )
}

#' Run multiple Fisher overlap tests
#'
#' @param query_sets Named list of query gene sets.
#' @param reference_sets Named list of reference gene sets.
#' @param universe Background gene universe.
#' @param p_adjust_method Multiple testing correction method.
#'
#' @return Data frame of overlap results.
#' @export
run_fisher_overlap_grid <- function(query_sets,
                                    reference_sets,
                                    universe,
                                    p_adjust_method = "BH") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  if (is.null(names(query_sets)) || any(names(query_sets) == "")) {
    stop("query_sets must be a named list.", call. = FALSE)
  }

  if (is.null(names(reference_sets)) || any(names(reference_sets) == "")) {
    stop("reference_sets must be a named list.", call. = FALSE)
  }

  out <- list()

  idx <- 1
  for (q in names(query_sets)) {
    for (r in names(reference_sets)) {
      out[[idx]] <- fisher_deg_overlap(
        deg_set1 = query_sets[[q]],
        deg_set2 = reference_sets[[r]],
        universe = universe
      ) |>
        dplyr::mutate(
          query_set = q,
          reference_set = r
        )

      idx <- idx + 1
    }
  }

  out <- dplyr::bind_rows(out) |>
    dplyr::mutate(
      q_BH = stats::p.adjust(.data$fisher_p, method = p_adjust_method)
    ) |>
    dplyr::select(
      .data$query_set,
      .data$reference_set,
      .data$overlap,
      .data$expected,
      .data$odds_ratio,
      .data$fisher_p,
      .data$q_BH,
      .data$direction,
      .data$set1_size,
      .data$set2_size,
      .data$universe_size
    )

  out
}

#' Run cellclass-vs-reference overlap workflow
#'
#' @param seu Seurat object.
#' @param celltype_col Metadata column containing cell classes.
#' @param reference_sets Named list of external gene sets.
#' @param assay Assay to use.
#' @param layer_counts Count layer for universe definition.
#' @param min_pct Minimum percent expression for FindMarkers.
#' @param min_cells_group Minimum cells per group.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param direction Direction of DE genes to test.
#' @param universe_genes Optional universe.
#' @param test_use Seurat FindMarkers test.
#'
#' @return Data frame of overlap results.
#' @export
run_cellclass_overlap <- function(seu,
                                  celltype_col = "cellclass",
                                  reference_sets,
                                  assay = "RNA",
                                  layer_counts = "counts",
                                  min_pct = 0.1,
                                  min_cells_group = 10,
                                  padj_cutoff = 0.05,
                                  direction = c("up", "down"),
                                  universe_genes = NULL,
                                  test_use = "wilcox") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  direction <- match.arg(direction)

  if (!celltype_col %in% colnames(seu@meta.data)) {
    stop("celltype_col not found in metadata: ", celltype_col, call. = FALSE)
  }

  if (is.null(universe_genes)) {
    Seurat::DefaultAssay(seu) <- assay
    counts <- Seurat::GetAssayData(seu, assay = assay, layer = layer_counts)
    universe_genes <- rownames(seu)[Matrix::rowSums(counts > 0) > 0]
  }

  celltypes <- sort(unique(as.character(seu@meta.data[[celltype_col]])))

  results <- purrr::map(celltypes, function(ct) {
    flag_col <- ".__tmp_flag__."
    seu[[flag_col]] <- ifelse(as.character(seu@meta.data[[celltype_col]]) == ct, ct, "Other")
    Seurat::Idents(seu) <- flag_col

    n_ct <- sum(Seurat::Idents(seu) == ct)
    n_ot <- sum(Seurat::Idents(seu) == "Other")

    if (n_ct < min_cells_group || n_ot < min_cells_group) {
      return(tibble::tibble(
        cellclass = ct,
        n_cellclass = n_ct,
        n_other = n_ot,
        n_sig = NA_integer_,
        n_set = NA_integer_,
        test = names(reference_sets),
        skipped = TRUE
      ))
    }

    markers <- Seurat::FindMarkers(
      seu,
      ident.1 = ct,
      ident.2 = "Other",
      test.use = test_use,
      logfc.threshold = 0,
      min.pct = min_pct,
      min.cells.group = min_cells_group
    ) |>
      tibble::rownames_to_column("gene") |>
      dplyr::filter(
        !is.na(.data$p_val_adj),
        .data$p_val_adj < padj_cutoff
      )

    gene_set <- if (direction == "up") {
      markers |>
        dplyr::filter(.data$avg_log2FC > 0) |>
        dplyr::pull(.data$gene) |>
        unique()
    } else {
      markers |>
        dplyr::filter(.data$avg_log2FC < 0) |>
        dplyr::pull(.data$gene) |>
        unique()
    }

    overlap_results <- lapply(names(reference_sets), function(ref_name) {
      fisher_deg_overlap(
        deg_set1 = gene_set,
        deg_set2 = reference_sets[[ref_name]],
        universe = universe_genes
      ) |>
        dplyr::mutate(
          cellclass = ct,
          n_cellclass = n_ct,
          n_other = n_ot,
          n_sig = nrow(markers),
          n_set = length(gene_set),
          test = paste0("Cellclass vs ", ref_name),
          skipped = FALSE
        )
    })

    dplyr::bind_rows(overlap_results)
  })

  dplyr::bind_rows(results) |>
    dplyr::mutate(q_BH = stats::p.adjust(.data$fisher_p, method = "BH"))
}
