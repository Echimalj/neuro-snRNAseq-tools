#' Atlas similarity utilities
#'
#' Functions for comparing Seurat clusters or subclusters to external marker
#' gene sets using cosine similarity and empirical null testing.
#'
#' @keywords internal
NULL

#' Compute cosine similarity
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#'
#' @return Numeric cosine similarity.
#' @export
cosine_similarity <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.", call. = FALSE)
  }

  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))

  if (denom == 0) {
    return(NA_real_)
  }

  sum(x * y) / denom
}

#' Compare clusters to a marker gene set using cosine similarity
#'
#' @param seurat_object Seurat object.
#' @param gene_list Character vector of external marker genes.
#' @param cluster_column Metadata column containing cluster labels.
#' @param assay Assay to use.
#' @param layer Expression layer. For Seurat v5, usually `"data"`.
#' @param logfc_threshold LogFC threshold for FindAllMarkers.
#' @param min_pct Minimum percent expression for FindAllMarkers.
#' @param top_n_genes Number of marker genes to keep per cluster. Use Inf for all.
#' @param n_null Number of null permutations.
#' @param seed Random seed.
#' @param p_adj_method P-value adjustment method.
#'
#' @return List with similarity table, top genes per cluster, and marker table.
#' @export
cluster_marker_cosine_similarity <- function(seurat_object,
                                             gene_list,
                                             cluster_column = "seurat_clusters",
                                             assay = "RNA",
                                             layer = "data",
                                             logfc_threshold = 0.25,
                                             min_pct = 0.1,
                                             top_n_genes = Inf,
                                             n_null = 1000,
                                             seed = 42,
                                             p_adj_method = "BH") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  set.seed(seed)

  if (!cluster_column %in% colnames(seurat_object@meta.data)) {
    stop("cluster_column not found in metadata: ", cluster_column, call. = FALSE)
  }

  gene_list <- unique(gene_list)
  gene_list <- gene_list[gene_list != ""]

  Seurat::Idents(seurat_object) <- seurat_object[[cluster_column]][, 1]

  all_markers <- Seurat::FindAllMarkers(
    object = seurat_object,
    assay = assay,
    slot = layer,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    only.pos = TRUE
  )

  filtered_markers <- all_markers |>
    dplyr::filter(.data$gene %in% gene_list)

  if (nrow(filtered_markers) == 0) {
    stop("None of the input genes were found among cluster markers.", call. = FALSE)
  }

  expr_data <- Seurat::GetAssayData(
    seurat_object,
    assay = assay,
    layer = layer
  )

  all_genes <- rownames(expr_data)

  top_markers <- filtered_markers |>
    dplyr::group_by(.data$cluster) |>
    dplyr::slice_max(
      order_by = .data$avg_log2FC,
      n = top_n_genes,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()

  clusters <- unique(top_markers$cluster)

  similarity_scores <- numeric(length(clusters))
  p_values <- numeric(length(clusters))
  names(similarity_scores) <- as.character(clusters)
  names(p_values) <- as.character(clusters)

  top_genes_per_cluster <- list()

  background_genes <- setdiff(all_genes, gene_list)

  for (cluster_id in clusters) {
    cluster_id_chr <- as.character(cluster_id)

    marker_genes <- top_markers |>
      dplyr::filter(.data$cluster == cluster_id) |>
      dplyr::pull(.data$gene) |>
      intersect(all_genes)

    if (length(marker_genes) == 0) {
      similarity_scores[cluster_id_chr] <- NA_real_
      p_values[cluster_id_chr] <- NA_real_
      next
    }

    cells_in_cluster <- colnames(seurat_object)[
      seurat_object@meta.data[[cluster_column]] == cluster_id
    ]

    if (length(cells_in_cluster) == 0) {
      similarity_scores[cluster_id_chr] <- NA_real_
      p_values[cluster_id_chr] <- NA_real_
      next
    }

    avg_expr <- Matrix::rowMeans(
      expr_data[marker_genes, cells_in_cluster, drop = FALSE]
    )

    observed_similarity <- cosine_similarity(
      x = as.numeric(avg_expr),
      y = rep(1, length(avg_expr))
    )

    null_similarities <- numeric(n_null)

    for (i in seq_len(n_null)) {
      random_genes <- sample(background_genes, length(marker_genes))
      random_expr <- Matrix::rowMeans(
        expr_data[random_genes, cells_in_cluster, drop = FALSE]
      )

      null_similarities[i] <- cosine_similarity(
        x = as.numeric(random_expr),
        y = rep(1, length(random_expr))
      )
    }

    p_val <- mean(null_similarities >= observed_similarity)

    similarity_scores[cluster_id_chr] <- observed_similarity
    p_values[cluster_id_chr] <- p_val
    top_genes_per_cluster[[cluster_id_chr]] <- sort(avg_expr, decreasing = TRUE)
  }

  adjusted_p_values <- stats::p.adjust(p_values, method = p_adj_method)

  results_df <- tibble::tibble(
    Cluster = names(similarity_scores),
    Similarity = as.numeric(similarity_scores),
    P_value = as.numeric(p_values),
    Adjusted_P_value = as.numeric(adjusted_p_values)
  ) |>
    dplyr::arrange(dplyr::desc(.data$Similarity))

  return(list(
    similarity_table = results_df,
    top_genes_per_cluster = top_genes_per_cluster,
    marker_table = top_markers
  ))
}

#' Plot cosine similarity results
#'
#' @param similarity_table Output similarity table.
#' @param significance_cutoff Adjusted p-value cutoff.
#'
#' @return ggplot object.
#' @export
plot_cosine_similarity <- function(similarity_table,
                                   significance_cutoff = 0.01) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }

  ggplot2::ggplot(
    similarity_table,
    ggplot2::aes(
      x = reorder(.data$Cluster, -.data$Similarity),
      y = .data$Similarity,
      fill = .data$Adjusted_P_value < significance_cutoff
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "darkgreen", "FALSE" = "grey70"),
      name = paste0("Adj. p < ", significance_cutoff)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Cosine Similarity by Cluster",
      x = "Cluster",
      y = "Cosine Similarity"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Run cosine similarity for multiple atlas gene sets
#'
#' @param seurat_object Seurat object.
#' @param gene_sets Named list of marker gene vectors.
#' @param cluster_column Metadata column containing cluster labels.
#' @param assay Assay to use.
#' @param layer Expression layer.
#' @param n_null Number of null permutations.
#'
#' @return Named list of similarity outputs.
#' @export
run_atlas_similarity_list <- function(seurat_object,
                                      gene_sets,
                                      cluster_column = "seurat_clusters",
                                      assay = "RNA",
                                      layer = "data",
                                      n_null = 1000) {
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("gene_sets must be a named list.", call. = FALSE)
  }

  results <- lapply(names(gene_sets), function(gs_name) {
    message("Running atlas similarity for: ", gs_name)

    cluster_marker_cosine_similarity(
      seurat_object = seurat_object,
      gene_list = gene_sets[[gs_name]],
      cluster_column = cluster_column,
      assay = assay,
      layer = layer,
      n_null = n_null
    )
  })

  names(results) <- names(gene_sets)
  return(results)
}

#' Save atlas similarity results
#'
#' @param similarity_results Output from cluster_marker_cosine_similarity() or run_atlas_similarity_list().
#' @param output_dir Output directory.
#' @param prefix File prefix.
#'
#' @return Invisibly returns file paths.
#' @export
save_atlas_similarity_results <- function(similarity_results,
                                          output_dir = "results/atlas_similarity",
                                          prefix = "atlas_similarity") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(similarity_results$similarity_table)) {
    file <- file.path(output_dir, paste0(prefix, "_similarity_table.csv"))
    write.csv(similarity_results$similarity_table, file = file, row.names = FALSE)
    return(invisible(file))
  }

  files <- vapply(names(similarity_results), function(nm) {
    file <- file.path(output_dir, paste0(prefix, "_", nm, "_similarity_table.csv"))
    write.csv(similarity_results[[nm]]$similarity_table, file = file, row.names = FALSE)
    file
  }, character(1))

  invisible(files)
}
