#' Module scoring utilities for Seurat objects
#'
#' Helper functions for disease-associated gene modules, module score plotting,
#' and percentage of cells expressing module genes.
#'
#' @keywords internal
NULL

#' Add a gene module score to a Seurat object
#'
#' @param seu Seurat object.
#' @param gene_set Character vector of genes.
#' @param module_name Character prefix for module score.
#' @param assay Assay to use.
#'
#' @return Seurat object with module score added.
#' @export
add_gene_module_score <- function(seu,
                                  gene_set,
                                  module_name = "Module_Score",
                                  assay = "RNA") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  Seurat::DefaultAssay(seu) <- assay

  gene_set <- unique(gene_set)
  gene_set <- gene_set[gene_set != ""]
  gene_set <- intersect(gene_set, rownames(seu))

  if (length(gene_set) == 0) {
    stop("None of the genes in gene_set were found in the Seurat object.", call. = FALSE)
  }

  seu <- Seurat::AddModuleScore(
    object = seu,
    features = list(gene_set),
    name = module_name,
    assay = assay
  )

  attr(seu, paste0(module_name, "_genes_used")) <- gene_set

  return(seu)
}

#' Prepare module score plotting data
#'
#' @param seu Seurat object.
#' @param score_col Module score metadata column, usually module_name + "1".
#' @param celltype_col Metadata column for cell type/subcluster.
#' @param group_col Metadata column for disease/control group.
#' @param celltype_levels Optional cell type order.
#' @param group_levels Optional group order.
#'
#' @return Data frame.
#' @export
make_module_score_df <- function(seu,
                                 score_col,
                                 celltype_col = "cellclass",
                                 group_col = "FDX",
                                 celltype_levels = NULL,
                                 group_levels = NULL) {
  required_cols <- c(score_col, celltype_col, group_col)
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df <- seu@meta.data

  if (!is.null(celltype_levels)) {
    df[[celltype_col]] <- factor(df[[celltype_col]], levels = celltype_levels)
  }

  if (!is.null(group_levels)) {
    df[[group_col]] <- factor(df[[group_col]], levels = group_levels)
  }

  df
}

#' Run Wilcoxon tests for module scores by cell type
#'
#' @param plot_df Data frame from make_module_score_df().
#' @param score_col Module score column.
#' @param celltype_col Cell type column.
#' @param group_col Group column.
#' @param group1 First group.
#' @param group2 Second group.
#'
#' @return Data frame of BH-adjusted p-values.
#' @export
wilcox_module_by_celltype <- function(plot_df,
                                      score_col,
                                      celltype_col = "cellclass",
                                      group_col = "FDX",
                                      group1 = "Control",
                                      group2 = "AD+CAA") {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  plot_df |>
    dplyr::filter(.data[[group_col]] %in% c(group1, group2)) |>
    dplyr::group_by(.data[[celltype_col]]) |>
    dplyr::group_modify(~{
      test <- stats::wilcox.test(
        stats::as.formula(paste(score_col, "~", group_col)),
        data = .x
      )

      tibble::tibble(
        group1 = group1,
        group2 = group2,
        p = test$p.value
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      p.adj = stats::p.adjust(p, method = "BH"),
      p.adj.label = paste0("BH p = ", signif(p.adj, 2))
    )
}

#' Calculate percent of cells expressing at least N module genes
#'
#' @param seu Seurat object.
#' @param gene_set Character vector of genes.
#' @param threshold Integer. Minimum number of genes expressed.
#' @param assay Assay to use.
#' @param layer Layer/slot to use.
#' @param celltype_col Cell type metadata column.
#' @param group_col Group metadata column.
#' @param celltype_levels Optional cell type order.
#' @param group_levels Optional group order.
#'
#' @return Data frame with percentages.
#' @export
percent_cells_expressing_module <- function(seu,
                                            gene_set,
                                            threshold,
                                            assay = "RNA",
                                            layer = "data",
                                            celltype_col = "cellclass",
                                            group_col = "FDX",
                                            celltype_levels = NULL,
                                            group_levels = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  gene_set <- unique(gene_set)
  gene_set <- gene_set[gene_set != ""]
  gene_set <- intersect(gene_set, rownames(seu))

  expr_mat <- Seurat::GetAssayData(seu, assay = assay, layer = layer)
  binary_mat <- expr_mat[gene_set, , drop = FALSE] > 0
  n_genes_per_cell <- Matrix::colSums(binary_mat)

  df <- data.frame(
    celltype = seu[[celltype_col]][, 1],
    group = seu[[group_col]][, 1],
    expressed_module = n_genes_per_cell > threshold
  )

  out <- df |>
    dplyr::group_by(celltype, group) |>
    dplyr::summarise(
      pct_expressing = 100 * mean(expressed_module),
      n_cells = dplyr::n(),
      .groups = "drop"
    )

  if (!is.null(celltype_levels)) {
    out$celltype <- factor(out$celltype, levels = celltype_levels)
  }

  if (!is.null(group_levels)) {
    out$group <- factor(out$group, levels = group_levels)
  }

  out
}
