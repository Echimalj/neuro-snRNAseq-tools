#' GSEA utilities
#'
#' Wrappers for running fgsea across cell types and generating
#' direction-specific enrichment heatmaps.
#'
#' @keywords internal
NULL


# ============================================================
# Gene set loading
# ============================================================

#' Load MSigDB pathways using msigdbr
#'
#' @param species Species name, e.g. "Homo sapiens" or "Mus musculus".
#' @param category MSigDB category.
#' @param subcategory MSigDB subcategory.
#' @param gene_col Gene column to use.
#'
#' @return Named list of pathways.
#' @export
load_msigdb_pathways <- function(species = "Homo sapiens",
                                 category = "C5",
                                 subcategory = "GO:BP",
                                 gene_col = "gene_symbol") {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required.", call. = FALSE)
  }

  msig <- msigdbr::msigdbr(
    species = species,
    category = category,
    subcategory = subcategory
  )

  split(msig[[gene_col]], msig$gs_name)
}


# ============================================================
# GSEA wrapper
# ============================================================

#' Run fgsea for one cell type
#'
#' @param df DEG/stat table for one cell type.
#' @param pathways Named pathway list.
#' @param gene_col Gene column.
#' @param stat_col Ranking statistic column.
#' @param min_size Minimum pathway size.
#' @param max_size Maximum pathway size.
#' @param remove_ensembl Logical. Remove Ensembl IDs when gene symbols are expected.
#'
#' @return fgsea result tibble.
#' @export
run_fgsea_one_celltype <- function(df,
                                   pathways,
                                   gene_col = "gene",
                                   stat_col = "stat",
                                   min_size = 15,
                                   max_size = 500,
                                   remove_ensembl = TRUE) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  required_cols <- c(gene_col, stat_col)
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df_ranked <- df |>
    dplyr::filter(
      !is.na(.data[[stat_col]]),
      !is.na(.data[[gene_col]])
    )

  if (remove_ensembl) {
    df_ranked <- df_ranked |>
      dplyr::filter(!grepl("^ENS", .data[[gene_col]]))
  }

  df_ranked <- df_ranked |>
    dplyr::group_by(.data[[gene_col]]) |>
    dplyr::slice_max(
      order_by = abs(.data[[stat_col]]),
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()

  ranks <- df_ranked[[stat_col]]
  names(ranks) <- df_ranked[[gene_col]]
  ranks <- sort(ranks, decreasing = TRUE)

  fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = min_size,
    maxSize = max_size
  ) |>
    tibble::as_tibble()
}


#' Run fgsea across cell types
#'
#' @param df Combined DEG/stat table.
#' @param pathways Named pathway list.
#' @param celltype_col Cell type column.
#' @param gene_col Gene column.
#' @param stat_col Ranking statistic column.
#' @param output_file Optional output CSV.
#' @param ... Additional arguments passed to run_fgsea_one_celltype().
#'
#' @return Combined fgsea result tibble.
#' @export
run_fgsea_by_celltype <- function(df,
                                  pathways,
                                  celltype_col = "cellclass",
                                  gene_col = "gene",
                                  stat_col = "stat",
                                  output_file = NULL,
                                  ...) {
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required.", call. = FALSE)
  }

  if (!celltype_col %in% colnames(df)) {
    stop("celltype_col not found: ", celltype_col, call. = FALSE)
  }

  gsea_list <- df |>
    split(df[[celltype_col]]) |>
    purrr::map(
      ~ run_fgsea_one_celltype(
        df = .x,
        pathways = pathways,
        gene_col = gene_col,
        stat_col = stat_col,
        ...
      )
    )

  gsea_all <- dplyr::bind_rows(
    lapply(names(gsea_list), function(ct) {
      gsea_list[[ct]] |>
        dplyr::mutate(cellclass = ct)
    })
  )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

    gsea_save <- gsea_all |>
      dplyr::mutate(
        leadingEdge = if ("leadingEdge" %in% colnames(gsea_all)) {
          sapply(.data$leadingEdge, paste, collapse = ";")
        } else {
          NA_character_
        }
      )

    write.csv(gsea_save, output_file, row.names = FALSE)
  }

  return(gsea_all)
}


# ============================================================
# GSEA cleaning/filtering
# ============================================================

#' Clean GO-style pathway names
#'
#' @param x Pathway names.
#'
#' @return Cleaned pathway names.
#' @export
clean_go_names <- function(x) {
  x |>
    gsub("^GOBP_", "", x = _) |>
    gsub("^GO_", "", x = _) |>
    gsub("_", " ", x = _) |>
    tools::toTitleCase()
}


#' Filter GSEA results
#'
#' @param gsea_df fgsea result dataframe.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param nes_cutoff Absolute NES cutoff.
#' @param min_size Minimum pathway size.
#' @param max_size Maximum pathway size.
#' @param positive_label Label for NES > 0.
#' @param negative_label Label for NES < 0.
#'
#' @return Filtered GSEA dataframe.
#' @export
filter_gsea_results <- function(gsea_df,
                                padj_cutoff = 0.05,
                                nes_cutoff = 1.5,
                                min_size = 15,
                                max_size = 500,
                                positive_label = "Condition1 enriched",
                                negative_label = "Condition2 enriched") {
  gsea_df |>
    dplyr::mutate(
      pathway_clean = clean_go_names(.data$pathway),
      direction = dplyr::case_when(
        .data$NES > 0 ~ positive_label,
        .data$NES < 0 ~ negative_label,
        TRUE ~ "Neutral"
      )
    ) |>
    dplyr::filter(
      !is.na(.data$padj),
      .data$padj < padj_cutoff,
      abs(.data$NES) >= nes_cutoff,
      .data$size >= min_size,
      .data$size <= max_size
    )
}


# ============================================================
# Heatmap matrix construction
# ============================================================

#' Build direction-specific GSEA heatmap matrix
#'
#' @param gsea_df Full GSEA dataframe.
#' @param filtered_df Filtered GSEA dataframe.
#' @param direction Either "positive" or "negative".
#' @param n_top Number of top pathways.
#' @param celltype_order Optional order for columns.
#' @param celltype_col Cell type column.
#'
#' @return Matrix for heatmap.
#' @export
build_gsea_direction_heatmap_matrix <- function(gsea_df,
                                                filtered_df,
                                                direction = c("positive", "negative"),
                                                n_top = 20,
                                                celltype_order = NULL,
                                                celltype_col = "cellclass") {
  direction <- match.arg(direction)

  gsea_clean <- gsea_df |>
    dplyr::mutate(pathway_clean = clean_go_names(.data$pathway))

  if (direction == "positive") {
    top_paths <- filtered_df |>
      dplyr::filter(.data$NES > 0) |>
      dplyr::group_by(.data$pathway, .data$pathway_clean) |>
      dplyr::summarise(
        max_NES = max(.data$NES, na.rm = TRUE),
        n_celltypes = dplyr::n_distinct(.data[[celltype_col]]),
        .groups = "drop"
      ) |>
      dplyr::arrange(desc(.data$n_celltypes), desc(.data$max_NES)) |>
      dplyr::slice_head(n = n_top) |>
      dplyr::pull(.data$pathway)

    heat <- gsea_clean |>
      dplyr::filter(.data$pathway %in% top_paths) |>
      dplyr::mutate(enrichment = ifelse(.data$NES > 0, .data$NES, 0))
  } else {
    top_paths <- filtered_df |>
      dplyr::filter(.data$NES < 0) |>
      dplyr::group_by(.data$pathway, .data$pathway_clean) |>
      dplyr::summarise(
        max_abs_NES = max(abs(.data$NES), na.rm = TRUE),
        n_celltypes = dplyr::n_distinct(.data[[celltype_col]]),
        .groups = "drop"
      ) |>
      dplyr::arrange(desc(.data$n_celltypes), desc(.data$max_abs_NES)) |>
      dplyr::slice_head(n = n_top) |>
      dplyr::pull(.data$pathway)

    heat <- gsea_clean |>
      dplyr::filter(.data$pathway %in% top_paths) |>
      dplyr::mutate(enrichment = ifelse(.data$NES < 0, abs(.data$NES), 0))
  }

  heat_wide <- heat |>
    dplyr::select(.data$pathway_clean, dplyr::all_of(celltype_col), .data$enrichment) |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(celltype_col),
      values_from = .data$enrichment,
      values_fill = 0
    )

  heat_mat <- as.data.frame(heat_wide)
  rownames(heat_mat) <- heat_mat$pathway_clean
  heat_mat <- heat_mat[, -1, drop = FALSE]

  if (!is.null(celltype_order)) {
    keep_cols <- intersect(celltype_order, colnames(heat_mat))
    heat_mat <- heat_mat[, keep_cols, drop = FALSE]
  }

  as.matrix(heat_mat)
}


# ============================================================
# Heatmap plotting
# ============================================================

#' Plot direction-specific GSEA heatmap
#'
#' @param mat Heatmap matrix.
#' @param direction "positive" or "negative".
#' @param title Plot title.
#' @param file Optional output file.
#' @param positive_color Color for positive heatmap.
#' @param negative_color Color for negative heatmap.
#' @param max_break Maximum heatmap scale.
#'
#' @return pheatmap object.
#' @export
plot_gsea_direction_heatmap <- function(mat,
                                        direction = c("positive", "negative"),
                                        title = NULL,
                                        file = NULL,
                                        positive_color = "#B2182B",
                                        negative_color = "#2166AC",
                                        max_break = 2.5) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.", call. = FALSE)
  }

  direction <- match.arg(direction)

  heat_color <- if (direction == "positive") positive_color else negative_color

  if (is.null(title)) {
    title <- if (direction == "positive") {
      "Positive NES-enriched pathways"
    } else {
      "Negative NES-enriched pathways"
    }
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  }

  pheatmap::pheatmap(
    mat,
    color = colorRampPalette(c("white", heat_color))(100),
    breaks = seq(0, max_break, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 9,
    fontsize_col = 10,
    border_color = NA,
    main = title,
    filename = file
  )
}


#' Run full GSEA heatmap workflow
#'
#' @param gsea_df Full GSEA result dataframe.
#' @param output_dir Output directory.
#' @param positive_label Label for positive NES condition.
#' @param negative_label Label for negative NES condition.
#' @param celltype_order Optional celltype order.
#' @param prefix File prefix.
#' @param n_top Number of pathways per heatmap.
#'
#' @return List containing filtered table and heatmap matrices.
#' @export
run_gsea_directional_heatmaps <- function(gsea_df,
                                          output_dir = "results/gsea",
                                          positive_label = "AD+CAA enriched",
                                          negative_label = "Control enriched",
                                          celltype_order = NULL,
                                          prefix = "GSEA_GO_BP",
                                          n_top = 20,
                                          padj_cutoff = 0.05,
                                          nes_cutoff = 1.5,
                                          max_break = 2.5) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  filtered <- filter_gsea_results(
    gsea_df = gsea_df,
    padj_cutoff = padj_cutoff,
    nes_cutoff = nes_cutoff,
    positive_label = positive_label,
    negative_label = negative_label
  )

  positive_mat <- build_gsea_direction_heatmap_matrix(
    gsea_df = gsea_df,
    filtered_df = filtered,
    direction = "positive",
    n_top = n_top,
    celltype_order = celltype_order
  )

  negative_mat <- build_gsea_direction_heatmap_matrix(
    gsea_df = gsea_df,
    filtered_df = filtered,
    direction = "negative",
    n_top = n_top,
    celltype_order = celltype_order
  )

  write.csv(
    positive_mat,
    file.path(output_dir, paste0(prefix, "_positive_enriched_heatmap_matrix.csv"))
  )

  write.csv(
    negative_mat,
    file.path(output_dir, paste0(prefix, "_negative_enriched_heatmap_matrix.csv"))
  )

  filtered_save <- filtered |>
    dplyr::mutate(
      leadingEdge = if ("leadingEdge" %in% colnames(filtered)) {
        sapply(.data$leadingEdge, paste, collapse = ";")
      } else {
        NA_character_
      }
    )

  write.csv(
    filtered_save,
    file.path(output_dir, paste0(prefix, "_filtered.csv")),
    row.names = FALSE
  )

  p_pos <- plot_gsea_direction_heatmap(
    mat = positive_mat,
    direction = "positive",
    title = paste0(positive_label, " pathways"),
    file = file.path(output_dir, paste0(prefix, "_positive_enriched_heatmap.pdf")),
    max_break = max_break
  )

  p_neg <- plot_gsea_direction_heatmap(
    mat = negative_mat,
    direction = "negative",
    title = paste0(negative_label, " pathways"),
    file = file.path(output_dir, paste0(prefix, "_negative_enriched_heatmap.pdf")),
    max_break = max_break
  )

  invisible(list(
    filtered = filtered,
    positive_matrix = positive_mat,
    negative_matrix = negative_mat,
    positive_plot = p_pos,
    negative_plot = p_neg
  ))
}
