#' volcano3D utilities for three-group single-cell comparisons
#'
#' Extended workflow adapted from the volcano3D package:
#' https://katrionagoldmann.github.io/volcano3D/index.html
#'
#' This implementation supports Seurat workflows, three-way genotype comparisons,
#' pairwise FindMarkers testing, gene-wise ANOVA, q-value correction, and export
#' of volcano3D significance categories.
#'
#' @keywords internal
NULL

#' Add combined cell type and group labels
#'
#' @param seu Seurat object.
#' @param group_col Metadata column containing condition/genotype.
#' @param output_col Name of combined metadata column.
#' @param sep Separator between identity and group.
#'
#' @return Seurat object.
#' @export
add_volcano3d_identity <- function(seu,
                                   group_col,
                                   output_col = "celltype.group",
                                   sep = "-",
                                   overwrite = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  if (!group_col %in% colnames(seu@meta.data)) {
    stop("group_col not found in metadata: ", group_col, call. = FALSE)
  }

  if (output_col %in% colnames(seu@meta.data) && !overwrite) {
    message(
      "Column '", output_col, "' already exists. ",
      "Keeping existing column. Use overwrite = TRUE to replace it."
    )
    return(seu)
  }

  seu[[output_col]] <- paste(
    as.character(Seurat::Idents(seu)),
    seu[[group_col]][, 1],
    sep = sep
  )

  message("Created column: ", output_col)

  return(seu)
}


#' Build volcano3D sample metadata for one cell type
#'
#' @param seu Seurat object.
#' @param celltype_label Cell type/cluster label.
#' @param groups Character vector of three groups.
#' @param combined_col Metadata column containing celltype-group labels.
#' @param group_col Metadata column containing group labels.
#' @param sample_col Metadata column containing sample IDs.
#' @param sep Separator used in combined labels.
#'
#' @return Filtered metadata data frame.
#' @export
build_volcano3d_sample_data <- function(seu,
                                        celltype_label,
                                        groups,
                                        combined_col = "celltype.group",
                                        group_col = "Genotype",
                                        sample_col = "orig.ident",
                                        sep = "-") {
  if (length(groups) != 3) {
    stop("groups must contain exactly three groups.", call. = FALSE)
  }

  required_cols <- c(combined_col, group_col, sample_col)
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  target_labels <- paste(celltype_label, groups, sep = sep)

  sample_data <- seu@meta.data
  sample_data$ID <- rownames(sample_data)

  sample_data <- sample_data[
    sample_data[[combined_col]] %in% target_labels,
  ]

  if (nrow(sample_data) == 0) {
    stop("No cells found for selected celltype_label and groups.", call. = FALSE)
  }

  sample_data[[group_col]] <- factor(sample_data[[group_col]], levels = groups)

  return(sample_data)
}


#' Extract expression matrix for volcano3D
#'
#' @param seu Seurat object.
#' @param sample_data Metadata from build_volcano3d_sample_data().
#' @param assay Assay to use.
#' @param layer Layer to extract.
#'
#' @return Gene x cell expression matrix.
#' @export
extract_volcano3d_expression <- function(seu,
                                         sample_data,
                                         assay = "RNA",
                                         layer = "scale.data") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  expr <- Seurat::GetAssayData(
    object = seu,
    assay = assay,
    layer = layer
  )

  expr <- as.matrix(expr)
  expr <- expr[, sample_data$ID, drop = FALSE]

  return(expr)
}


#' Run pairwise Seurat FindMarkers contrasts for volcano3D
#'
#' @param seu Seurat object.
#' @param celltype_label Cell type/cluster label.
#' @param groups Character vector of three groups.
#' @param combined_col Metadata column with celltype-group labels.
#' @param assay Assay to use.
#' @param sep Separator used in combined labels.
#'
#' @return Data frame of pairwise statistics.
#' @export
run_volcano3d_pairwise_markers <- function(seu,
                                           celltype_label,
                                           groups,
                                           combined_col = "celltype.group",
                                           assay = "RNA",
                                           sep = "-") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  if (length(groups) != 3) {
    stop("groups must contain exactly three groups.", call. = FALSE)
  }

  Seurat::DefaultAssay(seu) <- assay
  Seurat::Idents(seu) <- seu[[combined_col]][, 1]

  group_labels <- paste(celltype_label, groups, sep = sep)

  contrasts <- list(
    AB = c(group_labels[1], group_labels[2]),
    BC = c(group_labels[2], group_labels[3]),
    CA = c(group_labels[3], group_labels[1])
  )

  pvals <- lapply(names(contrasts), function(cn) {
    vars <- contrasts[[cn]]

    out <- Seurat::FindMarkers(
      object = seu,
      ident.1 = vars[1],
      ident.2 = vars[2],
      assay = assay,
      min.pct = 0,
      logfc.threshold = 0,
      verbose = FALSE
    )

    out <- out[, c("p_val", "p_val_adj", "avg_log2FC")]
    colnames(out) <- paste0(cn, "_", c("pvalue", "padj", "logFC"))

    out
  })

  pvals_df <- do.call(cbind, pvals)

  return(pvals_df)
}


#' Run gene-wise ANOVA for volcano3D
#'
#' @param expr Gene x cell expression matrix.
#' @param sample_data Metadata matching expression columns.
#' @param group_col Metadata column containing genotype/condition.
#' @param sample_col Optional metadata column for sample covariate.
#' @param use_sample_covariate Logical.
#'
#' @return Data frame with overall p-value and q-value.
#' @export
run_volcano3d_anova <- function(expr,
                                sample_data,
                                group_col = "Genotype",
                                sample_col = "orig.ident",
                                use_sample_covariate = TRUE) {
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package 'car' is required.", call. = FALSE)
  }
  if (!requireNamespace("qvalue", quietly = TRUE)) {
    stop("Package 'qvalue' is required.", call. = FALSE)
  }

  if (!all(colnames(expr) == sample_data$ID)) {
    stop("Expression matrix columns must match sample_data$ID.", call. = FALSE)
  }

  aov_stats <- apply(expr, 1, function(x) {
    df <- data.frame(
      Gene = as.numeric(x),
      Group = sample_data[[group_col]],
      Sample = sample_data[[sample_col]]
    )

    if (stats::var(df$Gene, na.rm = TRUE) == 0) {
      return(1)
    }

    fit <- try(
      if (use_sample_covariate) {
        stats::lm(Gene ~ Group + Sample, data = df)
      } else {
        stats::lm(Gene ~ Group, data = df)
      },
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      return(1)
    }

    stat_aov <- try(car::Anova(fit), silent = TRUE)

    if (inherits(stat_aov, "try-error")) {
      return(1)
    }

    stat_aov[1, "Pr(>F)"]
  })

  qvals <- qvalue::qvalue(aov_stats)$qvalue

  out <- data.frame(
    overall_pvalue = aov_stats,
    overall_padj = qvals
  )

  rownames(out) <- rownames(expr)

  return(out)
}


#' Build volcano3D statistics table
#'
#' @param pairwise_stats Output from run_volcano3d_pairwise_markers().
#' @param anova_stats Output from run_volcano3d_anova().
#'
#' @return Combined statistics table.
#' @export
build_volcano3d_stats_table <- function(pairwise_stats,
                                        anova_stats) {
  shared_genes <- intersect(rownames(pairwise_stats), rownames(anova_stats))

  combined <- cbind(
    pairwise_stats[shared_genes, , drop = FALSE],
    anova_stats[shared_genes, , drop = FALSE]
  )

  return(combined)
}


#' Run volcano3D polar coordinate classification
#'
#' @param expr Gene x cell expression matrix.
#' @param sample_data Metadata matching expression columns.
#' @param stats_table Combined pairwise and ANOVA statistics table.
#' @param group_col Metadata column containing group labels.
#' @param pcutoff Adjusted p-value cutoff.
#' @param scheme Color scheme for volcano3D categories.
#' @param labs Labels for volcano3D categories.
#' @param pval_cols Ordered p-value columns. Default order: overall, AB, BC, CA.
#' @param padj_cols Ordered adjusted p-value columns. Default order: overall, AB, BC, CA.
#'
#' @return volcano3D polar object.
#' @export
run_volcano3d_polar_coords <- function(expr,
                                       sample_data,
                                       stats_table,
                                       group_col = "Genotype",
                                       pcutoff = 0.05,
                                       scheme = c(
                                         "grey60", "red", "gold2", "green3",
                                         "cyan", "blue", "purple"
                                       ),
                                       labs = c(
                                         "ns", "GroupC", "GroupC+GroupA",
                                         "GroupA", "GroupB+GroupA",
                                         "GroupB", "GroupB+GroupC"
                                       ),
                                       pval_cols = c(
                                         "overall_pvalue",
                                         "AB_pvalue",
                                         "BC_pvalue",
                                         "CA_pvalue"
                                       ),
                                       padj_cols = c(
                                         "overall_padj",
                                         "AB_padj",
                                         "BC_padj",
                                         "CA_padj"
                                       )) {
  if (!requireNamespace("volcano3D", quietly = TRUE)) {
    stop("Package 'volcano3D' is required.", call. = FALSE)
  }

  if (!group_col %in% colnames(sample_data)) {
    stop("group_col not found in sample_data: ", group_col, call. = FALSE)
  }

  missing_cols <- setdiff(c(pval_cols, padj_cols), colnames(stats_table))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required volcano3D statistics columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!all(rownames(stats_table) %in% rownames(expr))) {
    stop(
      "Not all genes in stats_table are present in expr rownames.",
      call. = FALSE
    )
  }

  expr <- expr[rownames(stats_table), , drop = FALSE]

  if (!is.null(sample_data$ID)) {
    if (!all(colnames(expr) == sample_data$ID)) {
      stop(
        "Expression matrix columns must match sample_data$ID in the same order.",
        call. = FALSE
      )
    }
  }

  p_matrix <- stats_table[, pval_cols, drop = FALSE]
  padj_matrix <- stats_table[, padj_cols, drop = FALSE]

  polar <- volcano3D::polar_coords(
    outcome = sample_data[[group_col]],
    data = t(expr),
    pvals = as.matrix(p_matrix),
    padj = as.matrix(padj_matrix),
    pcutoff = pcutoff,
    scheme = scheme,
    labs = labs
  )

  return(polar)
}

#' Save volcano3D significance category gene lists
#'
#' @param polar volcano3D polar object.
#' @param categories Character vector of volcano3D categories.
#' @param file Output CSV file.
#'
#' @return Invisibly returns file path.
#' @export
save_volcano3d_gene_categories <- function(polar,
                                           categories,
                                           file) {
  if (!requireNamespace("volcano3D", quietly = TRUE)) {
    stop("Package 'volcano3D' is required.", call. = FALSE)
  }

  colnames_list <- list()

  for (cat in categories) {
    tryCatch({
      subset <- volcano3D::significance_subset(
        polar,
        cat,
        "data"
      )

      colnames_list[[cat]] <- colnames(subset)
    }, error = function(e) {
      message("Skipping ", cat, ": ", e$message)
    })
  }

  if (length(colnames_list) == 0) {
    warning("No categories returned genes. CSV not created.")
    return(invisible(NULL))
  }

  max_len <- max(vapply(colnames_list, length, numeric(1)))

  out <- as.data.frame(lapply(colnames_list, function(x) {
    length(x) <- max_len
    x
  }))

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, file = file, row.names = FALSE)

  invisible(file)
}


#' Run full volcano3D workflow for one cell type
#'
#' @param seu Seurat object.
#' @param celltype_label Cell type/cluster label.
#' @param groups Character vector of three groups in desired order.
#' @param group_col Metadata group column.
#' @param sample_col Metadata sample column.
#' @param combined_col Combined celltype-group metadata column.
#' @param assay Assay to use.
#' @param expr_layer Expression layer to extract.
#' @param sep Separator used in combined labels.
#' @param use_sample_covariate Logical.
#' @param pcutoff Adjusted p-value cutoff.
#' @param labs volcano3D category labels.
#' @param scheme volcano3D color scheme.
#'
#' @return List containing sample data, expression matrix, stats table, and polar object.
#' @export
run_volcano3d_workflow <- function(seu,
                                   celltype_label,
                                   groups,
                                   group_col = "Genotype",
                                   sample_col = "orig.ident",
                                   combined_col = "celltype.group",
                                   assay = "RNA",
                                   expr_layer = "scale.data",
                                   sep = "-",
                                   use_sample_covariate = TRUE,
                                   pcutoff = 0.05,
                                   labs = c(
                                     "ns", "GroupC", "GroupC+GroupA",
                                     "GroupA", "GroupB+GroupA",
                                     "GroupB", "GroupB+GroupC"
                                   ),
                                   scheme = c(
                                     "grey60", "red", "gold2", "green3",
                                     "cyan", "blue", "purple"
                                   )) {
  if (!combined_col %in% colnames(seu@meta.data)) {
    seu <- add_volcano3d_identity(
      seu = seu,
      group_col = group_col,
      output_col = combined_col,
      sep = sep
    )
  }

  sample_data <- build_volcano3d_sample_data(
    seu = seu,
    celltype_label = celltype_label,
    groups = groups,
    combined_col = combined_col,
    group_col = group_col,
    sample_col = sample_col,
    sep = sep
  )

  expr <- extract_volcano3d_expression(
    seu = seu,
    sample_data = sample_data,
    assay = assay,
    layer = expr_layer
  )

  pairwise_stats <- run_volcano3d_pairwise_markers(
    seu = seu,
    celltype_label = celltype_label,
    groups = groups,
    combined_col = combined_col,
    assay = assay,
    sep = sep
  )

  anova_stats <- run_volcano3d_anova(
    expr = expr,
    sample_data = sample_data,
    group_col = group_col,
    sample_col = sample_col,
    use_sample_covariate = use_sample_covariate
  )

  stats_table <- build_volcano3d_stats_table(
    pairwise_stats = pairwise_stats,
    anova_stats = anova_stats
  )

  expr <- expr[rownames(stats_table), , drop = FALSE]

  polar <- run_volcano3d_polar_coords(
    expr = expr,
    sample_data = sample_data,
    stats_table = stats_table,
    group_col = group_col,
    pcutoff = pcutoff,
    scheme = scheme,
    labs = labs
  )

  return(list(
    sample_data = sample_data,
    expression = expr,
    pairwise_stats = pairwise_stats,
    anova_stats = anova_stats,
    stats_table = stats_table,
    polar = polar
  ))
}
