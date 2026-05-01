#' Pseudobulk DESeq2 utilities
#'
#' Helper functions for cell-type-specific pseudobulk differential expression.
#'
#' @keywords internal
NULL

#' Run pseudobulk DESeq2 for one cell class
#'
#' @param seu Seurat object.
#' @param cellclass Cell class to subset.
#' @param cellclass_col Metadata column containing cell class labels.
#' @param donor_col Metadata column containing biological replicate/sample ID.
#' @param group_col Metadata column containing group/genotype.
#' @param group_levels Character vector of two groups; contrast is group_levels[2] vs group_levels[1].
#' @param assay Assay to use.
#' @param layer Count layer.
#' @param min_cells_per_donor Minimum cells required per donor.
#' @param min_gene_total Minimum total counts required per gene.
#' @param verbose Logical.
#'
#' @return List with dds, results, pseudobulk counts, and coldata.
#' @export
run_pb_deseq2 <- function(seu,
                          cellclass,
                          cellclass_col = "cellclass",
                          donor_col = "orig.ident",
                          group_col = "Genotype",
                          group_levels = c("WT", "TKO"),
                          assay = "RNA",
                          layer = "counts",
                          min_cells_per_donor = 20,
                          min_gene_total = 10,
                          verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  required_cols <- c(cellclass_col, donor_col, group_col)
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  cells_keep <- colnames(seu)[seu[[cellclass_col]][, 1] == cellclass]
  sub <- subset(seu, cells = cells_keep)

  if (ncol(sub) == 0) {
    if (verbose) message("Skipping ", cellclass, ": no cells found.")
    return(NULL)
  }

  donor_vec <- sub[[donor_col]][, 1]
  donor_n <- table(donor_vec)
  keep_donors <- names(donor_n)[donor_n >= min_cells_per_donor]

  sub <- subset(
    sub,
    cells = colnames(sub)[donor_vec %in% keep_donors]
  )

  if (ncol(sub) == 0) {
    if (verbose) message("Skipping ", cellclass, ": no donors passed min_cells_per_donor.")
    return(NULL)
  }

  Seurat::DefaultAssay(sub) <- assay

  donor_map <- sub@meta.data |>
    dplyr::select(
      donor = dplyr::all_of(donor_col),
      group = dplyr::all_of(group_col)
    ) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(.data$group), .data$group %in% group_levels)

  if (length(unique(donor_map$group)) < 2) {
    if (verbose) message("Skipping ", cellclass, ": only one group present after filtering.")
    return(NULL)
  }

  ae_formals <- names(formals(Seurat::AggregateExpression))

  ae_args <- list(
    object = sub,
    assays = assay,
    group.by = donor_col,
    return.seurat = FALSE
  )

  if ("layer" %in% ae_formals) ae_args$layer <- layer
  if ("slot" %in% ae_formals) ae_args$slot <- layer
  if ("normalization.method" %in% ae_formals) ae_args$normalization.method <- "none"
  if ("scale.factor" %in% ae_formals) ae_args$scale.factor <- 1

  pb <- do.call(Seurat::AggregateExpression, ae_args)[[assay]]

  if (verbose) {
    mx <- suppressWarnings(max(pb))
    message(cellclass, " pseudobulk max value: ", mx)

    if (is.finite(mx) && mx < 5) {
      message(
        "Warning: pseudobulk values look small and may represent averages. ",
        "Consider using a manual count-sum fallback."
      )
    }
  }

  donor_map$donor_sanitized <- gsub("_", "-", donor_map$donor)

  coldata <- donor_map |>
    dplyr::filter(.data$donor_sanitized %in% colnames(pb)) |>
    as.data.frame()

  rownames(coldata) <- coldata$donor_sanitized
  coldata$donor <- NULL
  coldata$donor_sanitized <- NULL

  pb <- pb[, rownames(coldata), drop = FALSE]

  count_mat <- as.matrix(pb)
  count_mat[is.na(count_mat)] <- 0
  storage.mode(count_mat) <- "numeric"
  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"

  keep_genes <- rowSums(count_mat) >= min_gene_total
  count_mat <- count_mat[keep_genes, , drop = FALSE]

  if (nrow(count_mat) == 0 || sum(count_mat) == 0) {
    if (verbose) message("Skipping ", cellclass, ": all-zero after gene filtering.")
    return(NULL)
  }

  coldata$group <- factor(coldata$group, levels = group_levels)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ group
  )

  dds <- DESeq2::DESeq(dds)

  res_default <- DESeq2::results(
  dds,
  contrast = c("group", group_levels[2], group_levels[1])
)

res_noIF <- DESeq2::results(
  dds,
  contrast = c("group", group_levels[2], group_levels[1]),
  independentFiltering = FALSE
)

res <- as.data.frame(res_default) |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(
    padj_noIF = res_noIF$padj,
    cellclass = cellclass,
    padj_safe = ifelse(is.na(padj), 1, padj)
  )
}

#' Run pseudobulk DESeq2 across multiple cell classes
#'
#' @param seu Seurat object.
#' @param cellclasses Character vector of cell classes.
#' @param ... Additional arguments passed to run_pb_deseq2().
#'
#' @return Named list of pseudobulk results.
#' @export
run_pb_deseq2_by_cellclass <- function(seu,
                                       cellclasses,
                                       ...) {
  results <- lapply(cellclasses, function(ct) {
    message("Running pseudobulk DESeq2 for: ", ct)
    run_pb_deseq2(seu = seu, cellclass = ct, ...)
  })

  names(results) <- cellclasses
  results
}

#' Extract significant upregulated genes
#'
#' @param res_df DESeq2 result data frame.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param lfc_cutoff Log2FC cutoff.
#' @param padj_col Column to use for adjusted p-values.
#'
#' @return Character vector of genes.
#' @export
get_sig_up <- function(res_df,
                       padj_cutoff = 0.05,
                       lfc_cutoff = 0,
                       padj_col = "padj") {
  res_df |>
    dplyr::filter(
      !is.na(.data[[padj_col]]),
      .data[[padj_col]] < padj_cutoff,
      .data$log2FoldChange > lfc_cutoff
    ) |>
    dplyr::pull(.data$gene) |>
    unique()
}

#' Extract significant downregulated genes
#'
#' @param res_df DESeq2 result data frame.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param lfc_cutoff Log2FC cutoff.
#' @param padj_col Column to use for adjusted p-values.
#'
#' @return Character vector of genes.
#' @export
get_sig_down <- function(res_df,
                         padj_cutoff = 0.05,
                         lfc_cutoff = 0,
                         padj_col = "padj") {
  res_df |>
    dplyr::filter(
      !is.na(.data[[padj_col]]),
      .data[[padj_col]] < padj_cutoff,
      .data$log2FoldChange < -lfc_cutoff
    ) |>
    dplyr::pull(.data$gene) |>
    unique()
}
