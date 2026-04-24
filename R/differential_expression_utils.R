#' Differential expression utilities
#'
#' Helper functions for single-cell DE, pseudo-bulk DE, and comparison between
#' single-cell and pseudo-bulk results.
#'
#' @keywords internal
NULL

#' Create combined cell type and group identity
#'
#' @param seu Seurat object.
#' @param celltype_col Metadata column or active identity source.
#' @param group_col Metadata column containing group labels.
#' @param output_col Name of new metadata column.
#' @param use_active_ident Logical. If TRUE, uses active identities as cell type.
#'
#' @return Seurat object with combined identity column.
#' @export
add_celltype_group_identity <- function(seu,
                                        celltype_col = NULL,
                                        group_col,
                                        output_col = "celltype.group",
                                        use_active_ident = TRUE) {
  if (!group_col %in% colnames(seu@meta.data)) {
    stop("group_col not found in metadata: ", group_col, call. = FALSE)
  }

  celltypes <- if (use_active_ident) {
    as.character(Seurat::Idents(seu))
  } else {
    if (is.null(celltype_col) || !celltype_col %in% colnames(seu@meta.data)) {
      stop("celltype_col not found in metadata.", call. = FALSE)
    }
    seu[[celltype_col]][, 1]
  }

  groups <- seu[[group_col]][, 1]

  seu[[output_col]] <- paste(celltypes, groups, sep = "_")

  return(seu)
}

#' Run single-cell DE for one cell type
#'
#' @param seu Seurat object.
#' @param celltype_label Cell type label.
#' @param group_1 First group.
#' @param group_2 Second group.
#' @param combined_col Metadata column with celltype_group labels.
#' @param assay Assay to use.
#' @param sep Separator used in combined identities.
#'
#' @return Single-cell DE table.
#' @export
run_singlecell_de_for_celltype <- function(seu,
                                           celltype_label,
                                           group_1,
                                           group_2,
                                           combined_col = "celltype.group",
                                           assay = "RNA",
                                           sep = "_",
                                           logfc_threshold = 0.25,
                                           min_pct = 0.1,
                                           test_use = "wilcox") {
  if (!combined_col %in% colnames(seu@meta.data)) {
    stop("combined_col not found in metadata: ", combined_col, call. = FALSE)
  }

  Seurat::DefaultAssay(seu) <- assay
  Seurat::Idents(seu) <- seu[[combined_col]][, 1]

  ident_1 <- paste(celltype_label, group_1, sep = sep)
  ident_2 <- paste(celltype_label, group_2, sep = sep)

  available_idents <- levels(Seurat::Idents(seu))

  if (!ident_1 %in% available_idents || !ident_2 %in% available_idents) {
    stop(
      "One or both identities not found: ",
      ident_1, ", ", ident_2,
      call. = FALSE
    )
  }

  sc_deg <- Seurat::FindMarkers(
    object = seu,
    ident.1 = ident_1,
    ident.2 = ident_2,
    assay = assay,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    test.use = test_use,
    verbose = FALSE
  )

  sc_deg$gene <- rownames(sc_deg)
  rownames(sc_deg) <- NULL

  names(sc_deg)[names(sc_deg) != "gene"] <- paste0(
    names(sc_deg)[names(sc_deg) != "gene"],
    ".sc"
  )

  return(sc_deg)
}

#' Run pseudo-bulk DE for one cell type
#'
#' @param seu Seurat object.
#' @param celltype_label Cell type label.
#' @param group_1 First group.
#' @param group_2 Second group.
#' @param group_by Columns used for pseudo-bulk aggregation.
#' @param combined_col Metadata column with celltype_group labels.
#' @param assay Assay to use.
#' @param sep Separator used in combined identities.
#'
#' @return Pseudo-bulk DE table.
#' @export
run_pseudobulk_de_for_celltype <- function(seu,
                                           celltype_label,
                                           group_1,
                                           group_2,
                                           group_by = c("Genotype", "orig.ident", "celltype.group"),
                                           combined_col = "celltype.group",
                                           assay = "RNA",
                                           sep = "_",
                                           test_use = "DESeq2") {
  if (!combined_col %in% colnames(seu@meta.data)) {
    stop("combined_col not found in metadata: ", combined_col, call. = FALSE)
  }

  missing_group_by <- setdiff(group_by, colnames(seu@meta.data))
  if (length(missing_group_by) > 0) {
    stop(
      "Missing group_by columns: ",
      paste(missing_group_by, collapse = ", "),
      call. = FALSE
    )
  }

  pseudo <- Seurat::AggregateExpression(
    seu,
    assays = assay,
    return.seurat = TRUE,
    group.by = group_by
  )

  Seurat::DefaultAssay(pseudo) <- assay
  Seurat::Idents(pseudo) <- pseudo[[combined_col]][, 1]

  # DESeq2 requires integer-like counts.
  if ("counts" %in% names(pseudo[[assay]]@layers)) {
    pseudo[[assay]]@layers$counts <- round(pseudo[[assay]]@layers$counts)
  }

  ident_1 <- paste(celltype_label, group_1, sep = sep)
  ident_2 <- paste(celltype_label, group_2, sep = sep)

  # AggregateExpression sometimes converts underscores to dashes.
  if (!ident_1 %in% levels(Seurat::Idents(pseudo))) {
    ident_1 <- gsub("_", "-", ident_1)
  }
  if (!ident_2 %in% levels(Seurat::Idents(pseudo))) {
    ident_2 <- gsub("_", "-", ident_2)
  }

  bulk_deg <- Seurat::FindMarkers(
    object = pseudo,
    ident.1 = ident_1,
    ident.2 = ident_2,
    assay = assay,
    test.use = test_use,
    verbose = FALSE
  )

  bulk_deg$gene <- rownames(bulk_deg)
  rownames(bulk_deg) <- NULL

  names(bulk_deg)[names(bulk_deg) != "gene"] <- paste0(
    names(bulk_deg)[names(bulk_deg) != "gene"],
    ".bulk"
  )

  return(bulk_deg)
}

#' Compare single-cell and pseudo-bulk DE results
#'
#' @param sc_deg Single-cell DE table.
#' @param bulk_deg Pseudo-bulk DE table.
#' @param p_col_sc Single-cell p-value column.
#' @param p_col_bulk Pseudo-bulk p-value column.
#' @param p_cutoff P-value cutoff.
#'
#' @return List containing merged table and DEG categories.
#' @export
compare_sc_and_pseudobulk_de <- function(sc_deg,
                                         bulk_deg,
                                         p_col_sc = "p_val.sc",
                                         p_col_bulk = "p_val.bulk",
                                         p_cutoff = 0.05) {
  merged <- merge(sc_deg, bulk_deg, by = "gene")
  merged <- merged[order(merged[[p_col_bulk]]), ]

  common <- merged$gene[
    merged[[p_col_bulk]] < p_cutoff & merged[[p_col_sc]] < p_cutoff
  ]

  only_sc <- merged$gene[
    merged[[p_col_bulk]] >= p_cutoff & merged[[p_col_sc]] < p_cutoff
  ]

  only_bulk <- merged$gene[
    merged[[p_col_bulk]] < p_cutoff & merged[[p_col_sc]] >= p_cutoff
  ]

  list(
    merged = merged,
    common = common,
    only_sc = only_sc,
    only_bulk = only_bulk,
    summary = data.frame(
      category = c("common", "only_single_cell", "only_pseudobulk"),
      n_genes = c(length(common), length(only_sc), length(only_bulk))
    )
  )
}

#' Run paired single-cell and pseudo-bulk DE workflow
#'
#' @param seu Seurat object.
#' @param celltype_label Cell type label to test.
#' @param group_col Group metadata column.
#' @param group_1 First group.
#' @param group_2 Second group.
#' @param sample_col Sample ID column.
#' @param combined_col Combined celltype-group column.
#' @param assay Assay to use.
#' @param output_dir Output directory.
#' @param prefix File prefix.
#'
#' @return List with single-cell, pseudo-bulk, merged, and summary results.
#' @export
run_sc_pseudobulk_de_workflow <- function(seu,
                                          celltype_label,
                                          group_col,
                                          group_1,
                                          group_2,
                                          sample_col = "orig.ident",
                                          combined_col = "celltype.group",
                                          assay = "RNA",
                                          output_dir = NULL,
                                          prefix = NULL) {
  if (!combined_col %in% colnames(seu@meta.data)) {
    seu <- add_celltype_group_identity(
      seu = seu,
      group_col = group_col,
      output_col = combined_col,
      use_active_ident = TRUE
    )
  }

  group_by <- c(group_col, sample_col, combined_col)

  sc_deg <- run_singlecell_de_for_celltype(
    seu = seu,
    celltype_label = celltype_label,
    group_1 = group_1,
    group_2 = group_2,
    combined_col = combined_col,
    assay = assay
  )

  bulk_deg <- run_pseudobulk_de_for_celltype(
    seu = seu,
    celltype_label = celltype_label,
    group_1 = group_1,
    group_2 = group_2,
    group_by = group_by,
    combined_col = combined_col,
    assay = assay
  )

  comparison <- compare_sc_and_pseudobulk_de(
    sc_deg = sc_deg,
    bulk_deg = bulk_deg
  )

  result <- list(
    single_cell = sc_deg,
    pseudobulk = bulk_deg,
    merged = comparison$merged,
    common = comparison$common,
    only_sc = comparison$only_sc,
    only_bulk = comparison$only_bulk,
    summary = comparison$summary
  )

  if (!is.null(output_dir)) {
    if (is.null(prefix)) {
      prefix <- paste(celltype_label, group_1, "vs", group_2, sep = "_")
    }

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    write.table(
      sc_deg,
      file = file.path(output_dir, paste0(prefix, "_singlecell_DEGs.txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    write.table(
      bulk_deg,
      file = file.path(output_dir, paste0(prefix, "_pseudobulk_DEGs.txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    write.table(
      comparison$merged,
      file = file.path(output_dir, paste0(prefix, "_merged_sc_pseudobulk_DEGs.txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    write.table(
      comparison$common,
      file = file.path(output_dir, paste0(prefix, "_common_DEGs.txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }

  return(result)
}
