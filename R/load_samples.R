check_required_packages(c(
  "Seurat", "SeuratObject", "dplyr", "tibble",
  "SoupX", "DoubletFinder", "sctransform", "harmony"
))


## load_sample_basic() expects a folder readable by Read10X()
## load_sample_with_soupx() expects a Cell Ranger output structure compatible with SoupX::load10X()


load_sample_with_soupx <- function(sample_path,
                                   sample_id,
                                   min_cells = 3,
                                   min_features = 200) {
  soup <- SoupX::load10X(sample_path)
  soup <- SoupX::autoEstCont(soup)
  adjusted <- SoupX::adjustCounts(soup)

  seu <- Seurat::CreateSeuratObject(
    counts = adjusted,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_features
  )

  rm(soup, adjusted)
  gc()
  return(seu)
}

load_sample_basic <- function(sample_path,
                              sample_id,
                              min_cells = 3,
                              min_features = 200) {
  counts <- Seurat::Read10X(data.dir = sample_path)

  seu <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_features
  )

  rm(counts)
  gc()
  return(seu)
}

load_sample <- function(sample_path,
                        sample_id,
                        min_cells = 3,
                        min_features = 200,
                        use_soupx = TRUE) {
  if (use_soupx) {
    load_sample_with_soupx(
      sample_path = sample_path,
      sample_id = sample_id,
      min_cells = min_cells,
      min_features = min_features
    )
  } else {
    load_sample_basic(
      sample_path = sample_path,
      sample_id = sample_id,
      min_cells = min_cells,
      min_features = min_features
    )
  }
}

load_samples_from_sheet <- function(sample_sheet,
                                    min_cells = 3,
                                    min_features = 200,
                                    use_soupx = NULL) {
  seurat_list <- lapply(seq_len(nrow(sample_sheet)), function(i) {
    sample_use_soupx <- if (!is.null(use_soupx)) {
      use_soupx
    } else if ("use_soupx" %in% colnames(sample_sheet)) {
      as.logical(sample_sheet$use_soupx[i])
    } else {
      TRUE
    }

    load_sample(
      sample_path = sample_sheet$sample_path[i],
      sample_id = sample_sheet$sample_id[i],
      min_cells = min_cells,
      min_features = min_features,
      use_soupx = sample_use_soupx
    )
  })

  names(seurat_list) <- sample_sheet$sample_id
  return(seurat_list)
}
