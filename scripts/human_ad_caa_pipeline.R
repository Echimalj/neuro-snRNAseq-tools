source("R/check_dependencies.R")
source("R/load_samples.R")
source("R/metadata_utils.R")

check_required_packages(c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tibble"
))

check_required_packages(c(
  "SoupX"
))

sample_sheet <- read.csv("human_sample_sheet.csv")

seurat_list <- load_samples_from_sheet(
  sample_sheet = sample_sheet,
  min_cells = 3,
  min_features = 200,
  use_soupx = TRUE
)

seurat_list <- annotate_samples_from_sheet(
  seurat_list = seurat_list,
  sample_sheet = sample_sheet
)

AD_CAA <- merge_seurat_samples(
  seurat_list = seurat_list,
  project_name = "AD_CAA_2025"
)

summarize_cells_by_group(AD_CAA, group_col = "orig.ident")
summarize_cells_by_group(AD_CAA, group_col = "condition")
