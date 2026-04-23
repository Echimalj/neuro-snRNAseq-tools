source("R/check_dependencies.R")
source("R/load_samples.R")

check_required_packages(c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tibble"
))

check_required_packages(c(
  "SoupX"
))

sample_sheet <- read.csv("inst/extdata/human_sample_sheet.csv")

seurat_list <- load_samples_from_sheet(
  sample_sheet = sample_sheet,
  min_cells = 3,
  min_features = 200,
  use_soupx = TRUE
)
