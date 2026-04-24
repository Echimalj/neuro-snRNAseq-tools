source("R/check_dependencies.R")
source("R/module_score_utils.R")
source("R/differential_expression_utils.R")
source("R/speckle_utils.R")

check_required_packages(c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "ggpubr"
))

AD_CAA <- readRDS("checkpoints/AD_CAA_final_processed.rds")


save_propeller_results(
  propeller_results = prop_sets,
  output_dir = "results",
  prefix = "CellTypeProportions_ADCAAvsCTRL"
)
