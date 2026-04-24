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

#Module Scores Utilities (DAA, DAM, DAO)
# Disease-associated astrocyte module
daa_genes <- c(
  "AMIGO2", "FBLN5", "FKBP5", "GBP2", "PSMB8", "SERPING1",
  "SRGN", "CD44", "HSPB1", "LTBP1", "TNC"
)

astro <- subset(
  AD_CAA,
  idents = c("Astrocytes1", "Astrocytes2", "Astrocytes3", "Astrocytes4", "Astrocytes5")
)

astro <- add_gene_module_score(
  seu = astro,
  gene_set = daa_genes,
  module_name = "DAA_Score",
  assay = "RNA"
)

daa_pct <- percent_cells_expressing_module(
  seu = astro,
  gene_set = daa_genes,
  threshold = 5,
  assay = "RNA",
  celltype_col = "cellclass",
  group_col = "FDX"
)


