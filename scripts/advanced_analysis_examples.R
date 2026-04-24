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

#Differential Expression Workflow

AD_CAA <- add_celltype_group_identity(
  seu = AD_CAA,
  group_col = "FDX",
  output_col = "celltype.FDX",
  use_active_ident = TRUE
)
#Define comparison
group_1 <- "AD+CAA"
group_2 <- "Control"

#Loop over all astrocyte populations: 
astro_clusters <- c(
  "Astrocytes1",
  "Astrocytes2",
  "Astrocytes3",
  "Astrocytes4",
  "Astrocytes5"
)

astro_results <- lapply(astro_clusters, function(ct) {
  message("Running DE for ", ct)

  run_sc_pseudobulk_de_workflow(
    seu = AD_CAA,
    celltype_label = ct,
    group_col = "FDX",
    group_1 = group_1,
    group_2 = group_2,
    sample_col = "orig.ident",
    combined_col = "celltype.FDX",
    assay = "RNA",
    output_dir = "results/de/astrocytes",
    prefix = paste0(ct, "_ADCAA_vs_Control")
  )
})

names(astro_results) <- astro_clusters

astro_results[["Astrocytes1"]]$summary

#common → robust DEGs (best for biology / figures)
#only_sc → sensitive, cell-level changes
#only_bulk → consistent across samples (stronger statistical power)

#Run for multiple cell types
celltypes_to_test <- c(
  "Astrocytes1", "Astrocytes2", "Astrocytes3",
  "Microglia1", "Microglia2", "Microglia3",
  "Oligodendrocytes1", "Oligodendrocytes2",
  "OPC1", "OPC2",
  "Endothelial", "Vascular"
)

de_results <- lapply(celltypes_to_test, function(ct) {
  run_sc_pseudobulk_de_workflow(
    seu = AD_CAA,
    celltype_label = ct,
    group_col = "FDX",
    group_1 = "AD+CAA",
    group_2 = "Control",
    sample_col = "orig.ident",
    combined_col = "celltype.FDX",
    assay = "RNA",
    output_dir = "results/de/all_celltypes",
    prefix = paste0(ct, "_ADCAA_vs_Control")
  )
})

names(de_results) <- celltypes_to_test




