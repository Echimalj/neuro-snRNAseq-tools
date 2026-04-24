source("R/check_dependencies.R")
source("R/module_score_utils.R")
source("R/differential_expression_utils.R")
source("R/speckle_utils.R")
source("R/subcluster_utils.R") #Advanced Functions
source("R/atlas_similarity_utils.R")

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

#Advanced Subclustering functions
check_required_packages(c(
  "Seurat",
  "sctransform",
  "speckle",
  "ggplot2"
))

#Subset by idents 
astro_idents <- c(
  "Astrocytes1",
  "Astrocytes2",
  "Astrocytes3",
  "Astrocytes4",
  "Astrocytes5"
)

AD_CAA_Astro <- save_subset_by_idents(
  seu = AD_CAA,
  idents = astro_idents,
  file = "checkpoints/AD_CAA_Astro_raw_subset.rds"
)

#Re-run subclustering for finer astrocyte states, and clean UMAP
AD_CAA_Astro <- run_standalone_subcluster_workflow(
  seu_sub = AD_CAA_Astro,
  assay = "RNA",
  resolution = 0.6,
  dims = 1:20,
  npcs = 20
)

save_umap_plot(
  seu = AD_CAA_Astro,
  file = "figures/UMAP_ADCAA_astrocytes.eps",
  width = 6,
  height = 6,
  label = FALSE
)

save_umap_plot(
  seu = AD_CAA_Astro,
  file = "figures/UMAP_ADCAA_astrocytes_split.eps",
  split_by = "FDX",
  width = 14,
  height = 6,
  label = FALSE
#Correlation Matrix to merge similar subclusters
cor_mat <- compute_identity_correlation(
  seu = AD_CAA_Astro,
  assay = "RNA",
  slot = "data"
)

save_correlation_matrix(
  cor_mat,
  "results/AD_CAA_Astro_correlation.txt"
)

plot_identity_correlation(cor_mat))

#Find Markers & Propeller
AD_CAA_Astro_markers <- find_subset_cluster_markers(
  seu_sub = AD_CAA_Astro,
  assay = "RNA",
  only_pos = TRUE,
  output_file = "results/AD_CAA_Astro.clusterDEGs.txt"
)

Prop_Astro <- run_subset_propeller(
  seu_sub = AD_CAA_Astro,
  sample_col = "orig.ident",
  group_col = "FDX",
  output_file = "results/CellTypeProportions_ADCAA_Astro.csv"
)

saveRDS(
  AD_CAA_Astro,
  file = "checkpoints/AD_CAA_Astro_2025.rds"
)

#Cosine Similarity Report
check_required_packages(c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tibble",
  "ggplot2"
))

AD_CAA_Astro <- readRDS("checkpoints/AD_CAA_Astro_2025.rds")

saddick_gene_sets <- list(
  Synapse = gene_list0,
  Oxidative_Stress_AB_Trafficking = gene_list1,
  ECM_Protective = gene_list2,
  Inflammatory = gene_list3,
  Synapse_Glutamate = gene_list4,
  ECM_Actin_Protective = gene_list5,
  Glutamate_Metallothioneins = gene_list6,
  Apoptosis_DNA_Damage = gene_list7,
  Synapse_2 = gene_list8
)

similarity_results <- run_atlas_similarity_list(
  seurat_object = AD_CAA_Astro,
  gene_sets = saddick_gene_sets,
  cluster_column = "seurat_clusters",
  assay = "RNA",
  layer = "data",
  n_null = 1000
)

save_atlas_similarity_results(
  similarity_results = similarity_results,
  output_dir = "results/atlas_similarity",
  prefix = "AD_CAA_Astro_Saddick"
)

#Single gene set:
res_synapse <- cluster_marker_cosine_similarity(
  seurat_object = AD_CAA_Astro,
  gene_list = gene_list0,
  cluster_column = "seurat_clusters",
  assay = "RNA",
  layer = "data",
  n_null = 1000
)

plot_cosine_similarity(res_synapse$similarity_table)


