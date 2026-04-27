source("R/check_dependencies.R")
source("R/module_score_utils.R")
source("R/differential_expression_utils.R")
source("R/speckle_utils.R")
source("R/subcluster_utils.R") #Advanced Functions
source("R/atlas_similarity_utils.R")
source("R/annotation_utils.R") #Advanced Functions
source("R/pseudobulk_utils.R")
source("R/overlap_utils.R")


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
#NOTE: This can also subset by Genotype or by FDX:
AD_CAA_AD <- save_subset_by_idents(
  seu = AD_CAA,
  idents = NULL,
  group_col = "FDX",
  groups = "AD+CAA",
  file = "checkpoints/AD_CAA_ADCAA_only.rds"
)

AD_CAA_Control <- save_subset_by_idents(
  seu = AD_CAA,
  idents = NULL,
  group_col = "FDX",
  groups = "Control",
  file = "checkpoints/AD_CAA_Control_only.rds"
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
#Propeller for more than 2 will generate an ANOVA, for Post hoc:
#For AD_CAA, this is only two groups, so regular run_propeller()
#is usually enough. But if you later have Control, AD, AD+CAA, or multiple pathology groups
posthoc_adcaa <- run_propeller_posthoc(
  seu = AD_CAA,
  sample_col = "orig.ident",
  group_col = "FDX",
  group_levels = c("Control", "AD", "AD+CAA")
)
filter_propeller_posthoc_cluster(posthoc_exneurons, "ExNeuron6")


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


#Collapse cellclass:
AD_CAA$cellclass <- as.character(Idents(AD_CAA))

AD_CAA <- add_collapsed_cellclass(
  seu = AD_CAA,
  input_col = "cellclass",
  output_col = "collapsed_cellclass",
  keep = character(),
  set_idents = FALSE
)

table(AD_CAA$collapsed_cellclass)
table(AD_CAA$collapsed_cellclass, AD_CAA$FDX)

#If want to keep any identity separate i.e. merge all ExNeuron except ExNeuron3
AD_CAA <- add_collapsed_cellclass(
  seu = AD_CAA,
  input_col = "cellclass",
  output_col = "collapsed_cellclass",
  keep = c("ExNeuron3"),
  set_idents = FALSE
)
#If you want to use broad classes as active identities
AD_CAA <- add_collapsed_cellclass(
  seu = AD_CAA,
  input_col = "cellclass",
  output_col = "collapsed_cellclass",
  keep = character(),
  set_idents = TRUE
)

#Pseudobulk Analysis by cell type
pb_ex <- run_pb_deseq2(
  seu = AD_CAA,
  cellclass = "ExNeuron",
  cellclass_col = "cellclass",
  donor_col = "orig.ident",
  group_col = "FDX",
  group_levels = c("Control", "AD+CAA")
)

pb_inh <- run_pb_deseq2(
  seu = AD_CAA,
  cellclass = "InhNeuron",
  cellclass_col = "cellclass",
  donor_col = "orig.ident",
  group_col = "FDX",
  group_levels = c("Control", "AD+CAA")
)

ex_up_pb <- get_sig_up(pb_ex$res)
ex_down_pb <- get_sig_down(pb_ex$res)

inh_up_pb <- get_sig_up(pb_inh$res)
inh_down_pb <- get_sig_down(pb_inh$res)

#Fischers Overlap Test (Need Comparison Gene List)
tests <- run_fisher_overlap_grid(
  query_sets = list(
    ExNeuron_ADCAA_up = ex_up_pb,
    ExNeuron_ADCAA_down = ex_down_pb,
    InhNeuron_ADCAA_up = inh_up_pb,
    InhNeuron_ADCAA_down = inh_down_pb
  ),
  reference_sets = list(
    Reference1_UP = Reference1_UP,
    Reference2_UP = Reference2_UP
  ),
  universe = unique(c(pb_ex$res$gene, pb_inh$res$gene))
)

tests


summary_all <- run_cellclass_overlap(
  seu = AD_CAA,
  celltype_col = "cellclass",
  reference_sets = list(
    Reference1_UP = Reference1_UP,
    Reference2_UP = Reference2_UP
  ),
  direction = "up",
  min_pct = 0.1,
  min_cells_group = 10,
  padj_cutoff = 0.05
)

write.table(
  summary_all,
  file = "OverlapFishers_Test.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


