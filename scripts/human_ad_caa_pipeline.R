source("R/check_dependencies.R")
source("R/load_samples.R")
source("R/metadata_utils.R")
source("R/qc_utils.R")
source("R/doublet_utils.R")
source("R/preprocess_utils.R")
source("R/integration_utils.R")
source("R/marker_utils.R")
source("R/checkpoint_utils.R")
source("R/annotation_utils.R")
source("R/subcluster_utils.R")
source("R/speckle_utils.R")

###Check Dependencies
check_required_packages(c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tibble"
))

check_required_packages(c(
  "SoupX"
))
#Check Point Directory
checkpoint_dir <- "checkpoints"

##Load Samples
sample_sheet <- read.csv("human_sample_sheet.csv")

seurat_list <- load_samples_from_sheet(
  sample_sheet = sample_sheet,
  min_cells = 3,
  min_features = 200,
  use_soupx = TRUE
)

##Annotate Metadata
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

##QC metrics 
AD_CAA <- add_qc_metrics(
  seu = AD_CAA,
  species = "human",
  add_ribo = FALSE,
  add_hb = FALSE
)

plot_basic_qc(AD_CAA)

AD_CAA_before_qc <- AD_CAA

AD_CAA <- filter_cells_basic(
  seu = AD_CAA,
  min_features = 200,
  max_mt = 1
)

summarize_filtering(
  before = AD_CAA_before_qc,
  after = AD_CAA,
  group_col = "orig.ident"
)

##Doublet Finder
check_required_packages(c(
  "DoubletFinder"
  ))

AD_CAA <- annotate_doublets_by_sample(
  seu = AD_CAA,
  split_by = "orig.ident",
  assay = "RNA",
  sct = FALSE,
  resolution = 0.1
)

summarize_doublets(AD_CAA, group_col = "orig.ident")

AD_CAA <- keep_singlets(AD_CAA)

##SCT Transform 
check_required_packages(c(
  "sctransform",
  "harmony",
  "clustree"
))

options(future.globals.maxSize = 5000 * 1024^2)

AD_CAA <- run_sct_pipeline(
  seu = AD_CAA,
  vst_flavor = "v2",
  npcs = 22,
  dims = 1:20,
  resolution = 0.5,
  assay = "RNA"
)

plot_clusters(AD_CAA)

AD_CAA <- checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_after_SCT.rds"),
  reload = TRUE
)
  
#Harmony Integration for Batch Correction
AD_CAA <- run_harmony_clustering(
  seu = AD_CAA,
  group_var = "condition",
  assay = "SCT",
  dims = 1:20,
  resolution = 0.5
)

plot_clusters(AD_CAA, reduction = "umap")

#Clustree to find adequate resolution
AD_CAA <- test_clustering_resolutions(
  seu = AD_CAA,
  reduction = "harmony",
  dims = 1:20,
  resolutions = c(0.025, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3)
)

plot_clustree(AD_CAA, prefix = "SCT_snn_res.")

AD_CAA <- run_harmony_clustering(
  seu = AD_CAA,
  group_var = "condition",
  assay = "SCT",
  dims = 1:20,
  resolution = 0.5
)

#Find all markers
AD_CAA <- prepare_rna_for_markers(
  seu = AD_CAA,
  assay = "RNA",
  run_join_layers = TRUE
)
AD_CAA <- checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_after_RNA_scaling.rds"),
  reload = TRUE
)

AD_CAA_markers <- find_cluster_markers(
  seu = AD_CAA,
  assay = "RNA",
  only_pos = TRUE
)

save_marker_table(
  markers = AD_CAA_markers,
  file = "AD_CAA_cluster_markers.txt"
)

save_checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_before-annotation.rds"),
  reload = TRUE
)

## Follow Intructions of the  `cluster_annotation_guide.md` then:

check_required_packages(c(
  "plyr",
  "pheatmap"
))

#Set labels for annotated clusters
new_cluster_ids <- c(
  "Oligodendrocytes1",
  "ExNeuron1",
  "Oligodendrocytes2",
  "Oligodendrocytes3",
  "Oligodendrocytes4",
  "Astrocytes1",
  "ExNeuron2",
  "ExNeuron3",
  "Astrocytes2",
  "ExNeuron4",
  "InhNeuron1",
  "Oligodendrocytes5",
  "InhNeuron2",
  "InhNeuron3",
  "Microglia1",
  "OPC1",
  "ExNeuron5",
  "ExNeuron6",
  "ExNeuron7",
  "Oligodendrocytes6",
  "ExNeuron8",
  "Oligodendrocytes7",
  "Microglia2",
  "ExNeuron9",
  "Oligodendrocytes8",
  "ExNeuron10",
  "OPC2",
  "InhNeuron4",
  "InhNeuron5",
  "Vascular",
  "Endothelial",
  "ExNeuron11",
  "ExNeuron12",
  "ExNeuron13",
  "ExNeuron14",
  "OPC3",
  "Astrocytes3"
)

new_cluster_ids <- make_cluster_label_vector(
  labels = new_cluster_ids,
  cluster_levels = levels(AD_CAA)
)

AD_CAA <- apply_cluster_labels(
  seu = AD_CAA,
  cluster_labels = new_cluster_ids,
  new_metadata_col = "celltype",
  set_idents = TRUE

save_checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_after-annotation.rds"),
  reload = TRUE
)

#OPTIONAL: Subclustering
  clusters_to_subcluster <- c(
  "Astrocytes1",
  "Astrocytes2",
  "Microglia1",
  "Microglia2",
  "Vascular"
)

AD_CAA <- subcluster_selected_idents(
  seu = AD_CAA,
  clusters_to_subcluster = clusters_to_subcluster,
  output_col = "subcluster",
  original_col = "SCT_original_cluster",
  resolution = 0.1,
  dims = 1:20
)

AD_CAA <- set_idents_from_metadata(AD_CAA, "subcluster")

cor_mat <- compute_identity_correlation(
  seu = AD_CAA,
  assay = "RNA",
  slot = "data",
  method = "pearson"
)

save_correlation_matrix(cor_mat, "Corr.mat.txt")
plot_identity_correlation(cor_mat)

#Merging similar subclusters
  merge_map <- c(
  "Astrocytes2" = "Astrocytes1",
  "Astrocytes3" = "Astrocytes1",
  "Astrocytes4" = "Astrocytes1",
  "Astrocytes5" = "Astrocytes1",
  "Astrocytes10" = "Astrocytes1",
  "Astrocytes12" = "Astrocytes1",
  "Astrocytes6" = "Astrocytes2",
  "Astrocytes7" = "Astrocytes3",
  "Astrocytes8" = "Astrocytes4",
  "Astrocytes9" = "Astrocytes5",
  "Astrocytes11" = "Astrocytes5",
  "Microglia2" = "Microglia1",
  "Microglia3" = "Microglia1",
  "Microglia6" = "Microglia1",
  "Microglia4" = "Microglia2",
  "Microglia5" = "Microglia3",
  "Microglia8" = "Microglia3",
  "Microglia7" = "Microglia4",
  "Microglia9" = "Microglia5"
)

AD_CAA <- merge_identity_labels(
  seu = AD_CAA,
  merge_map = merge_map,
  metadata_col = "final_celltype",
  set_idents = TRUE
)

save_checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_after-subcluster.rds"),
  reload = TRUE
)

#Cell Proportions Speckle 

check_required_packages(c("speckle", "limma"))

prop_all <- run_propeller(
  seu = AD_CAA,
  sample_col = "orig.ident",
  group_col = "FDX"
)

neuron_idents <- c(
  "ExNeuron1", "ExNeuron2", "ExNeuron3", "ExNeuron4", "ExNeuron5",
  "ExNeuron6", "ExNeuron7", "ExNeuron8", "ExNeuron9", "ExNeuron10",
  "ExNeuron11", "ExNeuron12", "ExNeuron13", "ExNeuron14",
  "InhNeuron1", "InhNeuron2", "InhNeuron3", "InhNeuron4", "InhNeuron5"
)

other_idents <- c(
  "Oligodendrocytes1", "Oligodendrocytes2", "Oligodendrocytes3",
  "Oligodendrocytes4", "Oligodendrocytes5", "Oligodendrocytes6",
  "Oligodendrocytes7", "Oligodendrocytes8",
  "OPC1", "OPC2", "OPC3",
  "SMC", "Astrocytes1", "Astrocytes2", "Astrocytes3",
  "Astrocytes4", "Astrocytes5",
  "Microglia1", "Microglia2", "Microglia3", "Microglia4", "Microglia5",
  "Endothelial", "VLMC1", "VLMC2", "Pericytes", "Fibroblast"
)

prop_sets <- run_propeller_sets(
  seu = AD_CAA,
  identity_sets = list(
    Neurons = neuron_idents,
    Other = other_idents
  ),
  sample_col = "orig.ident",
  group_col = "FDX"
)

save_propeller_results(
  propeller_results = prop_all,
  output_dir = "results",
  prefix = "CellTypeProportionsALL_ADCAAvsCTRL"
)

save_propeller_results(
  propeller_results = prop_sets,
  output_dir = "results",
  prefix = "CellTypeProportions_ADCAAvsCTRL"
)

save_checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_final_processed.rds"),
  reload = TRUE
)
