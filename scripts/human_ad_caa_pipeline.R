source("R/check_dependencies.R")
source("R/load_samples.R")
source("R/metadata_utils.R")
source("R/qc_utils.R")
source("R/doublet_utils.R")
source("R/preprocess_utils.R")
source("R/integration_utils.R")
source("R/marker_utils.R")
source("R/checkpoint_utils.R")


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
  file = file.path(checkpoint_dir, "AD_CAA_before-annotation.rds")
)
