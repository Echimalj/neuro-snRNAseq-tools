# ============================================================
# 3D Volcano Analysis – Example Workflow (ExNeuron6)
# ============================================================
#
# Description:
# This script demonstrates a three-way differential expression
# analysis workflow for single-nucleus RNA-seq data using a
# polar (3D volcano) visualization framework.
#
# The analysis compares three genotype conditions:
# WT vs FDD vs FDDTKO within a specific cell type (ExNeuron6).
#
# ------------------------------------------------------------
# Method Attribution:
# This workflow is adapted from the volcano3D R package:
 :contentReference[oaicite:0]{index=0}
 Developed by :contentReference[oaicite:1]{index=1}
#
# Original concept and discussion:
# https://github.com/KatrionaGoldmann/volcano3D/issues/15
#
# ------------------------------------------------------------
# Extensions implemented in this repository:
# - Integration with Seurat objects (v5)
# - Automated cell type × genotype identity construction
# - Pairwise differential expression using Seurat::FindMarkers
# - Gene-wise ANOVA with sample-level covariate (orig.ident)
#   to reduce pseudoreplication
# - q-value multiple testing correction for global significance
# - Unified classification of genes across all contrasts
# - Export of volcano3D gene categories for downstream analysis
#
# ------------------------------------------------------------
# Inputs:
# - Seurat object with:
#     - cell type identities
#     - genotype/condition metadata (Genotype)
#     - sample identifiers (orig.ident)
#
# ------------------------------------------------------------
# Outputs:
# - Interactive radial (3D volcano) plot
# - SVG figure for publication
# - CSV file with gene classification categories
#
# ------------------------------------------------------------
# Author:
# Enrique Chimal
# PhD Candidate – Medical Neuroscience
#
# ============================================================

# NOTE:
# This script is intended as an example. For reusable workflows,
# see R/volcano3d_utils.R


source("R/volcano3d_utils.R")
source("R/check_dependencies.R")

check_required_packages(c(
  "Seurat",
  "volcano3D",
  "car",
  "qvalue",
  "plotly"
))

WTFDDTKO <- readRDS("WTFDDTKO_2025.rds")

WTFDDTKO <- add_volcano3d_identity(
  seu = WTFDDTKO,
  group_col = "Genotype",
  output_col = "celltype.genotype",
  sep = "-"
)

res_ex6 <- run_volcano3d_workflow(
  seu = WTFDDTKO,
  celltype_label = "ExNeuron6",
  groups = c("WT", "FDD", "FDDTKO"),
  group_col = "Genotype",
  sample_col = "orig.ident",
  combined_col = "celltype.genotype",
  assay = "RNA",
  expr_layer = "scale.data",
  sep = "-",
  use_sample_covariate = TRUE,
  pcutoff = 0.05,
  labs = c(
    "ns",
    "FDDTKO",
    "FDDTKO+WT",
    "WT",
    "FDD+WT",
    "FDD",
    "FDD+FDDTKO"
  )
)

p <- volcano3D::radial_plotly(res_ex6$polar)
p

plotly::save_image(p, "figures/radial_ExNeuron6.svg")

save_volcano3d_gene_categories(
  polar = res_ex6$polar,
  categories = c(
    "FDDTKO",
    "FDDTKO+WT",
    "WT",
    "FDD+WT",
    "FDD",
    "FDD+FDDTKO"
  ),
  file = "results/ExNeuron6_DEG_WTvsFDDvsFDDTKO.csv"
)
