# neuro-snRNAseq-tools

Reusable functions for single-nucleus RNA-seq analysis:
- Seurat workflows
- pseudobulk analysis


🚧 Work in progress

![R](https://img.shields.io/badge/R-4.x-blue)
![Seurat](https://img.shields.io/badge/Seurat-v5-green)
![Status](https://img.shields.io/badge/status-in%20progress-orange)

Reusable, modular workflows for preprocessing and analyzing single-nucleus RNA-seq (snRNA-seq) data using Seurat.

## Overview

This repository provides a set of reusable functions to streamline common steps in single-nucleus RNA-seq analysis, 
from upstream Cell Ranger processing to downstream Seurat-based analysis

- HPC/Slurm-based Cell Ranger processing
- Data loading (10x Genomics, optional SoupX correction)
- Sample-sheet driven metadata integration
- Quality control and filtering
- Doublet detection (DoubletFinder)
- SCTransform normalization and dimensionality reduction
- Batch correction with Harmony
- Clustering and visualization
- Marker gene detection
- Cell-type proportion testing (Speckle / propeller)
- Pseudobulk differential expression (DESeq2)
- Gene set overlap testing (Fisher exact test)
- Module scoring (DAA, DAM, DAO, etc.)
- Multi-condition analysis (3D volcano framework)

The goal is to separate reusable computational workflows from project-specific biological analyses, enabling reproducibility and scalability across datasets.

---

## Key Features
**Modular design**  
Each step is implemented as reusable, independent functions.

**Sample-sheet driven workflows**  
No hard-coded metadata; pipelines are dataset-agnostic.

**Flexible preprocessing**  
Optional SoupX, SCTransform, Harmony integration.

**Memory-aware execution**  
Checkpoint utilities prevent crashes in large datasets.

**Extensible architecture**  
Designed to integrate with downstream tools such as CellChat, DESeq2, and custom gene set analyses.

---
## Example workflow
See `scripts/human_ad_caa_pipeline.R` for a full example using a human AD/CAA dataset.
See `scripts/advanced_analysis_examples.R` for a extended analysis in human AD/CAA dataset.
See `scripts/cellranger_slurm_template` for a full example of HPC Slurm cellranger submission.

## Expected input
A sample sheet CSV (.csv) with:
- `sample_id`
- `sample_path`
- `condition`
- `batch`
- `species`
- `SoupX ambient RNA removal required TRUE/FALSE`
- `orig.ident`
- `genotype`
- `FDX`
- any additional metadata fields (automatically added to Seurat object)
  
## Repository Structure

## Repository Structure

```text
neuro-snRNAseq-tools/
├── R/                                   # Core reusable R modules
│   ├── load_samples.R                   # 10x loading + optional SoupX correction
│   ├── metadata_utils.R                 # Sample sheet metadata integration
│   ├── qc_utils.R                       # QC metrics + filtering
│   ├── doublet_utils.R                  # DoubletFinder workflow (per-sample)
│   ├── preprocess_utils.R               # SCTransform + PCA + clustering
│   ├── integration_utils.R              # Harmony integration + clustering
│   ├── marker_utils.R                   # Marker detection + RNA assay prep
│   ├── annotation_utils.R               # Manual cluster relabeling
│   ├── subcluster_utils.R               # Targeted subclustering + correlation merging
│   ├── speckle_utils.R                  # Cell proportion testing (propeller)
│   ├── module_score_utils.R             # Gene module scoring (DAA, DAM, DAO, etc.)
│   ├── differential_expression_utils.R  # Two-group DE analysis (global + per cell type)
│   ├── atlas_similarity_utils.R         # Cosine Similarity Analysis with null distribution
│   ├── check_dependencies.R             # Check for required modules
│   ├── overlap_utils.R                  # Fisher overlap + enrichment testing
│   ├── pseudobulk_utils.R               # DESeq2 pseudobulk workflows
│   ├── volcano3d_utils.R                # Three group-condition DE + radial visualization
│   ├── human_sample_sheet.csv           # Example of Metadata required for load_samples utils
│   └── checkpoint_utils.R               # Memory-safe save/load checkpoints
│
├── scripts/                             # Example workflows
│   ├── human_ad_caa_pipeline.R          # End-to-end preprocessing pipeline
│   ├── advanced_analysis_examples.R     # Module scoring + DE + downstream analyses
│   ├── 3D_Volcano_Analysis_Example.R    # 3-group DE + radial visualization
│   └── cellranger_slurm_template.sh     # HPC/Slurm Cell Ranger pipeline
│
├── docs/                                # Documentation
│   ├── cellranger_pipeline.md           # HPC Cell Ranger workflow guide
│   ├── cluster_annotation_guide.md      # Manual cluster annotation strategy
│   └── volcano3d_guide.md               # Multi-condition DE framework
│
├── README.md                            # Project overview and usage
└── .gitignore                           # Ignore large files / outputs
```


## Documentation

- `docs/cellranger_pipeline.md`  
  HPC/Slurm workflow for running Cell Ranger on multiple samples using job arrays.

- `docs/cluster_annotation_guide.md`  
  Manual cluster annotation strategy using marker genes, enrichment resources, and reference atlases.

- `docs/volcano3d_guide.md`  
  Three-way differential expression and radial visualization framework.

## Advanced Analysis Workflows

In addition to preprocessing, this repository includes utilities for downstream biological analyses:

- Pseudobulk differential expression (DESeq2)
- Gene set overlap testing (Fisher exact test)
- Cell-type-specific enrichment workflows
- Module scoring and expression profiling
- Multi-condition differential expression (3-group comparisons)

See:

`scripts/advanced_analysis_examples.R`
`scripts/3D_Volcano_Analysis_Example.R`


## Design Philosophy

This toolkit is designed around three principles:

Reproducibility → eliminate hard-coded logic
Modularity → each step is independently reusable
Scalability → compatible with large snRNA-seq datasets and HPC environments

## Notes
This toolkit is designed to separate reusable preprocessing logic from project-specific biological analyses.
