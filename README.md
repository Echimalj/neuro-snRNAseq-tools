# neuro-snRNAseq-tools

Reusable functions for single-nucleus RNA-seq analysis:
- Seurat workflows
- pseudobulk analysis
- CellChat utilities

🚧 Work in progress

![R](https://img.shields.io/badge/R-4.x-blue)
![Seurat](https://img.shields.io/badge/Seurat-v5-green)
![Status](https://img.shields.io/badge/status-in%20progress-orange)

Reusable, modular workflows for preprocessing and analyzing single-nucleus RNA-seq (snRNA-seq) data using Seurat.

## Overview

This repository provides a set of reusable functions to streamline common steps in single-nucleus RNA-seq analysis, including:

- Data loading (10x Genomics, with optional SoupX correction)
- Metadata integration via sample sheets
- Quality control and filtering
- Doublet detection
- SCTransform-based normalization
- Batch correction with Harmony
- Clustering and visualization
- Marker gene identification

The goal is to separate reusable computational workflows from project-specific biological analyses, enabling reproducibility and scalability across datasets.

## Key Features
Modular design → each step is implemented as reusable functions
Sample-sheet driven workflows → no hard-coded sample metadata
Flexible preprocessing → optional SoupX, SCT, Harmony
Memory-aware pipeline → supports checkpointing for large datasets
Extensible → designed to integrate with downstream tools (CellChat, pseudobulk, etc.)

## Example workflow
See `scripts/human_ad_caa_pipeline.R` for a full example using a human AD/CAA dataset.

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
  
##Documentation

Cluster annotation workflow:

docs/cluster_annotation_guide.md

##Design Philosophy

This toolkit is designed around three principles:

Reproducibility → eliminate hard-coded logic
Modularity → each step is independently reusable
Scalability → compatible with large snRNA-seq datasets and HPC environments

## Notes
This toolkit is designed to separate reusable preprocessing logic from project-specific biological analyses.
