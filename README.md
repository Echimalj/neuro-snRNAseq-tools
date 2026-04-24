# neuro-snRNAseq-tools

Reusable functions for single-nucleus RNA-seq analysis:
- Seurat workflows
- pseudobulk analysis
- CellChat utilities

🚧 Work in progress


Reusable utilities for preprocessing and analyzing single-nucleus RNA-seq datasets in Seurat.

## Included functionality
- 10x sample loading with SoupX ambient RNA correction
- sample-sheet-based metadata annotation
- QC metric calculation and filtering
- DoubletFinder wrapper for per-sample doublet detection
- SCTransform preprocessing
- Harmony-based integration and clustering
- RNA assay preparation for marker detection

## Example workflow
See `scripts/human_ad_caa_pipeline.R` for a full example using a human AD/CAA dataset.

## Expected input
A sample sheet CSV with:
- `sample_id`
- `sample_path`
- `condition`
- `batch`
- `species`
- `SoupX ambient RNA removal required TRUE/FALSE`
- `orig.ident`
- `genotype`
- `FDX`

## Notes
This toolkit is designed to separate reusable preprocessing logic from project-specific biological analyses.
