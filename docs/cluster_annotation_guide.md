# Cluster Annotation Guide


This document describes the workflow used to annotate clusters in human snRNA-seq data using the AD/CAA dataset as an example.

---

## Overview

Cluster annotation was performed using:

- Differentially expressed genes (DEGs)
- Enrichment analysis (Enrichr)
- Reference atlases:
  - Azimuth Cell Types 2021 (primary)
  - Allen Brain Atlas scRNA-seq 2021
  - Tabula Sapiens (secondary validation for human datasets)
  - Tabula Muris (secondary validation for mouse datasets)

---

## Step 1: Prepare Cluster DEG Table

1. Start from the output file: `AD_CAA_cluster_markers.txt`
2. Filter DEGs using:
- ` Adjusted p-value (q-value) < 0.05`

3. For each cluster:
- Select marker genes
- Export into a CSV file where:
  - Each column = one cluster
  - Column header = cluster ID (e.g., C0, C1, ...)

4. Save as: `AD_CAA_ClusterDEGs.csv`
💡 Tip: This step is easiest in Excel for quick filtering and sorting.

---

## Step 2: Enrichment Analysis

For each cluster:

1. Upload marker genes to Enrichr
2. Use relevant libraries:
- Azimuth Cell Types 2021 (primary)
- Allen Brain Atlas (secondary)

3. Record:
- Top enriched cell types
- Combined scores
- Consistency across databases

---

## Step 3: Annotation Criteria

Cluster labels were assigned based on:

- Agreement between Azimuth and Allen Atlas
- Strength of enrichment signal
- Biological plausibility (Cannonical Markers - SLC17A7, GAD1, GAD2, GFAP, MBP, VCAN) 

### Interpretation Symbols

- `*` → Strong agreement between Azimuth and Allen Atlas  
- `-` → Allen Atlas disagrees with Azimuth (top 3 combined score)

---

## Step 4: Final Cluster Annotations

| Cluster | # DEGs | Cell Type            | Confidence |
|--------|--------|----------------------|------------|
| C0     | 6593   | Oligodendrocytes     | * |
| C1     | 8892   | Excitatory Neuron    | * |
| C2     | 8631   | Oligodendrocytes     | * |
| C3     | 5146   | Oligodendrocytes     | * |
| C4     | 5542   | Oligodendrocytes     | * |
| C5     | 7041   | Astrocytes           | * |
| C6     | 7411   | Excitatory Neuron    | * |
| C7     | 8908   | Excitatory Neuron    | * |
| C8     | 5576   | Astrocytes           | * |
| C9     | 8198   | Excitatory Neuron    | * |
| C10    | 6951   | Inhibitory Neuron    | * |
| C11    | 3919   | Oligodendrocytes     | * |
| C12    | 6201   | Inhibitory Neuron    | * |
| C13    | 5103   | Inhibitory Neuron    | - |
| C14    | 7699   | Microglia            | * |
| C15    | 5355   | OPC                  | * |
| C16    | 10591  | Excitatory Neuron    | * |
| C17    | 11296  | Excitatory Neuron    | * |
| C18    | 7304   | Excitatory Neuron    | * |
| C19    | 6712   | Oligodendrocytes     | * |
| C20    | 8021   | Excitatory Neuron    | * |
| C21    | 6177   | Oligodendrocytes     | * |
| C22    | 7028   | Microglia            | * |
| C23    | 6309   | Excitatory Neuron    | * |
| C24    | 4660   | Oligodendrocytes     | * |
| C25    | 8957   | Excitatory Neuron    | * |
| C26    | 3171   | OPC                  | * |
| C27    | 5703   | Inhibitory Neuron    | * |
| C28    | 6022   | Inhibitory Neuron    | - |
| C29    | 6370   | Vascular             | * |
| C30    | 6111   | Endothelial          | * |
| C31    | 7622   | Excitatory Neuron    | * |
| C32    | 11847  | Excitatory Neuron    | * |
| C33    | 10905  | Excitatory Neuron    | * |
| C34    | 3753   | Excitatory Neuron    | * |
| C35    | 4129   | OPC                  | * |
| C36    | 6081   | Astrocytes           | * |

---

## Step 5: Applying Annotations in Seurat

Example:

```r
cluster_ids <- c(
"0" = "Oligodendrocytes",
"1" = "Excitatory Neuron",
"5" = "Astrocytes",
"14" = "Microglia"
)

AD_CAA$celltype <- plyr::mapvalues(
AD_CAA$seurat_clusters,
from = names(cluster_ids),
to = cluster_ids
)
```

---

## Step 6: Annotation QC Plots

After applying manual annotations, generate UMAP and marker-expression plots to confirm that broad cell classes and canonical markers are biologically consistent.

### 6.1 UMAP of collapsed cell classes

```r
# Save UMAP without split
cairo_ps(filename = "UMAP-ADCAA-CTRL_Nosplit.eps", width = 14, height = 8)
DimPlot(
  AD_CAA,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5
)
dev.off()

# Save UMAP split by disease group
cairo_ps(filename = "UMAP-ADCAA-CTRL.eps", width = 18, height = 8)
DimPlot(
  AD_CAA,
  reduction = "umap",
  label = FALSE,
  split.by = "FDX",
  pt.size = 0.5
)
dev.off()
```
### 6.2 Cannonical Marker Plot 

```r
markers.to.plot <- c(
  "GRIN2A", "RBFOX3", "SLC17A7", "SYN3",        # Excitatory neurons
  "GAD1", "GAD2",                              # Inhibitory neurons
  "GJA1", "SLC1A3", "SLC1A2", "GFAP",         # Astrocytes
  "MBP", "PLP1", "MOBP", "MOG",               # Oligodendrocytes
  "VCAN", "NXPH1",                             # OPCs
  "INPP5D", "CSF1R", "TGFBR1", "TREM2",       # Microglia
  "VIM", "COL1A2",                             # Fibroblast-like cells
  "FLT1", "CLDN5",                             # Endothelial cells
  "PDGFRB", "KCNJ8",                           # Pericytes
  "ACTA2", "TAGLN",                            # Smooth muscle cells
  "PDGFRA", "DCN"                              # VLMC / meningeal-associated cells
)

markers.to.plot <- unique(markers.to.plot)

desired_order <- c(
  "ExNeuron",
  "InhNeuron",
  "Astrocytes",
  "Oligodendrocytes",
  "OPC",
  "Microglia",
  "Fibroblast",
  "Endothelial",
  "Pericytes",
  "SMC",
  "VLMC"
)
AD_CAA$cellclass <- as.character(Idents(AD_CAA))

AD_CAA$cellclass <- factor(
  AD_CAA$cellclass,
  levels = rev(desired_order)
)

Idents(AD_CAA) <- AD_CAA$cellclass

cairo_ps(filename = "ClusterFeatureMarkers-CollapsedAD_CAAProject.eps", width = 15, height = 4)
DotPlot(
  AD_CAA,
  scale = TRUE,
  features = markers.to.plot,
  dot.min = 0.05
) +
  RotatedAxis() +
  scale_colour_gradient2(
    low = "#298c8c",
    mid = "#b8b8b8",
    high = "#a00000"
  )
dev.off()
```
### 6.3 Interpretation

Use these plots to verify that:

Excitatory neurons express markers such as SLC17A7, GRIN2A, SYN3
Inhibitory neurons express GAD1 and GAD2
Astrocytes express GFAP, GJA1, SLC1A2, SLC1A3
Oligodendrocytes express MBP, PLP1, MOBP, MOG
OPCs express VCAN and NXPH1
Microglia express INPP5D, CSF1R, TREM2
Vascular-associated populations show expected endothelial, pericyte, SMC, and VLMC markers

These plots provide a final manual quality-control step before downstream analyses such as subclustering, differential expression, module scoring, or cell-cell communication analysis.
