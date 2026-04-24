# Cluster Annotation Guide

This document describes the workflow used to annotate clusters i.e. human snRNA-seq data (AD/CAA dataset).

---

## Overview

Cluster annotation was performed using:

- Differentially expressed genes (DEGs)
- Enrichment analysis (Enrichr)
- Reference atlases:
  - Azimuth Cell Types 2021 (primary)
  - Allen Brain Atlas scRNA-seq 2021
  - Tabula Sapiens (secondary validation) - Tabula Muris (For Mouse Experiments)

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
