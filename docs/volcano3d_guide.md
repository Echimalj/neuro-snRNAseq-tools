# Three-way Differential Expression and 3D Volcano Analysis

## Method Overview

To identify genotype-dependent transcriptional changes within specific cell types, we implemented a three-way differential expression framework adapted from the volcano3D package (v2.0.9), developed by Katriona Goldmann.

This implementation was adapted from discussions in the volcano3D repository:
https://github.com/KatrionaGoldmann/volcano3D/issues/15

The workflow has been extended to support modern Seurat objects and includes additional statistical controls to mitigate pseudoreplication.

---

## Key Modifications in This Repository

Compared to the original implementation:

* Updated for compatibility with **Seurat v5**
* Integrated directly with **single-cell cluster identities**
* Added **ANOVA-based global testing per gene**
* Included **sample-level covariate (orig.ident)** to reduce pseudobulking bias
* Implemented **q-value correction** for overall significance
* Standardized outputs for downstream pathway analysis

---

## Workflow Description

For each cell type of interest:

1. Cells were grouped by **cell type × genotype**
   (e.g., Astrocytes-WT, Astrocytes-Tg-FDD, Astrocytes-Tg-FDD;Mapt⁻/⁻)

2. Expression data (scaled RNA assay) were extracted for the selected cells

3. Three pairwise contrasts were computed:

   * A vs B
   * B vs C
   * C vs A

4. For each contrast:

   * Differential expression was computed using Seurat’s `FindMarkers`
   * Parameters:

     * `min.pct = 0`
     * `logfc.threshold = 0`
   * Outputs:

     * p-values
     * BH-adjusted p-values
     * log₂ fold changes

---

## Global Statistical Testing

In parallel, an ANOVA model was applied per gene:

* Model:

  ```text
  Expression ~ Genotype + Sample
  ```

* Purpose:

  * Test **overall expression differences across all groups**
  * Reduce **pseudoreplication bias** by accounting for sample-level variation

* Implementation:

  * Linear model (`lm`)
  * Type II ANOVA (`car::Anova`)
  * Multiple testing correction using **qvalue**

---

## Gene Classification Strategy

Genes were classified based on the **combined pattern of significance and directionality across all contrasts**, not on any single comparison.

Each gene was assigned to one of seven categories:

| Category | Interpretation          |
| -------- | ----------------------- |
| ns       | Not significant         |
| R        | B                       |
| Y        | B+C                     |
| G        | C                       |
| C        | A+C                     |
| DB       | A                       |
| P        | A+B                     |

---

## Visualization

Classification and visualization were performed using:

* `polar_coords()` from volcano3D
* `radial_plotly()`

These functions map genes into polar coordinate space, capturing:

* magnitude of expression change
* directionality across contrasts
* specificity to biological conditions

---

## Outputs

* Polar (radial) plots
* Integrated DEG tables (pairwise + ANOVA)
* Gene lists per classification category (for enrichment analysis)

---

## Notes

* This method is particularly useful for:

  * multi-condition experiments
  * genotype comparisons
  * identifying shared vs condition-specific transcriptional programs

* Interpretation should consider:

  * consistency across contrasts
  * biological plausibility
  * downstream pathway enrichment

---
