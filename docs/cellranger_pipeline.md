# Cell Ranger Pipeline (HPC / Slurm)

This document describes how to run 10x Genomics Cell Ranger on a high-performance computing (HPC) cluster using Slurm job arrays.

---

## Overview

This workflow processes multiple snRNA-seq samples in parallel using:

* **Cell Ranger v9**
* **Slurm job arrays**
* Shared reference transcriptome (GRCh38)

Each sample is processed independently using `cellranger count`.

---

## Key Features

* Parallel processing via Slurm arrays
* Scalable to large cohorts
* Centralized reference genome
* Per-sample output directories
* Resource-aware configuration (CPU, memory, runtime)

---

## Directory Structure

```text
project/
  fastq/
  reference/
    GRCh38/
  cellranger_outs/
  scripts/
    cellranger_slurm.sh
```

---

## Slurm Script

```bash
#!/bin/bash
#SBATCH --job-name=cellranger_ADCAA
#SBATCH -p general
#SBATCH --output=cellrangerAD_%A_%a.out
#SBATCH --error=cellrangerAD_%A_%a.err
#SBATCH --array=1-8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@domain.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-10:00:00
#SBATCH --mem=256G
#SBATCH -A your_allocation

module load cellranger/9.0.1

TRANSCRIPTOME="/path/to/GRCh38"
FASTQ_DIR="/path/to/fastqs"
BASE_OUTDIR="/path/to/output"

SAMPLES=(Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8)

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$INDEX]}
OUTDIR="$BASE_OUTDIR/${SAMPLE}"

mkdir -p "$OUTDIR"

cellranger count \
  --id="$SAMPLE" \
  --transcriptome="$TRANSCRIPTOME" \
  --fastqs="$FASTQ_DIR" \
  --sample="$SAMPLE" \
  --localcores=$SLURM_CPUS_PER_TASK \
  --create-bam=false \
  --output-dir="$OUTDIR"
```

---

## Running the Pipeline

Submit the job:

```bash
sbatch scripts/cellranger_slurm.sh
```

Each array task processes one sample:

```text
Task 1 → Sample 2052
Task 2 → Sample 2495
...
```

---

## Input Requirements

### FASTQ files

* Located in a single directory
* Must follow standard 10x naming:

  ```
  SAMPLE_S1_L001_R1_001.fastq.gz
  SAMPLE_S1_L001_R2_001.fastq.gz
  ```

### Reference

* Pre-built Cell Ranger reference (e.g., GRCh38)
* Generated using:

```bash
cellranger mkref
```

---

## Output

Each sample produces:

```text
cellranger_outs/
  SAMPLE/
    outs/
      filtered_feature_bc_matrix/
      raw_feature_bc_matrix/
      metrics_summary.csv
      web_summary.html
```

---

## Resource Guidelines

Typical requirements for snRNA-seq:

* CPUs: 8–16
* Memory: 128–256 GB
* Runtime: 1–4 days depending on dataset size

---

## Notes

* `--create-bam=false` reduces storage usage
* Ensure FASTQ directory contains all lanes for each sample
* Match `--sample` exactly to FASTQ prefixes
* Job arrays must match number of samples

---

## Design Rationale

This approach:

* Avoids manual per-sample submission
* Scales efficiently across HPC environments
* Standardizes preprocessing before downstream Seurat workflows

---

## Integration with This Repository

Outputs from Cell Ranger feed directly into:

```text
R/load_samples.R
```

using either:

* `filtered_feature_bc_matrix/` (standard Seurat input)
* full Cell Ranger output (for SoupX workflows)

---
