# Getting Started with scMMR

## Overview

**scMMR** (Single-Cell Multi-task Model in R) is a comprehensive toolkit
for single-cell RNA-seq analysis powered by PyTorch multi-task deep
neural networks. It provides: - Joint cell type classification and
embedding regression - Bulk RNA-seq deconvolution with cell type
proportion estimation - Gene importance attribution via Integrated
Gradients - Pathway enrichment analysis (GSEA, SCPA) - Perturbation
ranking with distribution-distance metrics - Pseudotime trajectory
analysis and cellular potency prediction - A rich visualization suite
(18+ publication-ready plot functions)

## Installation

### Step 1: Install scMMR from GitHub

``` r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("HUI950319/scMMR")
```

### Step 2: Set Up Python Environment

scMMR uses PyTorch via `reticulate`. You can either create a fresh conda
environment or point to an existing one.

**Option A: Automatic installation** (creates a new conda environment):

``` r
library(scMMR)
install_scMMR_python(envname = "scMMR-env", method = "conda")
```

This installs: `torch`, `scanpy`, `anndata`, `numpy`, `pandas`, `scipy`,
`h5py`.

**Option B: Use an existing conda environment**:

``` r
library(scMMR)
use_scMMR_python(condaenv = "/path/to/your/conda/env")
```

### Step 3: Verify Setup

``` r
library(scMMR)
use_scMMR_python(condaenv = "/path/to/your/conda/env")

# Check that Python modules are accessible
reticulate::py_module_available("torch")
reticulate::py_module_available("scanpy")
```

## Quick Start: Predict Cell Types

Here is a minimal example using the bundled test data and pre-trained
model:

``` r
library(scMMR)
library(Seurat)

# 1. Activate Python environment
use_scMMR_python(condaenv = "/path/to/your/conda/env")

# 2. Load test data and model
seu <- qs::qread(system.file("extdata", "toy_test.qs", package = "scMMR"))
model_path <- system.file("extdata", "model.pt", package = "scMMR")

# 3. Run prediction with gene importance
result <- DNN_predict(
  input     = seu,
  model_dir = dirname(model_path),
  explain   = TRUE,
  top_k_global = 20,
  device    = "cpu"
)

# 4. Merge predictions into Seurat object
seu <- AddMetaData(seu, result$predictions)

# 5. Visualize cell type composition
PlotAlluvia(seu, by = "group", fill = "cell_type_pred")
```

## What’s Next?

- [Cell Type
  Classification](https://hui950319.github.io/scMMR/articles/cell-type-classification.md)
  – Train and predict with multi-task DNN
- [Bulk
  Deconvolution](https://hui950319.github.io/scMMR/articles/deconvolution.md)
  – Estimate cell type proportions from bulk RNA-seq
- [Pathway
  Analysis](https://hui950319.github.io/scMMR/articles/pathway-analysis.md)
  – GSEA and SCPA enrichment
- [Perturbation
  Ranking](https://hui950319.github.io/scMMR/articles/perturbation-ranking.md)
  – Rank cell types by treatment response
- [Visualization
  Gallery](https://hui950319.github.io/scMMR/articles/visualization-gallery.md)
  – All 18+ plot functions
