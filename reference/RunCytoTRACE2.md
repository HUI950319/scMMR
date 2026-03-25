# Run CytoTRACE 2 Cellular Potency Prediction

Predict single-cell developmental potential (stemness / differentiation)
using CytoTRACE 2. Wraps `CytoTRACE2::cytotrace2()` with a cleaner
interface, automatic species handling, and Seurat v5 compatibility.

## Usage

``` r
RunCytoTRACE2(object, ...)

# S3 method for class 'Seurat'
RunCytoTRACE2(
  object,
  species = c("human", "mouse"),
  assay = NULL,
  layer = c("counts", "data"),
  batch_size = 10000,
  smooth_batch_size = 1000,
  cores = NULL,
  seed = 14,
  ...
)

# Default S3 method
RunCytoTRACE2(
  object,
  species = c("human", "mouse"),
  batch_size = 10000,
  smooth_batch_size = 1000,
  cores = NULL,
  seed = 14,
  ...
)
```

## Arguments

- object:

  A Seurat object, matrix (genes x cells), or data.frame.

- ...:

  Additional arguments passed to `CytoTRACE2::cytotrace2()`.

- species:

  Character: `"human"` or `"mouse"`. Default `"human"`.

- assay:

  Character. Seurat assay to use (default: active assay). Ignored for
  matrix input.

- layer:

  Character: `"counts"` (default) or `"data"`. Ignored for matrix input.

- batch_size:

  Integer or `NULL`. Number of cells per batch for KNN smoothing.
  Default 10000 (recommended for \>10K cells). Set `NULL` to disable
  subsampling.

- smooth_batch_size:

  Integer or `NULL`. Cells to subsample within each batch for diffusion
  smoothing. Default 1000.

- cores:

  Integer. Number of parallel cores. Default `NULL` (auto-detect, use
  half).

- seed:

  Integer. Random seed for reproducibility (default 14).

## Value

For Seurat input: the Seurat object with metadata columns added:

- CytoTRACE2_Score:

  Predicted potency score (0-1, higher = more potent)

- CytoTRACE2_Potency:

  Potency category (Differentiated, Unipotent, Oligopotent, Multipotent,
  Pluripotent, Totipotent)

- CytoTRACE2_Relative:

  Relative order normalised to 0-1

- preKNN_CytoTRACE2_Score:

  Score before KNN smoothing

- preKNN_CytoTRACE2_Potency:

  Category before KNN smoothing

For matrix/data.frame input: a data.frame with cell IDs as row names and
the same columns as above.

## Details

CytoTRACE 2 predicts cellular potency in three steps:

1.  Preprocess input data

2.  Predict potency (discrete categories + continuous score 0-1)

3.  Smooth predictions via diffusion + KNN rescaling

If CytoTRACE2 is not installed, the function will prompt you with
installation instructions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Seurat object
seu <- RunCytoTRACE2(seu, species = "human")
FeaturePlot(seu, features = "CytoTRACE2_Score")
VlnPlot(seu, features = "CytoTRACE2_Score", group.by = "cell_type")

# Raw matrix (genes x cells)
result <- RunCytoTRACE2(counts_matrix, species = "mouse")
} # }
```
