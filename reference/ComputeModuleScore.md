# Calculate gene module scores

Compute gene module activity scores using AUCell, Seurat
(AddModuleScore), or UCell methods. Supports multiple input formats for
gene sets: named list, data.frame, or GMT file.

## Usage

``` r
ComputeModuleScore(x, ...)

# Default S3 method
ComputeModuleScore(
  x,
  gene.sets,
  method = c("AUCell", "Seurat", "UCell"),
  min.size = 20,
  batch.size = 500,
  nbin = 24,
  ctrl = 100,
  cores = 1,
  seed = 11,
  ...
)

# S3 method for class 'Seurat'
ComputeModuleScore(
  x,
  gene.sets,
  method = c("AUCell", "Seurat", "UCell"),
  min.size = 20,
  batch.size = 500,
  nbin = 24,
  ctrl = 100,
  cores = 1,
  seed = 11,
  assay = Seurat::DefaultAssay(x),
  layer = NULL,
  store = c("assay", "metadata"),
  assay.name = NULL,
  prefix = NULL,
  ...
)
```

## Arguments

- x:

  A gene expression matrix (genes x cells) or a Seurat object.

- ...:

  Additional arguments passed to scoring methods.

- gene.sets:

  Gene sets in one of three formats:

  - A named list of character vectors (gene names).

  - A data.frame with at least two columns: term (gene set name) and
    gene (gene symbol).

  - A file path to a GMT file.

- method:

  Scoring method: `"AUCell"` (default), `"Seurat"`, or `"UCell"`.

- min.size:

  Minimum number of genes (after filtering) for a gene set to be scored.
  Default: 20.

- batch.size:

  Number of cells per batch for AUCell method to reduce memory. Default:
  500.

- nbin:

  Number of expression bins for Seurat method control gene selection.
  Default: 24.

- ctrl:

  Number of control genes per feature gene for Seurat method. Default:
  100.

- cores:

  Number of parallel cores. Default: 1.

- seed:

  Random seed for reproducibility. Default: 11.

- assay:

  Name of the Seurat assay to use. Default: `DefaultAssay(x)`.

- layer:

  Data layer to extract: `"counts"` or `"data"`. Defaults to `"counts"`
  for AUCell/UCell and `"data"` for Seurat method.

- store:

  Storage destination for scores: `"assay"` (default) stores as a new
  assay, `"metadata"` stores as columns in `meta.data`.

- assay.name:

  Name of the new assay when `store = "assay"`. Default: the method
  name.

- prefix:

  Column name prefix when `store = "metadata"`. Default: the method name
  followed by underscore (e.g. `"AUCell_setA"`).

## Value

A score matrix (gene sets x cells) or a Seurat object with scores stored
as an assay.
