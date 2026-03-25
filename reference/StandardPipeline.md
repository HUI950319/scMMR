# Standard single-cell analysis pipeline

Run a standard Seurat workflow: Normalize, FindVariableFeatures,
ScaleData, PCA, (optional Harmony batch correction), UMAP,
FindNeighbors, FindClusters.

## Usage

``` r
StandardPipeline(
  seu,
  nfeatures = 2000,
  npcs = 50,
  dims = 1:30,
  use.harmony = FALSE,
  harmony.group.by = NULL,
  harmony.dims = 1:50,
  resolution = 0.6,
  k.param = 20,
  vars.to.regress = NULL,
  verbose = FALSE,
  seed = 42
)
```

## Arguments

- seu:

  Seurat object

- nfeatures:

  Number of highly variable features. Default: 2000

- npcs:

  Number of PCA dimensions to compute. Default: 50

- dims:

  Dimensions to use for downstream steps (UMAP, neighbors). Default:
  1:30

- use.harmony:

  Whether to run Harmony batch correction. Default: FALSE

- harmony.group.by:

  Character vector of metadata columns for Harmony grouping. Required
  when `use.harmony = TRUE`. Default: NULL

- harmony.dims:

  Dimensions passed to `harmony::RunHarmony()`. Default: 1:50

- resolution:

  Clustering resolution for `FindClusters`. Default: 0.6

- k.param:

  Number of nearest neighbors for `FindNeighbors`. Default: 20

- vars.to.regress:

  Variables to regress out during `ScaleData`. Default: NULL

- verbose:

  Whether to print Seurat messages. Default: FALSE

- seed:

  Random seed. Default: 42

## Value

Seurat object with PCA, (Harmony), UMAP reductions and cluster
assignments.
