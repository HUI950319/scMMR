# Remove ambient RNA contamination via decontX

Estimate and remove ambient RNA contamination from a Seurat object using
celda::decontX. The corrected counts are stored in a new assay
"decontX".

## Usage

``` r
ComputeAmbientRNA(seu, split.by = NULL, cluster.name = NULL)
```

## Arguments

- seu:

  Seurat object

- split.by:

  The grouping variable used to split the input object (batch
  correction). Default: NULL.

- cluster.name:

  Name of cluster field in Seurat object metadata. Must be pre-defined.

## Value

Seurat object with an additional "decontX" assay (corrected counts) and
a \`decontX_contamination\` column in metadata.
