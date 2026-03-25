# Run GSEA for each cell type via FindMarkers + clusterProfiler::GSEA

For each cell type, calls `FindMarkers()` to compute log2 fold changes,
ranks genes by logFC, and runs `clusterProfiler::GSEA()`.

## Usage

``` r
.run_gsea_across_celltypes(
  seu,
  gene.sets,
  celltypes,
  ident.1,
  ident.2,
  pvalue.cutoff,
  logfc.threshold,
  eps,
  minGSSize,
  maxGSSize,
  verbose
)
```

## Arguments

- seu:

  Seurat object with Idents set to composite group.

- gene.sets:

  Named list of gene sets.

- celltypes:

  Character vector of cell types to process.

- ident.1, ident.2:

  Condition labels.

- pvalue.cutoff:

  Passed to `clusterProfiler::GSEA()`.

- logfc.threshold:

  Passed to `FindMarkers()`.

- eps:

  Passed to `clusterProfiler::GSEA()`.

- minGSSize, maxGSSize:

  Gene set size limits for GSEA.

- verbose:

  Logical; print progress.

## Value

A data.frame with standardized columns.
