# Run SCPA for each cell type via seurat_extract + compare_pathways

For each cell type, extracts expression matrices for the two conditions
and runs `SCPA::compare_pathways()` for distribution-based comparison.

## Usage

``` r
.run_scpa_across_celltypes(
  seu,
  gene.sets,
  celltypes,
  ident.1,
  ident.2,
  cores,
  verbose
)
```

## Arguments

- seu:

  Seurat object with composite group column `..pw_group..`.

- gene.sets:

  Named list of gene sets.

- celltypes:

  Character vector of cell types to process.

- ident.1, ident.2:

  Condition labels.

- cores:

  Number of parallel cores for SCPA.

- verbose:

  Logical; print progress.

## Value

A data.frame with standardized columns.
