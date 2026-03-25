# Mark doublets via DoubletFinder

Detect doublets in a Seurat object using DoubletFinder. Compatible with
Seurat V5 and DoubletFinder \>= 2.06.

## Usage

``` r
ComputeDoublets(seu, PCs = 1:10, split.by = NULL, num.cores = 1)
```

## Arguments

- seu:

  Seurat object

- PCs:

  Vectors indicating used principal components. Default: 1:10

- split.by:

  Name of a metadata column to split by (run DoubletFinder per group).
  Default: NULL

- num.cores:

  Threads for calculation. Default: 1

## Value

Seurat object with \`DF.classifications\` column added to metadata.
Values: "Doublet.hc" (high confidence), "Doublet.lc" (low confidence),
"Singlet".
