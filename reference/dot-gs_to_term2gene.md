# Convert named list of gene sets to TERM2GENE data.frame

Creates a two-column data.frame (term, gene) required by
`clusterProfiler::GSEA()`.

## Usage

``` r
.gs_to_term2gene(gene.sets)
```

## Arguments

- gene.sets:

  Named list of character vectors.

## Value

A 2-column data.frame with columns `term` and `gene`.
