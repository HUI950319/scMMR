# Convert named list of gene sets to SCPA pathway format

SCPA's `compare_pathways()` expects pathways as a list of data.frames,
each with columns `Pathway` and `Genes`.

## Usage

``` r
.gs_to_scpa_format(gene.sets)
```

## Arguments

- gene.sets:

  Named list of character vectors.

## Value

A list of data.frames suitable for `SCPA::compare_pathways()`.
