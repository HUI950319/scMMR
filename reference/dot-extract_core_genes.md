# Extract core_enrichment genes from RunGseaEnrich result

Extract core_enrichment genes from RunGseaEnrich result

## Usage

``` r
.extract_core_genes(res, type = "gsea", p.cutoff = 0.05)
```

## Arguments

- res:

  A `RunGseaEnrich` result list or a `compareClusterResult` object.

- type:

  Which result to extract: `"gsea"` or `"enricher"`.

- p.cutoff:

  Filter terms by p.adjust cutoff. Default: 0.05.

## Value

A named list: term name -\> character vector of core_enrichment genes.
