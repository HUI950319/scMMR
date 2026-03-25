# Parse geneset input for RunGseaEnrich

Converts various geneset input formats into a TERM2GENE data.frame with
columns `gs_name` and `gene_symbol`.

## Usage

``` r
.parse_geneset_to_term2gene(geneset)
```

## Arguments

- geneset:

  Gene set specification. One of:

  msigdbr args

  :   A list with msigdbr parameters (e.g.
      `list(species = "Homo sapiens", collection = "H")`). Recognized
      keys: `species`, `collection`, `subcollection`, `category`,
      `subcategory`.

  data.frame

  :   A data.frame with columns `gs_name` and `gene_symbol` (TERM2GENE
      format).

  named list

  :   A named list of character vectors, e.g.
      `list(PathA = c("TP53","BRCA1"))`.

  GMT file path

  :   A character string path to a GMT file.

## Value

A data.frame with columns `gs_name` and `gene_symbol`.
