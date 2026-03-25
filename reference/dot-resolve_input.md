# Resolve input to AnnData (internal)

Accepts either an h5ad file path or a Seurat object, returns a Python
AnnData object.

## Usage

``` r
.resolve_input(input, embedding_key = "umap")
```

## Arguments

- input:

  Character (h5ad path) or Seurat object.

- embedding_key:

  Embedding key (used for Seurat objects).

## Value

Python AnnData object.
