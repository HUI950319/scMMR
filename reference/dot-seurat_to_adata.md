# Convert Seurat object to Python AnnData (internal)

Extracts counts matrix, metadata, and embeddings from a Seurat object
and passes them to the Python helper.

## Usage

``` r
.seurat_to_adata(srt, embedding_key = "umap")
```

## Arguments

- srt:

  Seurat object.

- embedding_key:

  Which reduction to include (default "umap").

## Value

Python AnnData object.
