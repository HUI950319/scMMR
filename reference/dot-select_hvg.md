# Select highly variable genes with batch fallback (internal)

Wraps `mt_select_hvg` with a `tryCatch` that falls back to non-batch
mode when the data is too small for batch-aware selection.

## Usage

``` r
.select_hvg(adata, n_top = 6000L, batch_key = NULL)
```

## Arguments

- adata:

  Python AnnData object.

- n_top:

  Number of HVGs (default 6000).

- batch_key:

  obs column for batch-aware selection (NULL = no batch).

## Value

Character vector of HVG names.
