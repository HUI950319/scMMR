# Prepare bulk expression input for deconvolution (internal)

Handles multiple input formats (matrix, data.frame, CSV path, h5ad path,
Seurat object) and returns a standardised list.

## Usage

``` r
.prepare_bulk_input(bulk_expr)
```

## Arguments

- bulk_expr:

  Input bulk expression data.

## Value

List with: `matrix` (genes x samples, numeric), `sample_names`
(character), `gene_names` (character).
