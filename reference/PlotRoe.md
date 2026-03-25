# O/E Ratio Heatmap for Cell Type Composition

Create a heatmap displaying the Observed / Expected ratio of cell type
composition across groups, highlighting group-specific enrichment.
Optionally computes per-cell p-values via chi-squared standardised
residuals or Fisher's exact test, and marks significant cells with
stars.

## Usage

``` r
PlotRoe(
  cellmeta,
  by,
  fill,
  palette = "Blues",
  column_names_rot = 0,
  display = c("value", "symbol", "both"),
  show.pvalue = FALSE,
  p.method = c("chisq", "fisher"),
  p.adjust.method = "BH",
  sig.threshold = c(0.001, 0.01, 0.05),
  return.data = FALSE,
  ...
)
```

## Arguments

- cellmeta:

  A Seurat object or data.frame containing cell metadata (e.g.
  `seurat_obj` or `seurat_obj@meta.data`).

- by:

  Character string specifying the grouping variable (e.g. `"group"`,
  `"sample"`).

- fill:

  Character string specifying the cell-type variable (e.g.
  `"cell_type_pred"`).

- palette:

  RColorBrewer palette name. Default `"Blues"`.

- column_names_rot:

  Rotation angle (0-360) for column labels. Default 0.

- display:

  Character. What to show inside each cell:

  `"value"`

  :   (Default) Numeric O/E ratio (e.g. `1.58`).

  `"symbol"`

  :   Grade symbols based on O/E level: `+++` (\>1), `++` (0.8–1), `+`
      (0.2–0.8), `+/-` (0–0.2), `-` (0).

  `"both"`

  :   Both numeric value and grade symbol (value on top, symbol below).

- show.pvalue:

  Logical. Whether to compute and display significance stars on the
  heatmap. Default `FALSE`. When combined with `display = "symbol"`,
  stars are appended to the symbol (e.g. `+++ ***`).

- p.method:

  Character. Method for per-cell p-values: `"chisq"` (default,
  chi-squared standardised residuals) or `"fisher"` (Fisher's exact test
  per 2x2 cell).

- p.adjust.method:

  P-value adjustment method passed to
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`.

- sig.threshold:

  Numeric vector defining significance cut-offs for star labels. Default
  `c(0.001, 0.01, 0.05)` (`***`, `**`, `*`).

- return.data:

  Logical. If `TRUE`, return a list with the O/E matrix, p-value matrix
  and the heatmap object. Default `FALSE`.

- ...:

  Additional arguments passed to `ComplexHeatmap::Heatmap()`.

## Value

When `return.data = FALSE` (default): a ComplexHeatmap Heatmap object.

When `return.data = TRUE`: a list with components `oe` (O/E matrix),
`pvalue` (adjusted p-value matrix or `NULL`), and `heatmap`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic O/E heatmap (numeric values)
PlotRoe(seu, by = "group", fill = "cell_type_pred")

# Grade symbols (+++ / ++ / + / +/- / -)
PlotRoe(seu, by = "group", fill = "cell_type_pred",
        display = "symbol")

# Both value and symbol
PlotRoe(seu, by = "group", fill = "cell_type_pred",
        display = "both", show.pvalue = TRUE)

# Return data for downstream use
res <- PlotRoe(seu, by = "group", fill = "cell_type_pred",
               show.pvalue = TRUE, return.data = TRUE)
res$oe       # O/E ratio matrix
res$pvalue   # adjusted p-value matrix
} # }
```
