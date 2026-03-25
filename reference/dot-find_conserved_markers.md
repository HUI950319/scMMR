# Find conserved markers across grouping variable levels

Adapted from scop::FindConservedMarkers2. Splits cells by grouping.var,
runs FindMarkers within each level, intersects significant genes, and
combines p-values.

## Usage

``` r
.find_conserved_markers(
  object,
  assay,
  layer,
  cells.1,
  cells.2,
  grouping.var,
  features = NULL,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  base = 2,
  pseudocount.use = 1,
  mean.fxn = NULL,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  latent.vars = NULL,
  only.pos = FALSE,
  min.cells.group = 3,
  min.cells.feature = 3,
  meta.method = "maximump",
  norm.method = "LogNormalize",
  verbose = TRUE,
  ...
)
```
