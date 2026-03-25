# Plot Differential Expression Results

Unified visualization for differential expression analysis. Supports
volcano, MA, pct (diff_pct volcano), manhattan, ring, and PCA +
hierarchical clustering plots.

## Usage

``` r
PlotDE(
  x,
  type = c("volcano", "ma", "pct", "pct_mirror", "manhattan", "ring", "pca_hc"),
  fc_threshold = 1,
  p_threshold = 0.05,
  highlight_fc_threshold = 3,
  point_size = 1.2,
  highlight_size = 3,
  colors = c("darkred", "grey50", "royalblue4", "green2"),
  symmetry = FALSE,
  label_genes = NULL,
  nlabel = 5,
  jitter_width = 0.5,
  jitter_height = 0.4,
  label_size = 4,
  group.by = NULL,
  reduction = "pca",
  dims = 1:2,
  dist_method = "euclidean",
  hclust_method = "ward.D2",
  k = NULL,
  add_ellipses = TRUE,
  combine = TRUE,
  ncol = NULL,
  y_max = NULL
)
```

## Arguments

- x:

  For `type = "volcano"`, `"ma"`, `"pct"`, `"manhattan"`, `"ring"`: a
  data.frame of DE results (e.g., from `RunDE`). For `type = "pca_hc"`:
  a Seurat object.

- type:

  Plot type. One of `"volcano"`, `"ma"`, `"pct"`, `"pct_mirror"`,
  `"manhattan"`, `"ring"`, or `"pca_hc"`.

- fc_threshold:

  Log2 fold change threshold for significance lines (volcano/ma).
  Default: 1.

- p_threshold:

  Adjusted p-value threshold. Default: 0.05.

- highlight_fc_threshold:

  FC threshold for gene labels (volcano/ma). Default: 3.

- point_size:

  Point size. Default: 1.2.

- highlight_size:

  Size of highlighted points. Default: 3.

- colors:

  Color vector for Up/NS/Down/Highlight groups. Default:
  `c("darkred", "grey50", "royalblue4", "green2")`.

- symmetry:

  Logical; symmetric x-axis for volcano. Default: FALSE.

- label_genes:

  Character vector of gene names to label manually. Default: NULL
  (auto-select by FC).

- nlabel:

  Number of top genes to label per group (manhattan/ring). Default: 5.

- jitter_width:

  Jitter width for manhattan/ring. Default: 0.5.

- jitter_height:

  Jitter height for manhattan. Default: 0.4.

- label_size:

  Label text size. Default: 4.

- group.by:

  Column name for grouping (pca_hc). Default: NULL.

- reduction:

  Reduction to use for PCA (pca_hc). Default: "pca".

- dims:

  Dimensions to plot (pca_hc). Default: 1:2.

- dist_method:

  Distance method for dendrogram (pca_hc). Default: "euclidean".

- hclust_method:

  Clustering method (pca_hc). Default: "ward.D2".

- k:

  Number of clusters for dendrogram coloring. Default: NULL (auto =
  number of groups).

- add_ellipses:

  Logical; add confidence ellipses (pca_hc). Default: TRUE.

- combine:

  Logical; combine multi-group volcano plots. Default: TRUE.

- ncol:

  Number of columns for combined facets. Default: NULL.

- y_max:

  Numeric; cap -log10(p) y-axis at this value (volcano/pct). Default:
  NULL (auto: 99.5th percentile).

## Value

A ggplot object (or combined plot for pca_hc).

## See also

[`RunDE`](https://hui950319.github.io/scMMR/reference/RunDE.md)

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- RunDE(seu, group.by = "celltype")
PlotDE(markers, type = "volcano")
PlotDE(markers, type = "pct")
PlotDE(markers, type = "manhattan")
PlotDE(seu, type = "pca_hc", group.by = "celltype")
} # }
```
