# Expression distribution plots (violin / box / bar / dot)

Visualises the expression of one or more features stratified by a
grouping variable. Accepts a Seurat object (genes are fetched from the
specified assay / layer) or a plain data.frame.

## Usage

``` r
ExpressionStatPlot2(
  seu,
  features,
  group.by,
  split.by = NULL,
  assay = NULL,
  layer = "data",
  cells = NULL,
  plot.by = c("group", "feature"),
  type = c("violin", "box", "bar", "dot"),
  fill.by = c("group", "feature"),
  palette = NULL,
  alpha = 0.8,
  add_box = FALSE,
  box_color = "black",
  box_width = 0.1,
  add_point = FALSE,
  pt.color = "grey30",
  pt.size = NULL,
  pt.alpha = 1,
  jitter.width = 0.4,
  jitter.height = 0.1,
  add_trend = FALSE,
  trend_color = "black",
  trend_linewidth = 1,
  trend_ptsize = 2,
  comparisons = NULL,
  ref_group = NULL,
  pairwise_method = "wilcox.test",
  multiplegroup_comparisons = FALSE,
  multiple_method = "kruskal.test",
  sig_label = c("p.signif", "p.format"),
  sig_labelsize = 3.5,
  bg.by = NULL,
  bg_palette = NULL,
  bg_alpha = 0.15,
  same.y.lims = FALSE,
  y.min = NULL,
  y.max = NULL,
  y.nbreaks = 5,
  sort = FALSE,
  stack = FALSE,
  flip = FALSE,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  aspect.ratio = NULL,
  base_size = 14,
  facet_nrow = NULL,
  facet_ncol = NULL,
  combine = TRUE,
  force = FALSE,
  seed = 11
)
```

## Arguments

- seu:

  A Seurat object \*\*or\*\* a data.frame containing feature columns and
  the grouping column.

- features:

  Character vector of feature names (gene names or metadata / data.frame
  column names).

- group.by:

  Column name used for grouping (x-axis in `"group"` mode, fill colour
  in `"feature"` mode).

- split.by:

  Column name for splitting into separate sub-plots. Passed through to
  `plt_con()` in `"group"` mode; creates a patchwork row in `"feature"`
  mode. Default `NULL`.

- assay:

  Seurat assay to use. `NULL` = active assay.

- layer:

  Seurat layer to pull expression from. Default `"data"`.

- cells:

  Optional character vector of cell barcodes to subset. `NULL` = all
  cells.

- plot.by:

  Axis orientation. `"group"` (default) or `"feature"`.

- type:

  Plot geometry: `"violin"` (default), `"box"`, `"bar"`, or `"dot"`.

- fill.by:

  What drives fill colour: `"group"` (default) or `"feature"`.

- palette:

  Colour palette name (passed to `UtilsR::plt_con()`). `NULL` = package
  default.

- alpha:

  Fill transparency (0–1). Default `0.8`.

- add_box:

  Logical. Overlay a narrow box-plot (median + IQR) on violin or bar.
  Default `FALSE`.

- box_color:

  Colour of the overlay box-plot. Default `"black"`.

- box_width:

  Width of the overlay box-plot. Default `0.1`.

- add_point:

  Logical. Overlay jittered raw points. Default `FALSE`.

- pt.color:

  Point colour. Default `"grey30"`.

- pt.size:

  Point size. `NULL` = auto (scaled by cell count).

- pt.alpha:

  Point transparency. Default `1`.

- jitter.width:

  Horizontal jitter extent. Default `0.4`.

- jitter.height:

  Vertical jitter extent. Default `0.1`.

- add_trend:

  Logical. Connect group medians / means with a trend line. Default
  `FALSE`.

- trend_color:

  Trend line colour. Default `"black"`.

- trend_linewidth:

  Trend line width. Default `1`.

- trend_ptsize:

  Trend point size. Default `2`.

- comparisons:

  List of length-2 vectors for pairwise statistical tests. `NULL` =
  none.

- ref_group:

  Reference group for one-vs-all comparisons. `NULL` = no reference.

- pairwise_method:

  Statistical test for pairwise comparisons. Default `"wilcox.test"`.

- multiplegroup_comparisons:

  Logical. Add an overall group comparison label. Default `FALSE`.

- multiple_method:

  Test for overall comparison. Default `"kruskal.test"`.

- sig_label:

  Significance label format: `"p.signif"` (default) or `"p.format"`.

- sig_labelsize:

  Text size for significance labels. Default `3.5`.

- bg.by:

  Column name for background colour bands (must be a super-grouping of
  `group.by`). `NULL` = none.

- bg_palette:

  Palette for background bands. `NULL` = default.

- bg_alpha:

  Transparency of background bands. Default `0.15`.

- same.y.lims:

  Logical. Share y limits across all panels. Default `FALSE`.

- y.min:

  Minimum y value (numeric) or quantile string (e.g. `"q1"`). `NULL` =
  data minimum.

- y.max:

  Maximum y value or quantile string. `NULL` = data maximum.

- y.nbreaks:

  Number of y-axis tick breaks. Default `5`.

- sort:

  Logical or `"increasing"`. Sort groups by median of the first feature.
  Default `FALSE`.

- stack:

  Logical. Stack all features into a single faceted plot (`"group"` mode
  only). Default `FALSE`.

- flip:

  Logical. Flip coordinates (horizontal layout). Default `FALSE`.

- title:

  Plot title. `NULL` = none.

- subtitle:

  Plot subtitle. `NULL` = none.

- xlab:

  X-axis label. `NULL` = auto.

- ylab:

  Y-axis label. `NULL` = auto.

- legend.position:

  Legend position. Default `"right"`.

- legend.direction:

  Legend direction. Default `"vertical"`.

- aspect.ratio:

  Aspect ratio of each panel. `NULL` = free.

- base_size:

  Base font size. Default `14`.

- facet_nrow, facet_ncol:

  Number of rows / columns in the combined patchwork layout.

- combine:

  Logical. Combine panels with patchwork. Default `TRUE`.

- force:

  Logical. Skip the high-cardinality group check. Default `FALSE`.

- seed:

  Random seed for jitter reproducibility. Default `11`.

## Value

A ggplot2 or patchwork object.

## Details

Two axis orientations are available via `plot.by`:

- `"group"` (default) — one panel per feature; groups on the x-axis.
  Powered by `UtilsR::plt_con()`.

- `"feature"` — one combined panel; features on the x-axis, groups
  encoded by fill colour.

## Examples

``` r
if (FALSE) { # \dontrun{
# Seurat object — expression of two genes by cluster
ExpressionStatPlot2(seu, features = c("CD3D", "CD8A"),
                    group.by = "seurat_clusters")

# Feature-axis mode: features on x, cell types as colour
ExpressionStatPlot2(seu, features = c("CD3D", "CD8A", "MS4A1"),
                    group.by = "cell_type", plot.by = "feature",
                    type = "violin")

# Plain data.frame
df <- data.frame(GeneA = rnorm(100), GeneB = rnorm(100),
                 cluster = sample(letters[1:4], 100, replace = TRUE))
ExpressionStatPlot2(df, features = c("GeneA", "GeneB"),
                    group.by = "cluster", type = "box")
} # }
```
