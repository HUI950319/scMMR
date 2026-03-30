# Expression distribution plots for Seurat objects

Generates violin, box, bar, dot, or column plots for one or more
features (genes or metadata numeric columns) stratified by a grouping
variable. Full feature parity with `scop::FeatureStatPlot` with
optimised parameter names consistent with `DimPlot2` / `FeaturePlot2`.

## Usage

``` r
FeaturePlot3(
  seu,
  features,
  group.by = NULL,
  split.by = NULL,
  bg.by = NULL,
  plot.by = c("group", "feature"),
  fill.by = c("group", "feature", "expression"),
  cells = NULL,
  assay = NULL,
  layer = "data",
  keep_empty = FALSE,
  individual = FALSE,
  type = c("violin", "box", "bar", "dot", "col"),
  palette = "npg",
  palcolor = NULL,
  alpha = 1,
  bg_palette = "Set2",
  bg_palcolor = NULL,
  bg_alpha = 0.2,
  add_box = FALSE,
  box_color = "black",
  box_width = 0.1,
  box_ptsize = 2,
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
  add_stat = c("none", "mean", "median"),
  stat_color = "black",
  stat_size = 1,
  stat_stroke = 1,
  stat_shape = 25,
  add_line = NULL,
  line_color = "red",
  line_size = 1,
  line_type = 1,
  cells.highlight = NULL,
  cols.highlight = "red",
  sizes.highlight = 1,
  alpha.highlight = 1,
  calculate_coexp = FALSE,
  same.y.lims = FALSE,
  y.min = NULL,
  y.max = NULL,
  y.trans = "identity",
  y.nbreaks = 5,
  sort = FALSE,
  stack = FALSE,
  flip = FALSE,
  comparisons = NULL,
  ref_group = NULL,
  pairwise_method = "wilcox.test",
  multiplegroup_comparisons = FALSE,
  multiple_method = "kruskal.test",
  sig_label = c("p.signif", "p.format"),
  sig_labelsize = 3.5,
  aspect.ratio = NULL,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = "Expression level",
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = NULL,
  theme_use = NULL,
  theme_args = list(),
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  force = FALSE,
  seed = 11
)
```

## Arguments

- seu:

  A Seurat object or data.frame. When a data.frame is provided, `assay`
  and `layer` are ignored.

- features:

  Character vector of feature names (genes or numeric metadata columns).

- group.by:

  Column name for the x-axis grouping variable. `NULL` treats all cells
  as a single group.

- split.by:

  Column name for a secondary grouping (dodging within each x-axis
  position). Default `NULL`.

- bg.by:

  Column name whose levels define background colour bands. Must be a
  super-grouping of `group.by` (each `group.by` level belongs to exactly
  one `bg.by` level). Default `NULL`.

- plot.by:

  Axis orientation: `"group"` (default) places groups on the x-axis (one
  panel per feature); `"feature"` places features on the x-axis (one
  panel per group level).

- fill.by:

  What drives fill colour:

  - `"group"` (default) — group or split levels.

  - `"feature"` — feature name.

  - `"expression"` — median expression per group (continuous gradient).

- cells:

  Character vector of cell barcodes to include. `NULL` = all.

- assay:

  Seurat assay name. `NULL` = active assay.

- layer:

  Seurat layer. Default `"data"`.

- keep_empty:

  Logical. Keep empty factor levels on the x-axis. Default `FALSE`.

- individual:

  Logical. Create a separate panel for each group level × feature ×
  split level combination. Default `FALSE`.

- type:

  Plot geometry. One of `"violin"` (default), `"box"`, `"bar"`, `"dot"`,
  or `"col"` (per-cell column bar).

- palette:

  Colour palette name for fill colours. Default `"npg"`.

- palcolor:

  Custom colour vector. Overrides `palette`. Default `NULL`.

- alpha:

  Fill transparency (0–1). Default `1`.

- bg_palette:

  Colour palette for background bands. Default `"Set2"`.

- bg_palcolor:

  Custom colours for background bands. Default `NULL`.

- bg_alpha:

  Transparency of background bands. Default `0.2`.

- add_box:

  Logical. Overlay a narrow box-plot (violin/bar only). Default `FALSE`.

- box_color:

  Colour of overlay box-plot. Default `"black"`.

- box_width:

  Width of overlay box-plot. Default `0.1`.

- box_ptsize:

  Median point size in box-plot. Default `2`.

- add_point:

  Logical. Overlay jittered individual points. Default `FALSE`.

- pt.color:

  Point colour. Default `"grey30"`.

- pt.size:

  Point size. `NULL` = auto.

- pt.alpha:

  Point transparency. Default `1`.

- jitter.width:

  Horizontal jitter. Default `0.4`.

- jitter.height:

  Vertical jitter. Default `0.1`.

- add_trend:

  Logical. Connect group medians/means with a trend line (violin/box:
  median; bar: mean). Default `FALSE`.

- trend_color:

  Trend line colour. Default `"black"`.

- trend_linewidth:

  Trend line width. Default `1`.

- trend_ptsize:

  Trend point size. Default `2`.

- add_stat:

  Summary stat to overlay. One of `"none"` (default), `"mean"`, or
  `"median"`.

- stat_color:

  Stat point colour. Default `"black"`.

- stat_size:

  Stat point size. Default `1`.

- stat_stroke:

  Stat point stroke. Default `1`.

- stat_shape:

  Stat point shape (integer). Default `25` (filled downward triangle).

- add_line:

  Numeric y-intercept for a horizontal reference line. `NULL` = none.

- line_color:

  Reference line colour. Default `"red"`.

- line_size:

  Reference line width. Default `1`.

- line_type:

  Reference line type. Default `1` (solid).

- cells.highlight:

  Character vector of cell barcodes to highlight, or `TRUE` for all
  cells. Requires `add_point = TRUE`. Default `NULL`.

- cols.highlight:

  Highlight point colour. Default `"red"`.

- sizes.highlight:

  Highlight point size. Default `1`.

- alpha.highlight:

  Highlight point transparency. Default `1`.

- calculate_coexp:

  Logical. Compute geometric mean co-expression of all `features` and
  add it as an extra `"CoExp"` panel. Only valid when all features are
  genes (in the expression matrix). Default `FALSE`.

- same.y.lims:

  Logical. Share y-axis limits across all panels. Default `FALSE`.

- y.min:

  Minimum y limit: numeric value or quantile string (e.g. `"q1"`).
  `NULL` = data minimum.

- y.max:

  Maximum y limit: numeric value or quantile string (e.g. `"q99"`).
  `NULL` = data maximum.

- y.trans:

  Y-axis transformation. `"identity"` (default) or `"log2"`.

- y.nbreaks:

  Number of y-axis tick breaks. Default `5`.

- sort:

  Logical or `"increasing"` / `"decreasing"`. Sort groups by median of
  the first feature. Default `FALSE`.

- stack:

  Logical. Stack all feature panels vertically (or horizontally when
  `flip = TRUE`) with shared axes and a single shared legend. Default
  `FALSE`.

- flip:

  Logical. Flip coordinates (horizontal layout). Default `FALSE`.

- comparisons:

  List of length-2 character vectors specifying pairwise comparisons, or
  `TRUE` to auto-detect groups from `split.by`. `NULL` = none.

- ref_group:

  Reference group for one-vs-all comparisons. Default `NULL`.

- pairwise_method:

  Pairwise test. Default `"wilcox.test"`.

- multiplegroup_comparisons:

  Logical. Add an overall group test label. Default `FALSE`.

- multiple_method:

  Overall group test. Default `"kruskal.test"`.

- sig_label:

  Significance label format: `"p.signif"` (default) or `"p.format"`.

- sig_labelsize:

  Text size of significance labels. Default `3.5`.

- aspect.ratio:

  Panel aspect ratio. `NULL` = free.

- title:

  Plot title. `NULL` = none.

- subtitle:

  Plot subtitle. `NULL` = none.

- xlab:

  X-axis label. `NULL` = auto (`group.by`).

- ylab:

  Y-axis label. Default `"Expression level"`.

- legend.position:

  Legend position string or coordinates. Default `"right"`.

- legend.direction:

  Legend direction. Default `"vertical"`.

- legend.title:

  Custom legend title. `NULL` = auto.

- theme_use:

  Theme specification. `NULL` =
  [`ggplot2::theme_bw`](https://ggplot2.tidyverse.org/reference/ggtheme.html);
  accepts a character function name, a function, or a `theme` object.

- theme_args:

  List of extra arguments passed to `theme_use`. Default
  [`list()`](https://rdrr.io/r/base/list.html).

- combine:

  Logical. Combine panels into a single patchwork. Default `TRUE`.

- nrow:

  Number of rows in patchwork layout.

- ncol:

  Number of columns in patchwork layout.

- byrow:

  Logical. Fill patchwork by row. Default `TRUE`.

- force:

  Logical. Skip the \>100-level cardinality check. Default `FALSE`.

- seed:

  Random seed for jitter. Default `11`.

## Value

A ggplot2 / patchwork object (or list when `combine = FALSE`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic violin plot
FeaturePlot3(seu, features = c("CD3D", "CD8A"), group.by = "cell_type")

# Box plot with statistical comparisons
FeaturePlot3(seu, features = "MS4A1", group.by = "cluster",
             type = "box", comparisons = list(c("B", "T")))

# Features on x-axis (one panel per cell type)
FeaturePlot3(seu, features = c("CD3D", "CD8A", "MS4A1"),
             group.by = "cell_type", plot.by = "feature")

# Stacked panels with background bands
FeaturePlot3(seu,
  features = c("Gene1", "Gene2", "Gene3"),
  group.by = "SubCellType", bg.by = "CellType",
  stack = TRUE, add_box = TRUE)

# Col (per-cell bar) plot
FeaturePlot3(seu, features = "CD3D", group.by = "cell_type", type = "col")
} # }
```
