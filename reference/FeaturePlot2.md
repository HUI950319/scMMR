# Feature Expression on 2-D Reduction Plot

Visualise feature values (gene expression, metadata scores, reduction
embeddings) on a 2-D reduction plot. Accepts Seurat objects or plain
data.frames. Key options include:

- **Multiple features / split panels** with unified colour scale
  (`keep_scale`).

- **Co-expression** geometric mean across gene features
  (`calculate_coexp = TRUE`).

- **Compare mode** – blend all features per cell into a single colour
  (`compare_features = TRUE`).

## Usage

``` r
FeaturePlot2(
  seu,
  features,
  reduction = NULL,
  dims = c(1, 2),
  split.by = NULL,
  cells = NULL,
  layer = "data",
  assay = NULL,
  show_stat = NULL,
  palette = NULL,
  palcolor = NULL,
  pt.size = NULL,
  pt.alpha = 1,
  bg_cutoff = 0,
  bg_color = "grey80",
  keep_scale = "feature",
  lower_quantile = 0,
  upper_quantile = 0.99,
  lower_cutoff = NULL,
  upper_cutoff = NULL,
  add_density = FALSE,
  density_color = "grey80",
  density_filled = FALSE,
  density_filled_palette = "Greys",
  density_filled_palcolor = NULL,
  cells.highlight = NULL,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  calculate_coexp = FALSE,
  compare_features = FALSE,
  color_blend_mode = c("blend", "average", "screen", "multiply"),
  label = FALSE,
  label.size = 4,
  label.fg = "white",
  label.bg = "black",
  label.bg.r = 0.1,
  label_insitu = FALSE,
  label_repel = FALSE,
  label_repulsion = 20,
  label_point_size = 1,
  label_point_color = "black",
  label_segment_color = "black",
  lineages = NULL,
  lineages_trim = c(0.01, 0.99),
  lineages_span = 0.75,
  lineages_palette = "Dark2",
  lineages_palcolor = NULL,
  lineages_arrow = grid::arrow(length = grid::unit(0.1, "inches")),
  lineages_linewidth = 1,
  lineages_line_bg = "white",
  lineages_line_bg_stroke = 0.5,
  lineages_whiskers = FALSE,
  lineages_whiskers_linewidth = 0.5,
  lineages_whiskers_alpha = 0.5,
  raster = NULL,
  raster_method = c("scattermore", "rasterise"),
  raster.dpi = c(512, 512),
  aspect.ratio = 1,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
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

  A Seurat object or a `data.frame`. For data.frames the columns
  specified by `dims` are used as x/y coordinates; all feature columns
  must also be present in the data.frame.

- features:

  Character vector (or named character vector) of features to plot.
  Features can be gene names (Seurat assay rows), metadata columns,
  reduction embedding columns, or data.frame columns. Named features use
  the name as subtitle. A list is also accepted and is flattened.

- reduction:

  Character. Reduction to use for coordinates. `NULL` auto-detects
  UMAP/tSNE. Ignored for data.frame input.

- dims:

  Integer vector of length 2. Dimension indices within the reduction
  (Seurat) or column indices (data.frame). Default `c(1,2)`.

- split.by:

  Column name to split into facets. `NULL` = no split.

- cells:

  Character vector of cell/row names to include. `NULL` = all.

- layer:

  Seurat assay layer to use. Default `"data"`.

- assay:

  Seurat assay name. `NULL` = default assay.

- show_stat:

  Logical. Show positive-cell statistics in subtitle. `NULL` = auto
  (`FALSE` for blank theme).

- palette:

  Colour palette name for the gradient. Default `"Spectral"` (standard)
  or `"Set1"` (compare mode).

- palcolor:

  Custom colour vector. `NULL` = use palette.

- pt.size:

  Point size. `NULL` = auto.

- pt.alpha:

  Point transparency. Default `1`.

- bg_cutoff:

  Numeric. Feature values ≤ cutoff are treated as background (shown in
  `bg_color`). Default `0`.

- bg_color:

  Background colour for low/NA values. Default `"grey80"`.

- keep_scale:

  How to align colour scales across panels:

  - `NULL` – each panel scaled independently to its own range.

  - `"feature"` (default) – panels for the same feature share a scale.

  - `"all"` – all panels share a single global scale.

- lower_quantile, upper_quantile:

  Quantile bounds for colour scale. Default `0` / `0.99`.

- lower_cutoff, upper_cutoff:

  Explicit value bounds (override quantile). `NULL` = use quantiles.

- add_density:

  Logical. Overlay 2-D density contour. Default `FALSE`.

- density_color:

  Contour colour. Default `"grey80"`.

- density_filled:

  Logical. Filled density raster. Default `FALSE`.

- density_filled_palette:

  Palette for filled density. Default `"Greys"`.

- density_filled_palcolor:

  Custom colours for filled density. `NULL`.

- cells.highlight:

  Cells to highlight (names or `TRUE` for all expressed cells). `NULL` =
  none.

- cols.highlight:

  Border colour of highlighted cells. Default `"black"`.

- sizes.highlight:

  Size of highlighted points. Default `1`.

- alpha.highlight:

  Alpha of highlighted points. Default `1`.

- stroke.highlight:

  Border width around highlighted points. Default `0.5`.

- calculate_coexp:

  Logical. Replace gene features with their geometric mean co-expression
  value. Default `FALSE`.

- compare_features:

  Logical. Blend all features per cell into one colour. Default `FALSE`.

- color_blend_mode:

  Blend mode when `compare_features = TRUE`. One of `"blend"`,
  `"average"`, `"screen"`, `"multiply"`.

- label:

  Logical. Label the high-expression region. Default `FALSE`.

- label.size:

  Label text size. Default `4`.

- label.fg:

  Label foreground colour. Default `"white"`.

- label.bg:

  Label background colour. Default `"black"`.

- label.bg.r:

  Label background ratio. Default `0.1`.

- label_insitu:

  Logical. Use feature names as labels (insitu) instead of numbers.
  Default `FALSE`.

- label_repel:

  Logical. Use `ggrepel` for labels. Default `FALSE`.

- label_repulsion:

  Repulsion force. Default `20`.

- label_point_size:

  Center point size for repelled labels. Default `1`.

- label_point_color:

  Center point colour. Default `"black"`.

- label_segment_color:

  Segment colour. Default `"black"`.

- lineages:

  Character vector of metadata column names containing pseudotime values
  (e.g. from `RunSlingshot`). `NULL` = none.

- lineages_trim:

  Quantile range for LOESS fit. Default `c(0.01,0.99)`.

- lineages_span:

  LOESS span. Default `0.75`.

- lineages_palette:

  Palette for lineage colours. Default `"Dark2"`.

- lineages_palcolor:

  Custom colour vector. Default `NULL`.

- lineages_arrow:

  Arrow for path end. Default
  `grid::arrow(length = grid::unit(0.1, "inches"))`.

- lineages_linewidth:

  Foreground linewidth. Default `1`.

- lineages_line_bg:

  Outline colour. Default `"white"`.

- lineages_line_bg_stroke:

  Outline extra width. Default `0.5`.

- lineages_whiskers:

  Draw cell-to-curve segments. Default `FALSE`.

- lineages_whiskers_linewidth:

  Whisker linewidth. Default `0.5`.

- lineages_whiskers_alpha:

  Whisker alpha. Default `0.5`.

- raster:

  Logical. Rasterise points. `NULL` = auto (\> 100k).

- raster_method:

  Character. Rasterisation backend: `"rasterise"` (default, faithful
  colours via `rasterise_layer`) or `"scattermore"` (faster C-level
  rendering, may shift colours).

- raster.dpi:

  Numeric scalar or 2-length vector. Default `c(512,512)`.

- aspect.ratio:

  Aspect ratio. Default `1`.

- title:

  Plot title. `NULL` = none.

- subtitle:

  Subtitle(s). `NULL` = auto from feature names. Single value or vector
  of `length(features)`.

- xlab, ylab:

  Axis labels. `NULL` = auto.

- legend.position:

  Legend position. Default `"right"`.

- legend.direction:

  Legend direction. Default `"vertical"`.

- legend.title:

  Legend title. `NULL` = feature name.

- theme_use:

  Theme specification (function, character, or `theme` object). `NULL` →
  `UtilsR::theme_blank`.

- theme_args:

  Extra arguments for `theme_use`. Default
  [`list()`](https://rdrr.io/r/base/list.html).

- combine:

  Logical. Combine panels with patchwork. Default `TRUE`.

- nrow, ncol:

  Layout dimensions. `NULL` = auto.

- byrow:

  Logical. Fill by row. Default `TRUE`.

- force:

  Logical. Force plot when \> 50 features. Default `FALSE`.

- seed:

  Random seed. Default `11`.

## Value

A ggplot or patchwork object (or list when `combine = FALSE`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Seurat input
FeaturePlot2(seu, features = "CD3D")
FeaturePlot2(seu, features = c("CD3D", "CD8A"), keep_scale = "all")
FeaturePlot2(seu, features = c("CD3D", "CD8A"),
             compare_features = TRUE, color_blend_mode = "blend")

# data.frame input
df <- data.frame(UMAP1 = rnorm(500), UMAP2 = rnorm(500),
                 GeneA = abs(rnorm(500)), GeneB = abs(rnorm(500)))
FeaturePlot2(df, features = "GeneA", dims = c(1, 2))
} # }
```
