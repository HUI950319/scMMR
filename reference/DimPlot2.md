# Enhanced 2D Dimensional Reduction Plot

Visualize cell groups on a 2-dimensional reduction plot (UMAP, t-SNE,
etc.). Supports Seurat objects and plain data.frames. Provides label,
highlight, density, mark, and raster layers with flexible theming.

Ported from `scop::CellDimPlot` with code optimizations. Uses `UtilsR`
theme / grob utilities and
[`scMMR::palette_colors`](https://hui950319.github.io/scMMR/reference/palette_colors.md)
for color mapping.

## Usage

``` r
DimPlot2(
  seu,
  group.by,
  reduction = NULL,
  dims = c(1, 2),
  split.by = NULL,
  cells = NULL,
  show_na = FALSE,
  show_stat = NULL,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Paired",
  palcolor = NULL,
  bg_color = "grey80",
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
  cells.highlight = NULL,
  cols.highlight = "black",
  sizes.highlight = 1,
  alpha.highlight = 1,
  stroke.highlight = 0.5,
  add_density = FALSE,
  density_color = "grey80",
  density_filled = FALSE,
  density_filled_palette = "Greys",
  density_filled_palcolor = NULL,
  add_mark = FALSE,
  mark_type = c("hull", "ellipse", "rect", "circle"),
  mark_expand = grid::unit(3, "mm"),
  mark_alpha = 0.1,
  mark_linetype = 1,
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
  hex = FALSE,
  hex.count = TRUE,
  hex.bins = 50,
  hex.binwidth = NULL,
  hex.linewidth = 0.5,
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

  A Seurat object or a `data.frame`. When a data.frame is provided, the
  first two columns specified by `dims` are treated as embedding
  coordinates, and `group.by` / `split.by` columns must be present.

- group.by:

  Character. Name(s) of column(s) in metadata to group (color) cells by.

- reduction:

  Character. Name of the dimensionality reduction to use. Ignored when
  `seu` is a data.frame. Default auto-detects UMAP/tSNE.

- dims:

  Integer vector of length 2. Dimensions to plot. For Seurat input:
  component indices (default `c(1, 2)`). For data.frame input: column
  indices for x/y coordinates.

- split.by:

  Character. Column name to split into facets. Default `NULL`.

- cells:

  Character vector of cell/row names to include. Default `NULL` (all).

- show_na:

  Logical. Show NA group with `bg_color`. Default `FALSE`.

- show_stat:

  Logical. Show cell counts in labels/subtitles. Default `NULL` (auto:
  `FALSE` when theme is blank, `TRUE` otherwise).

- pt.size:

  Point size. Default auto-calculated.

- pt.alpha:

  Point transparency. Default `1`.

- palette:

  Character. Color palette name (passed to `palette_colors`). Default
  `"Paired"`.

- palcolor:

  Character vector. Custom colors. Default `NULL`.

- bg_color:

  Background color for NA points. Default `"grey80"`.

- label:

  Logical. Add group labels. Default `FALSE`.

- label.size:

  Label font size. Default `4`.

- label.fg:

  Label foreground color. Default `"white"`.

- label.bg:

  Label background color. Default `"black"`.

- label.bg.r:

  Label background ratio. Default `0.1`.

- label_insitu:

  Logical. Place raw group names at cell centers. Default `FALSE`.

- label_repel:

  Logical. Repel labels. Default `FALSE`.

- label_repulsion:

  Repulsion force. Default `20`.

- label_point_size:

  Center point size for labels. Default `1`.

- label_point_color:

  Center point color. Default `"black"`.

- label_segment_color:

  Segment color for repelled labels. Default `"black"`.

- cells.highlight:

  Character vector of cell names to highlight, or `TRUE` for all non-NA
  cells. Default `NULL`.

- cols.highlight:

  Highlight border color. Default `"black"`.

- sizes.highlight:

  Highlight point size. Default `1`.

- alpha.highlight:

  Highlight transparency. Default `1`.

- stroke.highlight:

  Highlight border width. Default `0.5`.

- add_density:

  Logical. Add density contours. Default `FALSE`.

- density_color:

  Contour line color. Default `"grey80"`.

- density_filled:

  Logical. Use filled contour bands. Default `FALSE`.

- density_filled_palette:

  Palette for filled contours. Default `"Greys"`.

- density_filled_palcolor:

  Custom colors for filled contours. Default `NULL`.

- add_mark:

  Logical. Add group boundary marks. Default `FALSE`.

- mark_type:

  Mark shape: `"hull"`, `"ellipse"`, `"rect"`, or `"circle"`. Default
  `"hull"`.

- mark_expand:

  Mark expansion. Default `grid::unit(3, "mm")`.

- mark_alpha:

  Mark transparency. Default `0.1`.

- mark_linetype:

  Mark border line type. Default `1`.

- lineages:

  Character vector of metadata column names containing pseudotime values
  (e.g. from `RunSlingshot`). Each column is one lineage; `NULL` = no
  lineage overlay. Default `NULL`.

- lineages_trim:

  Numeric pair. Quantile range of pseudotime cells to include in LOESS
  fit (removes extremes). Default `c(0.01, 0.99)`.

- lineages_span:

  LOESS smoother span. Larger = smoother. Default `0.75`.

- lineages_palette:

  Colour palette for lineage lines. Default `"Dark2"`.

- lineages_palcolor:

  Custom colour vector for lineages. Default `NULL`.

- lineages_arrow:

  Arrow specification for path end. Default
  `grid::arrow(length = grid::unit(0.1, "inches"))`.

- lineages_linewidth:

  Foreground path linewidth. Default `1`.

- lineages_line_bg:

  Colour of the outline stroke behind each path. Default `"white"`.

- lineages_line_bg_stroke:

  Extra width added to the outline. Default `0.5`.

- lineages_whiskers:

  Logical. Draw segments from each cell to the smooth curve. Default
  `FALSE`.

- lineages_whiskers_linewidth:

  Whisker segment linewidth. Default `0.5`.

- lineages_whiskers_alpha:

  Whisker transparency. Default `0.5`.

- hex:

  Logical. Use hexagonal binning. Default `FALSE`.

- hex.count:

  Logical. Map hex alpha to count. Default `TRUE`.

- hex.bins:

  Number of hex bins. Default `50`.

- hex.binwidth:

  Hex bin width. Default `NULL`.

- hex.linewidth:

  Hex border width. Default `0.5`.

- raster:

  Logical. Rasterize points (auto if \>100k cells). Default `NULL`.

- raster_method:

  Character. Rasterisation backend: `"rasterise"` (default, faithful
  colours via `rasterise_layer`) or `"scattermore"` (faster C-level
  rendering, may shift colours).

- raster.dpi:

  Numeric. Raster resolution in DPI. Accepts a scalar or a two-length
  vector (max is used). Default `c(512, 512)`.

- aspect.ratio:

  Aspect ratio. Default `1`.

- title:

  Plot title. Default `NULL`.

- subtitle:

  Plot subtitle. Default `NULL`.

- xlab:

  X-axis label. Default auto from reduction key.

- ylab:

  Y-axis label. Default auto from reduction key.

- legend.position:

  Legend position. Default `"right"`.

- legend.direction:

  Legend direction. Default `"vertical"`.

- legend.title:

  Legend title. Default `NULL` (uses group name).

- theme_use:

  Theme specification. Accepts a function (e.g. `UtilsR::theme_blank`),
  a character name (e.g. `"theme_sc"`), or a `theme` object (e.g.
  `theme_bw()`). Default `NULL` which uses `UtilsR::theme_blank`.

- theme_args:

  List of extra arguments to `theme_use`. Default
  [`list()`](https://rdrr.io/r/base/list.html).

- combine:

  Logical. Combine multiple panels via patchwork. Default `TRUE`.

- nrow, ncol:

  Layout dimensions for combined plot. Default `NULL`.

- byrow:

  Arrange panels by row. Default `TRUE`.

- force:

  Logical. Force plot even with \>100 group levels. Default `FALSE`.

- seed:

  Random seed. Default `11`.

## Value

A ggplot object (or list of ggplot objects if `combine = FALSE`).

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Seurat input ---
DimPlot2(seu, group.by = "celltype")
DimPlot2(seu, group.by = "celltype", label = TRUE, add_mark = TRUE)

# --- data.frame input ---
df <- data.frame(UMAP1 = rnorm(200), UMAP2 = rnorm(200),
                 cluster = sample(letters[1:5], 200, replace = TRUE))
DimPlot2(df, group.by = "cluster", dims = c(1, 2))
} # }
```
