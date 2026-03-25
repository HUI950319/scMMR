# Scatter Correlation Plot

Create a scatter plot showing the correlation between two variables
(genes and/or metadata columns) in single-cell data. Supports coloring
by groups, faceting, per-group or global correlation statistics, and
multiple correlation methods.

## Usage

``` r
PlotScatter(
  object,
  var1,
  var2,
  group.by = NULL,
  split.by = NULL,
  cells = NULL,
  assay = NULL,
  layer = "data",
  method = c("spearman", "pearson", "kendall"),
  smooth.method = "lm",
  show.cor = TRUE,
  show.smooth = TRUE,
  show.rug = FALSE,
  cor.digits = 3,
  cor.size = 4,
  point.size = 1,
  point.alpha = 0.6,
  smooth.size = 1,
  smooth.color = "#fdc086",
  show.se = TRUE,
  ncol = 3,
  palette = "Paired",
  palcolor = NULL,
  point.color = "#984ea3",
  rug.color = "#7fc97f",
  title = NULL,
  marginal = c("none", "density", "histogram", "boxplot", "violin", "densigram"),
  marginal.size = 5,
  raster = NULL,
  raster.dpi = 300,
  ...
)
```

## Arguments

- object:

  A Seurat object or a data.frame containing the variables to plot.

- var1:

  Character. First variable (x-axis): a gene name, metadata column, or
  data.frame column name.

- var2:

  Character. Second variable (y-axis): a gene name, metadata column, or
  data.frame column name.

- group.by:

  Character. Optional column to color points and compute per-group
  correlations. Default: `NULL` (single color).

- split.by:

  Character. Optional column for faceting (`facet_wrap`). Default:
  `NULL`.

- cells:

  Character vector. Cell barcodes to include (Seurat only). Default:
  `NULL` (all cells).

- assay:

  Character. Seurat assay to use for gene expression (Seurat only).
  Default: `DefaultAssay(object)`.

- layer:

  Character. Data layer to extract (Seurat only). Default: `"data"`
  (log-normalized).

- method:

  Correlation method for `ggpubr::stat_cor`: `"spearman"`, `"pearson"`,
  or `"kendall"`. Default: `"spearman"`.

- smooth.method:

  Smoothing method for `geom_smooth`. Default: `"lm"`. Set to `"loess"`,
  `"gam"`, etc. as needed.

- show.cor:

  Logical. Show correlation statistics on the plot. Default: `TRUE`.

- show.smooth:

  Logical. Show regression / smooth line. Default: `TRUE`.

- show.rug:

  Logical. Show rug plots on the axes. Default: `FALSE`.

- cor.digits:

  Integer. Number of decimal digits for correlation display. Default: 3.

- cor.size:

  Numeric. Font size for correlation text. Default: 4.

- point.size:

  Numeric. Size of scatter points. Default: 1.

- point.alpha:

  Numeric. Transparency of scatter points (0–1). Default: 0.6.

- smooth.size:

  Numeric. Width of the regression line. Default: 1.

- smooth.color:

  Character. Color of the smooth line when `group.by = NULL`. Default:
  `"#fdc086"`. Ignored when `group.by` is set (line color follows
  group).

- show.se:

  Logical. Show confidence interval around smooth line. Default: `TRUE`.

- ncol:

  Integer. Number of columns for faceting when `split.by` is set.
  Default: 3.

- palette:

  Character. Color palette name passed to
  [`palette_colors`](https://hui950319.github.io/scMMR/reference/palette_colors.md).
  Default: `"Paired"`.

- palcolor:

  Character vector. Custom colors overriding `palette`. Default: `NULL`.

- point.color:

  Character. Fixed point color when `group.by = NULL`. Default:
  `"#984ea3"`.

- rug.color:

  Character. Color for rug marks. Default: `"#7fc97f"`.

- title:

  Character. Plot title. Default: `NULL` (auto-generated).

- marginal:

  Character. Type of marginal distribution plot to add on the x- and
  y-axes: `"none"` (default), `"density"`, `"histogram"`, `"boxplot"`,
  `"violin"`, or `"densigram"` (histogram + density overlay). Requires
  patchwork. Ignored when `split.by` is used.

- marginal.size:

  Numeric. Relative size of marginal plots compared to the main scatter
  panel. Default: 5 (i.e. the main panel is 5 times larger than the
  marginal).

- raster:

  Logical. If `TRUE`, rasterize the point layer via
  `ggrastr::rasterise()` to reduce file size for large datasets.
  Default: `NULL` (auto: `TRUE` when \> 50 000 cells).

- raster.dpi:

  Integer. DPI for rasterized points. Default: 300.

- ...:

  Additional arguments passed to `geom_point`.

## Value

A `ggplot` object (or a `patchwork` object when `marginal != "none"`).

## Details

The first argument `object` can be either a **Seurat object** or a
**data.frame** (/ tibble).

**When a Seurat object is provided**, variables `var1` and `var2` can
be:

- Gene names (expression extracted from the specified assay/layer).

- Metadata column names in `object@meta.data`.

- Any combination of the two.

**When a data.frame is provided**, `var1`, `var2`, `group.by`, and
`split.by` must all be column names of the data.frame. The
Seurat-specific parameters (`assay`, `layer`, `cells`) are ignored.

When `group.by` is set, points are colored by the grouping variable and
correlation statistics (`ggpubr::stat_cor`) are shown per group. When
`split.by` is set, the plot is faceted by that variable.

When `marginal` is set to a plot type other than `"none"`, marginal
distribution plots are added to the x- and y-axes via `patchwork`
layout. Supported types include density, histogram, boxplot, violin, and
densigram (histogram + density overlay). When `group.by` is set, the
marginal plots are colored/filled by group. Note: marginal plots are
incompatible with faceting (`split.by`); if both are specified, marginal
plots are silently skipped.

**Note:** When marginal plots are added, the return value is a
`patchwork` object. You can still use `patchwork::&` or
`patchwork::plot_annotation()` for further modifications.

## See also

[`palette_colors`](https://hui950319.github.io/scMMR/reference/palette_colors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)

# --- Seurat object input ---
# Two genes
PlotScatter(seu, var1 = "METTL3", var2 = "SETD2")

# Gene vs metadata
PlotScatter(seu, var1 = "nFeature_RNA", var2 = "PTH",
            method = "pearson")

# Color by cell type
PlotScatter(seu, var1 = "TP53", var2 = "MDM2",
            group.by = "celltype", palette = "npg")

# With marginal boxplot
PlotScatter(seu, var1 = "TP53", var2 = "MDM2",
            group.by = "celltype", marginal = "boxplot")

# --- data.frame input ---
df <- data.frame(x = rnorm(200), y = rnorm(200),
                 grp = sample(c("A", "B"), 200, replace = TRUE))
PlotScatter(df, var1 = "x", var2 = "y")
PlotScatter(df, var1 = "x", var2 = "y", group.by = "grp",
            marginal = "density")
} # }
```
