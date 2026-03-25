# Plot Dynamic Features Along Pseudotime

Visualise gene expression or metadata features as a function of
pseudotime (or any continuous trajectory variable). Supports three
fitting methods with user-controllable smoothness:

- **GAM** (`mgcv`): control smoothness via `smooth_k` (smaller =
  smoother).

- **loess**: control smoothness via `loess_span` (larger = smoother).

- **B-spline**: control smoothness via `bspline_knot` (more knots = more
  flexible / less smooth).

## Usage

``` r
PlotDynamicFeatures(
  srt,
  pseudotime,
  features,
  group.by = NULL,
  assay = NULL,
  layer = "counts",
  fit_method = c("gam", "loess", "bspline"),
  smooth_k = 10,
  loess_span = 0.75,
  bspline_knot = 3,
  family = NULL,
  exp_method = c("log1p", "raw", "zscore"),
  lib_normalize = (layer == "counts"),
  add_point = TRUE,
  add_line = TRUE,
  add_interval = TRUE,
  add_rug = TRUE,
  pt.size = 1,
  line.size = 1,
  point_palette = "Paired",
  point_palcolor = NULL,
  line_palette = "Dark2",
  line_palcolor = NULL,
  compare_lineages = TRUE,
  compare_features = FALSE,
  ncol = NULL,
  nrow = NULL,
  reverse = FALSE,
  flip = FALSE,
  seed = 11,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- pseudotime:

  Character vector of column name(s) in `srt@meta.data` containing
  pseudotime values. Multiple names enable multi-lineage comparison.

- features:

  Character vector of gene names (looked up in the assay) and/or column
  names in `srt@meta.data`.

- group.by:

  Optional character string. A column in `srt@meta.data` used to colour
  scatter points.

- assay:

  Character. Seurat assay to use. Default: active assay.

- layer:

  Character. Data layer to use. Default `"counts"`.

- fit_method:

  Character, one of `"gam"`, `"loess"`, or `"bspline"`. Default `"gam"`.

- smooth_k:

  Integer. Number of basis dimensions for the GAM smooth term
  (`mgcv::s(k = ...)`). Smaller values produce smoother curves. Default
  10.

- loess_span:

  Numeric (0, 1\]. The `span` parameter for
  [`stats::loess()`](https://rdrr.io/r/stats/loess.html). Larger values
  produce smoother curves. Default 0.75.

- bspline_knot:

  Integer. Number of internal knots for B-spline fitting. More knots
  allow a more flexible (less smooth) curve. Default 3.

- family:

  Character or family object for GAM. If `NULL` (default), `"gaussian"`
  is used.

- exp_method:

  Character, one of `"log1p"`, `"raw"`, or `"zscore"`. Transformation
  applied to expression values before plotting. Default `"log1p"`.

- lib_normalize:

  Logical. Whether to normalise by library size. Default `TRUE` when
  `layer = "counts"`.

- add_point:

  Logical. Show scatter points. Default `TRUE`.

- add_line:

  Logical. Show fitted curve. Default `TRUE`.

- add_interval:

  Logical. Show 95% confidence interval ribbon. Default `TRUE`.

- add_rug:

  Logical. Show rug marks along x-axis. Default `TRUE`.

- pt.size:

  Numeric. Point size. Default 1.

- line.size:

  Numeric. Line width. Default 1.

- point_palette:

  Character. Palette name for scatter points / rug coloured by
  `group.by`. See
  [`show_palettes`](https://hui950319.github.io/scMMR/reference/show_palettes.md)
  for available names. Default `"Paired"`.

- point_palcolor:

  Optional character vector of custom colours for the `group.by`
  variable. Overrides `point_palette`.

- line_palette:

  Character. Palette name for fitted lines / ribbons. Default `"Dark2"`.

- line_palcolor:

  Optional character vector of custom colours for lineages or features.
  Overrides `line_palette`.

- compare_lineages:

  Logical. When multiple `pseudotime` columns are given, overlay them in
  a single panel per feature. Default `TRUE`.

- compare_features:

  Logical. When multiple `features` are given, overlay them in a single
  panel per lineage. Default `FALSE`.

- ncol:

  Integer. Number of columns in the facet layout.

- nrow:

  Integer. Number of rows in the facet layout.

- reverse:

  Logical. Reverse the pseudotime axis. Default `FALSE`.

- flip:

  Logical. Flip x and y axes. Default `FALSE`.

- seed:

  Integer. Random seed. Default 11.

- ...:

  Additional arguments passed to the fitting function
  ([`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html),
  [`stats::loess`](https://rdrr.io/r/stats/loess.html), or
  [`stats::lm.fit`](https://rdrr.io/r/stats/lmfit.html)).

## Value

A `ggplot` object (or list of ggplot objects if both `compare_lineages`
and `compare_features` are `FALSE`).

## Examples

``` r
if (FALSE) { # \dontrun{
# After running Slingshot or any trajectory method:
# Default GAM fitting
PlotDynamicFeatures(srt, pseudotime = "Lineage1",
                    features = c("Gene1", "Gene2"))

# Smoother curve (smaller k)
PlotDynamicFeatures(srt, pseudotime = "Lineage1",
                    features = "Gene1", smooth_k = 5)

# More flexible curve (larger k)
PlotDynamicFeatures(srt, pseudotime = "Lineage1",
                    features = "Gene1", smooth_k = 20)

# Use loess with custom span
PlotDynamicFeatures(srt, pseudotime = "Lineage1",
                    features = "Gene1",
                    fit_method = "loess", loess_span = 0.3)

# B-spline with 5 knots
PlotDynamicFeatures(srt, pseudotime = "Lineage1",
                    features = "Gene1",
                    fit_method = "bspline", bspline_knot = 5)

# Compare two lineages
PlotDynamicFeatures(srt,
                    pseudotime = c("Lineage1", "Lineage2"),
                    features = c("Gene1", "Gene2"),
                    compare_lineages = TRUE)
} # }
```
