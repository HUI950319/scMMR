# Alluvial Plot with Gaps Between Cell Types

Create an alluvial plot showing cell type composition changes across
conditions/groups, with visible gaps between strata (cell types). Based
on ggalluvial, identical to
[`PlotAlluvia`](https://hui950319.github.io/scMMR/reference/PlotAlluvia.md)
but with spacing between cell type strata for cleaner visual separation.
Gaps are implemented by inserting invisible (transparent) strata between
real cell types.

## Usage

``` r
PlotAlluvia2(
  cellmeta,
  by = "group",
  fill = "celltype",
  palette = "lancet",
  palcolor = NULL,
  flow.alpha = 0.4,
  bar.width = 0.5,
  bar.color = "gray50",
  gap = 0.01,
  show.label = FALSE,
  label.size = 3,
  label.color = "black",
  label.box = TRUE,
  label.fill = "white",
  label.alpha = 1,
  show.pct = TRUE,
  label.style = c("all", "default", "n", "pct", "name"),
  label.min.pct = 0,
  y.type = c("percent", "count"),
  legend.ncol = 1,
  base.size = 15,
  curve.type = "linear"
)
```

## Arguments

- cellmeta:

  A Seurat object or data.frame containing cell metadata.

- by:

  Character string specifying the grouping variable. Default `"group"`.

- fill:

  Character string specifying the cell type variable. Default
  `"celltype"`.

- palette:

  Character. Colour palette name passed to
  [`palette_colors`](https://hui950319.github.io/scMMR/reference/palette_colors.md).
  Default `"lancet"`.

- palcolor:

  Character vector. Custom colours overriding `palette`. Default `NULL`.

- flow.alpha:

  Numeric (0-1), transparency of alluvial flows. Default 0.4.

- bar.width:

  Numeric (0-1), width of stratum bars. Default 0.5.

- bar.color:

  Stratum bar border colour. Default `"gray50"`.

- gap:

  Numeric. Size of gap between strata as a fraction of total height.
  Default 0.005. Set 0 for no gaps (same as PlotAlluvia).

- show.label:

  Logical. Default `FALSE`.

- label.size:

  Numeric. Default 3.

- label.color:

  Character. Default `"black"`.

- label.box:

  Logical. If `TRUE`, boxed labels. Default `FALSE`.

- label.fill:

  Background fill for boxed labels. Default `"white"`.

- label.alpha:

  Transparency of boxed label background. Default `1`.

- show.pct:

  Logical. Default `FALSE`.

- label.style:

  Character. Label format: `"default"`, `"all"`, `"n"`, `"pct"`, or
  `"name"`.

- label.min.pct:

  Numeric (0-1). Minimum proportion to show label. Default 0.

- y.type:

  Character. `"percent"` (default) or `"count"`.

- legend.ncol:

  Integer. Default 1.

- base.size:

  Numeric. Default 15.

- curve.type:

  Character. Flow curve type. Default `"sigmoid"`. Options: `"linear"`,
  `"cubic"`, `"quintic"`, `"sine"`, `"arctangent"`, `"sigmoid"`,
  `"xspline"`.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ToyData)

# --- Basic usage (with default gaps) ---
PlotAlluvia2(toy_test)

# --- Adjust gap size ---
PlotAlluvia2(toy_test, gap = 0.01)   # larger gaps
PlotAlluvia2(toy_test, gap = 0)      # no gaps (same as PlotAlluvia)

# --- Labels ---
PlotAlluvia2(toy_test, label.style = "pct", label.min.pct = 0.05)
PlotAlluvia2(toy_test, label.style = "pct", label.min.pct = 0.05,
             label.box = TRUE)

# --- Custom colours ---
PlotAlluvia2(toy_test, palcolor = UtilsR::pal_paraSC)
PlotAlluvia2(toy_test, palette = "Paired")

# --- Raw counts ---
PlotAlluvia2(toy_test, y.type = "count", label.style = "n")
} # }
```
