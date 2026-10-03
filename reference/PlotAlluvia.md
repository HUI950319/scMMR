# Sankey Plot for Cell Population Composition

Create a sankey (alluvial flow) plot showing cell type composition
changes across conditions/groups, based on the ggalluvial package. Each
column represents a group, strata represent cell types, and flows
connect the same cell type across adjacent groups. Parameters mirror
`PlotArea` (internal), with an additional `curve.type` argument.

## Usage

``` r
PlotAlluvia(
  cellmeta,
  by = "group",
  fill = "celltype",
  palette = "lancet",
  palcolor = NULL,
  flow.alpha = 0.4,
  bar.width = 0.5,
  bar.color = "gray50",
  show.label = FALSE,
  label.size = 3,
  label.color = "black",
  label.box = FALSE,
  label.fill = "white",
  label.alpha = 1,
  show.pct = FALSE,
  label.style = c("default", "all", "n", "pct", "name"),
  label.min.pct = 0,
  y.type = c("percent", "count"),
  legend.ncol = 1,
  base.size = 15,
  curve.type = "sigmoid"
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

  Character string specifying the fill variable (e.g.
  `"cell_type_pred"`).

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

- show.label:

  Logical, whether to show cell type labels inside strata. Default
  `FALSE`.

- label.size:

  Numeric, font size for stratum labels. Default 3.

- label.color:

  Stratum label colour. Default `"black"`.

- label.box:

  Logical. If `TRUE`, draw labels with white background and border (like
  `geom_label`). Default `FALSE` (plain text).

- label.fill:

  Background fill for boxed labels. Default `"white"`.

- label.alpha:

  Transparency of boxed label background. Default `1`.

- show.pct:

  Logical, whether to append within-group percentage to labels. Default
  `FALSE`. Only used when `show.label = TRUE`.

- label.style:

  Character. Label format style:

  - `"default"` — uses `show.label`/`show.pct` flags.

  - `"all"` — name + count + percentage, e.g. `"T cells (50, 8.1%)"`.

  - `"n"` — name + count, e.g. `"T cells (50)"`.

  - `"pct"` — name + percentage, e.g. `"T cells 8.1%"`.

  - `"name"` — name only.

- label.min.pct:

  Numeric (0-1). Minimum within-group proportion to display a label.
  Cell types below this threshold are unlabelled. Default 0 (show all).
  E.g. `0.05` hides labels for cell types \< 5%.

- y.type:

  Character. Y-axis scale type:

  - `"percent"` (default) — all bars equal height, y-axis shows
    percentages.

  - `"count"` — bars reflect raw cell counts.

- legend.ncol:

  Integer, number of legend columns. Default 1.

- base.size:

  Numeric, base font size for the theme. Default 15.

- curve.type:

  Character, flow curve type. Default `"sigmoid"`. Other options:
  `"linear"`, `"cubic"`, etc.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ToyData)

# --- Basic usage (defaults: by="group", fill="celltype", lancet palette) ---
PlotAlluvia(toy_test)

# --- From data.frame ---
PlotAlluvia(toy_test@meta.data, by = "group", fill = "celltype")

# --- Label styles ---
# Name + percentage (hide small cell types < 5%)
PlotAlluvia(toy_test, label.style = "pct", label.min.pct = 0.05)

# Name + count + percentage
PlotAlluvia(toy_test, label.style = "all", label.min.pct = 0.05)

# Name only
PlotAlluvia(toy_test, label.style = "name", label.min.pct = 0.03)

# --- Boxed labels (white background + border) ---
PlotAlluvia(toy_test, label.style = "pct", label.min.pct = 0.05,
           label.box = TRUE)

# Semi-transparent boxed labels
PlotAlluvia(toy_test, label.style = "pct", label.min.pct = 0.05,
           label.box = TRUE, label.alpha = 0.7)

# --- Custom colours ---
# Use a different palette
PlotAlluvia(toy_test, palette = "Paired")

# Custom colour vector
PlotAlluvia(toy_test, palcolor = UtilsR::pal_paraSC)

# --- Appearance ---
# Sigmoid flow (default) vs linear
PlotAlluvia(toy_test, curve.type = "linear")

# Adjust flow transparency and bar width
PlotAlluvia(toy_test, flow.alpha = 0.2, bar.width = 0.3)

# Raw counts instead of percentages
PlotAlluvia(toy_test, y.type = "count", label.style = "n",
           label.min.pct = 0)

# Multi-column legend
PlotAlluvia(toy_test, legend.ncol = 2)
} # }
```
