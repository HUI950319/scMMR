# Sankey Plot for Cell Population Composition

Create a sankey (alluvial flow) plot showing cell type composition
changes across conditions/groups, based on the ggalluvial package. Each
column represents a group, strata represent cell types, and flows
connect the same cell type across adjacent groups. Parameters mirror
[`PlotAlluvia`](https://hui950319.github.io/scMMR/reference/PlotAlluvia.md),
with an additional `curve.type` argument.

## Usage

``` r
PlotSankey(
  cellmeta,
  by,
  fill,
  colors = NULL,
  palette = NULL,
  flow.alpha = 0.4,
  bar.width = 0.5,
  bar.color = "gray50",
  show.label = FALSE,
  label.size = 3,
  label.color = "black",
  show.pct = FALSE,
  y.percent = TRUE,
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

- colors:

  Optional character vector of colours. If `NULL` and `palette` is also
  `NULL`, default ggplot2 colours are used. A named vector maps specific
  levels to colours.

- palette:

  Optional RColorBrewer palette name (e.g. `"Paired"`, `"Set2"`) for
  auto-generating colours when `colors = NULL`. Default `NULL` (use
  ggplot2 defaults).

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

- show.pct:

  Logical, whether to append within-group percentage to labels. Default
  `FALSE`. Only used when `show.label = TRUE`.

- y.percent:

  Logical. If `TRUE` (default), y-axis shows percentages (all bars equal
  height). If `FALSE`, y-axis shows raw cell counts.

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
# From Seurat object directly
PlotSankey(seurat_obj, by = "group", fill = "celltype")

# From data.frame
PlotSankey(seurat_obj@meta.data, by = "group", fill = "cell_type_pred")
} # }
```
