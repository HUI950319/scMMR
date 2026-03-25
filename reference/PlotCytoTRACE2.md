# CytoTRACE 2 Potency Plot

Violin + box plot of CytoTRACE2 potency scores grouped by cell type,
with optional statistical comparisons.

## Usage

``` r
PlotCytoTRACE2(
  object,
  group.by = NULL,
  score = "CytoTRACE2_Score",
  order = TRUE,
  colors = NULL,
  base_size = 12,
  point.size = 0.3
)
```

## Arguments

- object:

  A Seurat object with CytoTRACE2 results (from `RunCytoTRACE2`).

- group.by:

  Character. Metadata column for grouping (e.g. `"cell_type"`). Default
  uses active idents.

- score:

  Character. Which score to plot: `"CytoTRACE2_Score"` (default) or
  `"CytoTRACE2_Relative"`.

- order:

  Logical. If `TRUE` (default), reorder groups by median score (most
  differentiated to most potent).

- colors:

  Character vector of colours, or `NULL` for default.

- base_size:

  Numeric. Base font size (default 12).

- point.size:

  Numeric. Jitter point size (default 0.3).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
seu <- RunCytoTRACE2(seu, species = "human")
PlotCytoTRACE2(seu, group.by = "cell_type")
} # }
```
