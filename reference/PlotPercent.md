# Differential Abundance Beeswarm Plot

Visualise the output of
[`RankPercent()`](https://hui950319.github.io/scMMR/reference/RankPercent.md):
a miloR-style beeswarm plot showing per-neighborhood logFC for each cell
type. Significant neighborhoods (FDR \< threshold) are coloured by their
logFC value on a diverging blue-white-red gradient; non-significant
neighborhoods are shown in grey.

## Usage

``` r
PlotPercent(
  da_result,
  fdr_threshold = 0.1,
  colors = c("blue", "white", "red"),
  na_color = "grey80",
  show_boxplot = FALSE,
  base_size = 12,
  point_size = 1.5
)
```

## Arguments

- da_result:

  Output list from
  [`RankPercent()`](https://hui950319.github.io/scMMR/reference/RankPercent.md).

- fdr_threshold:

  Numeric. FDR threshold for colouring significant neighborhoods
  (default 0.1). Neighborhoods with `p_adj >= fdr_threshold` are shown
  in `na_color`.

- colors:

  Character vector of length 3: colours for the diverging gradient
  `c(low, mid, high)`. Default `c("blue", "white", "red")`.

- na_color:

  Colour for non-significant neighborhoods. Default `"grey80"`.

- show_boxplot:

  Logical. If `TRUE`, overlay a transparent boxplot behind the beeswarm
  points to show the distribution summary for each cell type. Default
  `FALSE`.

- base_size:

  Base font size (default 12).

- point_size:

  Point size (default 1.5).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
da <- RankPercent(pred$shared_embedding, q1@meta.data,
                  conditions = c("PH", "SH"))
PlotPercent(da)
PlotPercent(da, fdr_threshold = 0.05)
} # }
```
