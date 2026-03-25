# Perturbation Ranking Plot

Visualise the output of
[`RankPerturbation()`](https://hui950319.github.io/scMMR/reference/RankPerturbation.md):
a ranked lollipop or bar chart of cell types ordered by perturbation
score, with significance markers. Reuses the same visual style as
[`PlotImportance()`](https://hui950319.github.io/scMMR/reference/PlotImportance.md).

## Usage

``` r
PlotPerturbation(
  rank_result,
  top_k = NULL,
  display = c("lollipop", "bar"),
  sig_threshold = 0.05,
  palette = "Reds 3",
  base_size = 12
)
```

## Arguments

- rank_result:

  Output list from
  [`RankPerturbation()`](https://hui950319.github.io/scMMR/reference/RankPerturbation.md).

- top_k:

  Integer. Show only top-k cell types (default `NULL`, show all).

- display:

  Character: `"lollipop"` (default) or `"bar"`.

- sig_threshold:

  Numeric. FDR threshold for significance stars (default 0.05).

- palette:

  HCL sequential palette name (default `"Reds 3"`).

- base_size:

  Base font size (default 12).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- RankPerturbation(pred$shared_embedding, q1@meta.data,
                        conditions = c("PH", "SH"))
PlotPerturbation(res)
PlotPerturbation(res, top_k = 10, display = "bar")
} # }
```
