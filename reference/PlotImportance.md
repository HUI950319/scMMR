# Importance Plot (Lollipop / Bar)

Visualise gene importance scores or pathway scores returned by
[`DNN_predict()`](https://hui950319.github.io/scMMR/reference/DNN_predict.md).

## Usage

``` r
PlotImportance(
  importance,
  top_k = 15L,
  display = c("lollipop", "bar"),
  palette = "Reds 3",
  ncol = 3,
  base_size = 12,
  group_by = NULL
)
```

## Arguments

- importance:

  Either:

  - A **data.frame** with gene importance scores:

    - **Global**: columns `gene`, `importance` (i.e. `pred$imp_global`).

    - **Per-class**: columns `cell_type`, `n_cells`, `rank`, `gene`,
      `importance` (i.e. `pred$imp_per_class`).

  - A **matrix** of pathway scores (i.e. `pred$pathway_scores`), with
    cells as rows and pathways as columns.

- top_k:

  Integer. Maximum number of top features to show. For per-class /
  per-group mode this is applied within each group. Default 15.

- display:

  Character, either `"lollipop"` (default, stem + dot) or `"bar"`
  (filled bar chart).

- palette:

  Character. An HCL sequential palette name recognised by
  `colorspace::sequential_hcl()`. Default `"Reds 3"`.

- ncol:

  Integer. Number of columns for `facet_wrap` in per-class / per-group
  mode. Default 3.

- base_size:

  Numeric. Base font size. Default 12.

- group_by:

  Optional character or factor vector of length `nrow(importance)` for
  per-group faceted pathway plots. Only used when `importance` is a
  matrix. If `NULL` (default), shows global mean pathway scores.

## Value

A ggplot object.

## Details

**Gene importance** (data.frame input): Automatically detects whether
the input is global importance or per-class importance based on the
presence of a `cell_type` column.

**Pathway scores** (matrix input): Accepts the `pathway_scores` matrix
from `DNN_predict(pathway_gmt = ...)`. Long pathway names are
automatically cleaned (database prefixes removed, underscores replaced,
sentence case applied).

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- DNN_predict(query, model_path, explain = TRUE,
                    pathway_gmt = "reactome.gmt")

# Gene importance (existing usage)
PlotImportance(pred$imp_global)
PlotImportance(pred$imp_per_class, top_k = 10, ncol = 4)

# Pathway scores: global (mean across all cells)
PlotImportance(pred$pathway_scores, top_k = 20)

# Pathway scores: per-group
PlotImportance(pred$pathway_scores, top_k = 10,
               group_by = q1$cell_type_pred, ncol = 4)

# Bar style
PlotImportance(pred$imp_global, display = "bar", palette = "Blues 3")
} # }
```
