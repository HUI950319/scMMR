# Rank Scatter Plot for Feature Scoring

A general-purpose scatter plot that ranks features (genes, regulons,
pathways, etc.) by their scores and highlights the top-ranked ones with
coloured points and text labels. Inspired by the SCENIC Regulon
Specificity Score (RSS) visualisation.

## Usage

``` r
PlotRankScatter(
  data,
  cell.type = NULL,
  topn = 5L,
  max.show = 200L,
  highlight.color = "#007D9B",
  base.color = "#BECEE3",
  label.size = 4,
  point.size = 3,
  title = NULL,
  ylab = "Score",
  base_size = 12,
  ncol = 3L,
  clean.names = TRUE
)
```

## Arguments

- data:

  Input data in one of the following formats:

  Named numeric vector

  :   Values are scores; names are feature labels. Produces a
      single-panel plot.

  Matrix

  :   Rows = features, columns = groups (e.g. cell types). Use
      `cell.type` to select one or several columns. If
      `cell.type = NULL` all columns are shown as a faceted plot.
      Typical inputs: RSS matrix from SCENIC, or
      `DNN_predict()$pathway_scores`.

  data.frame

  :   Must contain a name column (`gene`, `name`, `regulon`, `pathway`,
      or `feature`) and a score column (`importance`, `score`, `NES`, or
      `value`). If a grouping column (`cell_type`, `group`, or
      `cluster`) is present the plot is faceted.

- cell.type:

  Character vector. For matrix input: column name(s) to plot. `NULL`
  (default) = plot all columns. For vector/data.frame input: used only
  as the plot title.

- topn:

  Integer. Number of top-ranked features to highlight (default 5).

- max.show:

  Integer. Maximum number of features to display (default 200).

- highlight.color:

  Colour for the top-ranked points (default `"#007D9B"`).

- base.color:

  Colour for the remaining points (default `"#BECEE3"`).

- label.size:

  Numeric. Text label size (default 4).

- point.size:

  Numeric. Point size (default 3).

- title:

  Character. Plot title. If `NULL` (default) the group name is used.

- ylab:

  Character. Y-axis label (default `"Score"`).

- base_size:

  Numeric. Base font size (default 12).

- ncol:

  Integer. Number of columns for faceted layout (default 3).

- clean.names:

  Logical. If `TRUE` (default), strip common prefixes (`HALLMARK_`,
  `KEGG_`, etc.), SCENIC suffixes (`(+)`, `(-)`), and replace
  underscores with spaces.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# 1) SCENIC RSS matrix - single cell type
PlotRankScatter(rssMat, cell.type = "Parathyroid", topn = 10)

# 2) SCENIC RSS matrix - multiple cell types
PlotRankScatter(rssMat, cell.type = c("Parathyroid", "Stromal"),
                topn = 5, ncol = 2)

# 3) Named numeric vector
scores <- setNames(runif(100), paste0("Gene", 1:100))
PlotRankScatter(scores, topn = 10, ylab = "Activity")

# 4) DNN_predict importance (data.frame)
PlotRankScatter(pred$imp_global, topn = 15)

# 5) Per-class importance - faceted
PlotRankScatter(pred$imp_per_class, topn = 10, ncol = 4)
} # }
```
