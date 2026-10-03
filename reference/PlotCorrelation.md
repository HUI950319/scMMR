# Volcano-Style Correlation Plot

Visualise the output of
[`RunCorrelation`](https://hui950319.github.io/scMMR/reference/RunCorrelation.md)
or
[`RunTraceGene`](https://hui950319.github.io/scMMR/reference/RunTraceGene.md)
as a volcano-like scatter plot with `-log10(p)` on the x-axis and
correlation coefficient on the y-axis. Points are colored by
significance direction (positive / negative / non-significant), and
selected genes are labeled with `ggrepel`.

## Usage

``` r
PlotCorrelation(
  data,
  name.col = NULL,
  score.col = NULL,
  p.col = NULL,
  use.padj = TRUE,
  p.cutoff = 0.05,
  cor.cutoff = 0,
  label = "top",
  topn = 5L,
  col.pos = "#E64B35",
  col.neg = "#4DBBD5",
  col.ns = "grey70",
  size.by = "none",
  point.size = 2.5,
  size.range = c(0.5, 5),
  alpha.by = "none",
  point.alpha = 0.7,
  alpha.range = c(0.1, 1),
  label.size = 3.5,
  label.color = "black",
  box.padding = 0.5,
  max.overlaps = 20L,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  ncol = 3L
)
```

## Arguments

- data:

  A data.frame from
  [`RunCorrelation`](https://hui950319.github.io/scMMR/reference/RunCorrelation.md),
  [`RunTraceGene`](https://hui950319.github.io/scMMR/reference/RunTraceGene.md),
  or
  [`RunTraceGSEA`](https://hui950319.github.io/scMMR/reference/RunTraceGSEA.md).
  Must contain a name column (`gene` or `Pathway`), a score column
  (`score`, `rho`, `Score`/`NES`), and a p-value column
  (`pvalue`/`PValue`, `padj`/ `padjust`/`AdjPValue`). Column matching is
  case-insensitive. An optional grouping column (`group`/
  `lineage`/`Lineage`) enables faceting.

- name.col:

  Character. Column name to use as the label/name (e.g. gene or pathway
  names). Default: `NULL` (auto-detect from
  `gene`/`pathway`/`name`/`feature`, case-insensitive).

- score.col:

  Character. Column name to use as the score (y-axis). Default: `NULL`
  (auto-detect from `score`/`rho`/ `NES`).

- p.col:

  Character. Column name to use as the p-value (x-axis). Default: `NULL`
  (auto-detect from `padj`/`padjust`/ `pvalue`/`PValue`). Overrides
  `use.padj` when set.

- use.padj:

  Logical. Use adjusted p-value (`padj` or `padjust`) instead of raw
  `pvalue` for the x-axis and significance threshold. Ignored when
  `p.col` is specified. Default: `TRUE`.

- p.cutoff:

  Numeric. Significance threshold on the chosen p-value column. Default:
  0.05.

- cor.cutoff:

  Numeric. Minimum absolute correlation to be considered significant.
  Default: 0 (any direction).

- label:

  Character or character vector controlling which genes to label.

  - `"top"` (default) — label `topn` positive and `topn` negative genes
    among significant hits.

  - `"sig"` — label all significant genes.

  - `"all"` — label every point.

  - `"none"` — no labels.

  - A character vector of specific gene names.

- topn:

  Integer. Number of top genes to label per direction (positive and
  negative) when `label = "top"`. Default: 5.

- col.pos:

  Character. Color for significant positive correlations. Default:
  `"#E64B35"` (red).

- col.neg:

  Character. Color for significant negative correlations. Default:
  `"#4DBBD5"` (blue).

- col.ns:

  Character. Color for non-significant points. Default: `"grey70"`.

- size.by:

  Character. Variable to map to point size. Built-in shortcuts: `"none"`
  (default, fixed size), `"cor"` (absolute correlation), `"pvalue"`
  (\\-\log\_{10}(p)\\). Or pass any column name in `data`.

- point.size:

  Numeric. Fixed point size when `size.by = "none"`. Default: 2.5.

- size.range:

  Numeric vector of length 2. Range of point sizes when `size.by` is not
  `"none"`. Default: `c(0.5, 5)`.

- alpha.by:

  Character. Variable to map to point transparency. Built-in shortcuts:
  `"none"` (default, fixed alpha), `"cor"` (absolute correlation),
  `"pvalue"` (\\-\log\_{10}(p)\\). Or pass any column name in `data`.

- point.alpha:

  Numeric. Fixed point transparency when `alpha.by = "none"`. Default:
  0.7.

- alpha.range:

  Numeric vector of length 2. Range of alpha values when `alpha.by` is
  not `"none"`. Default: `c(0.1, 1)`.

- label.size:

  Numeric. Label text size. Default: 3.5.

- label.color:

  Character. Label text color. Default: `"black"`.

- box.padding:

  Numeric. Padding around label boxes. Default: 0.5.

- max.overlaps:

  Integer. Maximum overlapping labels. Default: 20.

- title:

  Character. Plot title. Default: `NULL` (auto).

- xlab:

  Character. X-axis label. Default: auto.

- ylab:

  Character. Y-axis label. Default: `NULL` (auto).

- ncol:

  Integer. Number of facet columns. Default: 3.

## Value

A `ggplot` object.

## See also

[`RunCorrelation`](https://hui950319.github.io/scMMR/reference/RunCorrelation.md),
[`PlotRankScatter`](https://hui950319.github.io/scMMR/reference/PlotRankScatter.md),
[`PlotScatter`](https://hui950319.github.io/scMMR/reference/PlotScatter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- RunCorrelation(seu, target = "PTH")
PlotCorrelation(res)                              # top 5 pos + 5 neg
PlotCorrelation(res, topn = 10)                   # top 10 pos + 10 neg
PlotCorrelation(res, label = c("GCM2", "CASR"))   # specific genes

# Custom score and p-value columns
PlotCorrelation(res, score.col = "prop_cor", p.col = "prop_padj")

# Map size and alpha to custom columns
PlotCorrelation(res, size.by = "prop_cor", alpha.by = "de_abs_logFC")

# Per-group
res <- RunCorrelation(seu, target = "PTH", group.by = "celltype")
PlotCorrelation(res, size.by = "cor", alpha.by = "cor")

# RunTraceGene output
tg <- RunTraceGene(seu, lineages = c("Lineage1", "Lineage2"))
PlotCorrelation(tg, size.by = "cor")

# RunTraceGSEA output
gsea <- RunTraceGSEA(seu, lineages = c("Lineage1"), gene.sets = gmt)
PlotCorrelation(gsea, size.by = "pvalue", alpha.by = "pvalue")
} # }
```
