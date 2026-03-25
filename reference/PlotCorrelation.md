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
  use.padj = TRUE,
  p.cutoff = 0.05,
  cor.cutoff = 0,
  label = "top",
  topn = 10L,
  col.pos = "#E64B35",
  col.neg = "#4DBBD5",
  col.ns = "grey70",
  size.by = c("none", "cor", "pvalue"),
  point.size = 2.5,
  size.range = c(0.5, 5),
  point.alpha = 0.7,
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

- use.padj:

  Logical. Use adjusted p-value (`padj` or `padjust`) instead of raw
  `pvalue` for the x-axis and significance threshold. Default: `TRUE`.

- p.cutoff:

  Numeric. Significance threshold on the chosen p-value column. Default:
  0.05.

- cor.cutoff:

  Numeric. Minimum absolute correlation to be considered significant.
  Default: 0 (any direction).

- label:

  Character or character vector controlling which genes to label.

  - `"top"` (default) — label `topn` genes with largest absolute
    correlation among significant hits.

  - `"sig"` — label all significant genes.

  - `"all"` — label every point.

  - `"none"` — no labels.

  - A character vector of specific gene names.

- topn:

  Integer. Number of top genes to label when `label = "top"`. Default:
  10.

- col.pos:

  Character. Color for significant positive correlations. Default:
  `"#E64B35"` (red).

- col.neg:

  Character. Color for significant negative correlations. Default:
  `"#4DBBD5"` (blue).

- col.ns:

  Character. Color for non-significant points. Default: `"grey70"`.

- size.by:

  Character. Variable to map to point size.

  - `"none"` (default) — all points use `point.size`.

  - `"cor"` — size mapped to absolute correlation value.

  - `"pvalue"` — size mapped to \\-\log\_{10}(p)\\.

- point.size:

  Numeric. Fixed point size when `size.by = "none"`, ignored otherwise.
  Default: 2.5.

- size.range:

  Numeric vector of length 2. Range of point sizes when `size.by` is not
  `"none"`. Default: `c(0.5, 5)`.

- point.alpha:

  Numeric. Point transparency. Default: 0.7.

- label.size:

  Numeric. Label text size. Default: 3.5.

- label.color:

  Character. Label text color. Default: `"black"`.

- box.padding:

  Numeric. Padding around label boxes (passed to
  `ggrepel::geom_text_repel`). Default: 0.5.

- max.overlaps:

  Integer. Maximum overlapping labels. Default: 20.

- title:

  Character. Plot title. Default: `NULL` (auto-generated from the
  `target` attribute).

- xlab:

  Character. X-axis label. Default: auto.

- ylab:

  Character. Y-axis label. Default: `NULL` (auto: `"Correlation"` or
  `"NES"` depending on input).

- ncol:

  Integer. Number of columns when a `group` column is present and
  faceting is used. Default: 3.

## Value

A `ggplot` object.

## See also

[`RunCorrelation`](https://hui950319.github.io/scMMR/reference/RunCorrelation.md),
[`PlotRankScatter`](https://hui950319.github.io/scMMR/reference/PlotRankScatter.md),
[`PlotScatter`](https://hui950319.github.io/scMMR/reference/PlotScatter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# --- RunCorrelation output ---
res <- RunCorrelation(seu, target = "PTH")
PlotCorrelation(res)
PlotCorrelation(res, label = "top", topn = 15)
PlotCorrelation(res, label = c("GCM2", "CASR", "VDR"))

# Map point size to |correlation|
PlotCorrelation(res, size.by = "cor")

# Map point size to -log10(p)
PlotCorrelation(res, size.by = "pvalue", size.range = c(1, 6))

# Per-group (RunCorrelation)
res <- RunCorrelation(seu, target = "PTH", group.by = "celltype")
PlotCorrelation(res, size.by = "cor")

# --- RunTraceGene output (rho + padjust + lineage) ---
tg <- RunTraceGene(seu, lineages = c("Lineage1", "Lineage2"))
PlotCorrelation(tg)                   # faceted by lineage
PlotCorrelation(tg, size.by = "cor")  # size = |rho|

# --- RunTraceGSEA output (Score/NES + AdjPValue + Lineage) ---
gsea <- RunTraceGSEA(seu, lineages = c("Lineage1", "Lineage2"),
                     gene.sets = gmt)
PlotCorrelation(gsea)                 # ylab auto = "NES"
PlotCorrelation(gsea, size.by = "pvalue")
} # }
```
