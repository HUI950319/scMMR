# Cell-Type Proportion vs Gene Expression Correlation

For each sample, compute the proportion of a given cell type and the
average expression of one or more genes, then draw a scatter-correlation
plot. Internally delegates to
[`PlotScatter`](https://hui950319.github.io/scMMR/reference/PlotScatter.md),
so all its features (grouping, marginal plots, regression lines, etc.)
are available.

## Usage

``` r
PlotPropCorrelation(
  seu,
  gene,
  celltype,
  celltype.col,
  sample.col,
  expr.in = NULL,
  group.by = NULL,
  assay = NULL,
  layer = "data",
  method = c("spearman", "pearson", "kendall"),
  show.cor = TRUE,
  show.smooth = TRUE,
  smooth.method = "lm",
  point.size = 3,
  point.alpha = 0.7,
  point.color = "#984ea3",
  cor.size = 4,
  marginal = "none",
  marginal.size = 5,
  palette = "Paired",
  palcolor = NULL,
  title = NULL,
  ncol = 3L,
  return.data = FALSE
)
```

## Arguments

- seu:

  A Seurat object.

- gene:

  Character vector. Gene(s) whose mean expression per sample is plotted
  on the y-axis. When length \> 1 the plot is faceted.

- celltype:

  Character vector. Cell type(s) whose proportion per sample is plotted
  on the x-axis. When length \> 1 the plot is faceted.

- celltype.col:

  Character. Column in `meta.data` containing cell type annotations.

- sample.col:

  Character. Column in `meta.data` containing sample / patient
  identifiers.

- expr.in:

  Character vector or `NULL`. Cell types in which to average gene
  expression. `NULL` (default) = all cells in each sample. Set to e.g.
  `"Parathyroid"` to average only within that cell type.

- group.by:

  Character or `NULL`. Metadata column for coloring points (e.g.
  `"disease"`). Each sample must have a unique group label. Default:
  `NULL`.

- assay:

  Character. Seurat assay. Default: `DefaultAssay(seu)`.

- layer:

  Character. Data layer. Default: `"data"`.

- method:

  Correlation method: `"spearman"`, `"pearson"`, or `"kendall"`.
  Default: `"spearman"`.

- show.cor:

  Logical. Show correlation statistics. Default: `TRUE`.

- show.smooth:

  Logical. Show regression line. Default: `TRUE`.

- smooth.method:

  Smoothing method for the trend line. Default: `"lm"`.

- point.size:

  Numeric. Point size. Default: 3.

- point.alpha:

  Numeric. Point transparency. Default: 0.7.

- point.color:

  Character. Point color when `group.by = NULL`. Default: `"#984ea3"`.

- cor.size:

  Numeric. Font size for correlation text. Default: 4.

- marginal:

  Character. Marginal plot type passed to
  [`PlotScatter`](https://hui950319.github.io/scMMR/reference/PlotScatter.md):
  `"none"`, `"density"`, `"histogram"`, etc. Default: `"none"`.

- marginal.size:

  Numeric. Relative size of marginal plots. Default: 5.

- palette:

  Character. Palette name. Default: `"Paired"`.

- palcolor:

  Character vector or `NULL`. Custom colors.

- title:

  Character or `NULL`. Plot title.

- ncol:

  Integer. Facet columns. Default: 3.

- return.data:

  Logical. If `TRUE`, return the sample-level data.frame instead of the
  plot. Default: `FALSE`.

## Value

A `ggplot` / `patchwork` object, or a data.frame if
`return.data = TRUE`.

## See also

[`PlotScatter`](https://hui950319.github.io/scMMR/reference/PlotScatter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic: Parathyroid proportion vs PTH expression
PlotPropCorrelation(seu, gene = "PTH", celltype = "Parathyroid",
                    celltype.col = "celltype", sample.col = "sample_id")

# Color by disease group
PlotPropCorrelation(seu, gene = "PTH", celltype = "Parathyroid",
                    celltype.col = "celltype", sample.col = "sample_id",
                    group.by = "disease")

# Multiple cell types
PlotPropCorrelation(seu, gene = "PTH",
                    celltype = c("Parathyroid", "Stromal"),
                    celltype.col = "celltype", sample.col = "sample_id")

# Multiple genes
PlotPropCorrelation(seu, gene = c("PTH", "GCM2"),
                    celltype = "Parathyroid",
                    celltype.col = "celltype", sample.col = "sample_id",
                    marginal = "density")

# Return data for custom plotting
df <- PlotPropCorrelation(seu, gene = "PTH", celltype = "Parathyroid",
                          celltype.col = "celltype",
                          sample.col = "sample_id",
                          return.data = TRUE)
} # }
```
