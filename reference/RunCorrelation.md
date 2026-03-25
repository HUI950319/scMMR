# Correlation Ranking of One Variable Against Many

Compute correlation (and p-value) between a target variable and multiple
other variables. The result is a data.frame directly compatible with
[`PlotRankScatter`](https://hui950319.github.io/scMMR/reference/PlotRankScatter.md).

## Usage

``` r
RunCorrelation(
  object,
  target,
  features = NULL,
  group.by = NULL,
  method = c("spearman", "pearson", "kendall"),
  assay = NULL,
  layer = "data",
  n_features = 2000L,
  min.pct = 0.05,
  min.cells = 20L,
  cl = 1L,
  padjust.method = "fdr",
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object or a data.frame.

- target:

  Character. Name of the target variable (gene or column).

- features:

  Character vector. Variables to correlate with `target`. For Seurat
  objects, `NULL` (default) uses VariableFeatures (up to `n_features`).
  For data.frames, `NULL` uses all numeric columns except `target`.

- group.by:

  Character. Optional column for per-group correlations. When set, the
  output contains a `group` column and can be plotted as a faceted
  `PlotRankScatter`. Default: `NULL`.

- method:

  Correlation method: `"spearman"` (default), `"pearson"`, or
  `"kendall"`.

- assay:

  Character. Seurat assay to use (Seurat only). Default:
  `DefaultAssay(object)`.

- layer:

  Character. Data layer to extract (Seurat only). Default: `"data"`.

- n_features:

  Integer. Maximum number of variable features to use when
  `features = NULL` (Seurat only). Default: 2000.

- min.pct:

  Numeric (0–1). Minimum fraction of cells expressing a gene for it to
  be included (Seurat only). Default: 0.05.

- min.cells:

  Integer. Minimum number of cells per group for correlation to be
  computed. Default: 20.

- cl:

  Integer. Number of parallel workers for per-group computation (via
  `pbapply::pblapply`). Only used when `group.by` is set and `cl > 1`.
  Default: 1 (serial).

- padjust.method:

  Character. p-value adjustment method passed to
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default:
  `"fdr"`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

A data.frame with columns:

- gene:

  Feature name.

- score:

  Correlation coefficient (Spearman rho, Pearson r, or Kendall tau).

- pvalue:

  Raw p-value from
  [`cor.test()`](https://rdrr.io/r/stats/cor.test.html).

- padj:

  Adjusted p-value.

- abs_score:

  Absolute value of `score`, for ranking.

- n_cells:

  Number of cells / observations used.

- group:

  (Only when `group.by` is set) Group label.

Rows are sorted by `abs_score` (descending). The data.frame is directly
compatible with
[`PlotRankScatter`](https://hui950319.github.io/scMMR/reference/PlotRankScatter.md):
the `gene` column maps to the name column and `score` maps to the score
column.

## Details

The first argument `object` can be either a **Seurat object** or a
**data.frame**.

**Seurat input:**

- `target` may be a gene name or metadata column.

- `features`: if `NULL`, the top 2 000 highly-variable genes
  (VariableFeatures) are used; genes with expression in fewer than
  `min.pct` of cells are dropped.

- When `group.by` is set, correlations are computed separately within
  each group (e.g. per cell type).

**data.frame input:**

- `target` and `features` must be column names.

- When `group.by` is set, data is split by that column and correlations
  are computed per group.

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Seurat: correlate PTH with all HVG ---
res <- RunCorrelation(seu, target = "PTH")
PlotRankScatter(res)

# --- Seurat: per cell-type ---
res <- RunCorrelation(seu, target = "PTH", group.by = "celltype")
PlotRankScatter(res)

# --- Seurat: specific genes ---
res <- RunCorrelation(seu, target = "PTH",
               features = c("GCM2", "CASR", "VDR", "CYP27B1"))

# --- data.frame ---
df <- data.frame(x = rnorm(100), a = rnorm(100),
                 b = rnorm(100), c = rnorm(100))
res <- RunCorrelation(df, target = "x")
PlotRankScatter(res)
} # }
```
