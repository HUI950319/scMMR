# Find Pseudotime-Associated Genes Along Trajectories

Identify genes whose expression significantly changes along pseudotime
trajectories. Supports two testing methods: fast Spearman correlation
and more rigorous GAM (Generalized Additive Model) fitting.

## Usage

``` r
RunTraceGene(
  seu,
  lineages,
  features = NULL,
  n_candidates = 2000,
  min.pct = 0.05,
  test.method = c("spearman", "gam"),
  assay = NULL,
  layer = "data",
  smooth_k = 10,
  family = "gaussian",
  padjust.method = "fdr",
  verbose = TRUE
)
```

## Arguments

- seu:

  A Seurat object with pseudotime stored in `meta.data`.

- lineages:

  Character vector of column names in `seu@meta.data` containing
  pseudotime values. Cells not on a lineage should have `NA`. Multiple
  lineages are analyzed independently.

- features:

  Character vector of gene names to test. Default: `NULL` (use top
  `n_candidates` variable features). Set to `"all"` to test all genes.

- n_candidates:

  Number of variable features to select when `features = NULL`. Default:
  2000.

- min.pct:

  Minimum fraction of lineage cells expressing a gene (expression \> 0).
  Default: 0.05.

- test.method:

  Testing method: `"spearman"` (fast) or `"gam"` (GAM smooth, more
  rigorous). Default: `"spearman"`.

- assay:

  Seurat assay to use. Default: `DefaultAssay(seu)`.

- layer:

  Data layer to extract. Default: `"data"`.

- smooth_k:

  Number of basis dimensions for GAM smooth term (`mgcv::s(k = ...)`).
  Only used when `test.method = "gam"`. Smaller values = smoother
  curves. Default: 10.

- family:

  Distribution family for GAM. Default: `"gaussian"`. Use `"nb"` for raw
  count data (requires `layer = "counts"`).

- padjust.method:

  Method for p-value adjustment. Default: `"fdr"`.

- verbose:

  Logical; print progress messages. Default: `TRUE`.

## Value

A data.frame with one row per gene-lineage combination:

- `gene`:

  Gene name.

- `lineage`:

  Lineage name.

- `rho`:

  Spearman correlation coefficient (positive = expression increases
  along pseudotime).

- `pvalue`:

  Raw p-value (Spearman test or GAM F-test).

- `padjust`:

  Adjusted p-value.

- `r_sq`:

  R-squared (GAM only; NA for Spearman).

- `dev_expl`:

  Deviance explained (GAM only; NA for Spearman).

- `peaktime`:

  Pseudotime of peak expression (median of top 1% fitted/raw values).

- `valleytime`:

  Pseudotime of valley expression (median of bottom 1% fitted/raw
  values).

- `exp_ncells`:

  Number of cells expressing the gene.

- `direction`:

  `"up"` (rho \> 0) or `"down"` (rho \< 0).

- `test_method`:

  Testing method used.

## Details

For each lineage, the function:

1.  Subsets cells with finite pseudotime
    ([`is.finite()`](https://rdrr.io/r/base/is.finite.html)).

2.  Selects candidate genes (HVG by default, or all / user-specified).

3.  Filters genes by minimum expression percentage (`min.pct`).

4.  Tests each gene for association with pseudotime.

5.  Computes effect sizes, peak/valley times, and adjusted p-values.

Two test methods are available:

- `"spearman"`:

  Computes Spearman correlation between gene expression and pseudotime.
  Fast (vectorized rho + per-gene p-value). Output includes `rho` and
  `pvalue`.

- `"gam"`:

  Fits a GAM smooth spline via
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) and tests
  whether the smooth term is significant (F-test). Output includes
  `r_sq` (R-squared), `dev_expl` (deviance explained), and `pvalue`.
  Additionally computes Spearman rho for direction. Requires mgcv.

## See also

[`RunTraceGSEA`](https://hui950319.github.io/scMMR/reference/RunTraceGSEA.md),
[`PlotDynamicFeatures`](https://hui950319.github.io/scMMR/reference/PlotDynamicFeatures.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)

# Fast Spearman correlation (default)
res <- RunTraceGene(seu, lineages = c("Lineage1", "Lineage2"))
sig <- res[res$padjust < 0.05, ]
head(sig[order(sig$padjust), ], 20)

# GAM fitting (more rigorous)
res_gam <- RunTraceGene(seu, lineages = "Lineage1",
                         test.method = "gam")
sig_gam <- res_gam[res_gam$padjust < 0.05, ]

# Visualize top up-regulated genes
up_genes <- head(sig$gene[sig$direction == "up"], 10)
PlotDynamicFeatures(seu, pseudotime = "Lineage1", features = up_genes)

# Test all genes
res_all <- RunTraceGene(seu, lineages = "Lineage1", features = "all")
} # }
```
