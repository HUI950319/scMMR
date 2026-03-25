# Differential Expression Analysis

## Overview

scMMR provides a streamlined workflow for differential expression (DE)
analysis across cell types, followed by GSEA enrichment and correlation
analysis with embedding features.

## Step 1: Differential Expression (RunDE)

[`RunDE()`](https://hui950319.github.io/scMMR/reference/RunDE.md) wraps
Seurat’s `FindMarkers` to run DE across all cell types between two
conditions:

``` r
library(scMMR)
library(Seurat)

de_results <- RunDE(
  seurat_obj    = seu,
  celltype.col  = "cell_type",
  ident.1       = "Treatment",
  ident.2       = "Control",
  test          = "wilcox",
  logfc.threshold = 0.25,
  min.pct       = 0.1
)
```

### Output

``` r
head(de_results)
#>   gene    avg_log2FC  p_val_adj  pct.1  pct.2  cell_type
#> 1 CXCL2   2.15        1.2e-15    0.82   0.15   T_cell
#> 2 JUNB    1.89        3.4e-12    0.75   0.20   T_cell
```

### Supported Tests

| Test        | Description                                    |
|-------------|------------------------------------------------|
| `"wilcox"`  | Wilcoxon rank-sum (default, non-parametric)    |
| `"t"`       | Student’s t-test                               |
| `"DESeq2"`  | DESeq2 negative binomial (requires raw counts) |
| `"poisson"` | Poisson generalized linear model               |

## Step 2: Visualize DE Results (PlotDE)

### Dot Plot

Shows top genes per cell type with size = -log10(p-value) and color =
logFC:

``` r
PlotDE(de_results, type = "dot", top_n = 5)

PlotDE(de_results, type = "dot", top_n = 10,
       top_by = "logfc", logfc.cutoff = 0.5)
```

### Volcano Plot

Per-cell-type volcano plots:

``` r
PlotDE(de_results, type = "volcano",
       logfc.cutoff = 1.0, pval.cutoff = 0.05,
       label = TRUE)
```

## Step 3: GSEA on DE Results (RunGsea)

Perform GSEA using the logFC rankings from DE:

``` r
gsea_results <- RunGsea(
  de_results   = de_results,
  gene.sets    = read_gmt(system.file("extdata/gmt/h.all.v2022.1.Hs.symbols.gmt",
                                       package = "scMMR")),
  pvalue.cutoff = 0.05,
  clean.names  = TRUE
)
```

### Visualize GSEA Results

``` r
# Bubble plot (cell type x pathway)
PlotGsea(gsea_results, top.n = 8)

# Faceted by enrichment direction
PlotGsea(gsea_results, top.n = 10, facet_by = "sign")
```

## Step 4: Correlation Analysis

### RunCorrelation

Correlate gene expression or pathway scores with DNN embedding axes:

``` r
cor_results <- RunCorrelation(
  seurat_obj = seu,
  features   = c("CXCL2", "JUNB", "IL7R", "GZMB"),
  lineage.col = "cell_type",
  cores       = 4
)
```

### PlotCorrelation

Volcano-style plot of correlations:

``` r
PlotCorrelation(cor_results, p.cutoff = 0.05, cor.cutoff = 0.3,
                label = TRUE, label.n = 10)
```

### PlotPropCorrelation

Correlate embedding dimensions with cell type proportions or module
scores:

``` r
PlotPropCorrelation(
  seu,
  emb.cols  = paste0("DNN_", 1:10),
  prop.data = proportions_df,
  method    = "spearman"
)
```

## Complete Workflow Example

``` r
# 1. DE analysis
de <- RunDE(seu, celltype.col = "cell_type",
            ident.1 = "Treatment", ident.2 = "Control")

# 2. Visualize top DE genes
PlotDE(de, type = "dot", top_n = 8)

# 3. GSEA enrichment
gsea <- RunGsea(de, gene.sets = read_gmt(hallmark_gmt))
PlotGsea(gsea, top.n = 10)

# 4. Correlate key genes with embedding
cor_res <- RunCorrelation(seu, features = de$gene[1:50])
PlotCorrelation(cor_res, label = TRUE)
```

## Tips

- Use `test = "wilcox"` for most cases (robust, non-parametric)
- Switch to `DESeq2` for small cell numbers per group
- Filter DE results by both logFC and adjusted p-value before GSEA
- `RunGsea` requires at least ~50 DE genes per cell type for meaningful
  results
