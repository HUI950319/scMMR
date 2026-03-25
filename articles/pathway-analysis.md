# Pathway Enrichment Analysis

## Overview

scMMR supports two methods for per-cell-type pathway enrichment
analysis:

- **GSEA** (Gene Set Enrichment Analysis) via `clusterProfiler` – tests
  for coordinated differential expression of gene sets
- **SCPA** (Single Cell Pathway Analysis) – tests for changes in
  multivariate gene set distributions between conditions

Both are accessed through the unified
[`RunPathwayAnalysis()`](https://hui950319.github.io/scMMR/reference/RunPathwayAnalysis.md)
interface.

## GSEA Method

### Basic Usage

``` r
library(scMMR)
library(Seurat)

# Load data with cell types and treatment groups
seu <- qs::qread("path/to/seurat.qs")

# Prepare gene sets
gmt_path <- system.file("extdata/gmt/h.all.v2022.1.Hs.symbols.gmt",
                         package = "scMMR")
gene_sets <- read_gmt(gmt_path)

# Run GSEA across all cell types
gsea_results <- RunPathwayAnalysis(
  seurat_obj   = seu,
  gene.sets    = gene_sets,
  method       = "GSEA",
  celltype.col = "cell_type",
  group.col    = "group",
  ident.1      = "Treatment",
  ident.2      = "Control",
  pvalue.cutoff = 0.05,
  clean.names  = TRUE,
  cores        = 4
)
```

### GSEA Output

``` r
head(gsea_results)
#>   Celltype                Pathway   NES     PValue   p.adjust Sign
#> 1 T_cell   Tnfa Signaling Via Nfkb  2.31  0.0001   0.003     Up
#> 2 T_cell   Interferon Gamma        1.98  0.0003   0.008     Up
#> 3 Macro    Inflammatory Response    2.15  0.0002   0.005     Up
```

## SCPA Method

SCPA tests whether the multivariate distribution of a gene set differs
between two conditions. It is particularly powerful for detecting
coordinated changes that may not show up in individual gene-level tests.

``` r
scpa_results <- RunPathwayAnalysis(
  seurat_obj   = seu,
  gene.sets    = gene_sets,
  method       = "SCPA",
  celltype.col = "cell_type",
  group.col    = "group",
  ident.1      = "Treatment",
  ident.2      = "Control",
  cores        = 4
)
```

### SCPA Output

``` r
head(scpa_results)
#>   Celltype          Pathway   FC     PValue  qvalue  Sign
#> 1 T_cell   Tnfa Signaling    3.21  0.001    0.02    Up
#> 2 T_cell   Il2 Stat5         2.85  0.003    0.04    Up
```

## Visualization with PlotPathwayBubble

``` r
# Default bubble plot
PlotPathwayBubble(gsea_results)

# Customized
PlotPathwayBubble(
  gsea_results,
  top.n        = 10,
  size.by      = "pvalue",
  pvalue.cutoff = 0.05,
  select.by    = "frequency",
  facet.by.sign = TRUE
)
```

## Advanced: GSEA on DE Results

For more control, you can first run differential expression, then
perform GSEA on the results:

``` r
# Step 1: Differential expression
de_results <- RunDE(
  seurat_obj   = seu,
  celltype.col = "cell_type",
  ident.1      = "Treatment",
  ident.2      = "Control"
)

# Step 2: GSEA on logFC rankings
gsea_on_de <- RunGsea(
  de_results   = de_results,
  gene.sets    = gene_sets,
  pvalue.cutoff = 0.05,
  clean.names  = TRUE
)

# Visualize
PlotGsea(gsea_on_de, top.n = 8, facet_by = "celltype")
```

## Cross-Enrichment Overlap

Find shared genes between TF regulons and enriched pathways:

``` r
overlap <- CrossEnrichOverlap(
  tf_results      = gsea_tf_results,
  pathway_results = gsea_pathway_results,
  min_overlap     = 3,
  method          = "fisher"
)
```

## Using Different GMT Databases

``` r
# Hallmark pathways (50 sets, broad biological processes)
gs_hallmark <- read_gmt(system.file("extdata/gmt/h.all.v2022.1.Hs.symbols.gmt",
                                     package = "scMMR"))

# KEGG pathways (186 sets, metabolic and signaling)
gs_kegg <- read_gmt(system.file("extdata/gmt/c2.cp.kegg.v2022.1.Hs.symbols.gmt",
                                 package = "scMMR"))

# PROGENy signaling (14 pathways, curated signaling)
gs_progeny <- read_gmt(system.file("extdata/gmt/progeny.human.top500.gmt",
                                    package = "scMMR"))

# CollecTRI regulons (TF-target gene sets)
gs_collectri <- read_gmt(system.file("extdata/gmt/collectri.human.gmt",
                                      package = "scMMR"))
```

## GSEA vs SCPA: When to Use Which?

| Feature          | GSEA                               | SCPA                                              |
|------------------|------------------------------------|---------------------------------------------------|
| Statistical test | Kolmogorov-Smirnov on ranked genes | Multivariate distribution test                    |
| Input            | Ranked gene list (logFC)           | Expression matrix per condition                   |
| Detects          | Coordinate up/down-regulation      | Distribution shifts (mean, variance, correlation) |
| Speed            | Fast                               | Moderate                                          |
| Best for         | Classic pathway enrichment         | Subtle regulatory changes                         |

## Tips

- Use `clean.names = TRUE` to remove MSigDB prefixes (e.g., “HALLMARK\_”
  -\> “Hallmark”)
- Set `cores > 1` for parallel processing across cell types
- Use Hallmark gene sets for initial exploration, then drill down with
  KEGG/Reactome
- SCPA requires sufficient cells per condition per cell type (~50+)
