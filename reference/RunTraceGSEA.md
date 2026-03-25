# Run Trajectory GSEA Along Pseudotime

Perform Gene Set Enrichment Analysis (GSEA) along pseudotime
trajectories. For each lineage, computes Spearman correlation between
gene expression and pseudotime, uses the correlation-based ranking to
run `clusterProfiler::GSEA()`, and identifies pathways activated or
repressed along the trajectory.

## Usage

``` r
RunTraceGSEA(
  seu,
  lineages,
  gene.sets,
  features = NULL,
  n_candidates = 2000,
  min.pct = 0.1,
  ranking.metric = c("rho", "signed_logp"),
  assay = NULL,
  layer = "data",
  pvalue.cutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  clean.names = TRUE,
  prefix = NULL,
  verbose = TRUE
)
```

## Arguments

- seu:

  A Seurat object with pseudotime stored in `meta.data`.

- lineages:

  Character vector of column names in `seu@meta.data` containing
  pseudotime values. Cells not on a lineage should have `NA` or `Inf`.
  Multiple lineages are analyzed independently.

- gene.sets:

  Gene sets in one of three formats:

  - A named list of character vectors (gene names).

  - A data.frame with columns term and gene.

  - A file path to a GMT file.

  Parsed internally by
  [`parse_gene_sets`](https://hui950319.github.io/scMMR/reference/parse_gene_sets.md).

- features:

  Character vector of gene names to include in the correlation analysis.
  Default: `NULL` (use all genes in the assay, which is recommended for
  GSEA).

- n_candidates:

  Not used when `features = NULL` (all genes are used by default). Only
  applies when you want to restrict to top variable features by setting
  `features = "HVG"`. Default: 2000.

- min.pct:

  Minimum fraction of lineage cells expressing a gene (expression \> 0)
  for the gene to be included. Default: 0.1.

- ranking.metric:

  Ranking method for GSEA: `"rho"` (Spearman rho, fast) or
  `"signed_logp"` (`-log10(p) * sign(rho)`, more discriminative).
  Default: `"rho"`.

- assay:

  Seurat assay to use. Default: `DefaultAssay(seu)`.

- layer:

  Data layer to extract. Default: `"data"` (log-normalized expression
  recommended for correlation analysis).

- pvalue.cutoff:

  P-value cutoff for `clusterProfiler::GSEA()`. Default: 1 (return all
  results; filter post-hoc).

- minGSSize:

  Minimum number of genes in a gene set. Default: 10.

- maxGSSize:

  Maximum number of genes in a gene set. Default: 500.

- eps:

  Boundary correction for GSEA p-value calculation. Default: 0 (exact
  p-values).

- clean.names:

  Logical; if `TRUE`, remove common database prefixes (e.g. `HALLMARK_`,
  `KEGG_`) and convert to Title Case. Default: `TRUE`.

- prefix:

  Custom regex prefix to remove when `clean.names = TRUE`. Default:
  `NULL` (auto-detect common prefixes).

- verbose:

  Logical; print progress messages. Default: `TRUE`.

## Value

A data.frame compatible with
[`PlotPathwayBubble`](https://hui950319.github.io/scMMR/reference/PlotPathwayBubble.md),
with columns:

- `Celltype`:

  Lineage name (for `PlotPathwayBubble` x-axis compatibility).

- `Lineage`:

  Lineage name.

- `Pathway`:

  Pathway name (cleaned if `clean.names = TRUE`).

- `Score`:

  NES (Normalized Enrichment Score). Positive = pathway activated along
  trajectory; negative = repressed.

- `PValue`:

  Raw p-value.

- `AdjPValue`:

  BH-adjusted p-value.

- `QValue`:

  Q-value.

- `Sign`:

  `"Activated"` (NES \> 0) or `"Repressed"` (NES \< 0).

- `SetSize`:

  Number of genes in set.

- `EnrichmentScore`:

  Raw enrichment score.

- `CoreEnrichment`:

  Core enrichment genes (slash-separated).

- `RankingMetric`:

  Ranking method used (`"rho"` or `"signed_logp"`).

- `NCells`:

  Number of cells on this lineage.

- `Method`:

  `"GSEA"`.

## Details

The workflow for each lineage is:

1.  Subset cells with finite pseudotime values
    ([`is.finite()`](https://rdrr.io/r/base/is.finite.html)).

2.  Use all genes by default (recommended for GSEA); or use HVG if
    `features = "HVG"`.

3.  Filter genes by minimum expression percentage (`min.pct`).

4.  Compute Spearman correlation (`cor.test(method = "spearman")`)
    between each gene's expression and pseudotime.

5.  Build a gene ranking from the correlation statistics.

6.  Run `clusterProfiler::GSEA()` with the ranking.

Two ranking metrics are available:

- `"rho"`:

  Uses Spearman rho directly. Genes positively correlated with
  pseudotime rank high (expression increases along trajectory);
  negatively correlated rank low. This mode uses vectorized
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) for speed.

- `"signed_logp"`:

  Uses `-log10(p) * sign(rho)` for sharper separation of highly
  significant correlations. This mode calls
  [`cor.test()`](https://rdrr.io/r/stats/cor.test.html) per gene and is
  slower but provides stronger discrimination.

## See also

[`PlotPathwayBubble`](https://hui950319.github.io/scMMR/reference/PlotPathwayBubble.md),
[`RunPathwayAnalysis`](https://hui950319.github.io/scMMR/reference/RunPathwayAnalysis.md),
[`PlotDynamicFeatures`](https://hui950319.github.io/scMMR/reference/PlotDynamicFeatures.md),
[`parse_gene_sets`](https://hui950319.github.io/scMMR/reference/parse_gene_sets.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)

# Hallmark gene sets
gmt <- system.file("extdata", "gmt",
  "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")

# Single lineage
res <- RunTraceGSEA(seu, lineages = "Lineage1", gene.sets = gmt)
PlotPathwayBubble(res)

# Multiple lineages with signed_logp ranking
res <- RunTraceGSEA(
  seu,
  lineages = c("Lineage1", "Lineage2", "Lineage3"),
  gene.sets = gmt,
  ranking.metric = "signed_logp"
)
PlotPathwayBubble(res, top.n = 15, facet.by.sign = TRUE)

# KEGG pathways with custom gene list
kegg_gmt <- system.file("extdata", "gmt",
  "c2.cp.kegg.v2022.1.Hs.symbols.gmt", package = "scMMR")
res <- RunTraceGSEA(
  seu,
  lineages = "Lineage1",
  gene.sets = kegg_gmt,
  features = Seurat::VariableFeatures(seu),
  min.pct = 0.05,
  ranking.metric = "signed_logp"
)
} # }
```
