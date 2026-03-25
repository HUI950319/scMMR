# Find Upstream Regulon Drivers of Pathways

Identify transcription factor regulons that may drive pathway activity
by integrating two complementary metrics:

1.  **Correlation**: Pearson/Spearman correlation between regulon
    activity scores and pathway scores at single-cell level.

2.  **Gene overlap**: intersection of regulon target genes with pathway
    gene sets, plus fraction of hits (`n_overlap / regulon_size`).

## Usage

``` r
FindRegulonDrivers(
  seu,
  pathway.genes,
  regulons,
  pathways = NULL,
  score.method = c("AUCell", "Seurat", "UCell"),
  min.size = 5,
  cor.method = c("pearson", "spearman"),
  min.overlap = 1,
  label.n = 5,
  label.overlap = NULL,
  plot = TRUE,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- seu:

  A Seurat object.

- pathway.genes:

  Pathway gene sets. Accepts:

  - Named list of character vectors.

  - GMT file path.

  - data.frame with columns (term, gene).

- regulons:

  Regulon target gene sets. Same formats as `pathway.genes` (typically
  from SCENIC `regulons.rds`).

- pathways:

  Character vector of specific pathways to analyze. Default: `NULL` (all
  pathways).

- score.method:

  Scoring method passed to `ComputeModuleScore`: `"AUCell"`, `"Seurat"`,
  or `"UCell"`. Default: `"AUCell"`.

- min.size:

  Minimum gene set size for scoring. Default: 5.

- cor.method:

  Correlation method: `"pearson"` or `"spearman"`. Default: `"pearson"`.

- min.overlap:

  Minimum number of overlapping genes to report. Default: 1.

- label.n:

  Number of top/bottom regulons to label in PCC rank bar plot. Default:
  5.

- label.overlap:

  Overlap count threshold for labeling in scatter plot. Default: `NULL`
  (auto: top 10% quantile).

- plot:

  Logical; generate plots. Default: `TRUE`.

- cores:

  Number of parallel cores for scoring. Default: 1.

- verbose:

  Logical; print progress. Default: `TRUE`.

## Value

A list with:

- `result`:

  Data.frame with columns: `regulon`, `PCC`, `rank`, `n_overlap`,
  `regulon_size`, `frac_of_hits`, `pathway`. One row per regulon ×
  pathway.

- `cor.mat`:

  Correlation matrix (pathways × regulons).

- `plots`:

  (if `plot = TRUE`) List of ggplot objects: `rank_bar` (PCC rank) and
  `scatter` (PCC vs overlap), one set per pathway.

## Details

Internally uses
[`ComputeModuleScore`](https://hui950319.github.io/scMMR/reference/ComputeModuleScore.md)
to score both regulon and pathway gene sets in the Seurat object.

## See also

[`ComputeModuleScore`](https://hui950319.github.io/scMMR/reference/ComputeModuleScore.md),
[`CrossEnrichOverlap`](https://hui950319.github.io/scMMR/reference/CrossEnrichOverlap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
regulons <- readRDS("regulons.rds")
hallmark <- system.file("extdata", "gmt",
  "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")

res <- FindRegulonDrivers(seu, pathway.genes = hallmark, regulons = regulons)
head(res$result)
res$plots[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]$scatter

# Focus on specific pathways
res2 <- FindRegulonDrivers(seu,
  pathway.genes = hallmark, regulons = regulons,
  pathways = c("HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"))
} # }
```
