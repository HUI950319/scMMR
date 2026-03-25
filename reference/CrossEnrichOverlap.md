# Cross-comparison of Core Enrichment Genes between Two GSEA Results

Compare core enrichment genes between two `RunGseaEnrich` results (e.g.,
TF regulons vs pathways) to identify shared regulatory genes. Returns an
overlap summary table and optionally a heatmap or tile plot.

## Usage

``` r
CrossEnrichOverlap(
  res_tf,
  res_pathway,
  type = c("gsea", "enricher"),
  p.cutoff = 0.05,
  min.overlap = 1,
  show.genes = TRUE,
  plot = TRUE,
  top.tf = 20,
  top.pathway = 20
)
```

## Arguments

- res_tf:

  TF/GRN `RunGseaEnrich` result (or `compareClusterResult`).

- res_pathway:

  Pathway `RunGseaEnrich` result (or `compareClusterResult`).

- type:

  Which sub-result to use: `"gsea"` or `"enricher"`. Default: `"gsea"`.

- p.cutoff:

  P.adjust cutoff for filtering enriched terms. Default: 0.05.

- min.overlap:

  Minimum number of overlapping genes to report. Default: 1.

- show.genes:

  Logical; include overlapping gene names in result. Default: `TRUE`.

- plot:

  Logical; generate a tile heatmap of overlap counts. Default: `TRUE`.

- top.tf:

  Maximum number of TF terms to display (by overlap count). Default: 20.

- top.pathway:

  Maximum number of pathway terms to display. Default: 20.

## Value

A list with:

- `overlap`:

  Data.frame with columns: `TF`, `Pathway`, `n_overlap`,
  `overlap_genes`, `n_tf_genes`, `n_pathway_genes`, `jaccard`.

- `tf_genes`:

  Named list of core_enrichment genes per TF term.

- `pathway_genes`:

  Named list of core_enrichment genes per pathway term.

- `plot`:

  (if `plot = TRUE`) A ggplot tile plot of overlap counts.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare TF regulon enrichment with pathway enrichment
ov <- CrossEnrichOverlap(res_PH_TF, res_PH)
head(ov$overlap)
ov$plot

# Only enricher results, stricter cutoff
ov2 <- CrossEnrichOverlap(res_TF, res_pathway,
  type = "enricher", p.cutoff = 0.01, min.overlap = 3)
} # }
```
