# Run Pathway Analysis across Cell Types

Perform GSEA or SCPA pathway analysis for each cell type in a Seurat
object. Compares two conditions (`ident.1` vs `ident.2`) within each
cell type, and returns a standardized data.frame suitable for
visualization with
[`PlotPathwayBubble`](https://hui950319.github.io/scMMR/reference/PlotPathwayBubble.md).

## Usage

``` r
RunPathwayAnalysis(
  seu,
  gene.sets,
  method = c("GSEA", "SCPA"),
  celltype.col = "celltype",
  group.col = "group",
  ident.1 = NULL,
  ident.2 = NULL,
  celltypes = NULL,
  pvalue.cutoff = 0.05,
  logfc.threshold = 0,
  eps = 0,
  minGSSize = 10,
  maxGSSize = 500,
  clean.names = TRUE,
  prefix = NULL,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- seu:

  A Seurat object with cell type and condition labels in `meta.data`.

- gene.sets:

  Gene sets in one of three formats:

  - A named list of character vectors (gene names).

  - A data.frame with at least two columns: term and gene.

  - A file path to a GMT file.

  Internally parsed by
  [`parse_gene_sets`](https://hui950319.github.io/scMMR/reference/parse_gene_sets.md).

- method:

  Pathway analysis method: `"GSEA"` or `"SCPA"`.

  `"GSEA"`

  :   Uses `FindMarkers()` to obtain log2 fold changes, ranks genes, and
      runs `clusterProfiler::GSEA()`. Requires clusterProfiler.

  `"SCPA"`

  :   Uses `SCPA::seurat_extract()` and `SCPA::compare_pathways()` for
      distribution-based pathway comparison. Requires SCPA.

- celltype.col:

  Column name in `meta.data` containing cell type labels. Default:
  `"celltype"`.

- group.col:

  Column name in `meta.data` containing condition/group labels. Default:
  `"group"`.

- ident.1:

  Treatment/test condition label (a value in `group.col`). Fold change
  direction: `ident.1 / ident.2`.

- ident.2:

  Control/reference condition label (a value in `group.col`).

- celltypes:

  Character vector specifying which cell types to analyze. Default:
  `NULL` (all cell types in `celltype.col`).

- pvalue.cutoff:

  P-value cutoff for GSEA enrichment results. Default: 0.05.

- logfc.threshold:

  Minimum log2 fold-change threshold for `FindMarkers()` (GSEA only).
  Default: 0 (include all genes).

- eps:

  Boundary correction for GSEA p-value calculation. Default: 0 (exact
  p-values).

- minGSSize:

  Minimum number of genes in a gene set for GSEA. Default: 10.

- maxGSSize:

  Maximum number of genes in a gene set for GSEA. Default: 500.

- clean.names:

  Logical; if `TRUE`, remove common database prefixes (e.g. `HALLMARK_`,
  `KEGG_`) and convert to Title Case. Default: `TRUE`.

- prefix:

  Custom regex prefix to remove when `clean.names = TRUE`. Default:
  `NULL` (auto-detect common prefixes).

- cores:

  Number of parallel cores for SCPA. Default: 1.

- verbose:

  Logical; print progress messages. Default: `TRUE`.

## Value

A data.frame with standardized columns:

- `Celltype`:

  Cell type name.

- `Pathway`:

  Pathway name (cleaned if `clean.names = TRUE`).

- `Score`:

  Enrichment score: NES (GSEA) or FC (SCPA).

- `PValue`:

  Raw p-value.

- `AdjPValue`:

  Adjusted p-value (BH).

- `QValue`:

  Q-value.

- `Sign`:

  `"Activated"` (Score \> 0) or `"Repressed"` (Score \< 0).

- `Method`:

  `"GSEA"` or `"SCPA"`.

GSEA results additionally include: `SetSize`, `EnrichmentScore`,
`CoreEnrichment`.

## See also

[`PlotPathwayBubble`](https://hui950319.github.io/scMMR/reference/PlotPathwayBubble.md),
[`parse_gene_sets`](https://hui950319.github.io/scMMR/reference/parse_gene_sets.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)
library(Seurat)

# Gene sets from GMT file
gmt_file <- system.file("extdata", "gmt",
  "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")

# --- GSEA ---
gsea_res <- RunPathwayAnalysis(
  seu, gene.sets = gmt_file, method = "GSEA",
  celltype.col = "celltype", group.col = "group",
  ident.1 = "PH", ident.2 = "PT"
)
PlotPathwayBubble(gsea_res)

# --- SCPA (specific cell types) ---
scpa_res <- RunPathwayAnalysis(
  seu, gene.sets = gmt_file, method = "SCPA",
  celltype.col = "celltype", group.col = "group",
  ident.1 = "PH", ident.2 = "PT",
  celltypes = c("Parathyroid cells", "Fibroblasts")
)
PlotPathwayBubble(scpa_res, size.by = "qvalue")
} # }
```
