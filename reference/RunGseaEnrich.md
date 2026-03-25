# Run compareCluster-based Enricher and GSEA Analysis

Perform enricher (ORA) and GSEA enrichment analysis across cell types
and groups using `clusterProfiler::compareCluster()`. Requires
`FindAllMarkers()` results grouped by cell type and condition.

## Usage

``` r
RunGseaEnrich(
  scobj,
  celltype = "celltype",
  group = "group",
  allMarkers.args = list(only.pos = FALSE, logfc.threshold = 0.1, min.pct = 0.25,
    thresh.use = 0.25),
  geneset = list(species = "Homo sapiens", collection = "H", subcollection = NULL),
  padj_cutoff = 0.001,
  plot = TRUE
)
```

## Arguments

- scobj:

  A Seurat object.

- celltype:

  Column name in `meta.data` for cell type labels. Default:
  `"celltype"`.

- group:

  Column name in `meta.data` for condition/group labels. Default:
  `"group"`.

- allMarkers.args:

  List of arguments passed to `Seurat::FindAllMarkers()`.

- geneset:

  Gene set specification. Supports multiple formats:

  msigdbr args

  :   A list with msigdbr parameters, e.g.
      `list(species = "Homo sapiens", collection = "H")`.

  data.frame

  :   With columns `gs_name` and `gene_symbol`.

  named list

  :   E.g. `list(PathA = c("TP53", "BRCA1"))`.

  GMT file

  :   Character path to a `.gmt` file.

- padj_cutoff:

  P-value adjusted cutoff for enricher ORA analysis. Default: 0.001.

- plot:

  Logical; whether to generate dotplots. Default: `TRUE`.

## Value

A list with components:

- `markers`:

  Data.frame from `FindAllMarkers()`, with added `celltype` and `group`
  columns.

- `enricher`:

  `compareClusterResult` from enricher.

- `gsea`:

  `compareClusterResult` from GSEA.

- `plots`:

  (if `plot = TRUE`) List of dotplot lists.

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)

# Using msigdbr Hallmark gene sets (default)
res <- RunGseaEnrich(seu)

# Using custom column names
res <- RunGseaEnrich(seu, celltype = "cell_type", group = "condition")

# Using a data.frame as geneset
my_gs <- data.frame(gs_name = c("PathA","PathA","PathB"),
                    gene_symbol = c("TP53","BRCA1","EGFR"))
res <- RunGseaEnrich(seu, geneset = my_gs)

# Using a named list
res <- RunGseaEnrich(seu,
  geneset = list(PathA = c("TP53","BRCA1"), PathB = c("EGFR","KRAS")))

# Using a GMT file
gmt <- system.file("extdata", "gmt",
  "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")
res <- RunGseaEnrich(seu, geneset = gmt)
} # }
```
