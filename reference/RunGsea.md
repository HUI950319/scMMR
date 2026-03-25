# Run Gene Set Enrichment Analysis on DE Results

Perform GSEA for each group in a
[`RunDE`](https://hui950319.github.io/scMMR/reference/RunDE.md) result
data.frame. Each group is analyzed independently using
`clusterProfiler::GSEA()`, and the individual `gseaResult` objects are
preserved for downstream visualization (e.g., classic GSEA running score
plots via `enrichplot::gseaplot2()`).

## Usage

``` r
RunGsea(
  de_result,
  geneset = list(species = "Homo sapiens", collection = "H"),
  group.by = "group1",
  score.by = "avg_log2FC",
  de_threshold = NULL,
  scoreType = c("std", "pos", "neg"),
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  p.adjust.method = "BH",
  clean.names = TRUE,
  prefix = NULL,
  cores = 1,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- de_result:

  A data.frame from
  [`RunDE`](https://hui950319.github.io/scMMR/reference/RunDE.md),
  containing at least `gene` and `score.by` columns.

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

- group.by:

  Column in `de_result` to split groups. Default: `"group1"`.

- score.by:

  Column used as gene ranking score. Default: `"avg_log2FC"`.

- de_threshold:

  Optional character string expression to filter `de_result` before GSEA
  (e.g. `"p_val_adj < 0.05"`). Default: `NULL` (use all genes, standard
  GSEA practice).

- scoreType:

  GSEA score type: `"std"`, `"pos"`, or `"neg"`. Default: `"std"`.

- minGSSize:

  Minimum gene set size. Default: 10.

- maxGSSize:

  Maximum gene set size. Default: 500.

- pvalueCutoff:

  P-value cutoff for GSEA results. Default: 1 (no filtering).

- p.adjust.method:

  P-value adjustment method. Default: `"BH"`.

- clean.names:

  Logical; clean pathway names by removing common prefixes (e.g.
  HALLMARK\_). Default: TRUE.

- prefix:

  Optional prefix to remove from pathway names.

- cores:

  Integer; number of parallel cores. Default: 1.

- seed:

  Random seed. Default: 11.

- verbose:

  Logical; print progress. Default: TRUE.

## Value

A list with components:

- `enrichment`:

  Data.frame of all GSEA results across groups.

- `results`:

  Named list of individual `gseaResult` objects (one per group), usable
  with `enrichplot::gseaplot2()`.

- `input`:

  Data.frame of the (filtered) input gene lists.

## See also

[`RunDE`](https://hui950319.github.io/scMMR/reference/RunDE.md),
[`RunGseaEnrich`](https://hui950319.github.io/scMMR/reference/RunGseaEnrich.md)

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- RunDE(seu, group.by = "celltype")
res <- RunGsea(markers,
               geneset = list(species = "Homo sapiens", collection = "H"))

# Summary table
head(res$enrichment)

# Classic GSEA running score plot
enrichplot::gseaplot2(res$results[["T_cell"]], geneSetID = 1)

# Custom GMT file
res2 <- RunGsea(markers, geneset = "path/to/custom.gmt")
} # }
```
