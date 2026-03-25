# Run Differential Expression Test

Perform differential expression (DE) testing across cell groups in a
Seurat object. Supports four marker types: one-vs-rest (`"all"`),
pairwise (`"paired"`), conserved across conditions (`"conserved"`), and
condition-specific (`"disturbed"`). Optionally downsamples cells before
testing for balanced comparisons.

## Usage

``` r
RunDE(
  seu,
  group.by = NULL,
  split.by = NULL,
  group1 = NULL,
  group2 = NULL,
  cells1 = NULL,
  cells2 = NULL,
  features = NULL,
  markers_type = c("all", "paired", "conserved", "disturbed"),
  grouping.var = NULL,
  meta.method = c("maximump", "minimump", "wilkinsonp", "meanp", "sump", "votep"),
  test.use = "wilcox",
  only.pos = TRUE,
  fc.threshold = 1.5,
  p_val_cutoff = NULL,
  logfc_cutoff = NULL,
  base = 2,
  pseudocount.use = 1,
  mean.fxn = NULL,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  max.cells.per.ident = Inf,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  norm.method = "LogNormalize",
  p.adjust.method = "bonferroni",
  downsample = NULL,
  layer = "data",
  assay = NULL,
  cores = max(1L, parallel::detectCores()%/%2),
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- seu:

  A Seurat object.

- group.by:

  Column name in `meta.data` for grouping. If `NULL`, uses
  `Idents(seu)`.

- split.by:

  Optional column name for a second grouping variable (e.g. condition).
  When provided, cells are grouped by `group.by x split.by` composite
  identity for DE testing. Result columns `celltype` and `condition` are
  added. Default: NULL.

- group1:

  Character vector of group labels for the first group.

- group2:

  Character vector of group labels for the second group.

- cells1:

  Character vector of cell barcodes for the first group. Overrides
  `group1`.

- cells2:

  Character vector of cell barcodes for the second group.

- features:

  Character vector of genes to test. Default: all genes.

- markers_type:

  Type of marker analysis. One of `"all"` (one-vs-rest), `"paired"` (all
  pairwise), `"conserved"` (across conditions), or `"disturbed"`
  (condition-specific).

- grouping.var:

  Column for condition variable (required for conserved/disturbed).

- meta.method:

  Method for combining p-values in conserved markers. One of
  `"maximump"`, `"minimump"`, `"wilkinsonp"`, `"meanp"`, `"sump"`,
  `"votep"`.

- test.use:

  DE test method. Passed to `Seurat::FindMarkers`.

- only.pos:

  Logical; return only positive markers. Default: TRUE.

- fc.threshold:

  Fold change threshold (linear scale, \>= 1). Pre-filter before
  testing: genes below this threshold are skipped. Default: 1.5.

- p_val_cutoff:

  Numeric or NULL. Post-filter: keep only results with
  `p_val_adj < p_val_cutoff`. Default: NULL (no filtering).

- logfc_cutoff:

  Numeric or NULL. Post-filter: keep only results with
  `|avg_log2FC| > logfc_cutoff`. Default: NULL (no filtering).

- base:

  Log base for fold change. Default: 2.

- pseudocount.use:

  Pseudocount for log fold change. Default: 1.

- mean.fxn:

  Custom mean function for fold change.

- min.pct:

  Minimum fraction of cells expressing the gene. Default: 0.1.

- min.diff.pct:

  Minimum difference in fraction between groups.

- max.cells.per.ident:

  Maximum cells per identity for testing.

- latent.vars:

  Latent variables for certain tests (MAST, LR, etc.).

- min.cells.feature:

  Minimum cells expressing feature. Default: 3.

- min.cells.group:

  Minimum cells per group. Default: 3.

- norm.method:

  Normalization method. Default: "LogNormalize".

- p.adjust.method:

  P-value adjustment method. Default: "bonferroni".

- downsample:

  Integer or NULL. If specified, downsample each identity class to at
  most this many cells before DE testing. Uses
  `subset(seu, downsample = downsample)`. Default: NULL.

- layer:

  Data layer to use. Default: "data".

- assay:

  Assay to use. Default: `DefaultAssay(seu)`.

- cores:

  Integer. Number of parallel cores. Default:
  `parallel::detectCores() %/% 2` (half of available cores). Set to 1 to
  disable parallelism.

- seed:

  Random seed. Default: 11.

- verbose:

  Logical; print progress messages. Default: TRUE.

## Value

A data.frame of DE results with columns: p_val, avg_log2FC, pct.1,
pct.2, p_val_adj, gene, group1, group2, test_group_number, test_group.

## See also

`FindMarkers`,
[`RunPathwayAnalysis`](https://hui950319.github.io/scMMR/reference/RunPathwayAnalysis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# One-vs-rest DE
markers <- RunDE(seu, group.by = "celltype")

# Pairwise with downsampling
markers <- RunDE(seu, group.by = "celltype",
                 markers_type = "paired", downsample = 200)

# Conserved markers across conditions
markers <- RunDE(seu, group.by = "celltype",
                 grouping.var = "group",
                 markers_type = "conserved")

# Joint celltype x condition DE (for downstream GSEA)
markers <- RunDE(seu, group.by = "celltype", split.by = "condition")

# With post-filtering: padj < 0.05 and |log2FC| > 1
markers <- RunDE(seu, group.by = "celltype",
                 p_val_cutoff = 0.05, logfc_cutoff = 1)
} # }
```
