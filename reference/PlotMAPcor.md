# Cell Type Correlation Heatmap

Visualise the similarity between cell types in query and reference
Seurat objects using a ComplexHeatmap-based heatmap. Cell types are
defined by collapsing (averaging) cells within each group, then
computing pairwise similarity (cosine, Pearson, or Spearman).

When `ref = NULL`, the function computes self-correlation within the
query object.

## Usage

``` r
PlotMAPcor(
  query,
  ref = NULL,
  query_group = NULL,
  ref_group = NULL,
  query_assay = NULL,
  ref_assay = NULL,
  features = NULL,
  nfeatures = 2000,
  method = c("cosine", "pearson", "spearman"),
  query_annotation = NULL,
  ref_annotation = NULL,
  heatmap_palette = "RdBu",
  heatmap_palcolor = NULL,
  query_palette = "Paired",
  query_palcolor = NULL,
  ref_palette = "Set1",
  ref_palcolor = NULL,
  query_annotation_palette = "Set2",
  query_annotation_palcolor = NULL,
  ref_annotation_palette = "Set3",
  ref_annotation_palcolor = NULL,
  border = TRUE,
  limits = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "top",
  nlabel = 0,
  label_cutoff = 0,
  label_by = "row",
  label_size = 10,
  title = NULL,
  width = NULL,
  height = NULL,
  seed = 11
)
```

## Arguments

- query:

  A Seurat object (query dataset).

- ref:

  A Seurat object (reference dataset). `NULL` = self-vs-self.

- query_group:

  Column name in query metadata for cell type grouping.

- ref_group:

  Column name in ref metadata for cell type grouping.

- query_assay:

  Assay name for query. `NULL` = DefaultAssay.

- ref_assay:

  Assay name for ref. `NULL` = DefaultAssay.

- features:

  Character vector of features. `NULL` = auto HVF.

- nfeatures:

  Max number of HVF features. Default 2000.

- method:

  Similarity method: `"cosine"`, `"pearson"`, or `"spearman"`. Default
  `"cosine"`.

- query_annotation:

  Character vector of categorical metadata columns in query to show as
  annotation tracks. Default `NULL`.

- ref_annotation:

  Character vector of categorical metadata columns in ref to show as
  annotation tracks. Default `NULL`.

- heatmap_palette:

  Palette for heatmap body. Default `"RdBu"`.

- heatmap_palcolor:

  Custom colours for heatmap body. Default `NULL`.

- query_palette:

  Palette for query group block. Default `"Paired"`.

- query_palcolor:

  Custom colours for query group. Default `NULL`.

- ref_palette:

  Palette for ref group block. Default `"Set1"`.

- ref_palcolor:

  Custom colours for ref group. Default `NULL`.

- query_annotation_palette:

  Palette(s) for query annotation tracks. Default `"Set2"`.

- query_annotation_palcolor:

  Custom colour(s) for query annotations.

- ref_annotation_palette:

  Palette(s) for ref annotation tracks. Default `"Set3"`.

- ref_annotation_palcolor:

  Custom colour(s) for ref annotations.

- border:

  Logical. Draw cell borders. Default `TRUE`.

- limits:

  Numeric vector of length 2 for colour scale limits. `NULL` = auto.

- cluster_rows:

  Logical. Cluster rows. Default `FALSE`.

- cluster_columns:

  Logical. Cluster columns. Default `FALSE`.

- show_row_names:

  Logical. Default `TRUE`.

- show_column_names:

  Logical. Default `TRUE`.

- row_names_side:

  Side for row names. Default `"left"`.

- column_names_side:

  Side for column names. Default `"top"`.

- nlabel:

  Integer. Top-N similarity values to label per row/column. Default `0`.

- label_cutoff:

  Numeric. Min similarity to show label. Default `0`.

- label_by:

  Label selection dimension: `"row"`, `"column"`, or `"both"`. Default
  `"row"`.

- label_size:

  Numeric. Label font size. Default `10`.

- title:

  Character. Heatmap title. Default `NULL` = auto.

- width:

  Numeric. Heatmap body width in inches. `NULL` = auto.

- height:

  Numeric. Heatmap body height in inches. `NULL` = auto.

- seed:

  Integer. Random seed. Default `11`.

## Value

A list with elements:

- plot:

  A patchwork/ggplot object.

- simil_matrix:

  Similarity matrix (query groups x ref groups).

- features:

  Features used for computation.

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)

# ---- simulate test data ----
set.seed(42)
ng <- 200; nc <- 100
counts <- matrix(0, nrow = ng, ncol = nc)
rownames(counts) <- paste0("Gene", 1:ng)
colnames(counts) <- paste0("Cell", 1:nc)
ct <- rep(c("A", "B", "C", "D", "E"), each = 20)
for (i in 1:nc) {
  base <- rpois(ng, lambda = 3)
  gi <- match(ct[i], c("A", "B", "C", "D", "E"))
  ms <- (gi - 1) * 20 + 1
  base[ms:(ms + 19)] <- base[ms:(ms + 19)] + rpois(20, 10)
  counts[, i] <- base
}
srt <- CreateSeuratObject(counts = counts)
srt <- NormalizeData(srt, verbose = FALSE)
srt <- AddMetaData(srt,
  factor(ct, levels = c("A", "B", "C", "D", "E")), "celltype")
srt <- AddMetaData(srt,
  factor(sample(c("G1", "S", "G2M"), nc, replace = TRUE)), "Phase")
srt <- AddMetaData(srt,
  factor(sample(c("B1", "B2"), nc, replace = TRUE)), "batch")

# ---- 1. Self-correlation (cosine, default) ----
res1 <- PlotMAPcor(query = srt, query_group = "celltype")
res1$plot

# ---- 2. Switch similarity method ----
res2 <- PlotMAPcor(query = srt, query_group = "celltype",
                    method = "pearson")

# ---- 3. Add categorical annotation tracks ----
res3 <- PlotMAPcor(query = srt, query_group = "celltype",
                    query_annotation = c("Phase", "batch"))
res3$plot

# ---- 4. Show top-N similarity labels ----
res4 <- PlotMAPcor(query = srt, query_group = "celltype",
                    nlabel = 2, label_cutoff = 0.9)

# ---- 5. Hierarchical clustering ----
res5 <- PlotMAPcor(query = srt, query_group = "celltype",
                    cluster_rows = TRUE, cluster_columns = TRUE)

# ---- 6. Custom colours and fixed limits ----
res6 <- PlotMAPcor(query = srt, query_group = "celltype",
                    heatmap_palette = "Spectral",
                    query_palette = "npg",
                    limits = c(0, 1))

# ---- 7. Custom feature set ----
res7 <- PlotMAPcor(query = srt, query_group = "celltype",
                    features = paste0("Gene", 1:50))

# ---- 8. Query vs Reference ----
counts2 <- matrix(0, nrow = ng, ncol = 60)
rownames(counts2) <- paste0("Gene", 1:ng)
colnames(counts2) <- paste0("R", 1:60)
ct2 <- rep(c("CT1", "CT2", "CT3"), each = 20)
for (i in 1:60) {
  base <- rpois(ng, lambda = 3)
  gi <- match(ct2[i], c("CT1", "CT2", "CT3"))
  ms <- (gi - 1) * 30 + 1
  base[ms:(ms + 29)] <- base[ms:(ms + 29)] + rpois(30, 8)
  counts2[, i] <- base
}
srt2 <- CreateSeuratObject(counts = counts2)
srt2 <- NormalizeData(srt2, verbose = FALSE)
srt2 <- AddMetaData(srt2, factor(ct2), "celltype")
srt2 <- AddMetaData(srt2,
  factor(sample(c("10X", "SS2"), 60, replace = TRUE)), "tech")

res8 <- PlotMAPcor(
  query = srt, ref = srt2,
  query_group = "celltype", ref_group = "celltype",
  query_annotation = "Phase",
  ref_annotation = "tech",
  method = "spearman",
  nlabel = 1
)
res8$plot

# ---- 9. Save output ----
ggplot2::ggsave("heatmap.pdf", res8$plot, width = 9, height = 7)

# ---- 10. Access similarity matrix ----
res1$simil_matrix   # query_groups x ref_groups
res1$features       # features used
} # }
```
