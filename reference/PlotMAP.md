# UMAP Projection Plot

Overlay query cells (with predicted UMAP coordinates) onto the reference
UMAP embedding. The reference cells are shown as a coloured scatter
plot; query cells are plotted as black points with red density contours.

## Usage

``` r
PlotMAP(
  ref,
  query,
  ref_emb = c("umap_1", "umap_2"),
  query_emb = c("umap_1_pred", "umap_2_pred"),
  color_by = "cell_type",
  facet_by = NULL,
  colors = NULL,
  point.size = 0.2,
  point.alpha = 0.1,
  query.size = 0.1,
  query.alpha = 1,
  density.color = "red",
  expand_mult = 0.05
)
```

## Arguments

- ref:

  A data.frame (or Seurat object) providing reference UMAP coordinates
  and cell-type labels.

  - If a **data.frame**: must contain columns specified by `ref_emb` and
    `color_by`.

  - If a **Seurat** object: columns are extracted via
    `Seurat::FetchData(ref, vars = c(ref_emb, color_by))`.

- query:

  A data.frame with query cell metadata, e.g. `q1@meta.data` after
  `AddMetaData(obj, pred$predictions)`. Must contain columns specified
  by `query_emb` (and optionally `facet_by`).

- ref_emb:

  Character vector of length 2: column names for reference UMAP axes.
  Default `c("umap_1", "umap_2")`.

- query_emb:

  Character vector of length 2: column names for predicted query UMAP
  axes. Default `c("umap_1_pred", "umap_2_pred")`.

- color_by:

  Column name in `ref` used for colouring reference cells. Default
  `"cell_type"`.

- facet_by:

  Optional column name in `query` for faceting (e.g. `"group"`). `NULL`
  = no faceting.

- colors:

  Optional named colour vector. If `NULL`, default ggplot2 colours are
  used.

- point.size:

  Size of reference points. Default 0.2.

- point.alpha:

  Alpha of reference points. Default 0.1.

- query.size:

  Size of query points. Default 0.1.

- query.alpha:

  Alpha of query points. Default 1.

- density.color:

  Colour of density contour lines. Default `"red"`.

- expand_mult:

  Fractional expansion applied to both axes so that density contour
  lines are not clipped at panel boundaries. Default 0.05 (5% on each
  side).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
ref_data <- Seurat::FetchData(seu,
  vars = c("umap_1", "umap_2", "cell_type"))
PlotMAP(ref = ref_data, query = q1@meta.data,
  facet_by = "group")
} # }
```
