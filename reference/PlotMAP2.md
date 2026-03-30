# Single-Cell Atlas Projection Plot

Overlay query cells onto the reference 2-D embedding (UMAP / tSNE / any
2-D space). Supports **three coordinate input modes** for both query and
reference — making the function independent of any specific mapping
method:

- **Mode 1 – Seurat reduction**:

  Pass a Seurat object and specify `query_reduction` / `ref_reduction`.
  Coordinates are read from `@reductions[[name]]@cell.embeddings`. If
  omitted, the function auto-detects (prefers UMAP/tSNE).

- **Mode 2 – metadata columns**:

  Pass a Seurat object or data.frame and specify `query_emb` / `ref_emb`
  as a character vector of two column names. Coordinates are read
  directly from `@meta.data` (Seurat) or the data.frame itself. Useful
  when predicted coordinates are stored in metadata (e.g. from
  `DNN_predict`).

- **Mode 3 – data.frame column indices**:

  Pass a plain data.frame; columns selected by `query_dims` / `ref_dims`
  (default `c(1, 2)`) are used as x/y.

Priority when multiple modes apply: **`*_emb` \> `*_reduction` \>
auto-detect reduction \> `*_dims` fallback**.

## Usage

``` r
PlotMAP2(
  ref,
  query,
  ref_emb = NULL,
  query_emb = NULL,
  ref_group = NULL,
  query_group = NULL,
  ref_reduction = NULL,
  ref_dims = c(1, 2),
  query_reduction = NULL,
  query_dims = c(1, 2),
  query_palette = "Set1",
  query_palcolor = NULL,
  query_color = "firebrick",
  ref_palette = "Paired",
  ref_palcolor = NULL,
  ref_color = "grey70",
  query.size = 0.8,
  query.alpha = 1,
  stroke = 0.5,
  stroke.color = "black",
  point.size = NULL,
  point.alpha = 0.4,
  show_density = FALSE,
  density.color = "red",
  density_filled = FALSE,
  density_filled_palette = "Greys",
  density_filled_palcolor = NULL,
  xlim = NULL,
  ylim = NULL,
  legend.position = "right",
  legend.direction = "vertical",
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  theme_use = NULL,
  theme_args = list(),
  aspect.ratio = 1,
  raster = NULL,
  raster_method = c("scattermore", "rasterise"),
  raster.dpi = c(512, 512),
  seed = 11
)
```

## Arguments

- ref:

  A Seurat object or `data.frame` for reference cells.

- query:

  A Seurat object or `data.frame` for query cells.

- ref_emb:

  See section above. Default `NULL`.

- query_emb:

  See section above. Default `NULL`.

- ref_group:

  Character. Column in ref metadata for colouring reference cells.
  `NULL` = single colour (`ref_color`).

- query_group:

  Character. Column in query metadata for colouring query cells. `NULL`
  = single colour (`query_color`).

- ref_reduction:

  See section above. Default `NULL`.

- ref_dims:

  See section above. Default `c(1, 2)`.

- query_reduction:

  See section above. Default `NULL`.

- query_dims:

  See section above. Default `c(1, 2)`.

- query_palette:

  Character. Palette for query group colours. Default `"Set1"`.

- query_palcolor:

  Character vector. Custom query colours. Default `NULL`.

- query_color:

  Character. Single colour when `query_group = NULL`. Default
  `"firebrick"`.

- ref_palette:

  Character. Palette for ref group colours. Default `"Paired"`.

- ref_palcolor:

  Character vector. Custom ref colours. Default `NULL`.

- ref_color:

  Character. Single colour when `ref_group = NULL`. Default `"grey70"`.

- query.size:

  Numeric. Query point size. Default `0.8`.

- query.alpha:

  Numeric. Query point alpha (0–1). Default `1`.

- stroke:

  Numeric. Border stroke width around query points. Set `0` to disable.
  Default `0.5`.

- stroke.color:

  Character. Stroke colour. Default `"black"`.

- point.size:

  Numeric. Reference point size. `NULL` = auto.

- point.alpha:

  Numeric. Reference point alpha. Default `0.4`.

- show_density:

  Logical. Add 2-D density contour for query cells. Default `FALSE`.

- density.color:

  Character. Contour line colour. Default `"red"`.

- density_filled:

  Logical. Filled contour bands. Default `FALSE`.

- density_filled_palette:

  Character. Palette for filled bands. Default `"Greys"`.

- density_filled_palcolor:

  Character vector. Custom fill colours. Default `NULL`.

- xlim, ylim:

  Numeric vectors of length 2 for axis limits. `NULL` = auto from
  combined data range.

- legend.position, legend.direction:

  Legend position/direction. Defaults `"right"` / `"vertical"`.

- title, subtitle, xlab, ylab:

  Plot labels. `NULL` = auto or empty.

- theme_use:

  Theme: function, character name, or `theme` object. Default `NULL` →
  `UtilsR::theme_blank` if available.

- theme_args:

  List of extra args for `theme_use`.

- aspect.ratio:

  Numeric. Panel aspect ratio. Default `1`.

- raster:

  Logical. Rasterise points. `NULL` = auto (\> 100k cells).

- raster_method:

  Character. Rasterisation backend: `"rasterise"` (default, faithful
  colours) or `"scattermore"` (faster, may shift colours).

- raster.dpi:

  Numeric scalar or 2-length vector. Default `c(512,512)`.

- seed:

  Integer. Random seed. Default `11`.

## Value

A `ggplot` object.

## Query coordinate specification (choose one)

- `query_emb`:

  Character vector of length 2: column names in `@meta.data` (Seurat) or
  in the data.frame to use as x/y. Takes priority over
  `query_reduction`.

- `query_reduction`:

  Character. Reduction name in the Seurat object. `NULL` = auto-detect
  (prefers UMAP/tSNE, then last reduction).

- `query_dims`:

  Integer vector of length 2. Dimension indices within the chosen
  reduction (Seurat) or column indices (data.frame). Default `c(1, 2)`.

## Reference coordinate specification (choose one)

- `ref_emb`:

  Character vector of length 2: column names to use as x/y. Takes
  priority over `ref_reduction`.

- `ref_reduction`:

  Character. Reduction name in the Seurat object. `NULL` = auto-detect.

- `ref_dims`:

  Integer vector of length 2. Default `c(1, 2)`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Mode 1: Seurat + reduction (after RunKNNMap)
q <- RunKNNMap(q, ref, ref_umap = "UMAP")
PlotMAP2(query = q, ref = ref_data,
         query_group = "celltype", ref_group = "celltype",
         query_reduction = "ref.embeddings")

# Mode 2: Seurat + meta.data column names (after DNN_predict)
PlotMAP2(query = q1, ref = ref_data,
         query_emb   = c("umap_1_pred", "umap_2_pred"),
         ref_emb     = c("umap_1", "umap_2"),
         query_group = "cell_type_pred",
         ref_group   = "celltype")

# Mode 3: plain data.frames
PlotMAP2(query = query_df, ref = ref_df,
         query_dims  = c(1, 2),
         ref_dims    = c(1, 2),
         query_group = "cluster",
         ref_group   = "celltype")
} # }
```
