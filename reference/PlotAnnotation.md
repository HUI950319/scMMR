# Annotation Heatmap for Single-Cell Metadata

Create a ComplexHeatmap-style annotation heatmap showing cell metadata
as colored annotation bars. Optionally includes a gene expression
heatmap matrix below the annotations.

## Usage

``` r
PlotAnnotation(
  seu,
  columns,
  sort.by = NULL,
  features = NULL,
  anno_point = NULL,
  downsample = NULL,
  palette = "Paired",
  palcolor = NULL,
  use_raster = TRUE,
  show_column_names = FALSE,
  assay = NULL,
  layer = "data",
  scale_rows = TRUE,
  ...
)
```

## Arguments

- seu:

  A Seurat object.

- columns:

  Character vector of metadata column names to display as annotation
  bars (e.g., `c("celltype", "sample", "Phase")`).

- sort.by:

  Column name to sort cells by. Default: first element of `columns`.

- features:

  Optional character vector of gene names. If provided, a scaled
  expression heatmap is drawn below the annotation bars. Default: NULL
  (annotation-only mode).

- anno_point:

  Optional metadata column name to display as a point annotation at the
  top (e.g., `"nCount_RNA"`).

- downsample:

  Integer; maximum number of cells per group (defined by `sort.by`) to
  keep. Default: NULL (no downsampling).

- palette:

  Palette name for discrete variables. Default: "Paired".

- palcolor:

  Optional custom color vector (overrides palette for the first discrete
  variable).

- use_raster:

  Logical; rasterize the expression heatmap for speed. Default: TRUE.

- show_column_names:

  Logical; show cell barcode labels. Default: FALSE.

- assay:

  Assay to use for expression data. Default: NULL (DefaultAssay).

- layer:

  Layer to pull expression from. Default: "data".

- scale_rows:

  Logical; z-score scale each gene across cells. Default: TRUE.

- ...:

  Additional arguments passed to `Heatmap`.

## Value

A `ComplexHeatmap::Heatmap` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Annotation-only (like clinical correlation heatmap)
PlotAnnotation(seu, columns = c("celltype", "sample", "Phase"))

# With gene expression heatmap
PlotAnnotation(seu, columns = c("celltype", "sample"),
               features = c("CD3D", "CD14", "MS4A1"),
               sort.by = "celltype")

# With point annotation and downsampling
PlotAnnotation(seu, columns = c("celltype", "condition"),
               anno_point = "nCount_RNA", downsample = 200)
} # }
```
