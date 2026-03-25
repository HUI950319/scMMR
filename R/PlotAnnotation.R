# ============================================================================
# PlotAnnotation — Annotation Heatmap for Single-Cell Metadata
# ============================================================================

#' @keywords internal
.auto_anno_colors <- function(meta, columns, palette, palcolor) {
  col_list <- list()
  n_discrete <- sum(!vapply(meta[columns], is.numeric, logical(1)))
  disc_idx <- 0L

  for (col in columns) {
    vals <- meta[[col]]
    if (is.numeric(vals)) {
      rng <- range(vals, na.rm = TRUE)
      col_list[[col]] <- circlize::colorRamp2(rng, c("white", "red"))
    } else {
      disc_idx <- disc_idx + 1L
      lvls <- if (is.factor(vals)) levels(vals) else sort(unique(as.character(vals)))
      if (n_discrete > 1) {
        # rotate palette so each variable gets visually distinct colors
        pal_choices <- c("Paired", "Set1", "Set2", "Set3", "Dark2",
                         "npg", "nejm", "lancet", "jama", "aaas", "d3")
        use_pal <- if (disc_idx == 1) palette else {
          pal_choices[((disc_idx - 1) %% length(pal_choices)) + 1]
        }
      } else {
        use_pal <- palette
      }
      col_list[[col]] <- palette_colors(lvls, palette = use_pal,
                                        palcolor = if (disc_idx == 1) palcolor else NULL)
    }
  }
  col_list
}

#' Annotation Heatmap for Single-Cell Metadata
#'
#' Create a ComplexHeatmap-style annotation heatmap showing cell metadata as
#' colored annotation bars. Optionally includes a gene expression heatmap
#' matrix below the annotations.
#'
#' @param seu A Seurat object.
#' @param columns Character vector of metadata column names to display as
#'   annotation bars (e.g., \code{c("celltype", "sample", "Phase")}).
#' @param sort.by Column name to sort cells by. Default: first element
#'   of \code{columns}.
#' @param features Optional character vector of gene names. If provided,
#'   a scaled expression heatmap is drawn below the annotation bars.
#'   Default: NULL (annotation-only mode).
#' @param anno_point Optional metadata column name to display as a point
#'   annotation at the top (e.g., \code{"nCount_RNA"}).
#' @param downsample Integer; maximum number of cells per group (defined
#'   by \code{sort.by}) to keep. Default: NULL (no downsampling).
#' @param palette Palette name for discrete variables. Default: "Paired".
#' @param palcolor Optional custom color vector (overrides palette for the
#'   first discrete variable).
#' @param use_raster Logical; rasterize the expression heatmap for speed.
#'   Default: TRUE.
#' @param show_column_names Logical; show cell barcode labels. Default: FALSE.
#' @param assay Assay to use for expression data. Default: NULL
#'   (DefaultAssay).
#' @param layer Layer to pull expression from. Default: "data".
#' @param scale_rows Logical; z-score scale each gene across cells.
#'   Default: TRUE.
#' @param ... Additional arguments passed to
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object.
#'
#' @examples
#' \dontrun{
#' # Annotation-only (like clinical correlation heatmap)
#' PlotAnnotation(seu, columns = c("celltype", "sample", "Phase"))
#'
#' # With gene expression heatmap
#' PlotAnnotation(seu, columns = c("celltype", "sample"),
#'                features = c("CD3D", "CD14", "MS4A1"),
#'                sort.by = "celltype")
#'
#' # With point annotation and downsampling
#' PlotAnnotation(seu, columns = c("celltype", "condition"),
#'                anno_point = "nCount_RNA", downsample = 200)
#' }
#'
#' @export
PlotAnnotation <- function(seu,
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
                           ...) {

  # --- dependency checks ---
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Package 'ComplexHeatmap' is required. Install with:\n",
         "  BiocManager::install('ComplexHeatmap')")
  if (!inherits(seu, "Seurat"))
    stop("'seu' must be a Seurat object.", call. = FALSE)

  meta <- seu@meta.data

  # validate columns
  missing_cols <- setdiff(columns, colnames(meta))
  if (length(missing_cols) > 0)
    stop("Columns not found in metadata: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  if (!is.null(anno_point) && !anno_point %in% colnames(meta))
    stop("anno_point '", anno_point, "' not found in metadata.", call. = FALSE)

  # --- sort.by default ---
  if (is.null(sort.by)) sort.by <- columns[1]
  if (!sort.by %in% colnames(meta))
    stop("sort.by '", sort.by, "' not found in metadata.", call. = FALSE)

  # --- downsample ---
  if (!is.null(downsample)) {
    grp <- as.character(meta[[sort.by]])
    keep <- unlist(lapply(split(rownames(meta), grp), function(ids) {
      if (length(ids) <= downsample) ids
      else sample(ids, downsample)
    }))
    meta <- meta[keep, , drop = FALSE]
  }

  # --- sort ---
  sort_vals <- meta[[sort.by]]
  if (is.factor(sort_vals)) {
    ord <- order(as.integer(sort_vals))
  } else {
    ord <- order(sort_vals)
  }
  meta <- meta[ord, , drop = FALSE]
  cells <- rownames(meta)

  # --- build annotation colors ---
  col_list <- .auto_anno_colors(meta, columns, palette, palcolor)

  # --- build annotation args ---
  anno_args <- list()
  # point annotation on top
  if (!is.null(anno_point)) {
    anno_args[[anno_point]] <- ComplexHeatmap::anno_points(
      meta[[anno_point]],
      size = grid::unit(0.5, "mm"),
      gp = grid::gpar(col = "grey40")
    )
  }
  # data frame of annotation bars
  anno_args[["df"]] <- meta[, columns, drop = FALSE]
  anno_args[["col"]] <- col_list

  ha <- do.call(ComplexHeatmap::HeatmapAnnotation, anno_args)

  # --- build heatmap ---
  if (is.null(features)) {
    # annotation-only: empty matrix
    ht <- ComplexHeatmap::Heatmap(
      matrix(nrow = 0, ncol = length(cells)),
      top_annotation = ha,
      show_column_names = show_column_names,
      ...
    )
  } else {
    # expression heatmap
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seu)
    expr_mat <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
    avail <- intersect(features, rownames(expr_mat))
    if (length(avail) == 0)
      stop("None of the requested features found in assay '", assay, "'.",
           call. = FALSE)
    mat <- as.matrix(expr_mat[avail, cells, drop = FALSE])
    if (scale_rows) {
      mat <- t(scale(t(mat)))
      mat[is.nan(mat)] <- 0
    }
    # clamp for better visualization
    cap <- quantile(abs(mat), 0.99, na.rm = TRUE)
    mat[mat > cap] <- cap
    mat[mat < -cap] <- -cap

    hm_col <- circlize::colorRamp2(
      c(-cap, 0, cap), c("blue", "white", "red")
    )
    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = "Expression",
      top_annotation = ha,
      col = hm_col,
      cluster_columns = FALSE,
      cluster_rows = length(avail) > 1,
      show_column_names = show_column_names,
      show_row_names = TRUE,
      use_raster = use_raster,
      ...
    )
  }

  ht
}
