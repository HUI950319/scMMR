# ============================================================================
# PlotMAP.R -- Universal single-cell atlas projection plot (renamed from PlotMAP2)
# ============================================================================

#' @title Single-Cell Atlas Projection Plot
#'
#' @description
#' Overlay query cells onto the reference 2-D embedding (UMAP / tSNE / any
#' 2-D space).  Supports **three coordinate input modes** for both query and
#' reference — making the function independent of any specific mapping method:
#'
#' \describe{
#'   \item{**Mode 1 – Seurat reduction**}{
#'     Pass a Seurat object and specify `query_reduction` / `ref_reduction`.
#'     Coordinates are read from `@@reductions[[name]]@@cell.embeddings`.
#'     If omitted, the function auto-detects (prefers UMAP/tSNE).}
#'   \item{**Mode 2 – metadata columns**}{
#'     Pass a Seurat object or data.frame and specify `query_emb` / `ref_emb`
#'     as a character vector of two column names.  Coordinates are read
#'     directly from `@@meta.data` (Seurat) or the data.frame itself.
#'     Useful when predicted coordinates are stored in metadata (e.g. from
#'     \code{DNN_predict}).}
#'   \item{**Mode 3 – data.frame column indices**}{
#'     Pass a plain data.frame; columns selected by `query_dims` / `ref_dims`
#'     (default \code{c(1, 2)}) are used as x/y.}}
#'
#' Priority when multiple modes apply: **`*_emb` > `*_reduction` > auto-detect
#' reduction > `*_dims` fallback**.
#'
#' @md
#' @param ref A Seurat object or \code{data.frame} for reference cells.
#' @param query A Seurat object or \code{data.frame} for query cells.
#' @param ref_group Character. Column in ref metadata for colouring reference
#'   cells.  \code{NULL} = single colour (\code{ref_color}).
#' @param query_group Character. Column in query metadata for colouring query
#'   cells.  \code{NULL} = single colour (\code{query_color}).
#'
#' @section Query coordinate specification (choose one):
#' \describe{
#'   \item{`query_emb`}{Character vector of length 2: column names in
#'     \code{@@meta.data} (Seurat) or in the data.frame to use as x/y.
#'     Takes priority over `query_reduction`.}
#'   \item{`query_reduction`}{Character. Reduction name in the Seurat object.
#'     \code{NULL} = auto-detect (prefers UMAP/tSNE, then last reduction).}
#'   \item{`query_dims`}{Integer vector of length 2. Dimension indices within
#'     the chosen reduction (Seurat) or column indices (data.frame).
#'     Default \code{c(1, 2)}.}
#' }
#'
#' @section Reference coordinate specification (choose one):
#' \describe{
#'   \item{`ref_emb`}{Character vector of length 2: column names to use as
#'     x/y.  Takes priority over `ref_reduction`.}
#'   \item{`ref_reduction`}{Character. Reduction name in the Seurat object.
#'     \code{NULL} = auto-detect.}
#'   \item{`ref_dims`}{Integer vector of length 2.  Default \code{c(1, 2)}.}
#' }
#'
#' @param query_emb See section above. Default \code{NULL}.
#' @param query_reduction See section above. Default \code{NULL}.
#' @param query_dims See section above. Default \code{c(1, 2)}.
#' @param ref_emb See section above. Default \code{NULL}.
#' @param ref_reduction See section above. Default \code{NULL}.
#' @param ref_dims See section above. Default \code{c(1, 2)}.
#'
#' @param query_palette Character. Palette for query group colours.
#'   Default \code{"Set1"}.
#' @param query_palcolor Character vector. Custom query colours.
#'   Default \code{NULL}.
#' @param query_color Character. Single colour when \code{query_group = NULL}.
#'   Default \code{"firebrick"}.
#' @param ref_palette Character. Palette for ref group colours.
#'   Default \code{"Paired"}.
#' @param ref_palcolor Character vector. Custom ref colours. Default \code{NULL}.
#' @param ref_color Character. Single colour when \code{ref_group = NULL}.
#'   Default \code{"grey70"}.
#' @param query.size Numeric. Query point size. Default \code{0.8}.
#' @param query.alpha Numeric. Query point alpha (0–1). Default \code{1}.
#' @param stroke Numeric. Border stroke width around query points.
#'   Set \code{0} to disable. Default \code{0.5}.
#' @param stroke.color Character. Stroke colour. Default \code{"black"}.
#' @param point.size Numeric. Reference point size. \code{NULL} = auto.
#' @param point.alpha Numeric. Reference point alpha. Default \code{0.4}.
#' @param show_density Logical. Add 2-D density contour for query cells.
#'   Default \code{FALSE}.
#' @param density.color Character. Contour line colour. Default \code{"red"}.
#' @param density_filled Logical. Filled contour bands. Default \code{FALSE}.
#' @param density_filled_palette Character. Palette for filled bands.
#'   Default \code{"Greys"}.
#' @param density_filled_palcolor Character vector. Custom fill colours.
#'   Default \code{NULL}.
#' @param xlim,ylim Numeric vectors of length 2 for axis limits.
#'   \code{NULL} = auto from combined data range.
#' @param legend.position,legend.direction Legend position/direction.
#'   Defaults \code{"right"} / \code{"vertical"}.
#' @param title,subtitle,xlab,ylab Plot labels. \code{NULL} = auto or empty.
#' @param theme_use Theme: function, character name, or \code{theme} object.
#'   Default \code{NULL} → \code{UtilsR::theme_blank} if available.
#' @param theme_args List of extra args for \code{theme_use}.
#' @param aspect.ratio Numeric. Panel aspect ratio. Default \code{1}.
#' @param raster Logical. Rasterise points. \code{NULL} = auto (> 100k cells).
#' @param raster_method Character. Rasterisation backend:
#'   \code{"rasterise"} (default, faithful colours) or
#'   \code{"scattermore"} (faster, may shift colours).
#' @param raster.dpi Numeric scalar or 2-length vector. Default \code{c(512,512)}.
#' @param seed Integer. Random seed. Default \code{11}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # Mode 1: Seurat + reduction (after RunKNNMap)
#' q <- RunKNNMap(q, ref, ref_umap = "UMAP")
#' PlotMAP(query = q, ref = ref_data,
#'          query_group = "celltype", ref_group = "celltype",
#'          query_reduction = "ref.embeddings")
#'
#' # Mode 2: Seurat + meta.data column names (after DNN_predict)
#' PlotMAP(query = q1, ref = ref_data,
#'          query_emb   = c("umap_1_pred", "umap_2_pred"),
#'          ref_emb     = c("umap_1", "umap_2"),
#'          query_group = "cell_type_pred",
#'          ref_group   = "celltype")
#'
#' # Mode 3: plain data.frames
#' PlotMAP(query = query_df, ref = ref_df,
#'          query_dims  = c(1, 2),
#'          ref_dims    = c(1, 2),
#'          query_group = "cluster",
#'          ref_group   = "celltype")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d stat_density_2d
#'   scale_color_manual scale_fill_gradientn scale_x_continuous
#'   scale_y_continuous labs theme guide_legend
#' @export
PlotMAP <- function(
    ref,
    query,
    ref_emb                 = NULL,
    query_emb               = NULL,
    ref_group               = NULL,
    query_group             = NULL,
    # --- ref coordinate specification ---
    ref_reduction           = NULL,
    ref_dims                = c(1, 2),
    # --- query coordinate specification ---
    query_reduction         = NULL,
    query_dims              = c(1, 2),
    # --- aesthetics ---
    query_palette           = "Set1",
    query_palcolor          = NULL,
    query_color             = "firebrick",
    ref_palette             = "Paired",
    ref_palcolor            = NULL,
    ref_color               = "grey70",
    query.size              = 0.8,
    query.alpha             = 1,
    stroke                  = 0.5,
    stroke.color            = "black",
    point.size              = NULL,
    point.alpha             = 0.4,
    show_density            = FALSE,
    density.color           = "red",
    density_filled          = FALSE,
    density_filled_palette  = "Greys",
    density_filled_palcolor = NULL,
    xlim                    = NULL,
    ylim                    = NULL,
    legend.position         = "right",
    legend.direction        = "vertical",
    title                   = NULL,
    subtitle                = NULL,
    xlab                    = NULL,
    ylab                    = NULL,
    theme_use               = NULL,
    theme_args              = list(),
    aspect.ratio            = 1,
    raster                  = NULL,
    raster_method           = c("scattermore", "rasterise"),
    raster.dpi              = c(512, 512),
    seed                    = 11) {

  set.seed(seed)

  # =========================================================================
  # Internal: extract 2-column embedding matrix from various input types
  # Priority: emb_cols (column names) > reduction name > auto-detect > dims
  # =========================================================================
  .extract_emb <- function(obj, emb_cols, reduction, dims, label) {
    # Helper to pull columns from a data.frame / matrix
    .cols_from_df <- function(df, cols, dims_fb, label) {
      if (!is.null(cols)) {
        if (!all(cols %in% colnames(df))) {
          missing <- cols[!cols %in% colnames(df)]
          stop(label, ": column(s) not found: ",
               paste(missing, collapse = ", "), call. = FALSE)
        }
        return(as.matrix(df[, cols, drop = FALSE]))
      }
      # fallback: column indices
      if (any(dims_fb > ncol(df))) {
        stop(label, ": 'dims' indices out of range (ncol = ", ncol(df), ").",
             call. = FALSE)
      }
      as.matrix(df[, dims_fb, drop = FALSE])
    }

    if (inherits(obj, "Seurat")) {
      if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("Package 'SeuratObject' required for Seurat input.", call. = FALSE)
      }
      # Mode 2: column names in meta.data
      if (!is.null(emb_cols)) {
        return(list(
          mat     = .cols_from_df(obj@meta.data, emb_cols, dims, label),
          meta    = obj@meta.data,
          red_key = NULL
        ))
      }
      # Mode 1: named reduction
      red_names <- names(obj@reductions)
      if (!is.null(reduction)) {
        if (!reduction %in% red_names) {
          stop(label, ": reduction '", reduction, "' not found. Available: ",
               paste(red_names, collapse = ", "), call. = FALSE)
        }
      } else {
        # Auto-detect
        umap_idx <- grep("umap|tsne", tolower(red_names))
        reduction <- if (length(umap_idx) > 0) red_names[umap_idx[1]] else tail(red_names, 1)
        message(label, ": using reduction '", reduction, "'")
      }
      red_obj <- obj@reductions[[reduction]]
      mat <- red_obj@cell.embeddings[, dims, drop = FALSE]
      return(list(
        mat     = mat,
        meta    = obj@meta.data,
        red_key = red_obj@key
      ))

    } else if (is.data.frame(obj)) {
      # Mode 2 or 3: column names or indices from data.frame
      mat <- .cols_from_df(obj, emb_cols, dims, label)
      return(list(
        mat     = mat,
        meta    = obj,
        red_key = NULL
      ))
    } else {
      stop(label, ": must be a Seurat object or data.frame.", call. = FALSE)
    }
  }

  # =========================================================================
  # 1. Extract coordinates
  # =========================================================================
  ref_parsed <- .extract_emb(ref, ref_emb, ref_reduction, ref_dims, "ref")
  ref_mat    <- ref_parsed$mat
  ref_meta   <- ref_parsed$meta
  red_key    <- ref_parsed$red_key

  query_parsed  <- .extract_emb(query, query_emb, query_reduction, query_dims, "query")
  query_mat     <- query_parsed$mat
  query_meta    <- query_parsed$meta

  # =========================================================================
  # 2. Axis labels (auto from reduction key or column names)
  # =========================================================================
  if (is.null(xlab)) {
    xlab <- if (!is.null(red_key)) paste0(red_key, ref_dims[1]) else colnames(ref_mat)[1]
  }
  if (is.null(ylab)) {
    ylab <- if (!is.null(red_key)) paste0(red_key, ref_dims[2]) else colnames(ref_mat)[2]
  }

  # =========================================================================
  # 3. Build tidy data frames
  # =========================================================================
  ref_df <- data.frame(x = ref_mat[, 1], y = ref_mat[, 2],
                       row.names = rownames(ref_mat))
  if (!is.null(ref_group)) {
    if (!ref_group %in% colnames(ref_meta)) {
      stop("ref_group '", ref_group, "' not found in reference metadata.",
           call. = FALSE)
    }
    ref_df$group.by <- factor(ref_meta[rownames(ref_df), ref_group])
  }

  query_df <- data.frame(x = query_mat[, 1], y = query_mat[, 2],
                         row.names = rownames(query_mat))
  if (!is.null(query_group)) {
    if (!query_group %in% colnames(query_meta)) {
      stop("query_group '", query_group, "' not found in query metadata.",
           call. = FALSE)
    }
    query_df$group.by <- factor(query_meta[rownames(query_df), query_group])
  }

  # =========================================================================
  # 4. Colour palettes
  # =========================================================================
  if (!is.null(ref_group)) {
    ref_colors <- palette_colors(levels(ref_df$group.by),
                                 palette = ref_palette, palcolor = ref_palcolor,
                                 NA_keep = TRUE)
  }
  if (!is.null(query_group)) {
    query_colors <- palette_colors(levels(query_df$group.by),
                                   palette = query_palette, palcolor = query_palcolor,
                                   NA_keep = TRUE)
  }

  # =========================================================================
  # 5. Axis limits, point size, raster
  # =========================================================================
  if (is.null(xlim)) xlim <- range(c(ref_df$x, query_df$x), na.rm = TRUE)
  if (is.null(ylim)) ylim <- range(c(ref_df$y, query_df$y), na.rm = TRUE)

  if (is.null(point.size)) point.size <- min(3000 / nrow(ref_df), 0.5)

  raster        <- raster %||% ((nrow(ref_df) + nrow(query_df)) > 1e5)
  raster_method <- match.arg(raster_method)
  if (!is.null(raster.dpi) && length(raster.dpi) == 2) raster.dpi <- max(raster.dpi)

  # Point geom helper
  # raster_method = "rasterise"  : faithful colours (default)
  # raster_method = "scattermore": faster but may shift colours
  .pt <- function(data, mapping, size_val, alpha_val, color = NULL, ...) {
    if (isTRUE(raster) && raster_method == "scattermore") {
      if (!requireNamespace("scattermore", quietly = TRUE))
        stop("Package 'scattermore' is required for raster_method = 'scattermore'.",
             call. = FALSE)
      args <- list(data = data, mapping = mapping,
                   pointsize = ceiling(size_val), alpha = alpha_val,
                   pixels = c(raster.dpi, raster.dpi), ...)
      if (!is.null(color)) args$color <- color
      return(do.call(scattermore::geom_scattermore, args))
    }
    args <- list(data = data, mapping = mapping,
                 size = size_val, alpha = alpha_val, ...)
    if (!is.null(color)) args$color <- color
    layer <- do.call(ggplot2::geom_point, args)
    if (isTRUE(raster)) rasterise_layer(layer, dpi = raster.dpi) else layer
  }

  # =========================================================================
  # 6. Theme
  # =========================================================================
  if (is.null(theme_use)) {
    if (requireNamespace("UtilsR", quietly = TRUE)) {
      theme_fn <- UtilsR::theme_blank; is_blank <- TRUE
    } else {
      theme_fn <- ggplot2::theme_bw;   is_blank <- FALSE
    }
  } else if (is.character(theme_use)) {
    theme_fn <- tryCatch(match.fun(theme_use), error = function(e) {
      if (requireNamespace("UtilsR", quietly = TRUE))
        tryCatch(utils::getFromNamespace(theme_use, "UtilsR"),
                 error = function(e2) ggplot2::theme_bw)
      else ggplot2::theme_bw
    })
    is_blank <- identical(theme_use, "theme_blank")
  } else if (is.function(theme_use)) {
    theme_fn <- theme_use
    is_blank <- isTRUE(attr(theme_use, "blank")) ||
      (requireNamespace("UtilsR", quietly = TRUE) &&
         identical(theme_fn, UtilsR::theme_blank))
  } else if (inherits(theme_use, "theme")) {
    theme_fn <- function(...) theme_use; is_blank <- FALSE
  } else {
    stop("'theme_use' must be NULL, a character, a function, or a theme object.",
         call. = FALSE)
  }
  if (isTRUE(is_blank)) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  # =========================================================================
  # 7. Shuffle query for unbiased foreground rendering
  # =========================================================================
  query_df <- query_df[sample(nrow(query_df)), , drop = FALSE]

  # =========================================================================
  # 8. Build plot
  # =========================================================================
  p <- ggplot2::ggplot()

  # --- Reference layer ------------------------------------------------------
  if (!is.null(ref_group)) {
    p <- p +
      .pt(ref_df,
          ggplot2::aes(x = .data$x, y = .data$y, color = .data$group.by),
          size_val = point.size, alpha_val = point.alpha) +
      ggplot2::scale_color_manual(
        name     = paste0(ref_group, " (Ref):"),
        values   = ref_colors, na.value = "grey80",
        guide    = ggplot2::guide_legend(
          title.hjust = 0, order = 1,
          override.aes = list(size = 4, alpha = 1))
      )
  } else {
    p <- p +
      .pt(ref_df, ggplot2::aes(x = .data$x, y = .data$y),
          size_val = point.size, alpha_val = point.alpha, color = ref_color)
  }

  # --- Density contour for query (optional) ---------------------------------
  if (isTRUE(show_density)) {
    if (isTRUE(density_filled)) {
      if (!requireNamespace("ggnewscale", quietly = TRUE))
        stop("Package 'ggnewscale' required for density_filled.", call. = FALSE)
      filled_cols <- palette_colors(palette  = density_filled_palette,
                                    palcolor = density_filled_palcolor)
      p <- p +
        ggplot2::stat_density_2d(
          data = query_df,
          mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                 fill = ggplot2::after_stat(density)),
          geom = "raster", contour = FALSE,
          inherit.aes = FALSE, show.legend = FALSE) +
        ggplot2::scale_fill_gradientn(colours = filled_cols) +
        ggnewscale::new_scale_fill()
    } else {
      p <- p +
        ggplot2::geom_density_2d(
          data = query_df,
          mapping = ggplot2::aes(x = .data$x, y = .data$y),
          color = density.color, inherit.aes = FALSE, show.legend = FALSE)
    }
  }

  # --- Separate colour scale for query overlay ------------------------------
  if (!is.null(ref_group)) p <- p + ggnewscale::new_scale_color()

  # Border stroke
  if (stroke > 0) {
    p <- p +
      .pt(query_df, ggplot2::aes(x = .data$x, y = .data$y),
          size_val = query.size + stroke, alpha_val = query.alpha,
          color = stroke.color, inherit.aes = FALSE)
  }

  # Query points
  if (!is.null(query_group)) {
    p <- p +
      .pt(query_df,
          ggplot2::aes(x = .data$x, y = .data$y, color = .data$group.by),
          size_val = query.size, alpha_val = query.alpha,
          inherit.aes = FALSE) +
      ggplot2::scale_color_manual(
        name   = paste0(query_group, " (Query):"),
        values = query_colors, na.value = "grey80",
        guide  = ggplot2::guide_legend(
          title.hjust = 0, order = 2,
          override.aes = list(
            size  = 4, alpha = 1, shape = 21, color = "black",
            fill  = query_colors[seq_along(levels(query_df$group.by))]
          ))
      )
  } else {
    p <- p +
      .pt(query_df, ggplot2::aes(x = .data$x, y = .data$y),
          size_val = query.size, alpha_val = query.alpha,
          color = query_color, inherit.aes = FALSE)
  }

  # --- Final scales + theme -------------------------------------------------
  p <- p +
    ggplot2::scale_x_continuous(limits = xlim) +
    ggplot2::scale_y_continuous(limits = ylim) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    do.call(theme_fn, theme_args) +
    ggplot2::theme(
      aspect.ratio     = aspect.ratio,
      legend.position  = legend.position,
      legend.direction = legend.direction
    )

  p
}
