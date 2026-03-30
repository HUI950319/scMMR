# ============================================================================
# DimPlot2.R -- Enhanced 2D dimensional reduction plot
# ============================================================================

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Validate and prepare input data for DimPlot2
#' @noRd
.dimplot2_prepare <- function(seu, group.by, split.by, reduction, dims, cells,
                              show_na, force) {
  # --- Accept Seurat or data.frame ---
  if (inherits(seu, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Package 'SeuratObject' is required when 'seu' is a Seurat object.",
           call. = FALSE)
    }
    # Resolve reduction
    if (is.null(reduction)) {
      red_names <- names(seu@reductions)
      reduction <- red_names[grep("umap|tsne", tolower(red_names))][1]
      if (is.na(reduction)) reduction <- red_names[length(red_names)]
    }
    if (!reduction %in% names(seu@reductions)) {
      stop("Reduction '", reduction, "' not found. Available: ",
           paste(names(seu@reductions), collapse = ", "), call. = FALSE)
    }
    reduction_key <- seu@reductions[[reduction]]@key
    dat_dim <- seu@reductions[[reduction]]@cell.embeddings
    colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
    if (is.null(rownames(dat_dim))) {
      rownames(dat_dim) <- colnames(seu)
    }
    dat_meta <- seu@meta.data
  } else if (is.data.frame(seu)) {
    # Expect columns: dim1, dim2 (or user-specified), plus group.by columns
    reduction_key <- "Dim"
    reduction <- "custom"
    dim_cols <- colnames(seu)[dims]
    if (any(is.na(dim_cols))) {
      stop("'dims' indices out of range for the input data.frame columns.",
           call. = FALSE)
    }
    dat_dim <- as.matrix(seu[, dims, drop = FALSE])
    colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
    rownames(dat_dim) <- rownames(seu)
    dat_meta <- seu
  } else {
    stop("'seu' must be a Seurat object or a data.frame.", call. = FALSE)
  }

  # --- Ensure factors ---
  all_groups <- unique(c(group.by, split.by))
  all_groups <- all_groups[!is.null(all_groups)]
  for (col in all_groups) {
    if (!col %in% colnames(dat_meta)) {
      stop("Column '", col, "' not found in metadata.", call. = FALSE)
    }
    if (!is.factor(dat_meta[[col]])) {
      dat_meta[[col]] <- factor(dat_meta[[col]],
                                 levels = unique(dat_meta[[col]]))
    }
    if (isTRUE(show_na) && any(is.na(dat_meta[[col]]))) {
      raw_levels <- unique(c(levels(dat_meta[[col]]), "NA"))
      dat_meta[[col]] <- as.character(dat_meta[[col]])
      dat_meta[[col]][is.na(dat_meta[[col]])] <- "NA"
      dat_meta[[col]] <- factor(dat_meta[[col]], levels = raw_levels)
    }
  }

  # --- Add split.by dummy if NULL ---
  if (is.null(split.by)) {
    split.by <- ".split_all"
    dat_meta[[split.by]] <- factor("")
  }

  # --- Combine dim + meta ---
  shared <- intersect(rownames(dat_dim), rownames(dat_meta))
  dat_use <- cbind(
    as.data.frame(dat_dim[shared, , drop = FALSE]),
    dat_meta[shared, unique(c(group.by, split.by)), drop = FALSE]
  )
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }

  # --- Force check ---
  nlev <- sapply(dat_use[, group.by, drop = FALSE], nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && isFALSE(force)) {
    warning("Variables with >100 levels: ",
            paste(names(nlev), collapse = ", "),
            ". Use force = TRUE to proceed.", call. = FALSE)
    return(NULL)
  }

  list(
    dat_use       = dat_use,
    dat_dim       = as.data.frame(dat_dim[shared, , drop = FALSE]),
    reduction     = reduction,
    reduction_key = reduction_key,
    split.by      = split.by
  )
}


#' Build label text for groups
#' @noRd
.build_label_text <- function(labels_tb, label, label_insitu, show_stat) {
  prefix <- if (!isTRUE(label_insitu) && isTRUE(label)) {
    paste0(seq_along(labels_tb), ": ")
  } else {
    ""
  }
  suffix <- if (isTRUE(show_stat)) paste0("(", labels_tb, ")") else ""
  paste0(prefix, names(labels_tb), suffix)
}


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

#' @title Enhanced 2D Dimensional Reduction Plot
#'
#' @description
#' Visualize cell groups on a 2-dimensional reduction plot (UMAP, t-SNE, etc.).
#' Supports Seurat objects and plain data.frames. Provides label, highlight,
#' density, mark, and raster layers with flexible theming.
#'
#' Uses \code{UtilsR} theme / grob utilities and \code{scMMR::palette_colors}
#' for color mapping.
#'
#' @md
#' @param seu A Seurat object or a \code{data.frame}. When a data.frame is
#'   provided, the first two columns specified by \code{dims} are treated as
#'   embedding coordinates, and \code{group.by} / \code{split.by} columns must
#'   be present.
#' @param group.by Character. Name(s) of column(s) in metadata to group
#'   (color) cells by.
#' @param reduction Character. Name of the dimensionality reduction to use.
#'   Ignored when \code{seu} is a data.frame. Default auto-detects UMAP/tSNE.
#' @param dims Integer vector of length 2. Dimensions to plot.
#'   For Seurat input: component indices (default \code{c(1, 2)}).
#'   For data.frame input: column indices for x/y coordinates.
#' @param split.by Character. Column name to split into facets. Default
#'   \code{NULL}.
#' @param cells Character vector of cell/row names to include. Default
#'   \code{NULL} (all).
#' @param show_na Logical. Show NA group with \code{bg_color}. Default
#'   \code{FALSE}.
#' @param show_stat Logical. Show cell counts in labels/subtitles. Default
#'   \code{NULL} (auto: \code{FALSE} when theme is blank, \code{TRUE}
#'   otherwise).
#' @param pt.size Point size. Default auto-calculated.
#' @param pt.alpha Point transparency. Default \code{1}.
#' @param palette Character. Color palette name (passed to
#'   \code{palette_colors}). Default \code{"Paired"}.
#' @param palcolor Character vector. Custom colors. Default \code{NULL}.
#' @param bg_color Background color for NA points. Default \code{"grey80"}.
#' @param label Logical. Add group labels. Default \code{FALSE}.
#' @param label.size Label font size. Default \code{4}.
#' @param label.fg Label foreground color. Default \code{"white"}.
#' @param label.bg Label background color. Default \code{"black"}.
#' @param label.bg.r Label background ratio. Default \code{0.1}.
#' @param label_insitu Logical. Place raw group names at cell centers.
#'   Default \code{FALSE}.
#' @param label_repel Logical. Repel labels. Default \code{FALSE}.
#' @param label_repulsion Repulsion force. Default \code{20}.
#' @param label_point_size Center point size for labels. Default \code{1}.
#' @param label_point_color Center point color. Default \code{"black"}.
#' @param label_segment_color Segment color for repelled labels. Default
#'   \code{"black"}.
#' @param cells.highlight Character vector of cell names to highlight, or
#'   \code{TRUE} for all non-NA cells. Default \code{NULL}.
#' @param cols.highlight Highlight border color. Default \code{"black"}.
#' @param sizes.highlight Highlight point size. Default \code{1}.
#' @param alpha.highlight Highlight transparency. Default \code{1}.
#' @param stroke.highlight Highlight border width. Default \code{0.5}.
#' @param add_density Logical. Add density contours. Default \code{FALSE}.
#' @param density_color Contour line color. Default \code{"grey80"}.
#' @param density_filled Logical. Use filled contour bands. Default
#'   \code{FALSE}.
#' @param density_filled_palette Palette for filled contours. Default
#'   \code{"Greys"}.
#' @param density_filled_palcolor Custom colors for filled contours. Default
#'   \code{NULL}.
#' @param add_mark Logical. Add group boundary marks. Default \code{FALSE}.
#' @param mark_type Mark shape: \code{"hull"}, \code{"ellipse"},
#'   \code{"rect"}, or \code{"circle"}. Default \code{"hull"}.
#' @param mark_expand Mark expansion. Default \code{grid::unit(3, "mm")}.
#' @param mark_alpha Mark transparency. Default \code{0.1}.
#' @param mark_linetype Mark border line type. Default \code{1}.
#' @param lineages Character vector of metadata column names containing
#'   pseudotime values (e.g. from \code{RunSlingshot}). Each column is one
#'   lineage; \code{NULL} = no lineage overlay. Default \code{NULL}.
#' @param lineages_trim Numeric pair. Quantile range of pseudotime cells to
#'   include in LOESS fit (removes extremes). Default \code{c(0.01, 0.99)}.
#' @param lineages_span LOESS smoother span. Larger = smoother. Default
#'   \code{0.75}.
#' @param lineages_palette Colour palette for lineage lines. Default
#'   \code{"Dark2"}.
#' @param lineages_palcolor Custom colour vector for lineages. Default
#'   \code{NULL}.
#' @param lineages_arrow Arrow specification for path end. Default
#'   \code{grid::arrow(length = grid::unit(0.1, "inches"))}.
#' @param lineages_linewidth Foreground path linewidth. Default \code{1}.
#' @param lineages_line_bg Colour of the outline stroke behind each path.
#'   Default \code{"white"}.
#' @param lineages_line_bg_stroke Extra width added to the outline. Default
#'   \code{0.5}.
#' @param lineages_whiskers Logical. Draw segments from each cell to the
#'   smooth curve. Default \code{FALSE}.
#' @param lineages_whiskers_linewidth Whisker segment linewidth. Default
#'   \code{0.5}.
#' @param lineages_whiskers_alpha Whisker transparency. Default \code{0.5}.
#' @param hex Logical. Use hexagonal binning. Default \code{FALSE}.
#' @param hex.count Logical. Map hex alpha to count. Default \code{TRUE}.
#' @param hex.bins Number of hex bins. Default \code{50}.
#' @param hex.binwidth Hex bin width. Default \code{NULL}.
#' @param hex.linewidth Hex border width. Default \code{0.5}.
#' @param raster Logical. Rasterize points (auto if >100k cells). Default
#'   \code{NULL}.
#' @param raster_method Character. Rasterisation backend:
#'   \code{"rasterise"} (default, faithful colours via \code{rasterise_layer})
#'   or \code{"scattermore"} (faster C-level rendering, may shift colours).
#' @param raster.dpi Numeric. Raster resolution in DPI. Accepts a scalar or
#'   a two-length vector (max is used). Default \code{c(512, 512)}.
#' @param aspect.ratio Aspect ratio. Default \code{1}.
#' @param title Plot title. Default \code{NULL}.
#' @param subtitle Plot subtitle. Default \code{NULL}.
#' @param xlab X-axis label. Default auto from reduction key.
#' @param ylab Y-axis label. Default auto from reduction key.
#' @param legend.position Legend position. Default \code{"right"}.
#' @param legend.direction Legend direction. Default \code{"vertical"}.
#' @param legend.title Legend title. Default \code{NULL} (uses group name).
#' @param theme_use Theme specification. Accepts a function (e.g.
#'   \code{UtilsR::theme_blank}), a character name (e.g. \code{"theme_sc"}),
#'   or a \code{theme} object (e.g. \code{theme_bw()}). Default \code{NULL}
#'   which uses \code{UtilsR::theme_blank}.
#' @param theme_args List of extra arguments to \code{theme_use}. Default
#'   \code{list()}.
#' @param combine Logical. Combine multiple panels via patchwork. Default
#'   \code{TRUE}.
#' @param nrow,ncol Layout dimensions for combined plot. Default \code{NULL}.
#' @param byrow Arrange panels by row. Default \code{TRUE}.
#' @param force Logical. Force plot even with >100 group levels. Default
#'   \code{FALSE}.
#' @param seed Random seed. Default \code{11}.
#'
#' @return A ggplot object (or list of ggplot objects if \code{combine =
#'   FALSE}).
#'
#' @examples
#' \dontrun{
#' # --- Seurat input ---
#' DimPlot2(seu, group.by = "celltype")
#' DimPlot2(seu, group.by = "celltype", label = TRUE, add_mark = TRUE)
#'
#' # --- data.frame input ---
#' df <- data.frame(UMAP1 = rnorm(200), UMAP2 = rnorm(200),
#'                  cluster = sample(letters[1:5], 200, replace = TRUE))
#' DimPlot2(df, group.by = "cluster", dims = c(1, 2))
#' }
#'
#' @export
DimPlot2 <- function(
    seu,
    group.by,
    reduction = NULL,
    dims = c(1, 2),
    split.by = NULL,
    cells = NULL,
    show_na = FALSE,
    show_stat = NULL,
    pt.size = NULL,
    pt.alpha = 1,
    palette = "Paired",
    palcolor = NULL,
    bg_color = "grey80",
    label = FALSE,
    label.size = 4,
    label.fg = "white",
    label.bg = "black",
    label.bg.r = 0.1,
    label_insitu = FALSE,
    label_repel = FALSE,
    label_repulsion = 20,
    label_point_size = 1,
    label_point_color = "black",
    label_segment_color = "black",
    cells.highlight = NULL,
    cols.highlight = "black",
    sizes.highlight = 1,
    alpha.highlight = 1,
    stroke.highlight = 0.5,
    add_density = FALSE,
    density_color = "grey80",
    density_filled = FALSE,
    density_filled_palette = "Greys",
    density_filled_palcolor = NULL,
    add_mark = FALSE,
    mark_type = c("hull", "ellipse", "rect", "circle"),
    mark_expand = grid::unit(3, "mm"),
    mark_alpha = 0.1,
    mark_linetype = 1,
    lineages = NULL,
    lineages_trim = c(0.01, 0.99),
    lineages_span = 0.75,
    lineages_palette = "Dark2",
    lineages_palcolor = NULL,
    lineages_arrow = grid::arrow(length = grid::unit(0.1, "inches")),
    lineages_linewidth = 1,
    lineages_line_bg = "white",
    lineages_line_bg_stroke = 0.5,
    lineages_whiskers = FALSE,
    lineages_whiskers_linewidth = 0.5,
    lineages_whiskers_alpha = 0.5,
    hex = FALSE,
    hex.count = TRUE,
    hex.bins = 50,
    hex.binwidth = NULL,
    hex.linewidth = 0.5,
    raster        = NULL,
    raster_method = c("scattermore", "rasterise"),
    raster.dpi    = c(512, 512),
    aspect.ratio = 1,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = NULL,
    theme_use = NULL,
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    force = FALSE,
    seed = 11) {

  set.seed(seed)
  mark_type <- match.arg(mark_type)

  # ---- Prepare data ----
  prep <- .dimplot2_prepare(
    seu = seu, group.by = group.by, split.by = split.by,
    reduction = reduction, dims = dims, cells = cells,
    show_na = show_na, force = force
  )
  if (is.null(prep)) return(invisible(NULL))

  dat_use       <- prep$dat_use
  dat_dim       <- prep$dat_dim
  reduction_key <- prep$reduction_key
  split.by      <- prep$split.by
  dim_x <- paste0(reduction_key, dims[1])
  dim_y <- paste0(reduction_key, dims[2])

  # ---- Lineage layers (computed once, reused per panel) ----
  lineage_layers <- NULL
  if (!is.null(lineages)) {
    if (!requireNamespace("ggnewscale", quietly = TRUE)) {
      stop("Package 'ggnewscale' is required for lineages.", call. = FALSE)
    }
    raw_meta <- if (inherits(seu, "Seurat")) seu@meta.data else seu
    missing_lin <- setdiff(lineages, colnames(raw_meta))
    if (length(missing_lin) > 0) {
      stop("Lineage column(s) not found in metadata: ",
           paste(missing_lin, collapse = ", "), call. = FALSE)
    }
    lin_shared <- intersect(rownames(prep$dat_dim), rownames(raw_meta))
    lin_dat    <- as.data.frame(prep$dat_dim[lin_shared, , drop = FALSE])
    for (l in lineages) lin_dat[[l]] <- raw_meta[lin_shared, l]
    if (!is.null(cells)) {
      lin_dat <- lin_dat[intersect(rownames(lin_dat), cells), , drop = FALSE]
    }
    lineage_layers <- c(
      list(ggnewscale::new_scale_color()),
      .compute_lineage_layers(
        dat                = lin_dat,
        x_col              = dim_x,
        y_col              = dim_y,
        lineages           = lineages,
        trim               = lineages_trim,
        span               = lineages_span,
        palette            = lineages_palette,
        palcolor           = lineages_palcolor,
        lineages_arrow     = lineages_arrow,
        linewidth          = lineages_linewidth,
        line_bg            = lineages_line_bg,
        line_bg_stroke     = lineages_line_bg_stroke,
        whiskers           = lineages_whiskers,
        whiskers_linewidth = lineages_whiskers_linewidth,
        whiskers_alpha     = lineages_whiskers_alpha
      )
    )
  }

  # ---- Auto pt.size & raster ----
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  raster        <- raster %||% (nrow(dat_use) > 1e5)
  raster_method <- match.arg(raster_method)
  if (!is.null(raster.dpi)) {
    if (is.numeric(raster.dpi) && length(raster.dpi) == 2) {
      raster.dpi <- max(raster.dpi)
    }
    if (!is.numeric(raster.dpi) || length(raster.dpi) != 1) {
      stop("'raster.dpi' must be a numeric scalar or two-length vector.",
           call. = FALSE)
    }
  } else {
    raster.dpi <- 300
  }

  # ---- Resolve theme function ----
  # Default: theme_blank from UtilsR
  if (is.null(theme_use)) {
    if (requireNamespace("UtilsR", quietly = TRUE)) {
      theme_fn <- UtilsR::theme_blank
    } else {
      theme_fn <- ggplot2::theme_bw
    }
    is_blank <- TRUE
  } else if (is.character(theme_use)) {
    theme_fn <- tryCatch(match.fun(theme_use), error = function(e) {
      if (requireNamespace("UtilsR", quietly = TRUE)) {
        tryCatch(utils::getFromNamespace(theme_use, "UtilsR"),
                 error = function(e2) ggplot2::theme_bw)
      } else {
        ggplot2::theme_bw
      }
    })
    is_blank <- identical(theme_use, "theme_blank")
  } else if (is.function(theme_use)) {
    theme_fn <- theme_use
    is_blank <- identical(theme_fn, UtilsR::theme_blank) ||
                isTRUE(attr(theme_use, "blank"))
  } else if (inherits(theme_use, "theme")) {
    # Already a theme object (e.g. theme_bw())
    theme_fn <- function(...) theme_use
    is_blank <- FALSE
  } else {
    stop("'theme_use' must be NULL, a character string, a function, or a theme object.",
         call. = FALSE)
  }
  # show_stat default: FALSE for blank theme, TRUE otherwise
  if (is.null(show_stat)) show_stat <- !is_blank

  # ---- Axis labels ----
  xlab <- xlab %||% dim_x
  ylab <- ylab %||% dim_y
  if (is_blank) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  # ---- Highlight cells prep ----
  if (is.character(cells.highlight)) {
    cells.highlight <- intersect(cells.highlight, rownames(dat_use))
    if (length(cells.highlight) == 0) cells.highlight <- NULL
  }

  # ---- Check optional packages ----
  if (isTRUE(add_mark)) {
    if (!requireNamespace("ggforce", quietly = TRUE)) {
      stop("Package 'ggforce' is required for add_mark = TRUE.", call. = FALSE)
    }
    if (!requireNamespace("ggnewscale", quietly = TRUE)) {
      stop("Package 'ggnewscale' is required for add_mark = TRUE.", call. = FALSE)
    }
  }
  if (isTRUE(add_density) && isTRUE(density_filled)) {
    if (!requireNamespace("ggnewscale", quietly = TRUE)) {
      stop("Package 'ggnewscale' is required for density_filled.", call. = FALSE)
    }
  }

  # ---- Helper: choose point geom ----
  # raster_method = "rasterise"  : faithful colours (default)
  # raster_method = "scattermore": faster but may shift colours
  .point_geom <- function(data, mapping, size_val, alpha_val, color = NULL, ...) {
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

  # ---- Build per-panel plots ----
  comb <- expand.grid(
    split = levels(dat_use[[split.by]]),
    group = group.by,
    stringsAsFactors = FALSE
  )
  rownames(comb) <- paste0(comb[["split"]], ":", comb[["group"]])

  # Axis limits from full data
  xlim_range <- range(dat_dim[[dim_x]], na.rm = TRUE)
  ylim_range <- range(dat_dim[[dim_y]], na.rm = TRUE)

  plist <- lapply(
    stats::setNames(rownames(comb), rownames(comb)),
    function(i) {
      g <- comb[i, "group"]
      s <- comb[i, "split"]

      # -- Colors --
      colors <- palette_colors(
        levels(dat_use[[g]]),
        palette = palette,
        palcolor = palcolor,
        NA_keep = TRUE
      )

      # -- Subset by split --
      dat <- dat_use
      cells_mask <- dat[[split.by]] != s
      dat[[g]][cells_mask] <- NA
      if (isFALSE(show_na)) {
        dat <- dat[!is.na(dat[[g]]), , drop = FALSE]
      }

      labels_tb <- table(dat[[g]])
      labels_tb <- labels_tb[labels_tb != 0]
      label_use <- .build_label_text(labels_tb, label, label_insitu, show_stat)

      # -- Prepare plot columns --
      dat[["x"]] <- dat[[dim_x]]
      dat[["y"]] <- dat[[dim_y]]
      dat[["group.by"]] <- dat[[g]]
      dat[, "split.by"] <- s

      # Sort: NA first (background), then shuffle foreground
      dat <- dat[order(dat[["group.by"]], decreasing = FALSE, na.last = FALSE), ,
                 drop = FALSE]
      na_end <- max(c(1, which(is.na(dat[["group.by"]]))))
      fg_idx <- seq(min(na_end + 1, nrow(dat)), nrow(dat))
      if (length(fg_idx) > 1) fg_idx <- sample(fg_idx)
      dat <- dat[c(seq_len(na_end), fg_idx), , drop = FALSE]

      # -- Subtitle --
      subtitle_use <- if (isTRUE(show_stat)) {
        subtitle %||% paste0(s, " nCells:", sum(!is.na(dat[["group.by"]])))
      } else {
        subtitle
      }

      # == Mark layer ==
      mark_layer <- NULL
      if (isTRUE(add_mark)) {
        mark_fun <- switch(mark_type,
          "ellipse" = ggforce::geom_mark_ellipse,
          "hull"    = ggforce::geom_mark_hull,
          "rect"    = ggforce::geom_mark_rect,
          "circle"  = ggforce::geom_mark_circle
        )
        mark_layer <- list(
          mark_fun(
            data = dat[!is.na(dat[["group.by"]]), , drop = FALSE],
            mapping = ggplot2::aes(
              x = .data[["x"]], y = .data[["y"]],
              color = .data[["group.by"]], fill = .data[["group.by"]]
            ),
            expand = mark_expand, alpha = mark_alpha,
            linetype = mark_linetype, show.legend = FALSE,
            inherit.aes = FALSE
          ),
          ggplot2::scale_fill_manual(values = colors[names(labels_tb)]),
          ggplot2::scale_color_manual(values = colors[names(labels_tb)]),
          ggnewscale::new_scale_fill(),
          ggnewscale::new_scale_color()
        )
      }

      # == Density layer ==
      density_layer <- NULL
      if (isTRUE(add_density)) {
        if (isTRUE(density_filled)) {
          filled_cols <- palette_colors(
            palette = density_filled_palette,
            palcolor = density_filled_palcolor
          )
          density_layer <- list(
            ggplot2::stat_density_2d(
              geom = "raster",
              ggplot2::aes(
                x = .data[["x"]], y = .data[["y"]],
                fill = ggplot2::after_stat(density)
              ),
              contour = FALSE, inherit.aes = FALSE, show.legend = FALSE
            ),
            ggplot2::scale_fill_gradientn(name = "Density", colours = filled_cols),
            ggnewscale::new_scale_fill()
          )
        } else {
          density_layer <- ggplot2::geom_density_2d(
            ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
            color = density_color, inherit.aes = FALSE, show.legend = FALSE
          )
        }
      }

      # == Base plot ==
      p <- ggplot2::ggplot(dat) +
        mark_layer +
        density_layer +
        ggplot2::labs(title = title, subtitle = subtitle_use,
                      x = xlab, y = ylab) +
        ggplot2::scale_x_continuous(limits = xlim_range) +
        ggplot2::scale_y_continuous(limits = ylim_range) +
        do.call(theme_fn, theme_args) +
        ggplot2::theme(
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction
        )

      # == Facet ==
      if (split.by != ".split_all") {
        p <- p + ggplot2::facet_grid(. ~ split.by)
      }

      # == Point layer (raster / hex / standard) ==
      if (isTRUE(raster)) {
        # Background NA points
        p <- p +
          .point_geom(
            data = dat[is.na(dat[["group.by"]]), , drop = FALSE],
            mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
            size_val = pt.size, alpha_val = pt.alpha, color = bg_color
          ) +
          .point_geom(
            data = dat[!is.na(dat[["group.by"]]), , drop = FALSE],
            mapping = ggplot2::aes(
              x = .data[["x"]], y = .data[["y"]],
              color = .data[["group.by"]]
            ),
            size_val = pt.size, alpha_val = pt.alpha
          )
      } else if (isTRUE(hex)) {
        if (!requireNamespace("hexbin", quietly = TRUE)) {
          stop("Package 'hexbin' is required for hex = TRUE.", call. = FALSE)
        }
        hex_aes <- ggplot2::aes(
          x = .data[["x"]], y = .data[["y"]],
          fill = .data[["group.by"]], color = .data[["group.by"]]
        )
        if (isTRUE(hex.count)) {
          hex_aes$alpha <- ggplot2::aes(alpha = ggplot2::after_stat(count))$alpha
        }
        p <- p + ggplot2::geom_hex(
          mapping = hex_aes, linewidth = hex.linewidth,
          bins = hex.bins, binwidth = hex.binwidth
        )
      } else {
        p <- p + ggplot2::geom_point(
          mapping = ggplot2::aes(
            x = .data[["x"]], y = .data[["y"]],
            color = .data[["group.by"]]
          ),
          size = pt.size, alpha = pt.alpha
        )
      }

      # == Highlight layer ==
      cells_hl <- cells.highlight
      if (isTRUE(cells_hl)) {
        cells_hl <- rownames(dat)[!is.na(dat[["group.by"]])]
      }
      if (!is.null(cells_hl) && isFALSE(hex)) {
        cell_df <- dat[rownames(dat) %in% cells_hl, , drop = FALSE]
        if (nrow(cell_df) > 0) {
          p <- p +
            .point_geom(
              data = cell_df,
              mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
              size_val = sizes.highlight + stroke.highlight,
              alpha_val = alpha.highlight, color = cols.highlight
            ) +
            .point_geom(
              data = cell_df,
              mapping = ggplot2::aes(
                x = .data[["x"]], y = .data[["y"]],
                color = .data[["group.by"]]
              ),
              size_val = sizes.highlight, alpha_val = alpha.highlight
            )
        }
      }

      # == Color scale ==
      legend_title_use <- if (is.null(legend.title)) {
        paste0(g, ":")
      } else {
        legend.title
      }
      p <- p + ggplot2::scale_color_manual(
        name = legend_title_use,
        values = colors[names(labels_tb)],
        labels = label_use,
        na.value = bg_color,
        guide = ggplot2::guide_legend(
          title.hjust = 0, order = 1,
          override.aes = list(size = 4, alpha = 1)
        )
      )
      if (isTRUE(hex)) {
        p <- p + ggplot2::scale_fill_manual(
          name = legend_title_use,
          values = colors[names(labels_tb)],
          labels = label_use,
          na.value = bg_color,
          guide = ggplot2::guide_legend(title.hjust = 0, order = 1)
        )
      }

      # == Lineage trajectory overlay ==
      if (!is.null(lineage_layers)) {
        p <- p + lineage_layers
      }

      # == Label layer ==
      if (isTRUE(label)) {
        if (!requireNamespace("ggrepel", quietly = TRUE)) {
          stop("Package 'ggrepel' is required for label = TRUE.", call. = FALSE)
        }
        label_df <- stats::aggregate(
          p$data[, c("x", "y")],
          by = list(p$data[["group.by"]]),
          FUN = stats::median
        )
        colnames(label_df)[1] <- "label"
        label_df <- label_df[!is.na(label_df[["label"]]), , drop = FALSE]
        if (isFALSE(label_insitu)) {
          label_df[["label"]] <- seq_len(nrow(label_df))
        }
        if (isTRUE(label_repel)) {
          p <- p +
            ggplot2::geom_point(
              data = label_df,
              mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
              color = label_point_color, size = label_point_size
            ) +
            ggrepel::geom_text_repel(
              data = label_df,
              ggplot2::aes(
                x = .data[["x"]], y = .data[["y"]],
                label = .data[["label"]]
              ),
              fontface = "bold", min.segment.length = 0,
              segment.color = label_segment_color,
              point.size = label_point_size, max.overlaps = 100,
              force = label_repulsion,
              color = label.fg, bg.color = label.bg, bg.r = label.bg.r,
              size = label.size, inherit.aes = FALSE
            )
        } else {
          p <- p +
            ggrepel::geom_text_repel(
              data = label_df,
              ggplot2::aes(
                x = .data[["x"]], y = .data[["y"]],
                label = .data[["label"]]
              ),
              fontface = "bold", min.segment.length = 0,
              segment.color = label_segment_color,
              point.size = NA, max.overlaps = 100, force = 0,
              color = label.fg, bg.color = label.bg, bg.r = label.bg.r,
              size = label.size, inherit.aes = FALSE
            )
        }
      }

      return(p)
    }
  )

  # ---- Combine ----
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop("Package 'patchwork' is required for combining plots.",
             call. = FALSE)
      }
      plot <- patchwork::wrap_plots(
        plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow
      )
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}
