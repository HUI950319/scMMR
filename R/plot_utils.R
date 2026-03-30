# ============================================================================
# Shared internal helpers for visualization functions
# ============================================================================

#' @keywords internal
.extract_cellmeta <- function(cellmeta) {
  if (inherits(cellmeta, "Seurat")) return(cellmeta@meta.data)
  if (is.data.frame(cellmeta)) return(cellmeta)
  stop("cellmeta must be a Seurat object, data.frame, or tibble.")
}

#' @keywords internal
.percentage_stat <- function(cellmeta, by, fill) {
  cellmeta <- .extract_cellmeta(cellmeta)
  if (!is.character(by) || !is.character(fill))
    stop("by and fill must be character vectors.")
  if (!(by %in% names(cellmeta)) || !(fill %in% names(cellmeta)))
    stop("by and fill must be columns in cellmeta.")

  data.stat <- as.data.frame(table(cellmeta[[by]], cellmeta[[fill]]))
  colnames(data.stat)[1:2] <- c(by, fill)
  data.stat <- data.stat %>%
    dplyr::group_by_at(by) %>%
    dplyr::mutate(margin.freq = sum(.data$Freq)) %>%
    dplyr::mutate(proportion = .data$Freq / .data$margin.freq)
  data.stat
}

# --------------------------------------------------------------------------
# Internal helper: auto-generate palette colours from RColorBrewer
# --------------------------------------------------------------------------

#' @keywords internal
.auto_palette <- function(fill_levels, palette) {
  n <- length(fill_levels)
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    stop("Package 'RColorBrewer' is required for palette colours.")
  max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  base_cols <- RColorBrewer::brewer.pal(min(max_n, max(3, n)), palette)
  setNames(grDevices::colorRampPalette(base_cols)(n), fill_levels)
}

# --------------------------------------------------------------------------
# Internal: compute LOESS-smoothed lineage trajectory layers
# --------------------------------------------------------------------------

#' Build ggplot2 layers for pseudotime lineage trajectories
#'
#' Fits a LOESS curve for each lineage column (pseudotime values stored in
#' metadata) and returns a list of ggplot2 layers (whiskers, background path,
#' coloured foreground path with arrow, and scale_color_manual).
#'
#' @param dat data.frame with x_col, y_col, and all lineage columns.
#' @param x_col,y_col Column names for the 2-D coordinates.
#' @param lineages Character vector: column names of pseudotime values.
#' @param trim Numeric [0,1] pair: quantile range to include. Default c(0.01,0.99).
#' @param span LOESS span. Default 0.75.
#' @param palette Colour palette for lineage colours. Default "Dark2".
#' @param palcolor Custom colour vector. NULL = use palette.
#' @param lineages_arrow Arrow specification for path ends.
#' @param linewidth Foreground path linewidth. Default 1.
#' @param line_bg Colour of the background stroke. Default "white".
#' @param line_bg_stroke Extra width added to background stroke. Default 0.5.
#' @param whiskers Logical. Draw per-cell whiskers to the smooth curve. Default FALSE.
#' @param whiskers_linewidth Whisker linewidth. Default 0.5.
#' @param whiskers_alpha Whisker alpha. Default 0.5.
#' @return A list of ggplot2 layers (ready to add directly to a plot).
#' @noRd
.compute_lineage_layers <- function(
    dat,
    x_col, y_col,
    lineages,
    trim                = c(0.01, 0.99),
    span                = 0.75,
    palette             = "Dark2",
    palcolor            = NULL,
    lineages_arrow      = grid::arrow(length = grid::unit(0.1, "inches")),
    linewidth           = 1,
    line_bg             = "white",
    line_bg_stroke      = 0.5,
    whiskers            = FALSE,
    whiskers_linewidth  = 0.5,
    whiskers_alpha      = 0.5) {

  colors <- palette_colors(lineages, palette = palette, palcolor = palcolor)

  # Fit LOESS for each lineage
  fitted_list <- lapply(stats::setNames(lineages, lineages), function(l) {
    pstime <- dat[[l]]
    q_lo   <- stats::quantile(pstime, trim[1], na.rm = TRUE)
    q_hi   <- stats::quantile(pstime, trim[2], na.rm = TRUE)
    idx    <- which(!is.na(pstime) & pstime > q_lo & pstime < q_hi)
    if (length(idx) < 4) return(NULL)
    idx    <- idx[order(pstime[idx])]
    xs     <- dat[idx, x_col]
    ys     <- dat[idx, y_col]
    pt     <- pstime[idx]
    wt     <- rep(1, length(idx))
    fit_x  <- stats::loess(xs ~ pt, weights = wt, span = span, degree = 2)$fitted
    fit_y  <- stats::loess(ys ~ pt, weights = wt, span = span, degree = 2)$fitted
    data.frame(
      Axis_1      = fit_x,
      Axis_2      = fit_y,
      raw_Axis_1  = xs,
      raw_Axis_2  = ys,
      Lineages    = factor(l, levels = lineages),
      stringsAsFactors = FALSE
    )
  })

  # Build layer list per lineage
  curve_layers <- lapply(lineages, function(l) {
    dat_sm <- fitted_list[[l]]
    if (is.null(dat_sm)) return(list())
    dat_sm <- unique(stats::na.omit(dat_sm))
    layers <- list()

    if (isTRUE(whiskers)) {
      layers <- c(layers, list(
        ggplot2::geom_segment(
          data    = dat_sm,
          mapping = ggplot2::aes(
            x = .data[["Axis_1"]], y = .data[["Axis_2"]],
            xend = .data[["raw_Axis_1"]], yend = .data[["raw_Axis_2"]],
            color = .data[["Lineages"]]
          ),
          linewidth   = whiskers_linewidth,
          alpha       = whiskers_alpha,
          show.legend = FALSE,
          inherit.aes = FALSE
        )
      ))
    }

    layers <- c(layers, list(
      # White background stroke (outline effect)
      ggplot2::geom_path(
        data    = dat_sm,
        mapping = ggplot2::aes(x = .data[["Axis_1"]], y = .data[["Axis_2"]]),
        color       = line_bg,
        linewidth   = linewidth + line_bg_stroke,
        arrow       = lineages_arrow,
        show.legend = FALSE,
        inherit.aes = FALSE
      ),
      # Coloured foreground path
      ggplot2::geom_path(
        data    = dat_sm,
        mapping = ggplot2::aes(
          x = .data[["Axis_1"]], y = .data[["Axis_2"]],
          color = .data[["Lineages"]]
        ),
        linewidth   = linewidth,
        arrow       = lineages_arrow,
        show.legend = TRUE,
        inherit.aes = FALSE
      )
    ))
    layers
  })

  c(
    unlist(curve_layers, recursive = FALSE),
    list(
      ggplot2::scale_color_manual(
        name   = "Lineage:",
        values = colors,
        guide  = ggplot2::guide_legend(
          title.hjust = 0, order = 99,
          override.aes = list(linewidth = 1.5, alpha = 1)
        )
      )
    )
  )
}


# --------------------------------------------------------------------------
# Internal: build a marginal distribution ggplot
# --------------------------------------------------------------------------

#' @keywords internal
.build_marginal <- function(plot_df, var, type, group.by, cols,
                            point.color, bins = 30) {
  # Base aes
  if (!is.null(group.by)) {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[[var]],
                   fill = .data[[group.by]],
                   color = .data[[group.by]])
    )
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[var]]))
  }

  # Geom layers by type
  if (type == "density") {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_density(alpha = 0.35, linewidth = 0.4)
    } else {
      p <- p + ggplot2::geom_density(
        fill = point.color, color = point.color, alpha = 0.35, linewidth = 0.4
      )
    }
  } else if (type == "histogram") {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_histogram(
        alpha = 0.5, position = "identity", bins = bins, linewidth = 0.2
      )
    } else {
      p <- p + ggplot2::geom_histogram(
        fill = point.color, color = "white", alpha = 0.7,
        bins = bins, linewidth = 0.2
      )
    }
  } else if (type == "boxplot") {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_boxplot(
        ggplot2::aes(y = .data[[group.by]]),
        alpha = 0.4, outlier.size = 0.3, linewidth = 0.3
      )
    } else {
      p <- p + ggplot2::geom_boxplot(
        ggplot2::aes(y = 0),
        fill = point.color, alpha = 0.4,
        outlier.size = 0.3, linewidth = 0.3
      )
    }
  } else if (type == "violin") {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_violin(
        ggplot2::aes(y = .data[[group.by]]),
        alpha = 0.35, linewidth = 0.3, scale = "width"
      )
    } else {
      p <- p + ggplot2::geom_violin(
        ggplot2::aes(y = 0),
        fill = point.color, alpha = 0.35,
        linewidth = 0.3, scale = "width"
      )
    }
  } else if (type == "densigram") {
    # histogram + density overlay
    if (!is.null(group.by)) {
      p <- p +
        ggplot2::geom_histogram(
          ggplot2::aes(y = ggplot2::after_stat(density)),
          alpha = 0.35, position = "identity",
          bins = bins, linewidth = 0.2
        ) +
        ggplot2::geom_density(alpha = 0.2, linewidth = 0.5)
    } else {
      p <- p +
        ggplot2::geom_histogram(
          ggplot2::aes(y = ggplot2::after_stat(density)),
          fill = point.color, color = "white", alpha = 0.4,
          bins = bins, linewidth = 0.2
        ) +
        ggplot2::geom_density(
          fill = point.color, color = point.color,
          alpha = 0.2, linewidth = 0.5
        )
    }
  }

  # Color scales
  if (!is.null(group.by)) {
    p <- p +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_color_manual(values = cols)
  }

  # Minimal theme: no axis labels, no legend, no background
  p <- p +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin     = ggplot2::margin(0, 0, 0, 0)
    )

  p
}
