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
