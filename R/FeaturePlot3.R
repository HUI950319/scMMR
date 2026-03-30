# ============================================================================
# FeaturePlot3 — Statistical expression plots (violin / box / bar / dot / col)
# ============================================================================

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Resolve a y-axis limit from numeric, quantile string, or NULL
#' @noRd
.fp3_resolve_ylim <- function(val, values, side = c("max", "min")) {
  side <- match.arg(side)
  if (is.null(val)) {
    if (side == "max") max(values, na.rm = TRUE) else min(values, na.rm = TRUE)
  } else if (is.character(val)) {
    q <- as.numeric(sub("^q(\\d+)$", "\\1", val)) / 100
    stats::quantile(values, q, na.rm = TRUE)
  } else {
    val
  }
}

#' Build background rect layer (bg.by bands)
#'
#' @param dat data.frame with group.by and optionally bg.by columns.
#' @param g The active group.by column name.
#' @param bg_col Named colour vector (name = bg.by level).
#' @param bg_map Named char vector mapping group.by → bg.by.
#' @param bg_alpha Transparency.
#' @param col_type "box" (uses x as discrete index) or "col" (uses numeric
#'   cell index).
#' @noRd
.fp3_bg_layer <- function(dat, g, bg_col, bg_map, bg_alpha,
                           col_type = FALSE) {
  if (isTRUE(col_type)) {
    # col type: x-axis is numeric cell index
    x_index <- split(dat[["cell"]], dat[["group.by"]])
    bg_data <- as.data.frame(t(sapply(x_index, range)))
    colnames(bg_data) <- c("xmin", "xmax")
    bg_data[["group.by"]] <- names(x_index)
    bg_data[["xmin"]] <- ifelse(
      bg_data[["xmin"]] == min(bg_data[["xmin"]]), -Inf,
      bg_data[["xmin"]] - 0.5
    )
    bg_data[["xmax"]] <- ifelse(
      bg_data[["xmax"]] == max(bg_data[["xmax"]]), Inf,
      bg_data[["xmax"]] + 0.5
    )
    bg_data[["ymin"]] <- -Inf
    bg_data[["ymax"]] <- Inf
    bg_data[["fill"]] <- bg_col[bg_map[as.character(bg_data[["group.by"]])]]
  } else {
    bg_data <- unique(dat[, "group.by", drop = FALSE])
    bg_data[["x"]] <- as.numeric(bg_data[["group.by"]])
    bg_data[["xmin"]] <- ifelse(
      bg_data[["x"]] == min(bg_data[["x"]]), -Inf,
      bg_data[["x"]] - 0.5
    )
    bg_data[["xmax"]] <- ifelse(
      bg_data[["x"]] == max(bg_data[["x"]]), Inf,
      bg_data[["x"]] + 0.5
    )
    bg_data[["ymin"]] <- -Inf
    bg_data[["ymax"]] <- Inf
    bg_data[["fill"]] <- bg_col[bg_map[as.character(bg_data[["group.by"]])]]
  }
  ggplot2::geom_rect(
    data    = bg_data,
    mapping = ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill    = bg_data[["fill"]],
    alpha   = bg_alpha,
    inherit.aes = FALSE
  )
}

#' Build a single expression stat panel
#' @noRd
.fp3_build_panel <- function(
    dat,             # data.frame: group.by, split.by, group.unique, value,
    #                   fill.by, cell (for col), features
    f,               # feature name (character)
    g,               # active group.by column name
    type,            # plot type
    fill.by,         # "group" / "feature" / "expression"
    colors,          # named fill colour vector
    colors_limits,   # range for "expression" fill (or NULL)
    keynm,           # legend key title
    bg_layer,        # geom_rect layer (or NULL)
    alpha, add_box, box_color, box_width, box_ptsize,
    add_point, pt.color, pt.size, pt.alpha, jitter.width, jitter.height,
    cells.highlight, cols.highlight, sizes.highlight, alpha.highlight,
    add_trend, trend_color, trend_linewidth, trend_ptsize,
    add_stat, stat_color, stat_size, stat_stroke, stat_shape,
    add_line, line_color, line_size, line_type,
    comparisons, ref_group, pairwise_method,
    multiplegroup_comparisons, multiple_method,
    sig_label, sig_labelsize,
    y_min_use, y_max_use, y.trans, y.nbreaks,
    flip, stack, keep_empty,
    title, subtitle, xlab, ylab,
    legend.position, legend.direction, legend.title,
    theme_fn, theme_args, aspect.ratio,
    levels_order,
    seed) {

  # ---- Base ggplot -----------------------------------------------------------
  if (type == "col") {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x    = .data[["cell"]],
        y    = .data[["value"]],
        fill = .data[["fill.by"]]
      )
    )
  } else {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x    = .data[["group.by"]],
        y    = .data[["value"]],
        fill = .data[["fill.by"]]
      )
    )
  }

  # ---- Background bands ------------------------------------------------------
  if (!is.null(bg_layer)) p <- p + bg_layer

  # ---- Reference line --------------------------------------------------------
  if (type %in% c("bar", "col")) {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = 2)
  }

  # ---- Main geom -------------------------------------------------------------
  if (type == "violin") {
    p <- p + ggplot2::geom_violin(
      scale    = "width",
      trim     = TRUE,
      alpha    = alpha,
      position = ggplot2::position_dodge()
    )
  } else if (type == "box") {
    p <- p +
      ggplot2::geom_boxplot(
        mapping  = ggplot2::aes(group = .data[["group.unique"]]),
        position = ggplot2::position_dodge(width = 0.9),
        color    = "black",
        width    = 0.8,
        outlier.shape = NA
      ) +
      ggplot2::stat_summary(
        fun     = stats::median, geom = "point",
        mapping = ggplot2::aes(group = .data[["split.by"]]),
        position = ggplot2::position_dodge(width = 0.9),
        color = "black", fill = "white", size = 1.5, shape = 21
      )
  } else if (type == "bar") {
    p <- p +
      ggplot2::stat_summary(
        fun     = mean, geom = "col",
        mapping = ggplot2::aes(group = .data[["split.by"]]),
        position = ggplot2::position_dodge(width = 0.9),
        width = 0.8, color = "black"
      ) +
      ggplot2::stat_summary(
        fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1),
        geom = "errorbar",
        mapping = ggplot2::aes(group = .data[["split.by"]]),
        position = ggplot2::position_dodge(width = 0.9),
        width = 0.2, color = "black"
      )
    y_min_use <- min(0, y_min_use)
  } else if (type == "dot") {
    brks <- seq(min(dat[["value"]], na.rm = TRUE),
                max(dat[["value"]], na.rm = TRUE), length.out = 15)
    bins <- cut(dat[["value"]], breaks = brks, include.lowest = TRUE)
    bins_mid <- sapply(strsplit(levels(bins), ","), function(x)
      stats::median(as.numeric(gsub("\\(|\\)|\\[|\\]", "", x)), na.rm = TRUE)
    )
    names(bins_mid) <- levels(bins)
    dat[["bins"]] <- bins_mid[as.character(bins)]
    p <- p +
      ggplot2::geom_count(
        data    = dat,
        ggplot2::aes(y = .data[["bins"]]),
        shape   = 21,
        alpha   = alpha,
        position = ggplot2::position_dodge(width = 0.9)
      ) +
      ggplot2::scale_size_area(name = "Count", max_size = 6, n.breaks = 4) +
      ggplot2::guides(size = ggplot2::guide_legend(
        override.aes = list(fill = "grey30", shape = 21), order = 2
      ))
  } else if (type == "col") {
    p <- p + ggplot2::geom_col()
    if (isTRUE(flip)) {
      p <- p + ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    } else {
      p <- p + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    }
    # Group border lines
    if (nlevels(dat[["group.by"]]) > 1) {
      x_index <- split(dat[["cell"]], dat[["group.by"]])
      border_x <- sapply(x_index, min)[-1] - 0.5
      p <- p + ggplot2::geom_vline(
        xintercept = border_x, linetype = 2, alpha = 0.5
      )
    }
  }

  # ---- Comparisons -----------------------------------------------------------
  if (length(comparisons) > 0 || isTRUE(comparisons)) {
    if (!requireNamespace("ggpubr", quietly = TRUE))
      stop("Package 'ggpubr' is required for comparisons.", call. = FALSE)

    if (isTRUE(comparisons)) {
      # Auto: use split.by groups
      row_sums <- Matrix::rowSums(
        table(dat[["group.by"]], dat[["split.by"]]) >= 2
      )
      group_use <- names(which(row_sums >= 2))
      method <- if (any(row_sums >= 3)) multiple_method else pairwise_method
      p <- p + ggpubr::stat_compare_means(
        data    = dat[dat[["group.by"]] %in% group_use, , drop = FALSE],
        mapping = ggplot2::aes(
          x = .data[["group.by"]], y = .data[["value"]],
          group = .data[["group.unique"]]
        ),
        label = sig_label, label.y = y_max_use,
        size = sig_labelsize, step.increase = 0.1,
        tip.length = 0.03, vjust = 1, method = method
      )
      y_max_use <- max(y_max_use * 1.15, y_max_use + (y_max_use - y_min_use) * 0.15)
    } else {
      # Manual comparisons
      split_active <- nlevels(dat[["split.by"]]) > 1
      if (split_active) {
        valid <- Filter(function(cp) length(cp) == 2 &&
                          all(cp %in% levels(dat[["split.by"]])), comparisons)
        if (length(valid) > 0) {
          p <- p + ggpubr::stat_compare_means(
            mapping = ggplot2::aes(
              x = .data[["group.by"]], y = .data[["value"]],
              group = .data[["split.by"]]
            ),
            comparisons = valid, ref.group = ref_group,
            method = pairwise_method, label = sig_label,
            label.y = y_max_use, size = sig_labelsize,
            step.increase = 0.1, tip.length = 0.03, vjust = 0
          )
        }
      } else {
        p <- p + ggpubr::stat_compare_means(
          mapping = ggplot2::aes(
            x = .data[["group.by"]], y = .data[["value"]],
            group = .data[["group.unique"]]
          ),
          comparisons = comparisons, ref.group = ref_group,
          method = pairwise_method, label = sig_label,
          label.y = y_max_use, size = sig_labelsize,
          step.increase = 0.1, tip.length = 0.03, vjust = 0
        )
      }
      y_max_use <- y_max_use + (y_max_use - y_min_use) * 0.15 * length(comparisons)
    }
  }

  if (isTRUE(multiplegroup_comparisons)) {
    if (!requireNamespace("ggpubr", quietly = TRUE))
      stop("Package 'ggpubr' is required.", call. = FALSE)
    p <- p + ggpubr::stat_compare_means(
      ggplot2::aes(
        x = .data[["group.by"]], y = .data[["value"]],
        group = .data[["group.unique"]]
      ),
      method = multiple_method, label = sig_label,
      label.y = y_max_use, size = sig_labelsize,
      vjust = 1.2, hjust = 0
    )
    y_max_use <- y_max_use + (y_max_use - y_min_use) * 0.1
  }

  # ---- Overlays --------------------------------------------------------------
  if (isTRUE(add_point) && type != "col") {
    set.seed(seed)
    p <- suppressWarnings(p +
      ggplot2::geom_point(
        ggplot2::aes(
          x = .data[["group.by"]], y = .data[["value"]],
          group = .data[["group.unique"]]
        ),
        inherit.aes = FALSE,
        color   = pt.color, size = pt.size, alpha = pt.alpha,
        position = ggplot2::position_jitterdodge(
          jitter.width = jitter.width, jitter.height = jitter.height,
          dodge.width = 0.9, seed = seed
        ),
        show.legend = FALSE
      )
    )
    # Highlighted cells
    if (!is.null(cells.highlight)) {
      hcells <- if (isTRUE(cells.highlight)) rownames(dat) else
        intersect(cells.highlight, rownames(dat))
      if (length(hcells) > 0) {
        hdf <- dat[hcells, , drop = FALSE]
        p <- p +
          ggplot2::geom_point(
            data = hdf,
            ggplot2::aes(
              x = .data[["group.by"]], y = .data[["value"]],
              group = .data[["group.unique"]]
            ),
            inherit.aes = FALSE,
            color   = cols.highlight, size = sizes.highlight,
            alpha   = alpha.highlight,
            position = ggplot2::position_jitterdodge(
              jitter.width = jitter.width, jitter.height = jitter.height,
              dodge.width = 0.9, seed = seed
            ),
            show.legend = FALSE
          )
      }
    }
  }

  if (isTRUE(add_box) && type %in% c("violin", "bar")) {
    p <- p +
      ggplot2::geom_boxplot(
        ggplot2::aes(group = .data[["group.unique"]]),
        position = ggplot2::position_dodge(width = 0.9),
        color    = box_color, fill = box_color,
        width    = box_width, outlier.shape = NA, show.legend = FALSE
      ) +
      ggplot2::stat_summary(
        fun     = stats::median, geom = "point",
        mapping = ggplot2::aes(group = .data[["split.by"]]),
        position = ggplot2::position_dodge(width = 0.9),
        color = "black", fill = "white", size = box_ptsize, shape = 21
      )
  }

  if (isTRUE(add_trend) && type %in% c("violin", "box", "bar")) {
    trend_fun <- if (type == "bar") mean else stats::median
    multi_split <- nlevels(dat[["split.by"]]) > 1
    if (multi_split) {
      # Per-group trend across split levels: extract dodge positions from a
      # temporary layer, then draw geom_line at those positions
      pt_layer <- ggplot2::stat_summary(
        fun = trend_fun, geom = "point",
        mapping = ggplot2::aes(
          group = .data[["split.by"]],
          color = .data[["group.by"]]
        ),
        position = ggplot2::position_dodge(width = 0.9),
        fill = "white", size = trend_ptsize, shape = 21
      )
      p_tmp <- p + pt_layer
      ld <- ggplot2::layer_data(p_tmp, length(p_tmp$layers))
      p <- p +
        ggplot2::geom_line(
          data = ld,
          ggplot2::aes(x = x, y = y, group = colour),
          color = trend_color, linewidth = trend_linewidth,
          inherit.aes = FALSE
        ) +
        ggplot2::stat_summary(
          fun     = trend_fun, geom = "point",
          mapping = ggplot2::aes(group = .data[["split.by"]]),
          position = ggplot2::position_dodge(width = 0.9),
          color = "black", fill = "white", size = trend_ptsize, shape = 21
        )
    } else {
      p <- p +
        ggplot2::stat_summary(
          fun = trend_fun, geom = "line",
          mapping = ggplot2::aes(group = .data[["split.by"]]),
          position = ggplot2::position_dodge(width = 0.9),
          color = trend_color, linewidth = trend_linewidth
        ) +
        ggplot2::stat_summary(
          fun = trend_fun, geom = "point",
          mapping = ggplot2::aes(group = .data[["split.by"]]),
          position = ggplot2::position_dodge(width = 0.9),
          color = "black", fill = "white", size = trend_ptsize, shape = 21
        )
    }
  }

  if (add_stat != "none" && type != "col") {
    p <- p +
      ggplot2::stat_summary(
        fun = add_stat, geom = "point",
        mapping = ggplot2::aes(
          group = .data[["split.by"]],
          shape = stat_shape
        ),
        position = ggplot2::position_dodge(width = 0.9),
        color  = stat_color, fill = stat_color,
        size   = stat_size, stroke = stat_stroke
      ) +
      ggplot2::scale_shape_identity()
  }

  if (!is.null(add_line)) {
    p <- p + ggplot2::geom_hline(
      yintercept = add_line,
      color      = line_color,
      linetype   = line_type,
      linewidth  = line_size
    )
  }

  # ---- Facet (features strip) ------------------------------------------------
  if (nrow(dat) == 0) {
    p <- p + ggplot2::facet_null()
  } else {
    if (isTRUE(stack) && !isTRUE(flip)) {
      p <- p +
        ggplot2::facet_grid(features ~ .) +
        ggplot2::theme(
          strip.text.y = ggplot2::element_text(angle = 0)
        )
    } else {
      p <- p + ggplot2::facet_grid(. ~ features)
    }
  }

  # ---- Labels ----------------------------------------------------------------
  p <- p + ggplot2::labs(
    title = title, subtitle = subtitle, x = xlab, y = ylab
  )

  # ---- x scale ---------------------------------------------------------------
  if (nrow(dat) != 0 && type != "col") {
    p <- p + ggplot2::scale_x_discrete(drop = !keep_empty)
  }

  # ---- Coordinates -----------------------------------------------------------
  if (isTRUE(flip)) {
    p <- p + ggplot2::coord_flip(ylim = c(y_min_use, y_max_use))
  } else {
    p <- p + ggplot2::coord_cartesian(ylim = c(y_min_use, y_max_use))
  }

  # ---- y scale ---------------------------------------------------------------
  if (isTRUE(stack)) {
    p <- p + ggplot2::scale_y_continuous(
      trans  = y.trans,
      breaks = c(y_min_use, y_max_use),
      labels = c(round(y_min_use, 1), round(y_max_use, 1))
    )
  } else {
    p <- p + ggplot2::scale_y_continuous(
      trans   = y.trans,
      n.breaks = y.nbreaks
    )
  }

  # ---- Fill scale ------------------------------------------------------------
  legend_title_use <- legend.title %||% paste0(keynm, ":")
  if (fill.by != "expression") {
    scale_args <- list(
      name   = legend_title_use,
      values = colors,
      breaks = levels_order,
      drop   = FALSE
    )
    if (isTRUE(stack)) scale_args[["limits"]] <- levels_order
    p <- p +
      do.call(ggplot2::scale_fill_manual, scale_args) +
      do.call(ggplot2::scale_color_manual, scale_args) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          title.hjust = 0, order = 1,
          override.aes = list(size = 4, color = "black", alpha = 1)
        )
      )
  } else {
    p <- p +
      ggplot2::scale_fill_gradientn(
        name    = legend_title_use,
        colours = colors,
        limits  = colors_limits
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_colorbar(
          frame.colour = "black", ticks.colour = "black",
          title.hjust = 0, order = 1
        )
      )
  }

  # ---- Theme -----------------------------------------------------------------
  theme_obj <- do.call(theme_fn, theme_args)
  p <- p + theme_obj +
    ggplot2::theme(
      aspect.ratio = aspect.ratio,
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      strip.text.y = ggplot2::element_text(angle = 0),
      panel.grid.major.y = ggplot2::element_line(
        color = "grey", linetype = 2, linewidth = 0.3
      ),
      legend.position  = legend.position,
      legend.direction = legend.direction
    )

  if (isTRUE(flip)) {
    extra <- if (isTRUE(stack)) {
      list(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        strip.text.x = ggplot2::element_text(angle = 90),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(
          color = "grey", linetype = 2, linewidth = 0.3
        )
      )
    } else {
      list(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(
          color = "grey", linetype = 2, linewidth = 0.3
        )
      )
    }
    p <- p + do.call(ggplot2::theme, extra)
  }

  p
}


# ---------------------------------------------------------------------------
# Stack combining with patchwork
# ---------------------------------------------------------------------------

#' Assemble stacked panels into a single patchwork figure
#' @noRd
.fp3_stack_combine <- function(plist, flip, legend.position, ylab) {
  n <- length(plist)
  if (n == 1) return(plist[[1]])

  # Strip redundant axes from all but the anchor panel
  if (isTRUE(flip)) {
    # Side-by-side (columns), keep left panel's y-axis labels intact
    plist_mod <- lapply(seq_along(plist), function(i) {
      p <- plist[[i]]
      if (i != 1) {
        p <- p + ggplot2::theme(
          legend.position  = "none",
          plot.title       = ggplot2::element_blank(),
          plot.subtitle    = ggplot2::element_blank(),
          axis.title.x     = ggplot2::element_blank(),
          axis.title.y     = ggplot2::element_blank(),
          axis.text.y      = ggplot2::element_blank(),
          axis.ticks.y     = ggplot2::element_blank(),
          plot.margin      = ggplot2::margin(0, 0, 0, -2)
        )
      } else {
        p <- p + ggplot2::theme(
          legend.position = "none",
          axis.title.x    = ggplot2::element_blank()
        )
      }
      p
    })
    combined <- patchwork::wrap_plots(plist_mod, nrow = 1) +
      patchwork::plot_annotation(caption = ylab) &
      ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0.5))
  } else {
    # Stacked (rows), keep bottom panel's x-axis labels intact
    plist_mod <- lapply(seq_along(plist), function(i) {
      p <- plist[[i]]
      if (i != n) {
        p <- p + ggplot2::theme(
          legend.position  = "none",
          plot.title       = ggplot2::element_blank(),
          plot.subtitle    = ggplot2::element_blank(),
          axis.title.x     = ggplot2::element_blank(),
          axis.text.x      = ggplot2::element_blank(),
          axis.ticks.x     = ggplot2::element_blank(),
          plot.margin      = ggplot2::margin(-2, 0, 0, 0)
        )
      } else {
        p <- p + ggplot2::theme(
          legend.position = "none",
          axis.title.y    = ggplot2::element_blank()
        )
      }
      p
    })
    combined <- patchwork::wrap_plots(plist_mod, ncol = 1) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = legend.position)
  }
  combined
}


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Expression distribution plots for Seurat objects
#'
#' Generates violin, box, bar, dot, or column plots for one or more features
#' (genes or metadata numeric columns) stratified by a grouping variable.
#' Parameter names are consistent with \code{DimPlot2} / \code{FeaturePlot2}.
#'
#' @param seu A Seurat object or data.frame. When a data.frame is provided,
#'   \code{assay} and \code{layer} are ignored.
#' @param features Character vector of feature names (genes or numeric metadata
#'   columns).
#' @param group.by Column name for the x-axis grouping variable. \code{NULL}
#'   treats all cells as a single group.
#' @param split.by Column name for a secondary grouping (dodging within each
#'   x-axis position). Default \code{NULL}.
#' @param bg.by Column name whose levels define background colour bands.
#'   Must be a super-grouping of \code{group.by} (each \code{group.by} level
#'   belongs to exactly one \code{bg.by} level). Default \code{NULL}.
#' @param plot.by Axis orientation: \code{"group"} (default) places groups on
#'   the x-axis (one panel per feature); \code{"feature"} places features on
#'   the x-axis (one panel per group level).
#' @param fill.by What drives fill colour:
#'   \itemize{
#'     \item \code{"group"} (default) — group or split levels.
#'     \item \code{"feature"} — feature name.
#'     \item \code{"expression"} — median expression per group (continuous
#'       gradient).
#'   }
#' @param cells Character vector of cell barcodes to include. \code{NULL} = all.
#' @param assay Seurat assay name. \code{NULL} = active assay.
#' @param layer Seurat layer. Default \code{"data"}.
#' @param keep_empty Logical. Keep empty factor levels on the x-axis. Default
#'   \code{FALSE}.
#' @param individual Logical. Create a separate panel for each group level ×
#'   feature × split level combination. Default \code{FALSE}.
#' @param type Plot geometry. One of \code{"violin"} (default), \code{"box"},
#'   \code{"bar"}, \code{"dot"}, or \code{"col"} (per-cell column bar).
#' @param palette Colour palette name for fill colours. Default \code{"npg"}.
#' @param palcolor Custom colour vector. Overrides \code{palette}. Default
#'   \code{NULL}.
#' @param alpha Fill transparency (0–1). Default \code{1}.
#' @param bg_palette Colour palette for background bands. Default \code{"Set2"}.
#' @param bg_palcolor Custom colours for background bands. Default \code{NULL}.
#' @param bg_alpha Transparency of background bands. Default \code{0.2}.
#' @param add_box Logical. Overlay a narrow box-plot (violin/bar only). Default
#'   \code{FALSE}.
#' @param box_color Colour of overlay box-plot. Default \code{"black"}.
#' @param box_width Width of overlay box-plot. Default \code{0.1}.
#' @param box_ptsize Median point size in box-plot. Default \code{2}.
#' @param add_point Logical. Overlay jittered individual points. Default
#'   \code{FALSE}.
#' @param pt.color Point colour. Default \code{"grey30"}.
#' @param pt.size Point size. \code{NULL} = auto.
#' @param pt.alpha Point transparency. Default \code{1}.
#' @param jitter.width Horizontal jitter. Default \code{0.4}.
#' @param jitter.height Vertical jitter. Default \code{0.1}.
#' @param add_trend Logical. Connect group medians/means with a trend line
#'   (violin/box: median; bar: mean). Default \code{FALSE}.
#' @param trend_color Trend line colour. Default \code{"black"}.
#' @param trend_linewidth Trend line width. Default \code{1}.
#' @param trend_ptsize Trend point size. Default \code{2}.
#' @param add_stat Summary stat to overlay. One of \code{"none"} (default),
#'   \code{"mean"}, or \code{"median"}.
#' @param stat_color Stat point colour. Default \code{"black"}.
#' @param stat_size Stat point size. Default \code{1}.
#' @param stat_stroke Stat point stroke. Default \code{1}.
#' @param stat_shape Stat point shape (integer). Default \code{25} (filled
#'   downward triangle).
#' @param add_line Numeric y-intercept for a horizontal reference line.
#'   \code{NULL} = none.
#' @param line_color Reference line colour. Default \code{"red"}.
#' @param line_size Reference line width. Default \code{1}.
#' @param line_type Reference line type. Default \code{1} (solid).
#' @param cells.highlight Character vector of cell barcodes to highlight, or
#'   \code{TRUE} for all cells. Requires \code{add_point = TRUE}. Default
#'   \code{NULL}.
#' @param cols.highlight Highlight point colour. Default \code{"red"}.
#' @param sizes.highlight Highlight point size. Default \code{1}.
#' @param alpha.highlight Highlight point transparency. Default \code{1}.
#' @param calculate_coexp Logical. Compute geometric mean co-expression of all
#'   \code{features} and add it as an extra \code{"CoExp"} panel. Only valid
#'   when all features are genes (in the expression matrix). Default
#'   \code{FALSE}.
#' @param same.y.lims Logical. Share y-axis limits across all panels. Default
#'   \code{FALSE}.
#' @param y.min Minimum y limit: numeric value or quantile string (e.g.
#'   \code{"q1"}). \code{NULL} = data minimum.
#' @param y.max Maximum y limit: numeric value or quantile string (e.g.
#'   \code{"q99"}). \code{NULL} = data maximum.
#' @param y.trans Y-axis transformation. \code{"identity"} (default) or
#'   \code{"log2"}.
#' @param y.nbreaks Number of y-axis tick breaks. Default \code{5}.
#' @param sort Logical or \code{"increasing"} / \code{"decreasing"}. Sort
#'   groups by median of the first feature. Default \code{FALSE}.
#' @param stack Logical. Stack all feature panels vertically (or horizontally
#'   when \code{flip = TRUE}) with shared axes and a single shared legend.
#'   Default \code{FALSE}.
#' @param flip Logical. Flip coordinates (horizontal layout). Default
#'   \code{FALSE}.
#' @param comparisons List of length-2 character vectors specifying pairwise
#'   comparisons, or \code{TRUE} to auto-detect groups from \code{split.by}.
#'   \code{NULL} = none.
#' @param ref_group Reference group for one-vs-all comparisons. Default
#'   \code{NULL}.
#' @param pairwise_method Pairwise test. Default \code{"wilcox.test"}.
#' @param multiplegroup_comparisons Logical. Add an overall group test label.
#'   Default \code{FALSE}.
#' @param multiple_method Overall group test. Default \code{"kruskal.test"}.
#' @param sig_label Significance label format: \code{"p.signif"} (default) or
#'   \code{"p.format"}.
#' @param sig_labelsize Text size of significance labels. Default \code{3.5}.
#' @param aspect.ratio Panel aspect ratio. \code{NULL} = free.
#' @param title Plot title. \code{NULL} = none.
#' @param subtitle Plot subtitle. \code{NULL} = none.
#' @param xlab X-axis label. \code{NULL} = auto (\code{group.by}).
#' @param ylab Y-axis label. Default \code{"Expression level"}.
#' @param legend.position Legend position string or coordinates. Default
#'   \code{"right"}.
#' @param legend.direction Legend direction. Default \code{"vertical"}.
#' @param legend.title Custom legend title. \code{NULL} = auto.
#' @param theme_use Theme specification. \code{NULL} = \code{ggplot2::theme_bw};
#'   accepts a character function name, a function, or a \code{theme} object.
#' @param theme_args List of extra arguments passed to \code{theme_use}. Default
#'   \code{list()}.
#' @param combine Logical. Combine panels into a single patchwork. Default
#'   \code{TRUE}.
#' @param nrow Number of rows in patchwork layout.
#' @param ncol Number of columns in patchwork layout.
#' @param byrow Logical. Fill patchwork by row. Default \code{TRUE}.
#' @param force Logical. Skip the >100-level cardinality check. Default
#'   \code{FALSE}.
#' @param seed Random seed for jitter. Default \code{11}.
#'
#' @return A ggplot2 / patchwork object (or list when \code{combine = FALSE}).
#'
#' @examples
#' \dontrun{
#' # Basic violin plot
#' FeaturePlot3(seu, features = c("CD3D", "CD8A"), group.by = "cell_type")
#'
#' # Box plot with statistical comparisons
#' FeaturePlot3(seu, features = "MS4A1", group.by = "cluster",
#'              type = "box", comparisons = list(c("B", "T")))
#'
#' # Features on x-axis (one panel per cell type)
#' FeaturePlot3(seu, features = c("CD3D", "CD8A", "MS4A1"),
#'              group.by = "cell_type", plot.by = "feature")
#'
#' # Stacked panels with background bands
#' FeaturePlot3(seu,
#'   features = c("Gene1", "Gene2", "Gene3"),
#'   group.by = "SubCellType", bg.by = "CellType",
#'   stack = TRUE, add_box = TRUE)
#'
#' # Col (per-cell bar) plot
#' FeaturePlot3(seu, features = "CD3D", group.by = "cell_type", type = "col")
#' }
#'
#' @export
FeaturePlot3 <- function(
    seu,
    features,
    group.by                  = NULL,
    split.by                  = NULL,
    bg.by                     = NULL,
    plot.by                   = c("group", "feature"),
    fill.by                   = c("group", "feature", "expression"),
    cells                     = NULL,
    assay                     = NULL,
    layer                     = "data",
    keep_empty                = FALSE,
    individual                = FALSE,
    type                      = c("violin", "box", "bar", "dot", "col"),
    palette                   = "npg",
    palcolor                  = NULL,
    alpha                     = 1,
    bg_palette                = "Set2",
    bg_palcolor               = NULL,
    bg_alpha                  = 0.2,
    add_box                   = FALSE,
    box_color                 = "black",
    box_width                 = 0.1,
    box_ptsize                = 2,
    add_point                 = FALSE,
    pt.color                  = "grey30",
    pt.size                   = NULL,
    pt.alpha                  = 1,
    jitter.width              = 0.4,
    jitter.height             = 0.1,
    add_trend                 = FALSE,
    trend_color               = "black",
    trend_linewidth           = 1,
    trend_ptsize              = 2,
    add_stat                  = c("none", "mean", "median"),
    stat_color                = "black",
    stat_size                 = 1,
    stat_stroke               = 1,
    stat_shape                = 25,
    add_line                  = NULL,
    line_color                = "red",
    line_size                 = 1,
    line_type                 = 1,
    cells.highlight           = NULL,
    cols.highlight            = "red",
    sizes.highlight           = 1,
    alpha.highlight           = 1,
    calculate_coexp           = FALSE,
    same.y.lims               = FALSE,
    y.min                     = NULL,
    y.max                     = NULL,
    y.trans                   = "identity",
    y.nbreaks                 = 5,
    sort                      = FALSE,
    stack                     = FALSE,
    flip                      = FALSE,
    comparisons               = NULL,
    ref_group                 = NULL,
    pairwise_method           = "wilcox.test",
    multiplegroup_comparisons = FALSE,
    multiple_method           = "kruskal.test",
    sig_label                 = c("p.signif", "p.format"),
    sig_labelsize             = 3.5,
    aspect.ratio              = NULL,
    title                     = NULL,
    subtitle                  = NULL,
    xlab                      = NULL,
    ylab                      = "Expression level",
    legend.position           = "right",
    legend.direction          = "vertical",
    legend.title              = NULL,
    theme_use                 = NULL,
    theme_args                = list(),
    combine                   = TRUE,
    nrow                      = NULL,
    ncol                      = NULL,
    byrow                     = TRUE,
    force                     = FALSE,
    seed                      = 11
) {
  set.seed(seed)

  # ---- Argument matching -----------------------------------------------------
  plot.by   <- match.arg(plot.by)
  fill.by   <- match.arg(fill.by)
  type      <- match.arg(type)
  sig_label <- match.arg(sig_label)
  add_stat  <- match.arg(add_stat)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- Input validation -------------------------------------------------------
  if (!is.null(add_line) && !is.numeric(add_line))
    stop("'add_line' must be numeric or NULL.", call. = FALSE)

  if (type == "col" && (isTRUE(add_box) || isTRUE(add_point) ||
                         isTRUE(add_trend) || add_stat != "none")) {
    message("Overlays (add_box/add_point/add_trend/add_stat) ignored for type = 'col'.")
    add_box <- add_point <- add_trend <- FALSE
    add_stat <- "none"
  }
  if (type == "col" && (isTRUE(multiplegroup_comparisons) ||
                         length(comparisons) > 0)) {
    message("Comparisons not supported for type = 'col'. Ignored.")
    multiplegroup_comparisons <- FALSE
    comparisons <- NULL
  }
  if (isTRUE(stack) && isTRUE(sort)) {
    message("sort forced to FALSE when stack = TRUE.")
    sort <- FALSE
  }

  # ---- Theme resolution ------------------------------------------------------
  if (is.null(theme_use)) {
    theme_fn  <- ggplot2::theme_bw
  } else if (is.character(theme_use)) {
    theme_fn <- tryCatch(match.fun(theme_use), error = function(e)
      tryCatch(utils::getFromNamespace(theme_use, "UtilsR"),
               error = function(e2) ggplot2::theme_bw)
    )
  } else if (is.function(theme_use)) {
    theme_fn <- theme_use
  } else if (inherits(theme_use, "theme")) {
    theme_fn <- function(...) theme_use
  } else {
    stop("'theme_use' must be NULL, character, function, or theme object.",
         call. = FALSE)
  }

  # ---- Extract data ----------------------------------------------------------
  if (is.null(group.by)) {
    group.by <- "All.groups"
    xlab     <- ""
  }

  dat_full <- .fp3_extract(seu, features = features, group.by = group.by,
                            split.by = split.by, bg.by = bg.by,
                            assay = assay, layer = layer, cells = cells)

  # ---- Co-expression ---------------------------------------------------------
  if (isTRUE(calculate_coexp) && length(features) > 1) {
    coexp_mat <- dat_full[, features, drop = FALSE]
    coexp_mat[!is.finite(as.matrix(coexp_mat))] <- NA
    dat_full[["CoExp"]] <- apply(as.matrix(coexp_mat), 1, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) return(NA_real_)
      exp(mean(log(x)))
    })
    features <- c(features, "CoExp")
    message("Co-expression added as 'CoExp'.")
  }

  # ---- High-cardinality check ------------------------------------------------
  all_group_cols <- unique(c(group.by, split.by))
  all_group_cols <- all_group_cols[all_group_cols %in% names(dat_full)]
  nlev <- sapply(dat_full[, all_group_cols, drop = FALSE], function(x)
    nlevels(as.factor(x)))
  nlev_hi <- nlev[nlev > 100]
  if (length(nlev_hi) > 0 && !isTRUE(force)) {
    warning("Columns with >100 levels: ",
            paste(names(nlev_hi), collapse = ", "),
            ". Use force = TRUE to continue.", call. = FALSE)
    return(invisible(NULL))
  }

  # ---- Ensure factors --------------------------------------------------------
  for (col in all_group_cols) {
    if (!is.factor(dat_full[[col]]))
      dat_full[[col]] <- factor(dat_full[[col]], levels = unique(dat_full[[col]]))
  }
  if (!is.null(bg.by) && bg.by %in% names(dat_full) && !is.factor(dat_full[[bg.by]]))
    dat_full[[bg.by]] <- factor(dat_full[[bg.by]], levels = unique(dat_full[[bg.by]]))

  # ---- Dummy split.by --------------------------------------------------------
  has_split <- !is.null(split.by) && split.by %in% names(dat_full)
  if (!has_split) {
    split.by <- ".split_dummy"
    dat_full[[split.by]] <- factor("")
    has_split <- FALSE
  }

  # ---- Global y limits -------------------------------------------------------
  global_ymin <- NULL
  global_ymax <- NULL
  if (isTRUE(same.y.lims)) {
    all_vals <- unlist(lapply(features[features %in% names(dat_full)], function(f)
      dat_full[[f]][is.finite(dat_full[[f]])]))
    global_ymin <- .fp3_resolve_ylim(y.min, all_vals, "min")
    global_ymax <- .fp3_resolve_ylim(y.max, all_vals, "max")
  }

  # ---- Auto pt.size ----------------------------------------------------------
  if (is.null(pt.size)) pt.size <- min(3000 / nrow(dat_full), 0.5)

  # ---- plot.by = "feature" ---------------------------------------------------
  if (plot.by == "feature") {
    return(.fp3_feature_axis(
      dat_full = dat_full, features = features,
      group.by = group.by, split.by = split.by,
      has_split = has_split, type = type,
      fill.by = fill.by, palette = palette, palcolor = palcolor, alpha = alpha,
      bg_palette = bg_palette, bg_palcolor = bg_palcolor, bg_alpha = bg_alpha,
      add_box = add_box, box_color = box_color, box_width = box_width,
      box_ptsize = box_ptsize,
      add_point = add_point, pt.color = pt.color, pt.size = pt.size,
      pt.alpha = pt.alpha, jitter.width = jitter.width, jitter.height = jitter.height,
      cells.highlight = cells.highlight, cols.highlight = cols.highlight,
      sizes.highlight = sizes.highlight, alpha.highlight = alpha.highlight,
      add_trend = add_trend, trend_color = trend_color,
      trend_linewidth = trend_linewidth, trend_ptsize = trend_ptsize,
      add_stat = add_stat, stat_color = stat_color, stat_size = stat_size,
      stat_stroke = stat_stroke, stat_shape = stat_shape,
      add_line = add_line, line_color = line_color,
      line_size = line_size, line_type = line_type,
      comparisons = comparisons, ref_group = ref_group,
      pairwise_method = pairwise_method,
      multiplegroup_comparisons = multiplegroup_comparisons,
      multiple_method = multiple_method,
      sig_label = sig_label, sig_labelsize = sig_labelsize,
      global_ymin = global_ymin, global_ymax = global_ymax,
      y.min = y.min, y.max = y.max, y.trans = y.trans, y.nbreaks = y.nbreaks,
      same.y.lims = same.y.lims,
      sort = sort, stack = stack, flip = flip, keep_empty = keep_empty,
      individual = individual,
      title = title, subtitle = subtitle, xlab = xlab, ylab = ylab,
      legend.position = legend.position, legend.direction = legend.direction,
      legend.title = legend.title,
      theme_fn = theme_fn, theme_args = theme_args, aspect.ratio = aspect.ratio,
      combine = combine, nrow = nrow, ncol = ncol, byrow = byrow,
      force = force, seed = seed
    ))
  }

  # ---- plot.by = "group" (main path) -----------------------------------------
  # Background colour map
  bg_map <- .fp3_bg_map(dat_full, group.by, bg.by)
  bg_col <- if (!is.null(bg.by) && bg.by %in% names(dat_full)) {
    palette_colors(levels(dat_full[[bg.by]]),
                   palette = bg_palette, palcolor = bg_palcolor)
  } else {
    grp_levels <- levels(dat_full[[group.by]])
    cols <- rep(c("transparent", "grey85"), length.out = length(grp_levels))
    stats::setNames(cols, grp_levels)
  }

  # ---- Colour resolution for fill.by = "expression" -------------------------
  median_expr <- NULL
  colors_limits_global <- NULL
  if (fill.by == "expression") {
    feat_cols <- features[features %in% names(dat_full)]
    agg_df <- stats::aggregate(
      dat_full[, feat_cols, drop = FALSE],
      by = list(dat_full[[group.by]], dat_full[[split.by]]),
      FUN = function(x) stats::median(x, na.rm = TRUE)
    )
    rownames(agg_df) <- paste0(agg_df[[1]], "-", agg_df[[2]])
    median_expr <- agg_df
    all_medians <- unlist(agg_df[, feat_cols, drop = FALSE])
    colors_limits_global <- range(all_medians, na.rm = TRUE)
  }

  # ---- Expand combinations ---------------------------------------------------
  comb <- expand.grid(
    group_name = group.by,
    stat_name  = features,
    stringsAsFactors = FALSE
  )

  # individual mode: per-group-level panel
  if (isTRUE(individual)) {
    comb <- merge(
      comb,
      expand.grid(
        group_name    = group.by,
        group_element = levels(dat_full[[group.by]]),
        split_name    = levels(dat_full[[split.by]]),
        stringsAsFactors = FALSE
      ),
      by = "group_name"
    )
  } else {
    comb <- merge(
      comb,
      data.frame(
        group_name    = group.by,
        group_element = I(list(levels(dat_full[[group.by]]))),
        split_name    = I(list(levels(dat_full[[split.by]]))),
        stringsAsFactors = FALSE
      ),
      by = "group_name"
    )
  }

  rownames(comb) <- paste0(
    comb[["stat_name"]], ":", comb[["group_name"]], ":",
    sapply(comb[["group_element"]], paste, collapse = ","), ":",
    sapply(comb[["split_name"]], paste, collapse = ",")
  )

  # ---- Build panels ----------------------------------------------------------
  plist <- stats::setNames(
    lapply(rownames(comb), function(i) {
      g  <- comb[i, "group_name"]
      f  <- comb[i, "stat_name"]
      grp_els <- comb[[i, "group_element"]]
      sp_els  <- comb[[i, "split_name"]]

      xlab_use <- xlab %||% g
      ylab_use <- ylab %||% "Expression level"

      # -- Colour resolution --
      colors <- .fp3_resolve_colors(
        fill.by = fill.by, f = f, g = g,
        dat = dat_full, split.by = split.by,
        features = features,
        palette = palette, palcolor = palcolor,
        median_expr = median_expr
      )
      colors_limits <- colors_limits_global

      # -- Subset data --
      if (!f %in% names(dat_full)) {
        message("Feature '", f, "' not found. Skipping.")
        return(NULL)
      }
      dat <- dat_full[
        dat_full[[g]] %in% grp_els & dat_full[[split.by]] %in% sp_els,
        unique(c(g, split.by, bg.by, f)),
        drop = FALSE
      ]
      dat[[g]] <- factor(dat[[g]], levels = levels(dat_full[[g]])[
        levels(dat_full[[g]]) %in% dat[[g]]
      ])

      # -- Sort groups --
      if (!isFALSE(sort) && nrow(dat) > 0) {
        df_sort <- stats::aggregate(
          dat[[f]], by = list(dat[[g]]), FUN = stats::median, na.rm = TRUE
        )
        decr <- !(is.character(sort) && sort == "increasing")
        sorted_lvl <- as.character(
          df_sort[order(df_sort[[2]], decreasing = decr), 1]
        )
        dat[[g]] <- factor(dat[[g]], levels = sorted_lvl)
      }

      # -- fill.by column --
      keynm <- .fp3_keynm(fill.by, f, g, split.by, has_split)
      if (fill.by == "feature") {
        dat[["fill.by"]] <- factor(f, levels = features)
        levels_order <- features
      } else if (fill.by == "group") {
        dat[["fill.by"]] <- if (has_split) dat[[split.by]] else dat[[g]]
        levels_order <- if (has_split) levels(dat_full[[split.by]]) else
          levels(dat[[g]])
      } else {
        # expression
        key_col <- paste0(dat[[g]], "-", dat[[split.by]])
        dat[["fill.by"]] <- median_expr[key_col, f]
        levels_order <- NULL
      }

      # -- group.unique (for dodging) --
      group_comb <- expand.grid(
        x = levels(dat[[split.by]]),
        y = levels(dat[[g]])
      )
      dat[["group.unique"]] <- factor(
        paste0("sp-", dat[[split.by]], "-gp-", dat[[g]]),
        levels = paste0("sp-", group_comb[[1]], "-gp-", group_comb[[2]])
      )
      dat[["split.by"]] <- dat[[split.by]]
      dat[["group.by"]] <- dat[[g]]
      dat[["features"]] <- f
      dat[["value"]] <- as.numeric(dat[[f]])
      dat <- dat[order(dat[["group.unique"]]), , drop = FALSE]

      # -- Violin: drop groups with < 2 obs --
      if (type == "violin") {
        grp_cnt <- table(dat[["group.unique"]])
        rm_grp  <- names(grp_cnt[grp_cnt < 2])
        if (length(rm_grp) > 0) {
          message("Dropped ", length(rm_grp), " group(s) with <2 obs in violin.")
          dat <- dat[!dat[["group.unique"]] %in% rm_grp, , drop = FALSE]
          dat[["group.unique"]] <- droplevels(dat[["group.unique"]])
        }
        if (nrow(dat) == 0) return(NULL)
      }

      # -- col type: assign cell index --
      if (type == "col") {
        if (isTRUE(flip)) {
          dat[["cell"]] <- rev(seq_len(nrow(dat)))
        } else {
          dat[["cell"]] <- seq_len(nrow(dat))
        }
      }

      # -- Finite check --
      fin_vals <- dat[["value"]][is.finite(dat[["value"]])]
      if (length(fin_vals) == 0) {
        message("No finite values for '", f, "'. Skipping.")
        return(NULL)
      }

      # -- flip: reverse group order --
      if (isTRUE(flip)) {
        dat[["group.by"]] <- factor(
          dat[["group.by"]],
          levels = rev(levels(dat[["group.by"]]))
        )
        ar <- if (!is.null(aspect.ratio)) 1 / aspect.ratio else NULL
      } else {
        ar <- aspect.ratio
      }

      # -- y limits --
      y_min_use <- global_ymin %||% .fp3_resolve_ylim(y.min, fin_vals, "min")
      y_max_use <- global_ymax %||% .fp3_resolve_ylim(y.max, fin_vals, "max")

      # -- Background layer --
      bg_layer_use <- .fp3_bg_layer(
        dat, g, bg_col, bg_map[[g]], bg_alpha,
        col_type = (type == "col")
      )

      # -- Build panel --
      .fp3_build_panel(
        dat           = dat, f = f, g = g,
        type          = type,
        fill.by       = fill.by,
        colors        = colors,
        colors_limits = colors_limits,
        keynm         = keynm,
        bg_layer      = bg_layer_use,
        alpha         = alpha,
        add_box       = add_box, box_color = box_color,
        box_width     = box_width, box_ptsize = box_ptsize,
        add_point     = add_point, pt.color = pt.color,
        pt.size       = pt.size, pt.alpha = pt.alpha,
        jitter.width  = jitter.width, jitter.height = jitter.height,
        cells.highlight = cells.highlight, cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight, alpha.highlight = alpha.highlight,
        add_trend     = add_trend, trend_color = trend_color,
        trend_linewidth = trend_linewidth, trend_ptsize = trend_ptsize,
        add_stat      = add_stat, stat_color = stat_color,
        stat_size     = stat_size, stat_stroke = stat_stroke,
        stat_shape    = stat_shape,
        add_line      = add_line, line_color = line_color,
        line_size     = line_size, line_type = line_type,
        comparisons   = comparisons, ref_group = ref_group,
        pairwise_method = pairwise_method,
        multiplegroup_comparisons = multiplegroup_comparisons,
        multiple_method = multiple_method,
        sig_label     = sig_label, sig_labelsize = sig_labelsize,
        y_min_use     = y_min_use, y_max_use = y_max_use,
        y.trans       = y.trans, y.nbreaks = y.nbreaks,
        flip          = flip, stack = stack, keep_empty = keep_empty,
        title         = title, subtitle = subtitle,
        xlab          = xlab_use, ylab = ylab_use,
        legend.position  = legend.position,
        legend.direction = legend.direction,
        legend.title  = legend.title,
        theme_fn      = theme_fn, theme_args = theme_args,
        aspect.ratio  = ar,
        levels_order  = levels_order,
        seed          = seed
      )
    }),
    rownames(comb)
  )
  plist <- Filter(Negate(is.null), plist)

  # ---- Stack layout ----------------------------------------------------------
  if (isTRUE(stack) && length(features) > 1 && !isTRUE(individual)) {
    plist_stacked <- lapply(group.by, function(g) {
      sub <- plist[sapply(strsplit(names(plist), ":"), function(x) x[2]) == g]
      if (length(sub) == 0) sub <- plist
      .fp3_stack_combine(sub, flip = flip,
                          legend.position = legend.position,
                          ylab = ylab %||% "Expression level")
    })
    names(plist_stacked) <- group.by
    plist <- plist_stacked
  }

  # ---- Combine ---------------------------------------------------------------
  if (isTRUE(combine)) {
    if (length(plist) == 1) return(plist[[1]])
    return(patchwork::wrap_plots(
      plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow
    ))
  }
  plist
}


# ============================================================================
# Supporting internals
# ============================================================================

#' Extract and validate feature data
#' @noRd
.fp3_extract <- function(seu, features, group.by, split.by, bg.by,
                          assay, layer, cells) {
  if (inherits(seu, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE))
      stop("SeuratObject is required for Seurat input.", call. = FALSE)
    if (!is.null(assay)) SeuratObject::DefaultAssay(seu) <- assay

    fetch_vars <- unique(c(features, group.by, split.by, bg.by))
    fetch_vars <- fetch_vars[!is.na(fetch_vars)]

    dat <- tryCatch(
      SeuratObject::FetchData(seu, vars = fetch_vars, layer = layer,
                              cells = cells),
      error = function(e) stop(conditionMessage(e), call. = FALSE)
    )
  } else if (is.data.frame(seu)) {
    dat <- seu
    if (!is.null(cells))
      dat <- dat[rownames(dat) %in% cells, , drop = FALSE]
    needed <- unique(c(features, group.by, split.by, bg.by))
    needed <- needed[!is.na(needed)]
    missing_cols <- setdiff(needed, names(dat))
    if (length(missing_cols) > 0)
      stop("Columns not found: ", paste(missing_cols, collapse = ", "),
           call. = FALSE)
  } else {
    stop("'seu' must be a Seurat object or a data.frame.", call. = FALSE)
  }
  dat
}

#' Build background colour map: group.by level → bg.by level
#' @noRd
.fp3_bg_map <- function(dat, group.by, bg.by) {
  # Returns a named list: name = group.by column, value = named char vector
  # mapping group level → bg.by level.
  if (is.null(bg.by) || !bg.by %in% names(dat)) {
    # Self-mapping: each group level maps to itself (alternating bg colours)
    bg_map <- lapply(group.by, function(g) {
      lvls <- levels(as.factor(dat[[g]]))
      stats::setNames(lvls, lvls)
    })
    return(stats::setNames(bg_map, group.by))
  }
  bg_map <- lapply(group.by, function(g) {
    tbl <- table(dat[[g]], dat[[bg.by]])
    if (max(Matrix::rowSums(tbl > 0)) > 1)
      warning("group.by must be nested within bg.by; background may be wrong.",
              call. = FALSE)
    stats::setNames(
      colnames(tbl)[apply(tbl, 1, function(x) which.max(x > 0))],
      rownames(tbl)
    )
  })
  stats::setNames(bg_map, group.by)
}

#' Resolve fill colours for one panel
#' @noRd
.fp3_resolve_colors <- function(fill.by, f, g, dat, split.by, features,
                                 palette, palcolor, median_expr) {
  if (fill.by == "feature") {
    palette_colors(features, palette = palette, palcolor = palcolor)
  } else if (fill.by == "group") {
    has_split <- nlevels(dat[[split.by]]) > 1
    lvls <- if (has_split) levels(dat[[split.by]]) else levels(dat[[g]])
    palette_colors(lvls, palette = palette, palcolor = palcolor)
  } else {
    # expression: continuous — median values determine gradient endpoints
    vals <- unlist(median_expr[, f, drop = FALSE])
    palette_colors(vals[is.finite(vals)], type = "continuous",
                   palette = palette, palcolor = palcolor)
  }
}

#' Legend key name based on fill.by / split.by status
#' @noRd
.fp3_keynm <- function(fill.by, f, g, split.by, has_split) {
  if (fill.by == "feature")   return("Features")
  if (fill.by == "expression") return("Median expression")
  if (has_split) return(split.by) else return(g)
}

#' feature-axis mode: one panel per group level, features on x-axis
#' @noRd
.fp3_feature_axis <- function(
    dat_full, features, group.by, split.by, has_split,
    type, fill.by, palette, palcolor, alpha,
    bg_palette, bg_palcolor, bg_alpha,
    add_box, box_color, box_width, box_ptsize,
    add_point, pt.color, pt.size, pt.alpha, jitter.width, jitter.height,
    cells.highlight, cols.highlight, sizes.highlight, alpha.highlight,
    add_trend, trend_color, trend_linewidth, trend_ptsize,
    add_stat, stat_color, stat_size, stat_stroke, stat_shape,
    add_line, line_color, line_size, line_type,
    comparisons, ref_group, pairwise_method,
    multiplegroup_comparisons, multiple_method, sig_label, sig_labelsize,
    global_ymin, global_ymax, y.min, y.max, y.trans, y.nbreaks,
    same.y.lims, sort, stack, flip, keep_empty, individual,
    title, subtitle, xlab, ylab,
    legend.position, legend.direction, legend.title,
    theme_fn, theme_args, aspect.ratio,
    combine, nrow, ncol, byrow, force, seed) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Resolve features
  features_avail <- features[features %in% names(dat_full)]
  if (length(features_avail) == 0) {
    stop("None of 'features' found in data.", call. = FALSE)
  }

  group_levels <- levels(dat_full[[group.by]])

  # Colours: by group level or by feature
  if (fill.by == "group") {
    colors <- palette_colors(group_levels, palette = palette, palcolor = palcolor)
    keynm <- group.by
    levels_order <- group_levels
  } else {
    colors <- palette_colors(features_avail, palette = palette, palcolor = palcolor)
    keynm <- "Features"
    levels_order <- features_avail
  }

  plist <- lapply(stats::setNames(group_levels, group_levels), function(g_lv) {
    sub_dat <- dat_full[dat_full[[group.by]] == g_lv, , drop = FALSE]
    if (nrow(sub_dat) == 0) return(NULL)

    # Pivot wide → long
    dat_long <- do.call(rbind, lapply(features_avail, function(f) {
      data.frame(
        group.by    = f,           # features become x-axis "groups"
        split.by    = if (has_split) sub_dat[[split.by]] else factor(""),
        value       = as.numeric(sub_dat[[f]]),
        fill.by_raw = if (fill.by == "group") g_lv else f,
        features    = f,
        stringsAsFactors = FALSE
      )
    }))
    dat_long[["group.by"]] <- factor(dat_long[["group.by"]], levels = features_avail)
    dat_long[["fill.by"]]  <- factor(dat_long[["fill.by_raw"]], levels = levels_order)
    dat_long[["group.unique"]] <- factor(
      paste0("sp-", dat_long[["split.by"]], "-gp-", dat_long[["group.by"]])
    )

    fin_vals <- dat_long[["value"]][is.finite(dat_long[["value"]])]
    if (length(fin_vals) == 0) return(NULL)

    y_min_use <- global_ymin %||% .fp3_resolve_ylim(y.min, fin_vals, "min")
    y_max_use <- global_ymax %||% .fp3_resolve_ylim(y.max, fin_vals, "max")
    ar <- if (isTRUE(flip) && !is.null(aspect.ratio)) 1 / aspect.ratio else aspect.ratio

    .fp3_build_panel(
      dat           = dat_long, f = features_avail[1], g = "group.by",
      type          = type, fill.by = fill.by,
      colors        = colors, colors_limits = NULL, keynm = keynm,
      bg_layer      = NULL,
      alpha         = alpha,
      add_box       = add_box, box_color = box_color,
      box_width     = box_width, box_ptsize = box_ptsize,
      add_point     = add_point, pt.color = pt.color,
      pt.size       = pt.size, pt.alpha = pt.alpha,
      jitter.width  = jitter.width, jitter.height = jitter.height,
      cells.highlight = cells.highlight, cols.highlight = cols.highlight,
      sizes.highlight = sizes.highlight, alpha.highlight = alpha.highlight,
      add_trend     = add_trend, trend_color = trend_color,
      trend_linewidth = trend_linewidth, trend_ptsize = trend_ptsize,
      add_stat      = add_stat, stat_color = stat_color,
      stat_size     = stat_size, stat_stroke = stat_stroke,
      stat_shape    = stat_shape,
      add_line      = add_line, line_color = line_color,
      line_size     = line_size, line_type = line_type,
      comparisons   = comparisons, ref_group = ref_group,
      pairwise_method = pairwise_method,
      multiplegroup_comparisons = multiplegroup_comparisons,
      multiple_method = multiple_method,
      sig_label     = sig_label, sig_labelsize = sig_labelsize,
      y_min_use     = y_min_use, y_max_use = y_max_use,
      y.trans       = y.trans, y.nbreaks = y.nbreaks,
      flip          = flip, stack = stack, keep_empty = keep_empty,
      title         = title %||% g_lv, subtitle = subtitle,
      xlab          = xlab %||% "Feature", ylab = ylab,
      legend.position  = legend.position,
      legend.direction = legend.direction,
      legend.title  = legend.title,
      theme_fn      = theme_fn, theme_args = theme_args,
      aspect.ratio  = ar,
      levels_order  = levels_order,
      seed          = seed
    )
  })

  plist <- Filter(Negate(is.null), plist)
  if (length(plist) == 0) return(invisible(NULL))

  # Stack
  if (isTRUE(stack) && length(plist) > 1) {
    return(.fp3_stack_combine(plist, flip = flip,
                               legend.position = legend.position,
                               ylab = ylab %||% "Expression level"))
  }
  if (isTRUE(combine)) {
    if (length(plist) == 1) return(plist[[1]])
    return(patchwork::wrap_plots(plotlist = plist, nrow = nrow,
                                  ncol = ncol, byrow = byrow))
  }
  plist
}
