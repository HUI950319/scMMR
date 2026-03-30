# ============================================================================
# ExpressionStatPlot2 — expression distribution plots for Seurat / data.frame
# ============================================================================
# Thin wrapper that:
#   1. Extracts feature data from a Seurat object or accepts a data.frame.
#   2. For plot.by = "group"  → delegates to UtilsR::plt_con().
#   3. For plot.by = "feature" → pivots to long format, builds a single ggplot
#      with features on the x-axis and groups as fill colour.
# ============================================================================

#' Expression distribution plots (violin / box / bar / dot)
#'
#' Visualises the expression of one or more features stratified by a grouping
#' variable.  Accepts a Seurat object (genes are fetched from the specified
#' assay / layer) or a plain data.frame.
#'
#' Two axis orientations are available via \code{plot.by}:
#' \itemize{
#'   \item \code{"group"} (default) — one panel per feature; groups on the
#'     x-axis.  Powered by \code{UtilsR::plt_con()}.
#'   \item \code{"feature"} — one combined panel; features on the x-axis,
#'     groups encoded by fill colour.
#' }
#'
#' @param seu A Seurat object **or** a data.frame containing feature columns
#'   and the grouping column.
#' @param features Character vector of feature names (gene names or metadata /
#'   data.frame column names).
#' @param group.by Column name used for grouping (x-axis in \code{"group"}
#'   mode, fill colour in \code{"feature"} mode).
#' @param split.by Column name for splitting into separate sub-plots.  Passed
#'   through to \code{plt_con()} in \code{"group"} mode; creates a patchwork
#'   row in \code{"feature"} mode.  Default \code{NULL}.
#' @param assay Seurat assay to use.  \code{NULL} = active assay.
#' @param layer Seurat layer to pull expression from.  Default \code{"data"}.
#' @param cells Optional character vector of cell barcodes to subset.
#'   \code{NULL} = all cells.
#' @param plot.by Axis orientation.  \code{"group"} (default) or
#'   \code{"feature"}.
#' @param type Plot geometry: \code{"violin"} (default), \code{"box"},
#'   \code{"bar"}, or \code{"dot"}.
#' @param fill.by What drives fill colour: \code{"group"} (default) or
#'   \code{"feature"}.
#' @param palette Colour palette name (passed to \code{UtilsR::plt_con()}).
#'   \code{NULL} = package default.
#' @param alpha Fill transparency (0–1).  Default \code{0.8}.
#' @param add_box Logical.  Overlay a narrow box-plot (median + IQR) on violin
#'   or bar.  Default \code{FALSE}.
#' @param box_color Colour of the overlay box-plot.  Default \code{"black"}.
#' @param box_width Width of the overlay box-plot.  Default \code{0.1}.
#' @param add_point Logical.  Overlay jittered raw points.  Default
#'   \code{FALSE}.
#' @param pt.color Point colour.  Default \code{"grey30"}.
#' @param pt.size Point size.  \code{NULL} = auto (scaled by cell count).
#' @param pt.alpha Point transparency.  Default \code{1}.
#' @param jitter.width Horizontal jitter extent.  Default \code{0.4}.
#' @param jitter.height Vertical jitter extent.  Default \code{0.1}.
#' @param add_trend Logical.  Connect group medians / means with a trend line.
#'   Default \code{FALSE}.
#' @param trend_color Trend line colour.  Default \code{"black"}.
#' @param trend_linewidth Trend line width.  Default \code{1}.
#' @param trend_ptsize Trend point size.  Default \code{2}.
#' @param comparisons List of length-2 vectors for pairwise statistical tests.
#'   \code{NULL} = none.
#' @param ref_group Reference group for one-vs-all comparisons.  \code{NULL} =
#'   no reference.
#' @param pairwise_method Statistical test for pairwise comparisons.  Default
#'   \code{"wilcox.test"}.
#' @param multiplegroup_comparisons Logical.  Add an overall group comparison
#'   label.  Default \code{FALSE}.
#' @param multiple_method Test for overall comparison.  Default
#'   \code{"kruskal.test"}.
#' @param sig_label Significance label format: \code{"p.signif"} (default) or
#'   \code{"p.format"}.
#' @param sig_labelsize Text size for significance labels.  Default \code{3.5}.
#' @param bg.by Column name for background colour bands (must be a
#'   super-grouping of \code{group.by}).  \code{NULL} = none.
#' @param bg_palette Palette for background bands.  \code{NULL} = default.
#' @param bg_alpha Transparency of background bands.  Default \code{0.15}.
#' @param same.y.lims Logical.  Share y limits across all panels.  Default
#'   \code{FALSE}.
#' @param y.min Minimum y value (numeric) or quantile string (e.g.
#'   \code{"q1"}).  \code{NULL} = data minimum.
#' @param y.max Maximum y value or quantile string.  \code{NULL} = data
#'   maximum.
#' @param y.nbreaks Number of y-axis tick breaks.  Default \code{5}.
#' @param sort Logical or \code{"increasing"}.  Sort groups by median of the
#'   first feature.  Default \code{FALSE}.
#' @param stack Logical.  Stack all features into a single faceted plot
#'   (\code{"group"} mode only).  Default \code{FALSE}.
#' @param flip Logical.  Flip coordinates (horizontal layout).  Default
#'   \code{FALSE}.
#' @param title Plot title.  \code{NULL} = none.
#' @param subtitle Plot subtitle.  \code{NULL} = none.
#' @param xlab X-axis label.  \code{NULL} = auto.
#' @param ylab Y-axis label.  \code{NULL} = auto.
#' @param legend.position Legend position.  Default \code{"right"}.
#' @param legend.direction Legend direction.  Default \code{"vertical"}.
#' @param aspect.ratio Aspect ratio of each panel.  \code{NULL} = free.
#' @param base_size Base font size.  Default \code{14}.
#' @param facet_nrow,facet_ncol Number of rows / columns in the combined
#'   patchwork layout.
#' @param combine Logical.  Combine panels with patchwork.  Default
#'   \code{TRUE}.
#' @param force Logical.  Skip the high-cardinality group check.  Default
#'   \code{FALSE}.
#' @param seed Random seed for jitter reproducibility.  Default \code{11}.
#'
#' @return A ggplot2 or patchwork object.
#'
#' @examples
#' \dontrun{
#' # Seurat object — expression of two genes by cluster
#' ExpressionStatPlot2(seu, features = c("CD3D", "CD8A"),
#'                     group.by = "seurat_clusters")
#'
#' # Feature-axis mode: features on x, cell types as colour
#' ExpressionStatPlot2(seu, features = c("CD3D", "CD8A", "MS4A1"),
#'                     group.by = "cell_type", plot.by = "feature",
#'                     type = "violin")
#'
#' # Plain data.frame
#' df <- data.frame(GeneA = rnorm(100), GeneB = rnorm(100),
#'                  cluster = sample(letters[1:4], 100, replace = TRUE))
#' ExpressionStatPlot2(df, features = c("GeneA", "GeneB"),
#'                     group.by = "cluster", type = "box")
#' }
#'
#' @export
ExpressionStatPlot2 <- function(
    seu,
    features,
    group.by,
    split.by              = NULL,
    assay                 = NULL,
    layer                 = "data",
    cells                 = NULL,
    plot.by               = c("group", "feature"),
    type                  = c("violin", "box", "bar", "dot"),
    fill.by               = c("group", "feature"),
    palette               = NULL,
    alpha                 = 0.8,
    add_box               = FALSE,
    box_color             = "black",
    box_width             = 0.1,
    add_point             = FALSE,
    pt.color              = "grey30",
    pt.size               = NULL,
    pt.alpha              = 1,
    jitter.width          = 0.4,
    jitter.height         = 0.1,
    add_trend             = FALSE,
    trend_color           = "black",
    trend_linewidth       = 1,
    trend_ptsize          = 2,
    comparisons           = NULL,
    ref_group             = NULL,
    pairwise_method       = "wilcox.test",
    multiplegroup_comparisons = FALSE,
    multiple_method       = "kruskal.test",
    sig_label             = c("p.signif", "p.format"),
    sig_labelsize         = 3.5,
    bg.by                 = NULL,
    bg_palette            = NULL,
    bg_alpha              = 0.15,
    same.y.lims           = FALSE,
    y.min                 = NULL,
    y.max                 = NULL,
    y.nbreaks             = 5,
    sort                  = FALSE,
    stack                 = FALSE,
    flip                  = FALSE,
    title                 = NULL,
    subtitle              = NULL,
    xlab                  = NULL,
    ylab                  = NULL,
    legend.position       = "right",
    legend.direction      = "vertical",
    aspect.ratio          = NULL,
    base_size             = 14,
    facet_nrow            = NULL,
    facet_ncol            = NULL,
    combine               = TRUE,
    force                 = FALSE,
    seed                  = 11
) {
  # ---- argument matching -------------------------------------------------------
  plot.by  <- match.arg(plot.by)
  type     <- match.arg(type)
  fill.by  <- match.arg(fill.by)
  sig_label <- match.arg(sig_label)

  # ---- dependency check --------------------------------------------------------
  if (!requireNamespace("UtilsR", quietly = TRUE))
    stop("Package 'UtilsR' is required. Install with: devtools::install('UtilsR')",
         call. = FALSE)

  # ---- extract data ------------------------------------------------------------
  dat <- .esp2_extract(seu, features = features, group.by = group.by,
                       split.by = split.by, bg.by = bg.by,
                       assay = assay, layer = layer, cells = cells)

  # ---- dispatch ----------------------------------------------------------------
  if (plot.by == "group") {
    # Delegate entirely to UtilsR::plt_con
    UtilsR::plt_con(
      data                      = dat,
      stat.by                   = features,
      group.by                  = group.by,
      split.by                  = split.by,
      bg.by                     = bg.by,
      type                      = type,
      fill.by                   = fill.by,
      palette                   = palette,
      alpha                     = alpha,
      add_box                   = add_box,
      box_color                 = box_color,
      box_width                 = box_width,
      add_point                 = add_point,
      pt.color                  = pt.color,
      pt.size                   = pt.size,
      pt.alpha                  = pt.alpha,
      jitter.width              = jitter.width,
      jitter.height             = jitter.height,
      add_trend                 = add_trend,
      trend_color               = trend_color,
      trend_linewidth           = trend_linewidth,
      trend_ptsize              = trend_ptsize,
      comparisons               = comparisons,
      ref_group                 = ref_group,
      pairwise_method           = pairwise_method,
      multiplegroup_comparisons = multiplegroup_comparisons,
      multiple_method           = multiple_method,
      sig_label                 = sig_label,
      sig_labelsize             = sig_labelsize,
      same.y.lims               = same.y.lims,
      y.min                     = y.min,
      y.max                     = y.max,
      y.nbreaks                 = y.nbreaks,
      sort                      = sort,
      stack                     = stack,
      flip                      = flip,
      title                     = title,
      subtitle                  = subtitle,
      xlab                      = xlab,
      ylab                      = ylab,
      legend.position           = legend.position,
      legend.direction          = legend.direction,
      aspect.ratio              = aspect.ratio,
      base_size                 = base_size,
      bg_palette                = bg_palette,
      bg_alpha                  = bg_alpha,
      facet_nrow                = facet_nrow,
      facet_ncol                = facet_ncol,
      combine                   = combine,
      force                     = force,
      seed                      = seed
    )
  } else {
    # plot.by == "feature": features on x-axis, groups as fill
    .esp2_feature_axis(
      dat          = dat,
      features     = features,
      group.by     = group.by,
      split.by     = split.by,
      bg.by        = bg.by,
      type         = type,
      fill.by      = fill.by,
      palette      = palette,
      alpha        = alpha,
      add_box      = add_box,
      box_color    = box_color,
      box_width    = box_width,
      add_point    = add_point,
      pt.color     = pt.color,
      pt.size      = pt.size,
      pt.alpha     = pt.alpha,
      jitter.width = jitter.width,
      jitter.height = jitter.height,
      add_trend    = add_trend,
      trend_color  = trend_color,
      trend_linewidth = trend_linewidth,
      trend_ptsize = trend_ptsize,
      comparisons  = comparisons,
      ref_group    = ref_group,
      pairwise_method = pairwise_method,
      multiplegroup_comparisons = multiplegroup_comparisons,
      multiple_method = multiple_method,
      sig_label    = sig_label,
      sig_labelsize = sig_labelsize,
      bg_palette   = bg_palette,
      bg_alpha     = bg_alpha,
      same.y.lims  = same.y.lims,
      y.min        = y.min,
      y.max        = y.max,
      y.nbreaks    = y.nbreaks,
      flip         = flip,
      title        = title,
      subtitle     = subtitle,
      xlab         = xlab %||% "Feature",
      ylab         = ylab %||% "Expression",
      legend.position  = legend.position,
      legend.direction = legend.direction,
      aspect.ratio = aspect.ratio,
      base_size    = base_size,
      facet_nrow   = facet_nrow,
      facet_ncol   = facet_ncol,
      combine      = combine,
      force        = force,
      seed         = seed
    )
  }
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Extract feature data from Seurat or data.frame
#' @noRd
.esp2_extract <- function(seu, features, group.by, split.by, bg.by,
                           assay, layer, cells) {
  if (inherits(seu, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE))
      stop("SeuratObject is required for Seurat input.", call. = FALSE)
    if (!is.null(assay)) SeuratObject::DefaultAssay(seu) <- assay

    # Columns to fetch: features + group.by + optional split.by / bg.by
    fetch_vars <- unique(c(features, group.by, split.by, bg.by))

    dat <- tryCatch(
      SeuratObject::FetchData(seu, vars = fetch_vars, layer = layer,
                              cells = cells),
      error = function(e) stop(conditionMessage(e), call. = FALSE)
    )
  } else if (is.data.frame(seu)) {
    dat <- seu
    if (!is.null(cells)) dat <- dat[rownames(dat) %in% cells, , drop = FALSE]
    # Check required columns
    needed <- unique(c(features, group.by, split.by, bg.by))
    missing_cols <- setdiff(needed, names(dat))
    if (length(missing_cols) > 0)
      stop("Columns not found in data.frame: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
  } else {
    stop("'seu' must be a Seurat object or a data.frame.", call. = FALSE)
  }
  dat
}


#' Build a "feature-axis" plot (features on x, groups as fill)
#'
#' Pivots wide → long, then builds a ggplot analogous to plt_con but with
#' Feature on the x-axis and group.by levels as fill colour.
#' @noRd
.esp2_feature_axis <- function(
    dat, features, group.by, split.by, bg.by,
    type, fill.by, palette, alpha,
    add_box, box_color, box_width,
    add_point, pt.color, pt.size, pt.alpha,
    jitter.width, jitter.height,
    add_trend, trend_color, trend_linewidth, trend_ptsize,
    comparisons, ref_group, pairwise_method,
    multiplegroup_comparisons, multiple_method,
    sig_label, sig_labelsize,
    bg_palette, bg_alpha,
    same.y.lims, y.min, y.max, y.nbreaks,
    flip, title, subtitle, xlab, ylab,
    legend.position, legend.direction, aspect.ratio, base_size,
    facet_nrow, facet_ncol, combine, force, seed) {

  # Helper for null-coalescing (mirrors UtilsR's %||%)
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Ensure group.by is a factor
  if (!is.factor(dat[[group.by]]))
    dat[[group.by]] <- factor(dat[[group.by]], levels = unique(dat[[group.by]]))

  # --- Resolve colours --------------------------------------------------------
  if (fill.by == "group") {
    color_levels <- levels(dat[[group.by]])
  } else {
    color_levels <- features
  }
  colors <- if (!is.null(palette)) {
    UtilsR::palette_colors(color_levels, palette = palette)
  } else {
    UtilsR::palette_colors(color_levels)
  }

  # --- Build one long-format data.frame per split level ----------------------
  build_one <- function(sub_dat, panel_title) {
    # Pivot features long
    dat_long <- do.call(rbind, lapply(features, function(f) {
      data.frame(
        Feature  = f,
        Value    = sub_dat[[f]],
        Group    = sub_dat[[group.by]],
        stringsAsFactors = FALSE
      )
    }))
    dat_long[["Feature"]] <- factor(dat_long[["Feature"]], levels = features)
    dat_long[["Group"]]   <- factor(dat_long[["Group"]],
                                    levels = levels(sub_dat[[group.by]]))

    # fill.by mapping
    if (fill.by == "group") {
      dat_long[["fill_var"]] <- dat_long[["Group"]]
    } else {
      dat_long[["fill_var"]] <- dat_long[["Feature"]]
    }

    # Remove non-finite values
    dat_long <- dat_long[is.finite(dat_long[["Value"]]), , drop = FALSE]
    if (nrow(dat_long) == 0) {
      warning("No finite values for this panel. Skipping.", call. = FALSE)
      return(NULL)
    }

    # Y limits
    vals <- dat_long[["Value"]]
    y_min_use <- if (!is.null(y.min)) {
      if (is.character(y.min)) {
        q <- as.numeric(sub("^q(\\d+)$", "\\1", y.min)) / 100
        stats::quantile(vals, q, na.rm = TRUE)
      } else y.min
    } else min(vals, na.rm = TRUE)

    y_max_use <- if (!is.null(y.max)) {
      if (is.character(y.max)) {
        q <- as.numeric(sub("^q(\\d+)$", "\\1", y.max)) / 100
        stats::quantile(vals, q, na.rm = TRUE)
      } else y.max
    } else max(vals, na.rm = TRUE)

    if (is.null(pt.size)) {
      pt_sz <- min(3000 / nrow(dat_long), 0.5)
    } else {
      pt_sz <- pt.size
    }

    # Base plot
    p <- ggplot2::ggplot(dat_long, ggplot2::aes(
      x    = .data[["Feature"]],
      y    = .data[["Value"]],
      fill = .data[["fill_var"]]
    ))

    # Type-specific geom
    if (type == "violin") {
      p <- p + ggplot2::geom_violin(
        scale = "width", trim = TRUE, alpha = alpha,
        position = ggplot2::position_dodge(width = 0.9)
      )
    } else if (type == "box") {
      p <- p +
        ggplot2::geom_boxplot(
          position = ggplot2::position_dodge(width = 0.9),
          color = "black", width = 0.8, alpha = alpha,
          outlier.shape = NA
        ) +
        ggplot2::stat_summary(
          fun = stats::median, geom = "point",
          position = ggplot2::position_dodge(width = 0.9),
          color = "black", fill = "white", size = 1.5, shape = 21
        )
    } else if (type == "bar") {
      p <- p +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::stat_summary(
          fun = mean, geom = "col",
          position = ggplot2::position_dodge(width = 0.9),
          width = 0.8, color = "black", alpha = alpha
        ) +
        ggplot2::stat_summary(
          fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1),
          geom = "errorbar",
          position = ggplot2::position_dodge(width = 0.9),
          width = 0.2, color = "black"
        )
      y_min_use <- min(0, y_min_use)
    } else if (type == "dot") {
      brks <- seq(min(dat_long[["Value"]], na.rm = TRUE),
                  max(dat_long[["Value"]], na.rm = TRUE), length.out = 15)
      bins  <- cut(dat_long[["Value"]], breaks = brks, include.lowest = TRUE)
      bin_lvls <- levels(bins)
      bins_mid <- vapply(bin_lvls, function(lbl) {
        nums <- regmatches(lbl, gregexpr("-?[0-9.]+[eE]?[+-]?[0-9]*", lbl))[[1]]
        mean(as.numeric(nums), na.rm = TRUE)
      }, numeric(1))
      dat_long[["bins"]] <- bins_mid[as.character(bins)]
      dot_df <- stats::aggregate(
        Value ~ Feature + fill_var + bins,
        data = dat_long, FUN = length
      )
      names(dot_df)[names(dot_df) == "Value"] <- "n"
      dot_df[["Feature"]] <- factor(dot_df[["Feature"]], levels = features)
      p <- p +
        ggplot2::geom_point(
          data = dot_df,
          ggplot2::aes(
            x    = .data[["Feature"]],
            y    = .data[["bins"]],
            fill = .data[["fill_var"]],
            size = .data[["n"]]
          ),
          shape = 21, alpha = alpha, inherit.aes = FALSE,
          position = ggplot2::position_dodge(width = 0.9)
        ) +
        ggplot2::scale_size_area(name = "Count", max_size = 6, n.breaks = 4) +
        ggplot2::guides(size = ggplot2::guide_legend(
          override.aes = list(fill = "grey30", shape = 21), order = 2
        ))
    }

    # Overlays
    if (isTRUE(add_box) && type %in% c("violin", "bar")) {
      p <- p +
        ggplot2::geom_boxplot(
          position = ggplot2::position_dodge(width = 0.9),
          color = box_color, fill = box_color,
          width = box_width, outlier.shape = NA, show.legend = FALSE
        ) +
        ggplot2::stat_summary(
          fun = stats::median, geom = "point",
          position = ggplot2::position_dodge(width = 0.9),
          color = "black", fill = "white", size = 1.5, shape = 21
        )
    }

    if (isTRUE(add_point) && type != "dot") {
      set.seed(seed)
      p <- p + ggplot2::geom_point(
        ggplot2::aes(x = .data[["Feature"]], y = .data[["Value"]]),
        inherit.aes = FALSE,
        color = pt.color, size = pt_sz, alpha = pt.alpha,
        position = ggplot2::position_jitterdodge(
          jitter.width = jitter.width, jitter.height = jitter.height,
          dodge.width = 0.9, seed = seed
        ),
        show.legend = FALSE
      )
    }

    if (isTRUE(add_trend) && type %in% c("violin", "box", "bar")) {
      trend_fun <- if (type == "bar") mean else stats::median
      p <- p +
        ggplot2::stat_summary(
          fun = trend_fun, geom = "line",
          ggplot2::aes(group = .data[["fill_var"]]),
          color = trend_color, linewidth = trend_linewidth,
          position = ggplot2::position_dodge(width = 0.9)
        ) +
        ggplot2::stat_summary(
          fun = trend_fun, geom = "point",
          ggplot2::aes(group = .data[["fill_var"]]),
          color = "black", fill = "white",
          size = trend_ptsize, shape = 21,
          position = ggplot2::position_dodge(width = 0.9)
        )
    }

    # Statistical comparisons
    if (!is.null(comparisons)) {
      p <- p + ggpubr::stat_compare_means(
        ggplot2::aes(x = .data[["Feature"]], y = .data[["Value"]]),
        comparisons = comparisons, ref.group = ref_group,
        method = pairwise_method, label = sig_label,
        label.y = y_max_use, size = sig_labelsize,
        step.increase = 0.1, tip.length = 0.03, vjust = 0
      )
      y_max_use <- y_max_use + (y_max_use - y_min_use) * 0.15 * length(comparisons)
    }

    if (isTRUE(multiplegroup_comparisons)) {
      p <- p + ggpubr::stat_compare_means(
        ggplot2::aes(x = .data[["Feature"]], y = .data[["Value"]]),
        method = multiple_method, label = sig_label,
        label.y = y_max_use, size = sig_labelsize,
        vjust = 1.2, hjust = 0
      )
      y_max_use <- y_max_use + (y_max_use - y_min_use) * 0.1
    }

    # Coordinates
    if (isTRUE(flip)) {
      p <- p + ggplot2::coord_flip(ylim = c(y_min_use, y_max_use))
    } else {
      p <- p + ggplot2::coord_cartesian(ylim = c(y_min_use, y_max_use))
    }

    # Scales
    p <- p +
      ggplot2::scale_y_continuous(n.breaks = y.nbreaks) +
      ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::guides(fill = ggplot2::guide_legend(
        title.hjust = 0, order = 1,
        override.aes = list(alpha = 1)
      ))

    # Theme
    p <- p +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        legend.position   = legend.position,
        legend.direction  = legend.direction,
        aspect.ratio      = aspect.ratio,
        panel.grid.major.y = ggplot2::element_line(
          color = "grey", linetype = 2, linewidth = 0.3
        )
      ) +
      ggplot2::labs(
        title = panel_title %||% title,
        subtitle = subtitle,
        x = xlab, y = ylab,
        fill = if (fill.by == "group") group.by else "Feature"
      )

    p
  }

  # --- Handle split.by --------------------------------------------------------
  if (!is.null(split.by)) {
    if (!is.factor(dat[[split.by]]))
      dat[[split.by]] <- factor(dat[[split.by]], levels = unique(dat[[split.by]]))
    split_levels <- levels(dat[[split.by]])
    plist <- lapply(split_levels, function(lv) {
      sub_dat <- dat[dat[[split.by]] == lv, , drop = FALSE]
      build_one(sub_dat, panel_title = lv)
    })
    plist <- Filter(Negate(is.null), plist)
    if (length(plist) == 0) return(invisible(NULL))
    return(
      patchwork::wrap_plots(plist, nrow = facet_nrow, ncol = facet_ncol) +
        patchwork::plot_layout(guides = "collect")
    )
  }

  # --- No split.by ------------------------------------------------------------
  p <- build_one(dat, panel_title = title)
  p
}
