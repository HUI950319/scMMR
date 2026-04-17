# ============================================================================
# PlotScatter — Scatter Correlation Plot
# ============================================================================

#' Scatter Correlation Plot
#'
#' Create a scatter plot showing the correlation between two variables (genes
#' and/or metadata columns) in single-cell data. Supports coloring by groups,
#' faceting, per-group or global correlation statistics, and multiple
#' correlation methods.
#'
#' @details
#' The first argument \code{object} can be either a \strong{Seurat object} or
#' a \strong{data.frame} (/ tibble).
#'
#' \strong{When a Seurat object is provided}, variables \code{var1} and
#' \code{var2} can be:
#' \itemize{
#'   \item Gene names (expression extracted from the specified assay/layer).
#'   \item Metadata column names in \code{object@@meta.data}.
#'   \item Any combination of the two.
#' }
#'
#' \strong{When a data.frame is provided}, \code{var1}, \code{var2},
#' \code{group.by}, and \code{split.by} must all be column names of the
#' data.frame.  The Seurat-specific parameters (\code{assay}, \code{layer},
#' \code{cells}) are ignored.
#'
#' When \code{group.by} is set, points are colored by the grouping variable
#' and correlation statistics (\code{ggpubr::stat_cor}) are shown per group.
#' When \code{split.by} is set, the plot is faceted by that variable.
#'
#' When \code{marginal} is set to a plot type other than \code{"none"},
#' marginal distribution plots are added to the x- and y-axes via
#' \code{patchwork} layout. Supported types include density, histogram,
#' boxplot, violin, and densigram (histogram + density overlay). When
#' \code{group.by} is set, the marginal plots are colored/filled by group.
#' Note: marginal plots are incompatible with faceting (\code{split.by});
#' if both are specified, marginal plots are silently skipped.
#'
#' \strong{Note:} When marginal plots are added, the return value is a
#' \code{patchwork} object. You can still use \code{patchwork::&} or
#' \code{patchwork::plot_annotation()} for further modifications.
#'
#' @param object A Seurat object or a data.frame containing the variables
#'   to plot.
#' @param var1 Character. First variable (x-axis): a gene name, metadata
#'   column, or data.frame column name.
#' @param var2 Character. Second variable (y-axis): a gene name, metadata
#'   column, or data.frame column name.
#' @param group.by Character. Optional column to color points and compute
#'   per-group correlations. Default: \code{NULL} (single color).
#' @param split.by Character. Optional column for faceting
#'   (\code{facet_wrap}). Default: \code{NULL}.
#' @param cells Character vector. Cell barcodes to include (Seurat only).
#'   Default: \code{NULL} (all cells).
#' @param assay Character. Seurat assay to use for gene expression (Seurat
#'   only). Default: \code{DefaultAssay(object)}.
#' @param layer Character. Data layer to extract (Seurat only). Default:
#'   \code{"data"} (log-normalized).
#' @param method Correlation method for \code{ggpubr::stat_cor}:
#'   \code{"spearman"}, \code{"pearson"}, or \code{"kendall"}.
#'   Default: \code{"spearman"}.
#' @param smooth.method Smoothing method for \code{geom_smooth}. Default:
#'   \code{"lm"}. Set to \code{"loess"}, \code{"gam"}, etc. as needed.
#' @param show.cor Logical. Show correlation statistics on the plot.
#'   Default: \code{TRUE}.
#' @param show.smooth Logical. Show regression / smooth line.
#'   Default: \code{TRUE}.
#' @param show.rug Logical. Show rug plots on the axes.
#'   Default: \code{FALSE}.
#' @param cor.digits Integer. Number of decimal digits for correlation display.
#'   Default: 3.
#' @param cor.size Numeric. Font size for correlation text. Default: 4.
#' @param point.size Numeric. Size of scatter points. Default: 1.
#' @param point.alpha Numeric. Transparency of scatter points (0--1).
#'   Default: 0.6.
#' @param smooth.size Numeric. Width of the regression line. Default: 1.
#' @param smooth.color Character. Color of the smooth line when
#'   \code{group.by = NULL}. Default: \code{"#fdc086"}.
#'   Ignored when \code{group.by} is set (line color follows group).
#' @param show.se Logical. Show confidence interval around smooth line.
#'   Default: \code{TRUE}.
#' @param se.fill Character. Fill color for the confidence interval ribbon
#'   when \code{group.by = NULL}. Default: \code{NULL} (follows
#'   \code{smooth.color}). When \code{group.by} is set, the fill follows
#'   group colors by default; set this to override with a fixed color.
#' @param se.alpha Numeric. Transparency of the confidence interval ribbon
#'   (0--1). Default: \code{0.2}.
#' @param ncol Integer. Number of columns for faceting when \code{split.by}
#'   is set. Default: 3.
#' @param palette Character. Color palette name passed to
#'   \code{\link{palette_colors}}. Default: \code{"Paired"}.
#' @param palcolor Character vector. Custom colors overriding \code{palette}.
#'   Default: \code{NULL}.
#' @param point.color Character. Fixed point color when \code{group.by = NULL}.
#'   Default: \code{"#984ea3"}.
#' @param rug.color Character. Color for rug marks. Default: \code{"#7fc97f"}.
#' @param title Character. Plot title. Default: \code{NULL} (auto-generated).
#' @param marginal Character. Type of marginal distribution plot to add on
#'   the x- and y-axes: \code{"none"} (default), \code{"density"},
#'   \code{"histogram"}, \code{"boxplot"}, \code{"violin"}, or
#'   \code{"densigram"} (histogram + density overlay). Requires \pkg{patchwork}.
#'   Ignored when \code{split.by} is used.
#' @param marginal.size Numeric. Relative size of marginal plots compared
#'   to the main scatter panel. Default: 5 (i.e. the main panel is 5 times
#'   larger than the marginal).
#' @param global.cor Logical. When \code{group.by} is set and
#'   \code{global.cor = TRUE}, show a single overall smooth line and
#'   correlation statistics (ignoring group), rather than per-group lines
#'   and statistics. The points are still colored by group. Default:
#'   \code{FALSE}.
#' @param raster Logical. If \code{TRUE}, rasterize the point layer via
#'   \code{ggrastr::rasterise()} to reduce file size for large datasets.
#'   Default: \code{NULL} (auto: \code{TRUE} when > 50 000 cells).
#' @param raster.dpi Integer. DPI for rasterized points. Default: 300.
#' @param ... Additional arguments passed to \code{geom_point}.
#'
#' @return A \code{ggplot} object (or a \code{patchwork} object when
#'   \code{marginal != "none"}).
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#'
#' # --- Seurat object input ---
#' # Two genes
#' PlotScatter(seu, var1 = "METTL3", var2 = "SETD2")
#'
#' # Gene vs metadata
#' PlotScatter(seu, var1 = "nFeature_RNA", var2 = "PTH",
#'             method = "pearson")
#'
#' # Color by cell type
#' PlotScatter(seu, var1 = "TP53", var2 = "MDM2",
#'             group.by = "celltype", palette = "npg")
#'
#' # With marginal boxplot
#' PlotScatter(seu, var1 = "TP53", var2 = "MDM2",
#'             group.by = "celltype", marginal = "boxplot")
#'
#' # --- data.frame input ---
#' df <- data.frame(x = rnorm(200), y = rnorm(200),
#'                  grp = sample(c("A", "B"), 200, replace = TRUE))
#' PlotScatter(df, var1 = "x", var2 = "y")
#' PlotScatter(df, var1 = "x", var2 = "y", group.by = "grp",
#'             marginal = "density")
#' }
#'
#' @seealso \code{\link{palette_colors}}
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_rug
#'   facet_wrap labs theme_bw theme element_text element_blank
#' @importFrom rlang .data
#' @export
PlotScatter <- function(object,
                        var1,
                        var2,
                        group.by      = NULL,
                        split.by      = NULL,
                        cells         = NULL,
                        assay         = NULL,
                        layer         = "data",
                        method        = c("spearman", "pearson", "kendall"),
                        smooth.method = "lm",
                        show.cor      = TRUE,
                        show.smooth   = TRUE,
                        show.rug      = FALSE,
                        cor.digits    = 3,
                        cor.size      = 4,
                        point.size    = 1,
                        point.alpha   = 0.6,
                        smooth.size   = 1,
                        smooth.color  = "#fdc086",
                        show.se       = TRUE,
                        se.fill       = NULL,
                        se.alpha      = 0.2,
                        ncol          = 3,
                        palette       = "Paired",
                        palcolor      = NULL,
                        point.color   = "#984ea3",
                        rug.color     = "#7fc97f",
                        title         = NULL,
                        marginal      = c("none", "density", "histogram",
                                          "boxplot", "violin", "densigram"),
                        marginal.size = 5,
                        global.cor    = FALSE,
                        raster        = NULL,
                        raster.dpi    = 300,
                        ...) {

  # ── Multi-var2 dispatch: wrap_plots when length(var2) > 1 ──
  if (length(var2) > 1) {
    p_list <- lapply(var2, function(v) {
      PlotScatter(
        object = object, var1 = var1, var2 = v,
        group.by = group.by, split.by = split.by, cells = cells,
        assay = assay, layer = layer, method = method,
        smooth.method = smooth.method, show.cor = show.cor,
        show.smooth = show.smooth, show.rug = show.rug,
        cor.digits = cor.digits, cor.size = cor.size,
        point.size = point.size, point.alpha = point.alpha,
        smooth.size = smooth.size, smooth.color = smooth.color,
        show.se = show.se, se.fill = se.fill, se.alpha = se.alpha,
        ncol = ncol, palette = palette, palcolor = palcolor,
        point.color = point.color, rug.color = rug.color,
        title = NULL, marginal = marginal, marginal.size = marginal.size,
        global.cor = global.cor, raster = raster, raster.dpi = raster.dpi,
        ...
      )
    })
    return(patchwork::wrap_plots(p_list, ncol = ncol))
  }

  # ── Input validation ──
  method   <- match.arg(method)
  marginal <- match.arg(marginal)

  is_seurat <- inherits(object, "Seurat")
  is_df     <- is.data.frame(object)

  if (!is_seurat && !is_df) {
    stop("'object' must be a Seurat object or a data.frame.", call. = FALSE)
  }

  # ── Build plot_df ──
  if (is_df) {
    # --- data.frame path ---
    vars_needed <- unique(c(var1, var2, group.by, split.by))
    missing_cols <- setdiff(vars_needed, colnames(object))
    if (length(missing_cols) > 0) {
      stop("Column(s) not found in data.frame: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    plot_df <- object[, vars_needed, drop = FALSE]

  } else {
    # --- Seurat path ---
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)

    # Determine which cells to use
    if (!is.null(cells)) {
      cells <- intersect(cells, colnames(object))
      if (length(cells) == 0) stop("No valid cells found.", call. = FALSE)
    } else {
      cells <- colnames(object)
    }

    # Extract variables
    meta_cols <- colnames(object@meta.data)
    all_genes <- rownames(object[[assay]])
    vars_needed <- unique(c(var1, var2, group.by, split.by))

    # Separate genes vs metadata
    is_gene <- vars_needed %in% all_genes
    is_meta <- vars_needed %in% meta_cols
    unknown <- vars_needed[!is_gene & !is_meta]

    if (length(unknown) > 0) {
      stop("Variable(s) not found in assay or meta.data: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }

    gene_vars <- vars_needed[is_gene & !is_meta]
    # If a variable is both gene and metadata, prefer gene
    # (unless it's group.by/split.by which should be metadata)
    force_meta <- c(group.by, split.by)
    gene_vars <- setdiff(gene_vars, force_meta)

    # Build data.frame
    plot_df <- data.frame(row.names = cells)

    # Add gene expression columns
    if (length(gene_vars) > 0) {
      expr_mat <- SeuratObject::GetAssayData(object, assay = assay,
                                             layer = layer)
      for (g in gene_vars) {
        plot_df[[g]] <- as.numeric(expr_mat[g, cells])
      }
    }

    # Add metadata columns
    for (v in setdiff(vars_needed, gene_vars)) {
      plot_df[[v]] <- object@meta.data[cells, v]
    }
  }

  # Check var1, var2 are numeric
  if (!is.numeric(plot_df[[var1]])) {
    stop("'var1' (", var1, ") must be numeric for correlation analysis.",
         call. = FALSE)
  }
  if (!is.numeric(plot_df[[var2]])) {
    stop("'var2' (", var2, ") must be numeric for correlation analysis.",
         call. = FALSE)
  }

  # Remove NA rows for the two main variables
  complete <- is.finite(plot_df[[var1]]) & is.finite(plot_df[[var2]])
  plot_df <- plot_df[complete, , drop = FALSE]
  if (nrow(plot_df) < 3) {
    stop("Insufficient data points for correlation analysis (n < 3).",
         call. = FALSE)
  }

  # ── Auto raster ──
  if (is.null(raster)) {
    raster <- nrow(plot_df) > 50000
  }

  # ── Build plot ──
  if (!is.null(group.by)) {
    # Convert to factor if not already
    if (!is.factor(plot_df[[group.by]])) {
      plot_df[[group.by]] <- factor(plot_df[[group.by]])
    }
    # Colors
    group_levels <- levels(plot_df[[group.by]])
    cols <- palette_colors(group_levels, palette = palette,
                           palcolor = palcolor)

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[[var1]], y = .data[[var2]],
                   color = .data[[group.by]])
    )
  } else {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[[var1]], y = .data[[var2]])
    )
  }

  # ── Point layer ──
  point_args <- list(size = point.size, alpha = point.alpha, ...)
  if (is.null(group.by)) {
    point_args$color <- point.color
  }
  point_layer <- do.call(ggplot2::geom_point, point_args)

  if (isTRUE(raster)) {
    if (requireNamespace("ggrastr", quietly = TRUE)) {
      point_layer <- ggrastr::rasterise(point_layer, dpi = raster.dpi)
    } else {
      message("Install 'ggrastr' for rasterized points. ",
              "Using default vector points.")
    }
  }
  p <- p + point_layer

  # ── Smooth line ──
  if (isTRUE(show.smooth)) {
    if (!is.null(group.by) && !isTRUE(global.cor)) {
      # Per-group smooth lines (color follows group)
      smooth_args <- list(
        method    = smooth.method,
        se        = show.se,
        na.rm     = TRUE,
        linewidth = smooth.size,
        alpha     = se.alpha
      )
      # Override SE fill with a fixed color if se.fill is specified
      if (!is.null(se.fill)) smooth_args$fill <- se.fill
      p <- p + do.call(ggplot2::geom_smooth, smooth_args)
    } else {
      # Global smooth line (single color, ignoring group)
      smooth_args <- list(
        method     = smooth.method,
        se         = show.se,
        na.rm      = TRUE,
        linewidth  = smooth.size,
        color      = smooth.color,
        fill       = if (!is.null(se.fill)) se.fill else smooth.color,
        alpha      = se.alpha,
        inherit.aes = if (!is.null(group.by) && isTRUE(global.cor)) FALSE else TRUE,
        mapping     = if (!is.null(group.by) && isTRUE(global.cor))
          ggplot2::aes(x = .data[[var1]], y = .data[[var2]]) else NULL
      )
      p <- p + do.call(ggplot2::geom_smooth, smooth_args)
    }
  }

  # ── Rug ──
  if (isTRUE(show.rug)) {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_rug(alpha = 0.4)
    } else {
      p <- p + ggplot2::geom_rug(color = rug.color, alpha = 0.4)
    }
  }

  # ── Correlation statistics ──
  if (isTRUE(show.cor)) {
    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      message("Install 'ggpubr' to display correlation statistics on the plot.")
    } else {
      # Use rho for spearman, R for pearson
      coef_name <- if (method == "spearman") "rho" else "R"

      if (!is.null(group.by) && isTRUE(global.cor)) {
        p <- p + ggpubr::stat_cor(
          mapping = ggplot2::aes(x = .data[[var1]], y = .data[[var2]]),
          inherit.aes = FALSE,
          method  = method,
          cor.coef.name = coef_name,
          digits  = cor.digits,
          size    = cor.size,
          color   = "black",
          label.x.npc = "left",
          label.y.npc = "top"
        )
      } else {
        p <- p + ggpubr::stat_cor(
          method  = method,
          cor.coef.name = coef_name,
          digits  = cor.digits,
          size    = cor.size,
          label.x.npc = "left",
          label.y.npc = "top"
        )
      }
    }
  }

  # ── Color scale ──
  if (!is.null(group.by)) {
    p <- p +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          title        = group.by,
          override.aes = list(size = 3, alpha = 1)
        )
      )
  }

  # ── Faceting ──
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", split.by)),
      ncol = ncol, scales = "free"
    )
  }

  # ── Title ──
  if (is.null(title)) {
    title <- paste0(var1, " vs ", var2)
  }

  p <- p +
    ggplot2::labs(x = var1, y = var2, title = title, color = group.by) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, size = 14,
                                               face = "bold"),
      axis.title       = ggplot2::element_text(size = 12),
      axis.text        = ggplot2::element_text(size = 10, color = "black"),
      strip.text       = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(face = "bold", size = 11),
      legend.text      = ggplot2::element_text(size = 10)
    )

  # ── Marginal distribution plots (patchwork layout) ──
  if (marginal != "none") {
    if (!is.null(split.by)) {
      message("Marginal plots are incompatible with faceting (split.by). ",
              "Skipping marginal plots.")
    } else if (!requireNamespace("patchwork", quietly = TRUE)) {
      message("Install 'patchwork' for marginal distribution plots: ",
              "install.packages('patchwork')")
    } else {
      # Resolve cols for ungrouped case
      cols_use <- if (!is.null(group.by)) cols else NULL

      # Top marginal: distribution of var1
      p_top <- .build_marginal(
        plot_df, var1, marginal, group.by,
        cols_use, point.color
      )

      # Right marginal: distribution of var2, flipped
      p_right <- .build_marginal(
        plot_df, var2, marginal, group.by,
        cols_use, point.color
      ) + ggplot2::coord_flip()

      # Move title to top annotation; remove from main
      # Override legend layout for bottom placement: horizontal, single row
      p <- p +
        ggplot2::ggtitle(NULL) +
        ggplot2::guides(
          color = ggplot2::guide_legend(
            title          = group.by,
            direction      = "horizontal",
            title.position = "left",
            title.hjust    = 0.5,
            nrow           = 1,
            override.aes   = list(size = 3, alpha = 1)
          )
        )

      # Assemble with patchwork:
      #   [top marginal] [spacer  ]
      #   [main scatter ] [right   ]
      # guides = "collect" pulls legend out of individual panels
      # plot_annotation(theme = ...) controls where the collected legend goes
      p <- p_top + patchwork::plot_spacer() +
        p + p_right +
        patchwork::plot_layout(
          ncol    = 2,
          nrow    = 2,
          widths  = c(marginal.size, 1),
          heights = c(1, marginal.size),
          guides  = "collect"
        ) +
        patchwork::plot_annotation(
          title = title,
          theme = ggplot2::theme(
            plot.title      = ggplot2::element_text(hjust = 0.5, size = 14,
                                                    face = "bold"),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box       = "horizontal",
            legend.title     = ggplot2::element_text(face = "bold", size = 11),
            legend.text      = ggplot2::element_text(size = 10),
            legend.key.size  = ggplot2::unit(0.4, "cm"),
            legend.margin    = ggplot2::margin(t = 2, b = 2)
          )
        )
    }
  }

  return(p)
}
