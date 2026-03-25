# ============================================================================
# PlotAlluvia — Alluvial (stacked area + bar) plot
# ============================================================================

#' Alluvial Plot for Cell Population Composition
#'
#' Create an alluvial (stacked area + bar) plot showing cell type composition
#' changes across conditions/groups.
#'
#' @param cellmeta A Seurat object or data.frame containing cell metadata
#'   (e.g. \code{seurat_obj} or \code{seurat_obj@@meta.data}).
#' @param by Character string specifying the grouping variable
#'   (e.g. \code{"group"}, \code{"sample"}).
#' @param fill Character string specifying the fill variable
#'   (e.g. \code{"cell_type_pred"}).
#' @param colors Optional character vector of colours. If \code{NULL} and
#'   \code{palette} is also \code{NULL}, default ggplot2 colours are used.
#'   A named vector maps specific levels to colours.
#' @param palette Optional RColorBrewer palette name (e.g. \code{"Paired"},
#'   \code{"Set2"}) for auto-generating colours when \code{colors = NULL}.
#'   Default \code{NULL} (use ggplot2 defaults).
#' @param flow.alpha Numeric (0-1), transparency of the area flow layer.
#'   Default 0.4.
#' @param bar.width Numeric (0-1) specifying bar width. Default 0.5.
#' @param bar.color Bar border colour. Default \code{"gray50"}.
#' @param show.label Logical, whether to show cell type labels inside bars.
#'   Default \code{FALSE}.
#' @param label.size Numeric, font size for bar labels. Default 3.
#' @param label.color Label text colour. Default \code{"black"}.
#' @param show.pct Logical, whether to append within-group percentage to
#'   labels. Default \code{FALSE}. Only used when \code{show.label = TRUE}.
#' @param y.percent Logical. If \code{TRUE} (default), y-axis shows
#'   percentages (all bars equal height). If \code{FALSE}, y-axis shows
#'   raw cell counts.
#' @param legend.ncol Integer, number of legend columns. Default 1.
#' @param base.size Numeric, base font size for the theme. Default 15.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # From Seurat object directly
#' PlotAlluvia(seurat_obj, by = "group", fill = "celltype")
#'
#' # From data.frame
#' PlotAlluvia(seurat_obj@@meta.data, by = "group", fill = "cell_type_pred")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_area geom_col scale_x_continuous
#'   scale_y_continuous guides guide_legend theme_classic theme element_blank
#'   element_text scale_fill_manual geom_text
#' @importFrom dplyr group_by arrange mutate ungroup summarise filter
#'   group_by_at %>%
#' @importFrom rlang .data
#' @export
PlotAlluvia <- function(cellmeta, by, fill, colors = NULL,
                        palette = NULL, flow.alpha = 0.4,
                        bar.width = 0.5, bar.color = "gray50",
                        show.label = FALSE, label.size = 3,
                        label.color = "black", show.pct = FALSE,
                        y.percent = TRUE, legend.ncol = 1,
                        base.size = 15) {

  cellmeta <- .extract_cellmeta(cellmeta)
  pop.stat <- .percentage_stat(cellmeta, by, fill)

  # ---- colours ----
  if (is.null(colors) && !is.null(palette)) {
    colors <- .auto_palette(levels(factor(cellmeta[[fill]])), palette)
  }

  # ---- y variable ----
  pos <- if (y.percent) "fill" else "stack"

  alluvia <- pop.stat %>%
    dplyr::group_by(get(by)) %>%
    dplyr::arrange(get(by), get(fill), by_group = TRUE) %>%
    dplyr::mutate(y = if (y.percent) .data$proportion else .data$Freq) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x = as.numeric(as.factor(get(by)))) %>%
    dplyr::group_by(get(by), get(fill)) %>%
    dplyr::reframe(
      x = c(.data$x - 0.25, .data$x, .data$x + 0.25),
      y = .data$y
    )
  colnames(alluvia)[1:2] <- c(by, fill)

  x.labels <- unique(alluvia[[by]])

  p <- ggplot2::ggplot(
    alluvia %>% dplyr::filter(.data$x %% 1 == 0),
    ggplot2::aes(.data$x, .data$y, fill = get(fill))
  ) +
    ggplot2::geom_area(
      data = alluvia, alpha = flow.alpha, position = pos
    ) +
    ggplot2::geom_col(
      width = bar.width, color = bar.color, position = pos
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(x.labels), labels = x.labels
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend.ncol))

  # y-axis label
  if (y.percent) {
    p <- p + ggplot2::scale_y_continuous(labels = scales::label_percent())
  }

  # stratum labels
  if (show.label || show.pct) {
    bar_data <- alluvia %>% dplyr::filter(.data$x %% 1 == 0)
    bar_data <- bar_data %>%
      dplyr::group_by(.data$x) %>%
      dplyr::mutate(pct_in_group = .data$y / sum(.data$y)) %>%
      dplyr::ungroup()

    bar_data$lbl <- ifelse(bar_data$y > 0, {
      nm  <- if (show.label) as.character(bar_data[[fill]]) else ""
      pct <- if (show.pct) scales::percent(bar_data$pct_in_group, accuracy = 0.1) else ""
      trimws(paste(nm, pct))
    }, "")

    pos_label <- if (y.percent) {
      ggplot2::position_fill(vjust = 0.5)
    } else {
      ggplot2::position_stack(vjust = 0.5)
    }

    p <- p + ggplot2::geom_text(
      data = bar_data,
      ggplot2::aes(.data$x, .data$y, group = get(fill), label = .data$lbl),
      position = pos_label, inherit.aes = FALSE,
      size = label.size, color = label.color
    )
  }

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }

  p <- p +
    ggplot2::theme_classic(base_size = base.size) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.text    = ggplot2::element_text(color = "black"),
      axis.title   = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y  = ggplot2::element_blank()
    )
  p
}
