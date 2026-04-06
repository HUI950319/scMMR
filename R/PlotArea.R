# ============================================================================
# PlotArea — Alluvial (stacked area + bar) plot
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
#' @param palette Character. Colour palette name passed to
#'   \code{\link{palette_colors}}. Default \code{"lancet"}.
#' @param palcolor Character vector. Custom colours overriding
#'   \code{palette}. Default \code{NULL}.
#' @param flow.alpha Numeric (0-1), transparency of the area flow layer.
#'   Default 0.4.
#' @param bar.width Numeric (0-1) specifying bar width. Default 0.5.
#' @param bar.color Bar border colour. Default \code{"gray50"}.
#' @param show.label Logical, whether to show cell type labels inside bars.
#'   Default \code{FALSE}.
#' @param label.size Numeric, font size for bar labels. Default 3.
#' @param label.color Label text colour. Default \code{"black"}.
#' @param label.box Logical. If \code{TRUE}, draw labels with white
#'   background and border (like \code{geom_label}). Default \code{FALSE}.
#' @param label.fill Background fill for boxed labels. Default \code{"white"}.
#' @param label.alpha Transparency of boxed label background. Default \code{1}.
#' @param show.pct Logical, whether to append within-group percentage to
#'   labels. Default \code{FALSE}. Only used when \code{show.label = TRUE}.
#' @param label.style Character. Label format style:
#'   \itemize{
#'     \item \code{"default"} — uses \code{show.label}/\code{show.pct} flags.
#'     \item \code{"all"} — name + count + percentage, e.g. \code{"T cells (50, 8.1\%)"}.
#'     \item \code{"n"} — name + count, e.g. \code{"T cells (50)"}.
#'     \item \code{"pct"} — name + percentage, e.g. \code{"T cells 8.1\%"}.
#'     \item \code{"name"} — name only.
#'   }
#' @param label.min.pct Numeric (0-1). Minimum within-group proportion to
#'   display a label. Cell types below this threshold are unlabelled.
#'   Default 0 (show all). E.g. \code{0.05} hides labels for cell types
#'   < 5\%.
#' @param y.type Character. Y-axis scale type:
#'   \itemize{
#'     \item \code{"percent"} (default) — all bars equal height, y-axis
#'       shows percentages.
#'     \item \code{"count"} — bars reflect raw cell counts.
#'   }
#' @param legend.ncol Integer, number of legend columns. Default 1.
#' @param base.size Numeric, base font size for the theme. Default 15.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' library(ToyData)
#'
#' # --- Basic usage (defaults: by="group", fill="celltype", lancet palette) ---
#' PlotArea(toy_test)
#'
#' # --- From data.frame ---
#' PlotArea(toy_test@@meta.data, by = "group", fill = "celltype")
#'
#' # --- Label styles ---
#' # Name + percentage (hide small cell types < 5%)
#' PlotArea(toy_test, label.style = "pct", label.min.pct = 0.05)
#'
#' # Name + count + percentage
#' PlotArea(toy_test, label.style = "all", label.min.pct = 0.05)
#'
#' # Name only
#' PlotArea(toy_test, label.style = "name", label.min.pct = 0.03)
#'
#' # --- Boxed labels (white background + border) ---
#' PlotArea(toy_test, label.style = "pct", label.min.pct = 0.05,
#'             label.box = TRUE)
#'
#' # Semi-transparent boxed labels
#' PlotArea(toy_test, label.style = "pct", label.min.pct = 0.05,
#'             label.box = TRUE, label.alpha = 0.7)
#'
#' # --- Custom colours ---
#' # Use a different palette
#' PlotArea(toy_test, palette = "Paired")
#'
#' # Custom colour vector (parathyroid cell-type colours)
#' PlotArea(toy_test, palcolor = UtilsR::pal_paraSC)
#'
#' # --- Appearance ---
#' # Adjust flow transparency and bar width
#' PlotArea(toy_test, flow.alpha = 0.2, bar.width = 0.3)
#'
#' # Raw counts instead of percentages
#' PlotArea(toy_test, y.type = "count", label.style = "n",
#'             label.min.pct = 0)
#'
#' # Multi-column legend
#' PlotArea(toy_test, legend.ncol = 2)
#'
#' # --- Legacy style (show.label + show.pct flags) ---
#' PlotArea(toy_test, show.label = TRUE, show.pct = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_area geom_col scale_x_continuous
#'   scale_y_continuous guides guide_legend theme_classic theme element_blank
#'   element_text scale_fill_manual geom_text
#' @importFrom dplyr group_by arrange mutate ungroup summarise filter
#'   group_by_at %>%
#' @importFrom rlang .data
#' @export
PlotArea <- function(cellmeta, by = "group", fill = "celltype",
                        palette = "lancet", palcolor = NULL,
                        flow.alpha = 0.4,
                        bar.width = 0.5, bar.color = "gray50",
                        show.label = FALSE, label.size = 3,
                        label.color = "black", label.box = FALSE,
                        label.fill = "white", label.alpha = 1,
                        show.pct = FALSE,
                        label.style = c("default", "all", "n", "pct", "name"),
                        label.min.pct = 0,
                        y.type = c("percent", "count"), legend.ncol = 1,
                        base.size = 15) {

  label.style <- match.arg(label.style)
  y.type <- match.arg(y.type)
  cellmeta <- .extract_cellmeta(cellmeta)
  pop.stat <- .percentage_stat(cellmeta, by, fill)

  # ---- colours ----
  colors <- palette_colors(levels(factor(cellmeta[[fill]])),
                           palette = palette, palcolor = palcolor)

  # ---- y variable ----
  pos <- if (y.type == "percent") "fill" else "stack"

  alluvia <- pop.stat %>%
    dplyr::group_by(get(by)) %>%
    dplyr::arrange(get(by), get(fill), by_group = TRUE) %>%
    dplyr::mutate(y = if (y.type == "percent") .data$proportion else .data$Freq) %>%
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
  if (y.type == "percent") {
    p <- p + ggplot2::scale_y_continuous(labels = scales::label_percent())
  }

  # stratum labels
  use_label <- if (label.style != "default") TRUE else (show.label || show.pct)
  if (use_label) {
    bar_data <- alluvia %>% dplyr::filter(.data$x %% 1 == 0)
    bar_data <- bar_data %>%
      dplyr::group_by(.data$x) %>%
      dplyr::mutate(
        pct_in_group = .data$y / sum(.data$y),
        grp_total = sum(.data$y)
      ) %>%
      dplyr::ungroup()

    # Compute count per cell type per group
    bar_data$.n <- round(bar_data$pct_in_group * nrow(cellmeta) /
                           length(unique(bar_data$x)))

    bar_data$lbl <- ifelse(
      bar_data$y > 0 & bar_data$pct_in_group >= label.min.pct, {
        nm  <- as.character(bar_data[[fill]])
        n   <- bar_data$.n
        pct <- scales::percent(bar_data$pct_in_group, accuracy = 0.1)
        switch(label.style,
          "all"  = sprintf("%s (%s, %s)", nm, n, pct),
          "n"    = sprintf("%s (%s)", nm, n),
          "pct"  = sprintf("%s %s", nm, pct),
          "name" = nm,
          "default" = {
            nm_part  <- if (show.label) nm else ""
            pct_part <- if (show.pct) pct else ""
            trimws(paste(nm_part, pct_part))
          }
        )
      }, NA_character_)

    pos_label <- if (y.type == "percent") {
      ggplot2::position_fill(vjust = 0.5)
    } else {
      ggplot2::position_stack(vjust = 0.5)
    }

    if (label.box) {
      p <- p + ggplot2::geom_label(
        data = bar_data,
        ggplot2::aes(.data$x, .data$y, group = get(fill), label = .data$lbl),
        position = pos_label, inherit.aes = FALSE,
        size = label.size, color = label.color,
        fill = label.fill, alpha = label.alpha
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = bar_data,
        ggplot2::aes(.data$x, .data$y, group = get(fill), label = .data$lbl),
        position = pos_label, inherit.aes = FALSE,
        size = label.size, color = label.color
      )
    }
  }

  p <- p + ggplot2::scale_fill_manual(values = colors)

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
