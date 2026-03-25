# ============================================================================
# PlotSankey — Sankey Plot for Cell Population Composition
# ============================================================================

#' Sankey Plot for Cell Population Composition
#'
#' Create a sankey (alluvial flow) plot showing cell type composition changes
#' across conditions/groups, based on the ggalluvial package. Each column
#' represents a group, strata represent cell types, and flows connect the
#' same cell type across adjacent groups. Parameters mirror
#' \code{\link{PlotAlluvia}}, with an additional \code{curve.type} argument.
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
#' @param flow.alpha Numeric (0-1), transparency of alluvial flows.
#'   Default 0.4.
#' @param bar.width Numeric (0-1), width of stratum bars. Default 0.5.
#' @param bar.color Stratum bar border colour. Default \code{"gray50"}.
#' @param show.label Logical, whether to show cell type labels inside
#'   strata. Default \code{FALSE}.
#' @param label.size Numeric, font size for stratum labels. Default 3.
#' @param label.color Stratum label colour. Default \code{"black"}.
#' @param show.pct Logical, whether to append within-group percentage to
#'   labels. Default \code{FALSE}. Only used when \code{show.label = TRUE}.
#' @param y.percent Logical. If \code{TRUE} (default), y-axis shows
#'   percentages (all bars equal height). If \code{FALSE}, y-axis shows
#'   raw cell counts.
#' @param legend.ncol Integer, number of legend columns. Default 1.
#' @param base.size Numeric, base font size for the theme. Default 15.
#' @param curve.type Character, flow curve type. Default \code{"sigmoid"}.
#'   Other options: \code{"linear"}, \code{"cubic"}, etc.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # From Seurat object directly
#' PlotSankey(seurat_obj, by = "group", fill = "celltype")
#'
#' # From data.frame
#' PlotSankey(seurat_obj@@meta.data, by = "group", fill = "cell_type_pred")
#' }
#'
#' @export
PlotSankey <- function(cellmeta, by, fill, colors = NULL,
                       palette = NULL, flow.alpha = 0.4,
                       bar.width = 0.5, bar.color = "gray50",
                       show.label = FALSE, label.size = 3,
                       label.color = "black", show.pct = FALSE,
                       y.percent = TRUE, legend.ncol = 1,
                       base.size = 15, curve.type = "sigmoid") {

  if (!requireNamespace("ggalluvial", quietly = TRUE))
    stop("Package 'ggalluvial' is required. Install with: ",
         "install.packages('ggalluvial')")

  cellmeta <- .extract_cellmeta(cellmeta)
  if (!is.character(by) || !is.character(fill))
    stop("by and fill must be character vectors.")
  if (!(by %in% names(cellmeta)) || !(fill %in% names(cellmeta)))
    stop("by and fill must be columns in cellmeta.")

  # ---- colours ----
  if (is.null(colors) && !is.null(palette)) {
    colors <- .auto_palette(levels(factor(cellmeta[[fill]])), palette)
  }

  # ---- frequency table ----
  stat <- as.data.frame(
    table(cellmeta[[by]], cellmeta[[fill]]),
    stringsAsFactors = FALSE
  )
  colnames(stat) <- c("Group", "Celltype", "Freq")

  if (is.factor(cellmeta[[by]])) {
    stat$Group <- factor(stat$Group, levels = levels(cellmeta[[by]]))
  } else {
    stat$Group <- factor(stat$Group)
  }
  if (is.factor(cellmeta[[fill]])) {
    stat$Celltype <- factor(stat$Celltype, levels = levels(cellmeta[[fill]]))
  } else {
    stat$Celltype <- factor(stat$Celltype)
  }

  # ---- drop empty levels & compute within-group proportions ----
  stat <- stat[stat$Freq > 0 | ave(stat$Freq, stat$Celltype, FUN = sum) > 0, ]
  stat <- droplevels(stat)

  stat <- stat %>%
    dplyr::group_by(.data$Group) %>%
    dplyr::mutate(
      group_total = sum(.data$Freq),
      pct = ifelse(.data$group_total > 0,
                   .data$Freq / .data$group_total, 0)
    ) %>%
    dplyr::ungroup()

  stat$y_val <- if (y.percent) stat$pct else stat$Freq

  # ---- labels ----
  stat$label <- ifelse(stat$Freq > 0, {
    nm  <- if (show.label) as.character(stat$Celltype) else ""
    pct <- if (show.pct) scales::percent(stat$pct, accuracy = 0.1) else ""
    trimws(paste(nm, pct))
  }, "")

  # ---- build plot ----
  p <- ggplot2::ggplot(
    stat,
    ggplot2::aes(
      x = .data$Group, y = .data$y_val,
      alluvium = .data$Celltype, stratum = .data$Celltype,
      fill = .data$Celltype
    )
  ) +
    ggalluvial::geom_alluvium(
      alpha = flow.alpha, width = bar.width, curve_type = curve.type
    ) +
    ggalluvial::geom_stratum(
      width = bar.width, color = bar.color, linewidth = 0.3
    )

  if (show.label || show.pct) {
    p <- p + ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      ggplot2::aes(label = .data$label),
      size = label.size, color = label.color
    )
  }

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }

  if (y.percent) {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::label_percent())
  }

  p <- p +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend.ncol)) +
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
