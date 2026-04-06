# ============================================================================
# PlotAlluvia — Sankey Plot for Cell Population Composition
# ============================================================================

#' Sankey Plot for Cell Population Composition
#'
#' Create a sankey (alluvial flow) plot showing cell type composition changes
#' across conditions/groups, based on the ggalluvial package. Each column
#' represents a group, strata represent cell types, and flows connect the
#' same cell type across adjacent groups. Parameters mirror
#' \code{PlotArea} (internal), with an additional \code{curve.type} argument.
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
#' @param flow.alpha Numeric (0-1), transparency of alluvial flows.
#'   Default 0.4.
#' @param bar.width Numeric (0-1), width of stratum bars. Default 0.5.
#' @param bar.color Stratum bar border colour. Default \code{"gray50"}.
#' @param show.label Logical, whether to show cell type labels inside
#'   strata. Default \code{FALSE}.
#' @param label.size Numeric, font size for stratum labels. Default 3.
#' @param label.color Stratum label colour. Default \code{"black"}.
#' @param label.box Logical. If \code{TRUE}, draw labels with white
#'   background and border (like \code{geom_label}). Default \code{FALSE}
#'   (plain text).
#' @param label.fill Background fill for boxed labels. Default \code{"white"}.
#' @param label.alpha Transparency of boxed label background.
#'   Default \code{1}.
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
#' @param curve.type Character, flow curve type. Default \code{"sigmoid"}.
#'   Other options: \code{"linear"}, \code{"cubic"}, etc.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' library(ToyData)
#'
#' # --- Basic usage (defaults: by="group", fill="celltype", lancet palette) ---
#' PlotAlluvia(toy_test)
#'
#' # --- From data.frame ---
#' PlotAlluvia(toy_test@@meta.data, by = "group", fill = "celltype")
#'
#' # --- Label styles ---
#' # Name + percentage (hide small cell types < 5%)
#' PlotAlluvia(toy_test, label.style = "pct", label.min.pct = 0.05)
#'
#' # Name + count + percentage
#' PlotAlluvia(toy_test, label.style = "all", label.min.pct = 0.05)
#'
#' # Name only
#' PlotAlluvia(toy_test, label.style = "name", label.min.pct = 0.03)
#'
#' # --- Boxed labels (white background + border) ---
#' PlotAlluvia(toy_test, label.style = "pct", label.min.pct = 0.05,
#'            label.box = TRUE)
#'
#' # Semi-transparent boxed labels
#' PlotAlluvia(toy_test, label.style = "pct", label.min.pct = 0.05,
#'            label.box = TRUE, label.alpha = 0.7)
#'
#' # --- Custom colours ---
#' # Use a different palette
#' PlotAlluvia(toy_test, palette = "Paired")
#'
#' # Custom colour vector
#' PlotAlluvia(toy_test, palcolor = UtilsR::pal_paraSC)
#'
#' # --- Appearance ---
#' # Sigmoid flow (default) vs linear
#' PlotAlluvia(toy_test, curve.type = "linear")
#'
#' # Adjust flow transparency and bar width
#' PlotAlluvia(toy_test, flow.alpha = 0.2, bar.width = 0.3)
#'
#' # Raw counts instead of percentages
#' PlotAlluvia(toy_test, y.type = "count", label.style = "n",
#'            label.min.pct = 0)
#'
#' # Multi-column legend
#' PlotAlluvia(toy_test, legend.ncol = 2)
#' }
#'
#' @export
PlotAlluvia <- function(cellmeta, by = "group", fill = "celltype",
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
                       base.size = 15, curve.type = "sigmoid") {

  label.style <- match.arg(label.style)
  y.type <- match.arg(y.type)

  if (!requireNamespace("ggalluvial", quietly = TRUE))
    stop("Package 'ggalluvial' is required. Install with: ",
         "install.packages('ggalluvial')")

  cellmeta <- .extract_cellmeta(cellmeta)
  if (!is.character(by) || !is.character(fill))
    stop("by and fill must be character vectors.")
  if (!(by %in% names(cellmeta)) || !(fill %in% names(cellmeta)))
    stop("by and fill must be columns in cellmeta.")

  # ---- colours ----
  colors <- palette_colors(levels(factor(cellmeta[[fill]])),
                           palette = palette, palcolor = palcolor)

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

  stat$y_val <- if (y.type == "percent") stat$pct else stat$Freq


  # ---- labels (NA for filtered-out, so geom_label/geom_text skip them) ----
  stat$label <- ifelse(stat$Freq > 0 & stat$pct >= label.min.pct, {
    nm  <- as.character(stat$Celltype)
    n   <- stat$Freq
    pct <- scales::percent(stat$pct, accuracy = 0.1)
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

  use_label <- if (label.style != "default") TRUE else (show.label || show.pct)
  if (use_label) {
    if (label.box) {
      p <- p + ggplot2::geom_label(
        stat = ggalluvial::StatStratum,
        ggplot2::aes(label = .data$label),
        size = label.size, color = label.color,
        fill = label.fill, alpha = label.alpha
      )
    } else {
      p <- p + ggplot2::geom_text(
        stat = ggalluvial::StatStratum,
        ggplot2::aes(label = .data$label),
        size = label.size, color = label.color
      )
    }
  }

  p <- p + ggplot2::scale_fill_manual(values = colors)

  if (y.type == "percent") {
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
