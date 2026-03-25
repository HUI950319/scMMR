# ============================================================================
# PlotPercent — Differential Abundance Beeswarm Plot
# ============================================================================

#' Differential Abundance Beeswarm Plot
#'
#' Visualise the output of \code{RankPercent()}: a miloR-style beeswarm plot
#' showing per-neighborhood logFC for each cell type.  Significant
#' neighborhoods (FDR < threshold) are coloured by their logFC value on a
#' diverging blue-white-red gradient; non-significant neighborhoods are
#' shown in grey.
#'
#' @param da_result Output list from \code{RankPercent()}.
#' @param fdr_threshold Numeric. FDR threshold for colouring significant
#'   neighborhoods (default 0.1).  Neighborhoods with
#'   \code{p_adj >= fdr_threshold} are shown in \code{na_color}.
#' @param colors Character vector of length 3: colours for the diverging
#'   gradient \code{c(low, mid, high)}.  Default
#'   \code{c("blue", "white", "red")}.
#' @param na_color Colour for non-significant neighborhoods.
#'   Default \code{"grey80"}.
#' @param show_boxplot Logical. If \code{TRUE}, overlay a transparent
#'   boxplot behind the beeswarm points to show the distribution summary
#'   for each cell type. Default \code{FALSE}.
#' @param base_size Base font size (default 12).
#' @param point_size Point size (default 1.5).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' da <- RankPercent(pred$shared_embedding, q1@@meta.data,
#'                   conditions = c("PH", "SH"))
#' PlotPercent(da)
#' PlotPercent(da, fdr_threshold = 0.05)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_boxplot
#'   scale_color_gradient2 coord_flip guides labs theme_bw theme
#'   element_text element_blank
#' @export
PlotPercent <- function(da_result,
                        fdr_threshold = 0.1,
                        colors        = c("blue", "white", "red"),
                        na_color      = "grey80",
                        show_boxplot  = FALSE,
                        base_size     = 12,
                        point_size    = 1.5) {

  if (!requireNamespace("ggbeeswarm", quietly = TRUE))
    stop("Package 'ggbeeswarm' is required for PlotPercent(). Install with:\n",
         "  install.packages('ggbeeswarm')")

  if (!is.list(da_result) || !"da_results" %in% names(da_result))
    stop("'da_result' must be the output of RankPercent().")

  df <- da_result$da_results

  # Significant neighborhoods: colour by logFC; non-significant: NA → grey
  df$logFC_color <- ifelse(df$p_adj < fdr_threshold, df$logFC, NA)

  # Order cell types by mean logFC
  ct_order <- stats::aggregate(logFC ~ cell_type_majority, data = df,
                               FUN = mean)
  ct_order <- ct_order[order(ct_order$logFC), ]
  df$cell_type_ordered <- factor(
    df$cell_type_majority,
    levels = ct_order$cell_type_majority
  )

  # Axis label from condition pair
  cond_labels <- da_result$params$conditions
  y_lab <- if (length(cond_labels) == 2) {
    paste0("log2FC (", cond_labels[2], " / ", cond_labels[1], ")")
  } else {
    "logFC"
  }

  n_sig <- sum(df$p_adj < fdr_threshold, na.rm = TRUE)
  sub_text <- sprintf("%d / %d neighborhoods significant (FDR < %.2f)",
                      n_sig, nrow(df), fdr_threshold)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x     = .data$cell_type_ordered,
      y     = .data$logFC,
      color = .data$logFC_color
    )
  )

  # Optional boxplot layer (behind beeswarm)
  if (show_boxplot) {
    p <- p +
      ggplot2::geom_boxplot(
        ggplot2::aes(group = .data$cell_type_ordered),
        outlier.shape = NA, fill = NA, color = "grey50",
        width = 0.6, linewidth = 0.3
      )
  }

  p <- p +
    ggbeeswarm::geom_quasirandom(size = point_size) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "grey40",
      linewidth = 0.5
    ) +
    ggplot2::scale_color_gradient2(
      midpoint = 0,
      low      = colors[1],
      mid      = colors[2],
      high     = colors[3],
      na.value = na_color
    ) +
    ggplot2::coord_flip() +
    ggplot2::guides(color = "none") +
    ggplot2::labs(
      x        = NULL,
      y        = y_lab,
      subtitle = sub_text
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text        = ggplot2::element_text(color = "black"),
      strip.text.y     = ggplot2::element_text(angle = 0),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}
