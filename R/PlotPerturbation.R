# ============================================================================
# PlotPerturbation — Perturbation Ranking Plot
# ============================================================================

#' Perturbation Ranking Plot
#'
#' Visualise the output of \code{RankPerturbation()}: a ranked lollipop or
#' bar chart of cell types ordered by perturbation score, with significance
#' markers.  Reuses the same visual style as \code{PlotImportance()}.
#'
#' @param rank_result Output list from \code{RankPerturbation()}.
#' @param top_k Integer. Show only top-k cell types (default \code{NULL},
#'   show all).
#' @param display Character: \code{"lollipop"} (default) or \code{"bar"}.
#' @param sig_threshold Numeric. FDR threshold for significance stars
#'   (default 0.05).
#' @param palette HCL sequential palette name (default \code{"Reds 3"}).
#' @param base_size Base font size (default 12).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' res <- RankPerturbation(pred$shared_embedding, q1@@meta.data,
#'                         conditions = c("PH", "SH"))
#' PlotPerturbation(res)
#' PlotPerturbation(res, top_k = 10, display = "bar")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_text labs
#' @export
PlotPerturbation <- function(rank_result,
                             top_k         = NULL,
                             display       = c("lollipop", "bar"),
                             sig_threshold = 0.05,
                             palette       = "Reds 3",
                             base_size     = 12) {

  display <- match.arg(display)

  if (!is.list(rank_result) || !"results" %in% names(rank_result))
    stop("'rank_result' must be the output of RankPerturbation().")

  df <- rank_result$results

  # top-k filter
  if (!is.null(top_k)) {
    df <- utils::head(df[order(-df$score), ], top_k)
  }

  # Prepare data for .plot_bar (expects: importance, var_ordered)
  plot_data <- data.frame(
    importance  = df$score,
    var_ordered = stats::reorder(df$cell_type, df$score),
    p_adj       = df$p_adj,
    stringsAsFactors = FALSE
  )

  # Significance stars
  plot_data$sig_label <- ifelse(
    plot_data$p_adj < 0.001, "***",
    ifelse(plot_data$p_adj < 0.01, "**",
           ifelse(plot_data$p_adj < sig_threshold, "*", ""))
  )

  x_lab <- paste0("Perturbation Score (",
                   rank_result$method, ")")
  subtitle <- paste(rank_result$conditions, collapse = " vs ")

  p <- .plot_bar(plot_data, display = display,
                 palette = palette, base_size = base_size,
                 x_label = x_lab)

  # Add significance stars
  p <- p +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$sig_label),
      hjust = -0.3, size = 5, color = "black"
    ) +
    ggplot2::labs(subtitle = subtitle)

  p
}
