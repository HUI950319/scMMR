# ============================================================================
# PlotMAP — UMAP Projection Plot
# ============================================================================

#' UMAP Projection Plot
#'
#' Overlay query cells (with predicted UMAP coordinates) onto the reference
#' UMAP embedding.  The reference cells are shown as a coloured scatter plot;
#' query cells are plotted as black points with red density contours.
#'
#' @param ref A data.frame (or Seurat object) providing reference UMAP
#'   coordinates and cell-type labels.
#'   \itemize{
#'     \item If a \strong{data.frame}: must contain columns specified by
#'       \code{ref_emb} and \code{color_by}.
#'     \item If a \strong{Seurat} object: columns are extracted via
#'       \code{Seurat::FetchData(ref, vars = c(ref_emb, color_by))}.
#'   }
#' @param query A data.frame with query cell metadata, e.g.
#'   \code{q1@@meta.data} after \code{AddMetaData(obj, pred$predictions)}.
#'   Must contain columns specified by \code{query_emb} (and optionally
#'   \code{facet_by}).
#' @param ref_emb Character vector of length 2: column names for reference
#'   UMAP axes.  Default \code{c("umap_1", "umap_2")}.
#' @param query_emb Character vector of length 2: column names for predicted
#'   query UMAP axes.  Default \code{c("umap_1_pred", "umap_2_pred")}.
#' @param color_by Column name in \code{ref} used for colouring reference
#'   cells.  Default \code{"cell_type"}.
#' @param facet_by Optional column name in \code{query} for faceting
#'   (e.g. \code{"group"}).  \code{NULL} = no faceting.
#' @param colors Optional named colour vector.  If \code{NULL}, default
#'   ggplot2 colours are used.
#' @param point.size Size of reference points. Default 0.2.
#' @param point.alpha Alpha of reference points. Default 0.1.
#' @param query.size Size of query points. Default 0.1.
#' @param query.alpha Alpha of query points. Default 1.
#' @param density.color Colour of density contour lines. Default \code{"red"}.
#' @param expand_mult Fractional expansion applied to both axes so that
#'   density contour lines are not clipped at panel boundaries.
#'   Default 0.05 (5\% on each side).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' ref_data <- Seurat::FetchData(seu,
#'   vars = c("umap_1", "umap_2", "cell_type"))
#' PlotMAP(ref = ref_data, query = q1@@meta.data,
#'   facet_by = "group")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d
#'   scale_color_manual guides guide_legend theme_void theme
#'   element_blank coord_equal facet_grid
#' @export
PlotMAP <- function(ref, query,
                    ref_emb       = c("umap_1", "umap_2"),
                    query_emb     = c("umap_1_pred", "umap_2_pred"),
                    color_by      = "cell_type",
                    facet_by      = NULL,
                    colors        = NULL,
                    point.size    = 0.2,
                    point.alpha   = 0.1,
                    query.size    = 0.1,
                    query.alpha   = 1,
                    density.color = "red",
                    expand_mult   = 0.05) {

  # --- resolve ref ---
  if (inherits(ref, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Package 'Seurat' is required when 'ref' is a Seurat object.")
    ref <- Seurat::FetchData(ref, vars = c(ref_emb, color_by))
  }
  stopifnot(is.data.frame(ref))
  stopifnot(all(c(ref_emb, color_by) %in% names(ref)))
  stopifnot(all(query_emb %in% names(query)))

  # --- build reference layer ---
  p <- ggplot2::ggplot(
    data = ref,
    mapping = ggplot2::aes(
      x = .data[[ref_emb[1]]],
      y = .data[[ref_emb[2]]],
      color = .data[[color_by]]
    )
  ) +
    ggplot2::geom_point(size = point.size, alpha = point.alpha)

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }

  p <- p +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = expand_mult)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = expand_mult)) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1), ncol = 2
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(10, 10, 10, 10)
    )

  # --- overlay query points + density ---
  p <- p +
    ggplot2::geom_point(
      data    = query,
      mapping = ggplot2::aes(
        x = .data[[query_emb[1]]],
        y = .data[[query_emb[2]]]
      ),
      color       = "black",
      size        = query.size,
      alpha       = query.alpha,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_density_2d(
      data    = query,
      mapping = ggplot2::aes(
        x = .data[[query_emb[1]]],
        y = .data[[query_emb[2]]]
      ),
      color       = density.color,
      contour_var = "ndensity",
      inherit.aes = FALSE
    )

  # --- facet ---
  if (!is.null(facet_by)) {
    stopifnot(facet_by %in% names(query))
    p <- p + ggplot2::facet_grid(
      stats::reformulate(facet_by)
    ) +
      ggplot2::theme(
        panel.spacing = ggplot2::unit(1.5, "cm")
      )
  }

  p
}
