# ============================================================================
# PlotRoe — Observed / Expected Ratio Heatmap
# ============================================================================

#' O/E Ratio Heatmap for Cell Type Composition
#'
#' Create a heatmap displaying the Observed / Expected ratio of cell type
#' composition across groups, highlighting group-specific enrichment.
#' Optionally computes per-cell p-values via chi-squared standardised
#' residuals or Fisher's exact test, and marks significant cells with stars.
#'
#' @param cellmeta A Seurat object or data.frame containing cell metadata
#'   (e.g. \code{seurat_obj} or \code{seurat_obj@@meta.data}).
#' @param by Character string specifying the grouping variable
#'   (e.g. \code{"group"}, \code{"sample"}).
#' @param fill Character string specifying the cell-type variable
#'   (e.g. \code{"cell_type_pred"}).
#' @param palette RColorBrewer palette name. Default \code{"Blues"}.
#' @param column_names_rot Rotation angle (0-360) for column labels.
#'   Default 0.
#' @param display Character. What to show inside each cell:
#'   \describe{
#'     \item{\code{"value"}}{(Default) Numeric O/E ratio (e.g. \code{1.58}).}
#'     \item{\code{"symbol"}}{Grade symbols based on O/E level:
#'       \code{+++} (>1), \code{++} (0.8--1), \code{+} (0.2--0.8),
#'       \code{+/-} (0--0.2), \code{-} (0).}
#'     \item{\code{"both"}}{Both numeric value and grade symbol
#'       (value on top, symbol below).}
#'   }
#' @param show.pvalue Logical. Whether to compute and display significance
#'   stars on the heatmap.  Default \code{FALSE}.  When combined with
#'   \code{display = "symbol"}, stars are appended to the symbol
#'   (e.g. \code{+++ ***}).
#' @param p.method Character. Method for per-cell p-values:
#'   \code{"chisq"} (default, chi-squared standardised residuals) or
#'   \code{"fisher"} (Fisher's exact test per 2x2 cell).
#' @param p.adjust.method P-value adjustment method passed to
#'   \code{p.adjust()}.  Default \code{"BH"}.
#' @param sig.threshold Numeric vector defining significance cut-offs for
#'   star labels.  Default \code{c(0.001, 0.01, 0.05)} (\code{***},
#'   \code{**}, \code{*}).
#' @param return.data Logical. If \code{TRUE}, return a list with the O/E
#'   matrix, p-value matrix and the heatmap object.  Default \code{FALSE}.
#' @param ... Additional arguments passed to
#'   \code{ComplexHeatmap::Heatmap()}.
#'
#' @return When \code{return.data = FALSE} (default): a ComplexHeatmap
#'   Heatmap object.
#'
#'   When \code{return.data = TRUE}: a list with components \code{oe}
#'   (O/E matrix), \code{pvalue} (adjusted p-value matrix or \code{NULL}),
#'   and \code{heatmap}.
#'
#' @examples
#' \dontrun{
#' # Basic O/E heatmap (numeric values)
#' PlotRoe(seu, by = "group", fill = "cell_type_pred")
#'
#' # Grade symbols (+++ / ++ / + / +/- / -)
#' PlotRoe(seu, by = "group", fill = "cell_type_pred",
#'         display = "symbol")
#'
#' # Both value and symbol
#' PlotRoe(seu, by = "group", fill = "cell_type_pred",
#'         display = "both", show.pvalue = TRUE)
#'
#' # Return data for downstream use
#' res <- PlotRoe(seu, by = "group", fill = "cell_type_pred",
#'                show.pvalue = TRUE, return.data = TRUE)
#' res$oe       # O/E ratio matrix
#' res$pvalue   # adjusted p-value matrix
#' }
#'
#' @export
PlotRoe <- function(cellmeta, by, fill,
                    palette          = "Blues",
                    column_names_rot = 0,
                    display          = c("value", "symbol", "both"),
                    show.pvalue      = FALSE,
                    p.method         = c("chisq", "fisher"),
                    p.adjust.method  = "BH",
                    sig.threshold    = c(0.001, 0.01, 0.05),
                    return.data      = FALSE,
                    ...) {

  display  <- match.arg(display)
  p.method <- match.arg(p.method)

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Package 'ComplexHeatmap' is required. Install with:\n",
         "  BiocManager::install('ComplexHeatmap')")
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    stop("Package 'RColorBrewer' is required. Install with:\n",
         "  install.packages('RColorBrewer')")

  cellmeta <- .extract_cellmeta(cellmeta)
  if (!(by %in% names(cellmeta)) || !(fill %in% names(cellmeta)))
    stop("'by' and 'fill' must be columns in cellmeta.")

  ## O/E ratio ---------------------------------------------------------------
  tab  <- table(cellmeta[[by]], cellmeta[[fill]])
  obs  <- apply(tab, 1, function(xx) xx / sum(xx)) %>% t()
  exp  <- colSums(tab) / sum(tab)
  exp_mat <- matrix(rep(exp, nrow(tab)), ncol = nrow(tab)) %>% t()
  oe   <- t(obs / exp_mat)
  oe[is.na(oe)] <- 1
  oe   <- as.matrix(oe)

  ## P-values ----------------------------------------------------------------
  p_mat <- NULL
  if (show.pvalue) {
    if (p.method == "chisq") {
      chi  <- suppressWarnings(stats::chisq.test(tab))
      sres <- as.matrix(chi[["stdres"]])
      sres[is.nan(sres)] <- 0
      raw_p <- 2 * stats::pnorm(abs(sres), lower.tail = FALSE)
    } else {
      nr <- nrow(tab); nc <- ncol(tab)
      raw_p <- matrix(1, nr, nc, dimnames = dimnames(tab))
      for (i in seq_len(nr)) for (j in seq_len(nc)) {
        a <- tab[i, j]; b <- sum(tab[i, ]) - a
        cc <- sum(tab[, j]) - a; d <- sum(tab) - a - b - cc
        raw_p[i, j] <- stats::fisher.test(matrix(c(a, b, cc, d), 2))$p.value
      }
    }
    adj_vec <- stats::p.adjust(as.vector(raw_p), method = p.adjust.method)
    p_mat <- t(matrix(adj_vec, nrow(raw_p), ncol(raw_p),
                       dimnames = dimnames(raw_p)))
    p_mat <- as.matrix(p_mat)
  }

  ## Ordering ----------------------------------------------------------------
  col_lvl <- levels(cellmeta[[by]])
  row_lvl <- rev(levels(cellmeta[[fill]]))

  if (is.null(row_lvl)) {
    ri <- order(oe[, 1], decreasing = TRUE)
  } else {
    ri <- match(row_lvl, rownames(oe))
  }
  oe <- oe[ri, , drop = FALSE]
  if (!is.null(p_mat)) p_mat <- p_mat[ri, , drop = FALSE]

  if (is.null(col_lvl)) {
    ci <- order(oe[1, ], decreasing = TRUE)
  } else {
    ci <- match(col_lvl, colnames(oe))
  }
  oe <- oe[, ci, drop = FALSE]
  if (!is.null(p_mat)) p_mat <- p_mat[, ci, drop = FALSE]

  ## Helper: O/E → grade symbol ----------------------------------------------
  .oe_symbol <- function(v) {
    if (v > 1)                  return("+++")
    if (v > 0.8 && v <= 1)     return("++")
    if (v > 0.2 && v <= 0.8)   return("+")
    if (v > 0   && v <= 0.2)   return("+/\u2212")
    return("\u2212")
  }

  ## Helper: p-value → stars -------------------------------------------------
  sig.threshold <- sort(sig.threshold)
  .star <- function(p) {
    if (is.na(p) || p >= max(sig.threshold)) return("")
    strrep("*", sum(p < sig.threshold))
  }

  ## Plot --------------------------------------------------------------------
  q95 <- stats::quantile(oe, .95, na.rm = TRUE)
  q90 <- round(stats::quantile(oe, .9, na.rm = TRUE), 1)
  oe_cap <- ifelse(oe > q95, q95, oe)

  hm_col <- tryCatch(
    RColorBrewer::brewer.pal(7, palette),
    error = function(e) {
      message(e$message, " -- using Blues"); RColorBrewer::brewer.pal(7, "Blues")
    }
  )

  ht <- ComplexHeatmap::Heatmap(
    oe_cap,
    name             = "Ratio (o/e)",
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    col              = hm_col,
    column_names_rot = column_names_rot,
    cell_fun = function(j, i, x, y, width, height, fill) {
      val     <- oe[i, j]
      txt_col <- if (val > q90) "white" else "black"
      txt_face <- if (val > q90) "bold" else "plain"
      star    <- if (!is.null(p_mat)) .star(p_mat[i, j]) else ""
      sym     <- .oe_symbol(val)

      if (display == "value") {
        ## Line 1: O/E value
        grid::grid.text(
          sprintf("%.2f", val), x, y,
          gp = grid::gpar(fontsize = 10, col = txt_col, fontface = txt_face)
        )
        ## Line 2: p-value stars (if any)
        if (nchar(star) > 0) {
          grid::grid.text(
            star, x, y - grid::unit(0.35, "lines"),
            gp = grid::gpar(fontsize = 8, col = txt_col)
          )
        }

      } else if (display == "symbol") {
        ## Single line: symbol + stars
        lbl <- if (nchar(star) > 0) paste(sym, star) else sym
        grid::grid.text(
          lbl, x, y,
          gp = grid::gpar(fontsize = 11, col = txt_col, fontface = txt_face)
        )

      } else {
        ## "both": line 1 = value, line 2 = symbol + stars
        grid::grid.text(
          sprintf("%.2f", val), x, y + grid::unit(0.55, "lines"),
          gp = grid::gpar(fontsize = 10, col = txt_col, fontface = txt_face)
        )
        lbl <- if (nchar(star) > 0) paste(sym, star) else sym
        grid::grid.text(
          lbl, x, y - grid::unit(0.45, "lines"),
          gp = grid::gpar(fontsize = 9, col = txt_col)
        )
      }
    },
    ...
  )

  if (return.data) return(list(oe = oe, pvalue = p_mat, heatmap = ht))
  ht
}
