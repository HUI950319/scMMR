# ============================================================================
# PlotRankScatter — Rank scatter plot for feature scoring
# ============================================================================

# --------------------------------------------------------------------------
# Internal: parse diverse inputs into a named list of named numeric vectors
# --------------------------------------------------------------------------

#' @keywords internal
.rank_scatter_parse <- function(data, cell.type) {

  ## -- Named numeric vector --
  if (is.numeric(data) && !is.null(names(data))) {
    lbl <- if (!is.null(cell.type)) cell.type[1] else ""
    return(setNames(list(data), lbl))
  }

  ## -- Matrix (RSS / pathway_scores / any feature x group matrix) --
  if (is.matrix(data)) {
    if (!is.null(cell.type)) {
      bad <- setdiff(cell.type, colnames(data))
      if (length(bad))
        stop("cell.type not found in matrix columns: ",
             paste(bad, collapse = ", "))
      if (length(cell.type) == 1L) {
        v <- data[, cell.type]
        names(v) <- rownames(data)
        return(setNames(list(v), cell.type))
      }
      return(lapply(setNames(cell.type, cell.type), function(ct) {
        v <- data[, ct]; names(v) <- rownames(data); v
      }))
    }
    # all columns
    return(lapply(
      setNames(colnames(data), colnames(data)),
      function(ct) { v <- data[, ct]; names(v) <- rownames(data); v }
    ))
  }

  ## -- data.frame --
  if (is.data.frame(data)) {
    # Case-insensitive column matching helper
    .match_col <- function(candidates, cols) {
      cols_lc <- tolower(cols)
      for (cand in candidates) {
        idx <- which(cols_lc == tolower(cand))[1]
        if (!is.na(idx)) return(cols[idx])
      }
      NA_character_
    }
    name_col  <- .match_col(
      c("gene", "name", "regulon", "pathway", "feature"), names(data)
    )
    score_col <- .match_col(
      c("importance", "score", "RSS", "NES", "value", "rho"), names(data)
    )

    if (is.na(name_col) || is.na(score_col))
      stop("data.frame must have a name column ",
           "(gene/name/regulon/pathway/feature) ",
           "and a score column (importance/score/NES/RSS/value/rho).")

    group_col <- .match_col(
      c("cell_type", "group", "cluster", "lineage", "celltype"), names(data)
    )

    if (!is.na(group_col) && is.null(cell.type)) {
      # faceted by group
      groups <- unique(data[[group_col]])
      return(lapply(setNames(groups, groups), function(g) {
        sub_df <- data[data[[group_col]] == g, ]
        v <- sub_df[[score_col]]; names(v) <- sub_df[[name_col]]; v
      }))
    }

    v   <- data[[score_col]]
    names(v) <- data[[name_col]]
    lbl <- if (!is.null(cell.type)) cell.type[1] else ""
    return(setNames(list(v), lbl))
  }

  stop("Unsupported input type. ",
       "Provide a named numeric vector, matrix, or data.frame.")
}


# --------------------------------------------------------------------------
# Internal: build rank data.frame from a named numeric vector
# --------------------------------------------------------------------------

#' @keywords internal
.rank_scatter_df <- function(scores, max.show, topn) {
  scores <- sort(scores, decreasing = TRUE)
  n      <- min(length(scores), max.show)
  data.frame(
    Rank  = seq_len(n),
    Score = scores[seq_len(n)],
    Label = names(scores)[seq_len(n)],
    IsTop = c(rep(TRUE,  min(topn, n)),
              rep(FALSE, max(0L, n - topn))),
    stringsAsFactors = FALSE
  )
}


# --------------------------------------------------------------------------
# Internal: single-panel rank scatter
# --------------------------------------------------------------------------

#' @keywords internal
.rank_scatter_single <- function(df, topn, hi_col, base_col,
                                  lab_size, pt_size, title,
                                  ylab, base_size) {

  df$pt_color <- ifelse(df$IsTop, hi_col, base_col)
  label_df    <- df[df$IsTop, ]

  ggplot2::ggplot(df, ggplot2::aes(.data$Rank, .data$Score)) +
    ggplot2::geom_point(size = pt_size, color = df$pt_color) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(.data$Rank, .data$Score, label = .data$Label),
      size = lab_size, inherit.aes = FALSE,
      max.overlaps = 20
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::ylab(ylab) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line  = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.text  = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}


# --------------------------------------------------------------------------
# Internal: faceted multi-panel rank scatter
# --------------------------------------------------------------------------

#' @keywords internal
.rank_scatter_facet <- function(df, topn, hi_col, base_col,
                                 lab_size, pt_size, title,
                                 ylab, base_size, ncol) {

  df$TopGroup <- ifelse(df$IsTop, "top", "other")
  label_df    <- df[df$IsTop, ]

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$Rank, .data$Score, color = .data$TopGroup)
  ) +
    ggplot2::geom_point(size = pt_size) +
    ggplot2::scale_color_manual(
      values = c("top" = hi_col, "other" = base_col),
      guide  = "none"
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(.data$Rank, .data$Score, label = .data$Label),
      size = lab_size, color = "black", inherit.aes = FALSE,
      max.overlaps = 20
    ) +
    ggplot2::facet_wrap(~ Group, ncol = ncol, scales = "free") +
    ggplot2::ylab(ylab) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line  = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.text  = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}


#' Rank Scatter Plot for Feature Scoring
#'
#' A general-purpose scatter plot that ranks features (genes, regulons,
#' pathways, etc.) by their scores and highlights the top-ranked ones
#' with coloured points and text labels.  Inspired by the SCENIC
#' Regulon Specificity Score (RSS) visualisation.
#'
#' @param data Input data in one of the following formats:
#'   \describe{
#'     \item{Named numeric vector}{Values are scores; names are feature
#'       labels.  Produces a single-panel plot.}
#'     \item{Matrix}{Rows = features, columns = groups (e.g. cell types).
#'       Use \code{cell.type} to select one or several columns.  If
#'       \code{cell.type = NULL} all columns are shown as a faceted plot.
#'       Typical inputs: RSS matrix from SCENIC, or
#'       \code{DNN_predict()$pathway_scores}.}
#'     \item{data.frame}{Must contain a name column
#'       (\code{gene}, \code{name}, \code{regulon}, \code{pathway}, or
#'       \code{feature}) and a score column
#'       (\code{importance}, \code{score}, \code{NES}, or \code{value}).
#'       If a grouping column (\code{cell_type}, \code{group}, or
#'       \code{cluster}) is present the plot is faceted.}
#'   }
#' @param cell.type Character vector.
#'   For matrix input: column name(s) to plot. \code{NULL} (default) =
#'   plot all columns.
#'   For vector/data.frame input: used only as the plot title.
#' @param topn Integer.
#'   Number of top-ranked features to highlight (default 5).
#' @param max.show Integer.
#'   Maximum number of features to display (default 200).
#' @param highlight.color Colour for the top-ranked points
#'   (default \code{"#007D9B"}).
#' @param base.color Colour for the remaining points
#'   (default \code{"#BECEE3"}).
#' @param label.size Numeric. Text label size (default 4).
#' @param point.size Numeric. Point size (default 3).
#' @param title Character. Plot title.
#'   If \code{NULL} (default) the group name is used.
#' @param ylab Character. Y-axis label (default \code{"Score"}).
#' @param base_size Numeric. Base font size (default 12).
#' @param ncol Integer. Number of columns for faceted layout (default 3).
#' @param clean.names Logical. If \code{TRUE} (default), strip common
#'   prefixes (\code{HALLMARK_}, \code{KEGG_}, etc.), SCENIC suffixes
#'   (\code{(+)}, \code{(-)}), and replace underscores with spaces.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # 1) SCENIC RSS matrix - single cell type
#' PlotRankScatter(rssMat, cell.type = "Parathyroid", topn = 10)
#'
#' # 2) SCENIC RSS matrix - multiple cell types
#' PlotRankScatter(rssMat, cell.type = c("Parathyroid", "Stromal"),
#'                 topn = 5, ncol = 2)
#'
#' # 3) Named numeric vector
#' scores <- setNames(runif(100), paste0("Gene", 1:100))
#' PlotRankScatter(scores, topn = 10, ylab = "Activity")
#'
#' # 4) DNN_predict importance (data.frame)
#' PlotRankScatter(pred$imp_global, topn = 15)
#'
#' # 5) Per-class importance - faceted
#' PlotRankScatter(pred$imp_per_class, topn = 10, ncol = 4)
#' }
#'
#' @export
PlotRankScatter <- function(data,
                            cell.type       = NULL,
                            topn            = 5L,
                            max.show        = 200L,
                            highlight.color = "#007D9B",
                            base.color      = "#BECEE3",
                            label.size      = 4,
                            point.size      = 3,
                            title           = NULL,
                            ylab            = "Score",
                            base_size       = 12,
                            ncol            = 3L,
                            clean.names     = TRUE) {

  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Package 'ggrepel' is required. ",
         "Install with:  install.packages('ggrepel')")

  # -- Parse input into a named list of named numeric vectors --
  score_list <- .rank_scatter_parse(data, cell.type)

  # -- Optionally clean feature names --
  if (clean.names) {
    score_list <- lapply(score_list, function(v) {
      nm <- names(v)
      # SCENIC suffixes
      nm <- sub("\\(\\+\\)$", "", nm)
      nm <- sub("\\(\\-\\)$", "", nm)
      nm <- trimws(nm)
      # Database prefixes
      nm <- sub("^HALLMARK_", "", nm)
      nm <- sub("^KEGG_",     "", nm)
      nm <- sub("^REACTOME_", "", nm)
      nm <- sub("^GOBP_",     "", nm)
      nm <- sub("^PROGENY_",  "", nm)
      nm <- sub("^WP_",       "", nm)
      nm <- sub("^PID_",      "", nm)
      nm <- sub("^BIOCARTA_", "", nm)
      # Underscores -> spaces
      nm <- gsub("_", " ", nm)
      names(v) <- nm
      v
    })
  }

  # -- Single panel vs multi-panel --
  if (length(score_list) == 1L) {
    df  <- .rank_scatter_df(score_list[[1]], max.show, topn)
    ttl <- if (!is.null(title)) title else names(score_list)
    p   <- .rank_scatter_single(df, topn, highlight.color, base.color,
                                label.size, point.size, ttl, ylab, base_size)
  } else {
    all_df <- lapply(names(score_list), function(nm) {
      d <- .rank_scatter_df(score_list[[nm]], max.show, topn)
      d$Group <- nm
      d
    })
    plot_df <- do.call(rbind, all_df)
    p <- .rank_scatter_facet(plot_df, topn, highlight.color, base.color,
                              label.size, point.size, title, ylab,
                              base_size, ncol)
  }

  p
}
