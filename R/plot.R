# ============================================================================
# Visualization functions for scMMR prediction results
# ============================================================================

# --------------------------------------------------------------------------
# Internal helper: calculate cell population percentage
# --------------------------------------------------------------------------

.percentage_stat <- function(cellmeta, by, fill) {
  if (!is.data.frame(cellmeta))
    stop("cellmeta must be a data.frame or tibble.")
  if (!is.character(by) || !is.character(fill))
    stop("by and fill must be character vectors.")
  if (!(by %in% names(cellmeta)) || !(fill %in% names(cellmeta)))
    stop("by and fill must be columns in cellmeta.")

  data.stat <- as.data.frame(table(cellmeta[[by]], cellmeta[[fill]]))
  colnames(data.stat)[1:2] <- c(by, fill)
  data.stat <- data.stat %>%
    dplyr::group_by_at(by) %>%
    dplyr::mutate(margin.freq = sum(.data$Freq)) %>%
    dplyr::mutate(proportion = .data$Freq / .data$margin.freq)
  data.stat
}

# ============================================================================
# 1. PlotAlluvia
# ============================================================================

#' Alluvial Plot for Cell Population Composition
#'
#' Create an alluvial (stacked area + bar) plot showing cell type composition
#' changes across conditions/groups.
#'
#' @param cellmeta A data.frame containing cell metadata (e.g.
#'   \code{seurat_obj@@meta.data}).
#' @param by Character string specifying the grouping variable
#'   (e.g. \code{"group"}, \code{"sample"}).
#' @param fill Character string specifying the fill variable
#'   (e.g. \code{"cell_type_pred"}).
#' @param colors Optional character vector of colours. If \code{NULL},
#'   default ggplot2 colours are used.
#' @param bar.width Numeric (0-1) specifying bar width. Default 0.5.
#' @param legend.ncol Integer, number of legend columns. Default 1.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' pred <- DNN_predict(query, model_path)
#' q1 <- Seurat::AddMetaData(toy_test, pred$predictions)
#' PlotAlluvia(q1@@meta.data, by = "group", fill = "cell_type_pred")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_area geom_col scale_x_continuous
#'   scale_y_continuous guides guide_legend theme_classic theme element_blank
#'   element_text scale_fill_manual
#' @importFrom dplyr group_by arrange mutate ungroup summarise filter
#'   group_by_at %>%
#' @importFrom rlang .data
#' @export
PlotAlluvia <- function(cellmeta, by, fill, colors = NULL,
                        bar.width = 0.5, legend.ncol = 1) {

  pop.stat <- .percentage_stat(cellmeta, by, fill)

  alluvia <- pop.stat %>%
    dplyr::group_by(get(by)) %>%
    dplyr::arrange(get(by), get(fill), by_group = TRUE) %>%
    dplyr::mutate(y = .data$proportion) %>%
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
    ggplot2::geom_area(data = alluvia, alpha = 0.4, position = "fill") +
    ggplot2::geom_col(width = bar.width, color = "gray50", position = "fill") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(x.labels), labels = x.labels
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend.ncol))

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }

  p <- p +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.text    = ggplot2::element_text(color = "black"),
      axis.title   = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y  = ggplot2::element_blank()
    )
  p
}


# ============================================================================
# 7. PlotRankScatter
# ============================================================================

# --------------------------------------------------------------------------
# Internal: parse diverse inputs into a named list of named numeric vectors
# --------------------------------------------------------------------------

.rank_scatter_parse <- function(data, cell.type) {

  ## ── Named numeric vector ──────────────────────────────────────────────────
  if (is.numeric(data) && !is.null(names(data))) {
    lbl <- if (!is.null(cell.type)) cell.type[1] else ""
    return(setNames(list(data), lbl))
  }

  ## ── Matrix (RSS / pathway_scores / any feature × group matrix) ────────────
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

  ## ── data.frame ────────────────────────────────────────────────────────────
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
#'       (\code{importance}, \code{score}, \code{RSS}, or \code{value}).
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
#' # 1) SCENIC RSS matrix ─ single cell type
#' PlotRankScatter(rssMat, cell.type = "Parathyroid", topn = 10)
#'
#' # 2) SCENIC RSS matrix ─ multiple cell types
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
#' # 5) Per-class importance ─ faceted
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

  # ── Parse input into a named list of named numeric vectors ─────────────
  score_list <- .rank_scatter_parse(data, cell.type)

  # ── Optionally clean feature names ─────────────────────────────────────
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
      # Underscores → spaces
      nm <- gsub("_", " ", nm)
      names(v) <- nm
      v
    })
  }

  # ── Single panel vs multi-panel ────────────────────────────────────────
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


# ============================================================================
# 4. PlotImportance
# ============================================================================

# --------------------------------------------------------------------------
# Internal helper: clean pathway names for display
# --------------------------------------------------------------------------

.clean_pathway_names <- function(x) {
  # Strip known database prefixes (MSigDB conventions)
  cleaned <- sub("^REACTOME_", "", x)
  cleaned <- sub("^GOBP_",     "", cleaned)
  cleaned <- sub("^GOCC_",     "", cleaned)
  cleaned <- sub("^GOMF_",     "", cleaned)
  cleaned <- sub("^KEGG_",     "", cleaned)
  cleaned <- sub("^HALLMARK_", "", cleaned)
  cleaned <- sub("^WP_",       "", cleaned)
  cleaned <- sub("^BIOCARTA_", "", cleaned)
  cleaned <- sub("^PID_",      "", cleaned)

  # Replace underscores with spaces
  cleaned <- gsub("_", " ", cleaned)

  # Sentence case: first letter uppercase, rest lowercase
  cleaned <- paste0(
    toupper(substr(cleaned, 1, 1)),
    tolower(substr(cleaned, 2, nchar(cleaned)))
  )
  cleaned
}

# --------------------------------------------------------------------------
# Internal helper: lollipop / bar plot for importance values
# --------------------------------------------------------------------------

.plot_bar <- function(plot_data, display = "lollipop",
                      palette = "Reds 3", base_size = 12,
                      x_label = "Importance") {

  # --- normalise importance to [0, 1] for colour scale ---
  rng <- range(plot_data$importance, na.rm = TRUE)
  plot_data$imp_scaled <- (plot_data$importance - rng[1]) /
    (rng[2] - rng[1] + 1e-10)

  # --- colour palette (colorspace ships with ggplot2) ---
  hcl_cols <- colorspace::sequential_hcl(
    palette = palette, n = 10, rev = TRUE
  )[3:10]

  color_guide <- ggplot2::guide_colorbar(
    barwidth = 1, barheight = 12, ticks = TRUE,
    title = "Scaled", title.position = "top",
    frame.colour = "black", frame.linewidth = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5
  )

  scale_fill <- ggplot2::scale_fill_gradientn(
    colors = hcl_cols, limits = c(0, 1),
    breaks = c(0, 0.5, 1), guide = color_guide
  )

  # --- build plot ---
  if (display == "lollipop") {
    scale_color <- ggplot2::scale_color_gradientn(
      colors = hcl_cols, limits = c(0, 1),
      breaks = c(0, 0.5, 1), guide = color_guide
    )

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$importance, y = .data$var_ordered)
    ) +
      ggplot2::geom_col(
        ggplot2::aes(fill = .data$imp_scaled), width = 0.1
      ) +
      ggplot2::geom_point(
        ggplot2::aes(size = .data$importance, color = .data$imp_scaled),
        shape = 16
      ) +
      ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
      scale_fill + scale_color
  } else {
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$importance, y = .data$var_ordered)
    ) +
      ggplot2::geom_col(ggplot2::aes(fill = .data$imp_scaled)) +
      scale_fill
  }

  p <- p +
    ggplot2::labs(x = x_label, y = NULL,
                  fill = "Scaled", color = "Scaled") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.005, 0.05))
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "none",
      axis.text        = ggplot2::element_text(color = "black"),
      panel.grid.major.x = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      ),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", hjust = 0)
    )

  p
}


#' Importance Plot (Lollipop / Bar)
#'
#' Visualise gene importance scores or pathway scores returned by
#' \code{DNN_predict()}.
#'
#' \strong{Gene importance} (data.frame input): Automatically detects whether
#' the input is global importance or per-class importance based on the
#' presence of a \code{cell_type} column.
#'
#' \strong{Pathway scores} (matrix input): Accepts the
#' \code{pathway_scores} matrix from \code{DNN_predict(pathway_gmt = ...)}.
#' Long pathway names are automatically cleaned (database prefixes removed,
#' underscores replaced, sentence case applied).
#'
#' @param importance Either:
#'   \itemize{
#'     \item A \strong{data.frame} with gene importance scores:
#'       \itemize{
#'         \item \strong{Global}: columns \code{gene}, \code{importance}
#'           (i.e. \code{pred$imp_global}).
#'         \item \strong{Per-class}: columns \code{cell_type}, \code{n_cells},
#'           \code{rank}, \code{gene}, \code{importance}
#'           (i.e. \code{pred$imp_per_class}).
#'       }
#'     \item A \strong{matrix} of pathway scores
#'       (i.e. \code{pred$pathway_scores}), with cells as rows and
#'       pathways as columns.
#'   }
#' @param top_k Integer. Maximum number of top features to show.
#'   For per-class / per-group mode this is applied within each group.
#'   Default 15.
#' @param display Character, either \code{"lollipop"} (default, stem + dot)
#'   or \code{"bar"} (filled bar chart).
#' @param palette Character. An HCL sequential palette name recognised by
#'   \code{colorspace::sequential_hcl()}.  Default \code{"Reds 3"}.
#' @param ncol Integer. Number of columns for \code{facet_wrap} in
#'   per-class / per-group mode.  Default 3.
#' @param base_size Numeric. Base font size. Default 12.
#' @param group_by Optional character or factor vector of length
#'   \code{nrow(importance)} for per-group faceted pathway plots.
#'   Only used when \code{importance} is a matrix.
#'   If \code{NULL} (default), shows global mean pathway scores.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' pred <- DNN_predict(query, model_path, explain = TRUE,
#'                     pathway_gmt = "reactome.gmt")
#'
#' # Gene importance (existing usage)
#' PlotImportance(pred$imp_global)
#' PlotImportance(pred$imp_per_class, top_k = 10, ncol = 4)
#'
#' # Pathway scores: global (mean across all cells)
#' PlotImportance(pred$pathway_scores, top_k = 20)
#'
#' # Pathway scores: per-group
#' PlotImportance(pred$pathway_scores, top_k = 10,
#'                group_by = q1$cell_type_pred, ncol = 4)
#'
#' # Bar style
#' PlotImportance(pred$imp_global, display = "bar", palette = "Blues 3")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_point labs
#'   scale_x_continuous scale_size_continuous scale_fill_gradientn
#'   scale_color_gradientn guide_colorbar theme_classic theme
#'   element_text element_line element_blank expansion facet_wrap
#'   scale_y_discrete
#' @importFrom dplyr group_by slice_max ungroup mutate
#' @importFrom rlang .data
#' @export
PlotImportance <- function(importance,
                           top_k     = 15L,
                           display   = c("lollipop", "bar"),
                           palette   = "Reds 3",
                           ncol      = 3,
                           base_size = 12,
                           group_by  = NULL) {

  display <- match.arg(display)

  # ========================================================================
  # PATHWAY MODE: matrix input (from DNN_predict()$pathway_scores)
  # ========================================================================
  if (is.matrix(importance)) {

    pw_mat <- importance
    if (is.null(colnames(pw_mat)))
      stop("Pathway scores matrix must have column names (pathway names).")

    x_lab <- "Pathway Score"

    if (is.null(group_by)) {
      ## ── Global pathway mode ────────────────────────────────────────────
      mean_scores <- colMeans(pw_mat, na.rm = TRUE)
      pw_df <- data.frame(
        pathway    = names(mean_scores),
        importance = as.numeric(mean_scores),
        stringsAsFactors = FALSE
      )

      if (!is.null(top_k)) {
        pw_df <- pw_df %>%
          dplyr::slice_max(.data$importance, n = top_k,
                           with_ties = FALSE)
      }

      pw_df$var_ordered <- stats::reorder(
        .clean_pathway_names(pw_df$pathway),
        pw_df$importance
      )

      p <- .plot_bar(pw_df, display = display,
                     palette = palette, base_size = base_size,
                     x_label = x_lab)

    } else {
      ## ── Per-group pathway mode ─────────────────────────────────────────
      if (length(group_by) != nrow(pw_mat))
        stop("length(group_by) must equal nrow(importance) ",
             "(one label per cell).")

      groups <- as.character(group_by)
      unique_groups <- sort(unique(groups))

      pw_list <- lapply(unique_groups, function(g) {
        idx <- which(groups == g)
        g_means <- colMeans(pw_mat[idx, , drop = FALSE], na.rm = TRUE)
        data.frame(
          group      = g,
          n_cells    = length(idx),
          pathway    = names(g_means),
          importance = as.numeric(g_means),
          stringsAsFactors = FALSE
        )
      })
      pw_df <- do.call(rbind, pw_list)

      if (!is.null(top_k)) {
        pw_df <- pw_df %>%
          dplyr::group_by(.data$group) %>%
          dplyr::slice_max(.data$importance, n = top_k,
                           with_ties = FALSE) %>%
          dplyr::ungroup()
      }

      pw_df$facet_label <- paste0(
        pw_df$group, " (n=", pw_df$n_cells, ")"
      )

      pw_df$var_ordered <- stats::reorder(
        paste(.clean_pathway_names(pw_df$pathway),
              pw_df$group, sep = "___"),
        pw_df$importance
      )

      p <- .plot_bar(pw_df, display = display,
                     palette = palette, base_size = base_size,
                     x_label = x_lab) +
        ggplot2::facet_wrap(~ facet_label, ncol = ncol,
                            scales = "free_y") +
        ggplot2::scale_y_discrete(
          labels = function(x) gsub("___.*$", "", x)
        )
    }

    return(p)
  }

  # ========================================================================
  # GENE MODE: data.frame input (existing behavior, unchanged)
  # ========================================================================
  if (!is.data.frame(importance))
    stop("importance must be a data.frame (from DNN_predict with explain=TRUE) ",
         "or a matrix (from DNN_predict()$pathway_scores).")
  if (!all(c("gene", "importance") %in% names(importance)))
    stop("importance must contain 'gene' and 'importance' columns")

  is_per_class <- "cell_type" %in% names(importance)

  ## ── Per-class mode ────────────────────────────────────────────────────────
  if (is_per_class) {

    # top-k within each cell type
    if (!is.null(top_k)) {
      importance <- importance %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::slice_max(.data$importance, n = top_k,
                         with_ties = FALSE) %>%
        dplyr::ungroup()
    }

    # facet label: "CellType (n = xxx)"
    if ("n_cells" %in% names(importance)) {
      importance$facet_label <- paste0(
        importance$cell_type, " (n=", importance$n_cells, ")"
      )
    } else {
      importance$facet_label <- importance$cell_type
    }

    # reorder gene within each cell type for free_y faceting
    importance$var_ordered <- stats::reorder(
      paste(importance$gene, importance$cell_type, sep = "___"),
      importance$importance
    )

    p <- .plot_bar(importance, display = display,
                   palette = palette, base_size = base_size) +
      ggplot2::facet_wrap(~ facet_label, ncol = ncol, scales = "free_y") +
      ggplot2::scale_y_discrete(
        labels = function(x) gsub("___.*$", "", x)
      )

  ## ── Global mode ───────────────────────────────────────────────────────────
  } else {

    if (!is.null(top_k)) {
      importance <- importance %>%
        dplyr::slice_max(.data$importance, n = top_k,
                         with_ties = FALSE)
    }

    importance$var_ordered <- stats::reorder(
      importance$gene, importance$importance
    )

    p <- .plot_bar(importance, display = display,
                   palette = palette, base_size = base_size)
  }

  p
}


# ============================================================================
# 2. PlotGroupPreference
# ============================================================================

#' Group Preference Heatmap (O/E Ratio)
#'
#' Create a heatmap displaying the Observed / Expected ratio of cell type
#' composition across groups, highlighting group-specific enrichment.
#'
#' @param cellmeta A data.frame containing cell metadata.
#' @param group.by Character string for the grouping column.
#' @param preference.on Character string for the cell-type column.
#' @param palette RColorBrewer palette name. Default \code{"Blues"}.
#' @param column_names_rot Rotation angle (0-360) for column labels.
#'   Default 0.
#' @param ... Additional arguments passed to
#'   \code{ComplexHeatmap::Heatmap()}.
#'
#' @return A ComplexHeatmap Heatmap object.
#'
#' @examples
#' \dontrun{
#' PlotGroupPreference(q1@@meta.data,
#'   group.by = "group",
#'   preference.on = "cell_type_pred",
#'   column_names_rot = 45
#' )
#' }
#'
#' @export
PlotGroupPreference <- function(cellmeta, group.by, preference.on,
                                palette = "Blues",
                                column_names_rot = 0, ...) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Package 'ComplexHeatmap' is required. Install with:\n",
         "  BiocManager::install('ComplexHeatmap')")
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    stop("Package 'RColorBrewer' is required. Install with:\n",
         "  install.packages('RColorBrewer')")

  if (!is.data.frame(cellmeta))
    stop("cellmeta must be a data frame")
  if (!(group.by %in% names(cellmeta)) ||
      !(preference.on %in% names(cellmeta)))
    stop("group.by and preference.on must be columns in cellmeta")

  ## O/E ratio ---------------------------------------------------------------
  pop.stat.mat   <- table(cellmeta[[group.by]], cellmeta[[preference.on]])
  pop.stat.mat.o <- apply(pop.stat.mat, 1, function(xx) xx / sum(xx)) %>% t()
  pop.stat.mat.e <- colSums(pop.stat.mat) / sum(pop.stat.mat)
  pop.stat.mat.e <- matrix(
    rep(pop.stat.mat.e, nrow(pop.stat.mat)),
    ncol = nrow(pop.stat.mat)
  ) %>% t()
  pop.stat.mat.oe <- t(pop.stat.mat.o / pop.stat.mat.e)
  pop.stat.mat.oe[is.na(pop.stat.mat.oe)] <- 1

  ## Ordering ----------------------------------------------------------------
  col.names.order <- levels(cellmeta[[group.by]])
  row.names.order <- rev(levels(cellmeta[[preference.on]]))
  if (is.null(row.names.order)) {
    pop.stat.mat.oe <- pop.stat.mat.oe[
      order(pop.stat.mat.oe[, 1], decreasing = TRUE), ]
  } else {
    pop.stat.mat.oe <- pop.stat.mat.oe[row.names.order, ]
  }
  if (is.null(col.names.order)) {
    pop.stat.mat.oe <- pop.stat.mat.oe[
      , order(pop.stat.mat.oe[1, ], decreasing = TRUE)]
  } else {
    pop.stat.mat.oe <- pop.stat.mat.oe[, col.names.order]
  }

  ## Plot --------------------------------------------------------------------
  cut.off  <- stats::quantile(pop.stat.mat.oe, .95, na.rm = TRUE)
  cut.off2 <- round(stats::quantile(pop.stat.mat.oe, .9, na.rm = TRUE), 1)
  pop.stat.mat.oe.tile <- ifelse(
    pop.stat.mat.oe > cut.off, cut.off, pop.stat.mat.oe
  )

  hm_colors <- tryCatch(
    RColorBrewer::brewer.pal(n = 7, name = palette),
    error = function(e) {
      message(e$message, " -- using default: Blues")
      RColorBrewer::brewer.pal(n = 7, name = "Blues")
    }
  )

  ComplexHeatmap::Heatmap(
    pop.stat.mat.oe.tile,
    name             = "Ratio (o/e)",
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    col              = hm_colors,
    column_names_rot = column_names_rot,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (pop.stat.mat.oe[i, j] > cut.off2) {
        grid::grid.text(
          sprintf("%.2f", pop.stat.mat.oe[i, j]), x, y,
          gp = grid::gpar(fontsize = 10, col = "white", fontface = "bold")
        )
      } else {
        grid::grid.text(
          sprintf("%.2f", pop.stat.mat.oe[i, j]), x, y,
          gp = grid::gpar(fontsize = 10, col = "black")
        )
      }
    },
    ...
  )
}


# ============================================================================
# 3. PlotMAP
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
#' @param query_meta A data.frame with query cell metadata, e.g.
#'   \code{q1@@meta.data} after \code{AddMetaData(obj, pred$predictions)}.
#'   Must contain columns specified by \code{query_emb} (and optionally
#'   \code{facet_by}).
#' @param ref_emb Character vector of length 2: column names for reference
#'   UMAP axes.  Default \code{c("umap_1", "umap_2")}.
#' @param query_emb Character vector of length 2: column names for predicted
#'   query UMAP axes.  Default \code{c("umap_1_pred", "umap_2_pred")}.
#' @param color_by Column name in \code{ref} used for colouring reference
#'   cells.  Default \code{"cell_type"}.
#' @param facet_by Optional column name in \code{query_meta} for faceting
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
#' PlotMAP(ref = ref_data, query_meta = q1@@meta.data,
#'   facet_by = "group")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d
#'   scale_color_manual guides guide_legend theme_void theme
#'   element_blank coord_equal facet_grid
#' @export
PlotMAP <- function(ref, query_meta,
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
  stopifnot(all(query_emb %in% names(query_meta)))

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
      data    = query_meta,
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
      data    = query_meta,
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
    stopifnot(facet_by %in% names(query_meta))
    p <- p + ggplot2::facet_grid(
      stats::reformulate(facet_by)
    ) +
      ggplot2::theme(
        panel.spacing = ggplot2::unit(1.5, "cm")
      )
  }

  p
}


# ============================================================================
# 5. PlotPerturbation
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


# ============================================================================
# 6. PlotPercent
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


# ============================================================================
# PlotScatter — Scatter correlation plot for single-cell data
# ============================================================================

# --------------------------------------------------------------------------
# Internal: build a marginal distribution ggplot
# --------------------------------------------------------------------------

#' @keywords internal
.build_marginal <- function(plot_df, var, type, group.by, cols,
                            point.color, bins = 30) {
  # Base aes
  if (!is.null(group.by)) {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[[var]],
                   fill = .data[[group.by]],
                   color = .data[[group.by]])
    )
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[var]]))
  }

  # Geom layers by type
  if (type == "density") {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_density(alpha = 0.35, linewidth = 0.4)
    } else {
      p <- p + ggplot2::geom_density(
        fill = point.color, color = point.color, alpha = 0.35, linewidth = 0.4
      )
    }
  } else if (type == "histogram") {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_histogram(
        alpha = 0.5, position = "identity", bins = bins, linewidth = 0.2
      )
    } else {
      p <- p + ggplot2::geom_histogram(
        fill = point.color, color = "white", alpha = 0.7,
        bins = bins, linewidth = 0.2
      )
    }
  } else if (type == "boxplot") {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_boxplot(
        ggplot2::aes(y = .data[[group.by]]),
        alpha = 0.4, outlier.size = 0.3, linewidth = 0.3
      )
    } else {
      p <- p + ggplot2::geom_boxplot(
        ggplot2::aes(y = 0),
        fill = point.color, alpha = 0.4,
        outlier.size = 0.3, linewidth = 0.3
      )
    }
  } else if (type == "violin") {
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_violin(
        ggplot2::aes(y = .data[[group.by]]),
        alpha = 0.35, linewidth = 0.3, scale = "width"
      )
    } else {
      p <- p + ggplot2::geom_violin(
        ggplot2::aes(y = 0),
        fill = point.color, alpha = 0.35,
        linewidth = 0.3, scale = "width"
      )
    }
  } else if (type == "densigram") {
    # histogram + density overlay
    if (!is.null(group.by)) {
      p <- p +
        ggplot2::geom_histogram(
          ggplot2::aes(y = ggplot2::after_stat(density)),
          alpha = 0.35, position = "identity",
          bins = bins, linewidth = 0.2
        ) +
        ggplot2::geom_density(alpha = 0.2, linewidth = 0.5)
    } else {
      p <- p +
        ggplot2::geom_histogram(
          ggplot2::aes(y = ggplot2::after_stat(density)),
          fill = point.color, color = "white", alpha = 0.4,
          bins = bins, linewidth = 0.2
        ) +
        ggplot2::geom_density(
          fill = point.color, color = point.color,
          alpha = 0.2, linewidth = 0.5
        )
    }
  }

  # Color scales
  if (!is.null(group.by)) {
    p <- p +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_color_manual(values = cols)
  }

  # Minimal theme: no axis labels, no legend, no background
  p <- p +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin     = ggplot2::margin(0, 0, 0, 0)
    )

  p
}


# ============================================================================
# PlotCorrelation
# ============================================================================

#' Volcano-Style Correlation Plot
#'
#' Visualise the output of \code{\link{RunCorrelation}} or
#' \code{\link{RunTraceGene}} as a volcano-like scatter plot with
#' \code{-log10(p)} on the x-axis and correlation coefficient on the y-axis.
#' Points are colored by significance direction
#' (positive / negative / non-significant), and selected genes are labeled
#' with \code{ggrepel}.
#'
#' @param data A data.frame from \code{\link{RunCorrelation}},
#'   \code{\link{RunTraceGene}}, or \code{\link{RunTraceGSEA}}.
#'   Must contain a name column (\code{gene} or \code{Pathway}),
#'   a score column (\code{score}, \code{rho}, \code{Score}/\code{NES}),
#'   and a p-value column (\code{pvalue}/\code{PValue}, \code{padj}/
#'   \code{padjust}/\code{AdjPValue}).  Column matching is
#'   case-insensitive.  An optional grouping column (\code{group}/
#'   \code{lineage}/\code{Lineage}) enables faceting.
#' @param use.padj Logical.  Use adjusted p-value (\code{padj} or
#'   \code{padjust}) instead of raw \code{pvalue} for the x-axis and
#'   significance threshold. Default: \code{TRUE}.
#' @param p.cutoff Numeric.  Significance threshold on the chosen p-value
#'   column. Default: 0.05.
#' @param cor.cutoff Numeric.  Minimum absolute correlation to be
#'   considered significant. Default: 0 (any direction).
#' @param label Character or character vector controlling which genes to
#'   label.
#'   \itemize{
#'     \item \code{"top"} (default) — label \code{topn} genes with largest
#'       absolute correlation among significant hits.
#'     \item \code{"sig"} — label all significant genes.
#'     \item \code{"all"} — label every point.
#'     \item \code{"none"} — no labels.
#'     \item A character vector of specific gene names.
#'   }
#' @param topn Integer.  Number of top genes to label when
#'   \code{label = "top"}. Default: 10.
#' @param col.pos Character. Color for significant positive correlations.
#'   Default: \code{"#E64B35"} (red).
#' @param col.neg Character. Color for significant negative correlations.
#'   Default: \code{"#4DBBD5"} (blue).
#' @param col.ns Character. Color for non-significant points.
#'   Default: \code{"grey70"}.
#' @param size.by Character. Variable to map to point size.
#'   \itemize{
#'     \item \code{"none"} (default) — all points use \code{point.size}.
#'     \item \code{"cor"} — size mapped to absolute correlation value.
#'     \item \code{"pvalue"} — size mapped to \eqn{-\log_{10}(p)}.
#'   }
#' @param point.size Numeric. Fixed point size when \code{size.by = "none"},
#'   ignored otherwise. Default: 2.5.
#' @param size.range Numeric vector of length 2. Range of point sizes when
#'   \code{size.by} is not \code{"none"}. Default: \code{c(0.5, 5)}.
#' @param point.alpha Numeric. Point transparency. Default: 0.7.
#' @param label.size Numeric. Label text size. Default: 3.5.
#' @param label.color Character. Label text color. Default: \code{"black"}.
#' @param box.padding Numeric. Padding around label boxes (passed to
#'   \code{ggrepel::geom_text_repel}). Default: 0.5.
#' @param max.overlaps Integer. Maximum overlapping labels.
#'   Default: 20.
#' @param title Character. Plot title. Default: \code{NULL} (auto-generated
#'   from the \code{target} attribute).
#' @param xlab Character. X-axis label.  Default: auto.
#' @param ylab Character. Y-axis label.  Default: \code{NULL} (auto:
#'   \code{"Correlation"} or \code{"NES"} depending on input).
#' @param ncol Integer.  Number of columns when a \code{group} column is
#'   present and faceting is used. Default: 3.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # --- RunCorrelation output ---
#' res <- RunCorrelation(seu, target = "PTH")
#' PlotCorrelation(res)
#' PlotCorrelation(res, label = "top", topn = 15)
#' PlotCorrelation(res, label = c("GCM2", "CASR", "VDR"))
#'
#' # Map point size to |correlation|
#' PlotCorrelation(res, size.by = "cor")
#'
#' # Map point size to -log10(p)
#' PlotCorrelation(res, size.by = "pvalue", size.range = c(1, 6))
#'
#' # Per-group (RunCorrelation)
#' res <- RunCorrelation(seu, target = "PTH", group.by = "celltype")
#' PlotCorrelation(res, size.by = "cor")
#'
#' # --- RunTraceGene output (rho + padjust + lineage) ---
#' tg <- RunTraceGene(seu, lineages = c("Lineage1", "Lineage2"))
#' PlotCorrelation(tg)                   # faceted by lineage
#' PlotCorrelation(tg, size.by = "cor")  # size = |rho|
#'
#' # --- RunTraceGSEA output (Score/NES + AdjPValue + Lineage) ---
#' gsea <- RunTraceGSEA(seu, lineages = c("Lineage1", "Lineage2"),
#'                      gene.sets = gmt)
#' PlotCorrelation(gsea)                 # ylab auto = "NES"
#' PlotCorrelation(gsea, size.by = "pvalue")
#' }
#'
#' @seealso \code{\link{RunCorrelation}}, \code{\link{PlotRankScatter}},
#'   \code{\link{PlotScatter}}
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline
#'   scale_color_manual scale_size_continuous facet_wrap labs theme_bw
#'   theme element_text element_blank guide_legend
#' @importFrom rlang .data
#' @export
PlotCorrelation <- function(data,
                            use.padj     = TRUE,
                            p.cutoff     = 0.05,
                            cor.cutoff   = 0,
                            label        = "top",
                            topn         = 10L,
                            col.pos      = "#E64B35",
                            col.neg      = "#4DBBD5",
                            col.ns       = "grey70",
                            size.by      = c("none", "cor", "pvalue"),
                            point.size   = 2.5,
                            size.range   = c(0.5, 5),
                            point.alpha  = 0.7,
                            label.size   = 3.5,
                            label.color  = "black",
                            box.padding  = 0.5,
                            max.overlaps = 20L,
                            title        = NULL,
                            xlab         = NULL,
                            ylab         = NULL,
                            ncol         = 3L) {

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame (from RunCorrelation / RunTraceGene / ",
         "RunTraceGSEA).", call. = FALSE)
  }

  # Case-insensitive column resolver
  .mc <- function(candidates, cols) {
    cols_lc <- tolower(cols)
    for (cand in candidates) {
      idx <- which(cols_lc == tolower(cand))[1]
      if (!is.na(idx)) return(cols[idx])
    }
    NA_character_
  }

  # Resolve name column: gene > pathway > name > feature
  name_col <- .mc(c("gene", "pathway", "name", "feature"), colnames(data))
  if (is.na(name_col)) {
    stop("Need a name column ('gene' or 'pathway').", call. = FALSE)
  }

  # Resolve score column: score > rho > NES
  score_col <- .mc(c("score", "rho", "NES"), colnames(data))
  if (is.na(score_col)) {
    stop("Need a score column ('score', 'rho', or 'NES').", call. = FALSE)
  }

  # Resolve p-value column: padj > padjust > AdjPValue > pvalue > PValue
  if (isTRUE(use.padj)) {
    p_col <- .mc(c("padj", "padjust", "AdjPValue"), colnames(data))
  } else {
    p_col <- NA_character_
  }
  if (is.na(p_col)) {
    p_col <- .mc(c("pvalue", "PValue"), colnames(data))
  }
  if (is.na(p_col)) {
    stop("No p-value column found (pvalue/PValue/padj/padjust/AdjPValue).",
         call. = FALSE)
  }

  # Resolve group column: group > lineage > Lineage
  group_col <- .mc(c("group", "lineage"), colnames(data))
  has_group <- !is.na(group_col)

  # ── Derived columns ───────────────────────────────────────────────────────
  df <- data
  df$gene  <- df[[name_col]]           # normalise to 'gene' for labels
  df$score <- df[[score_col]]          # normalise to 'score'
  if (has_group) df$group <- df[[group_col]]   # normalise to 'group'
  df$.p     <- df[[p_col]]
  df$.logp  <- -log10(pmax(df$.p, .Machine$double.xmin))

  df$.sig <- ifelse(
    df$.p <= p.cutoff & abs(df$score) >= cor.cutoff,
    ifelse(df$score > 0, "Positive", "Negative"),
    "NS"
  )
  df$.sig <- factor(df$.sig, levels = c("Positive", "Negative", "NS"))

  # ── Auto title / axis labels ──────────────────────────────────────────────
  target_name <- attr(data, "target")
  if (is.null(title)) {
    title <- if (!is.null(target_name)) {
      paste0("Correlation with ", target_name)
    } else {
      "Correlation Volcano"
    }
  }
  if (is.null(xlab)) {
    is_adj <- tolower(p_col) %in% c("padj", "padjust", "adjpvalue")
    xlab <- bquote(-log[10] ~ italic(.(ifelse(is_adj, "p.adj", "p"))))
  }
  if (is.null(ylab)) {
    # Auto-detect GSEA output: NES column or Method == "GSEA"
    is_gsea <- tolower(score_col) == "nes" ||
      ("Method" %in% colnames(data) && any(data[["Method"]] == "GSEA"))
    ylab <- if (is_gsea) "NES" else "Correlation"
  }

  # ── Resolve label subset ──────────────────────────────────────────────────
  sig_df <- df[df$.sig != "NS", , drop = FALSE]

  if (length(label) == 1 && label == "top") {
    if (has_group) {
      # Per-group top N
      label_df <- do.call(rbind, lapply(
        split(sig_df, sig_df$group), function(sub) {
          sub <- sub[order(-abs(sub$score)), ]
          head(sub, topn)
        }
      ))
    } else {
      sig_df <- sig_df[order(-abs(sig_df$score)), ]
      label_df <- head(sig_df, topn)
    }
  } else if (length(label) == 1 && label == "sig") {
    label_df <- sig_df
  } else if (length(label) == 1 && label == "all") {
    label_df <- df
  } else if (length(label) == 1 && label == "none") {
    label_df <- df[0, ]
  } else {
    # Specific gene names
    label_df <- df[df$gene %in% label, , drop = FALSE]
  }

  # ── size.by ──────────────────────────────────────────────────────────────
  size.by <- match.arg(size.by)

  if (size.by == "cor") {
    df$.size_var <- abs(df$score)
    size_lab     <- "|Correlation|"
  } else if (size.by == "pvalue") {
    df$.size_var <- df$.logp
    size_lab     <- bquote(-log[10] ~ italic(p))
  }

  # ── Build plot ────────────────────────────────────────────────────────────
  sig_colors <- c("Positive" = col.pos, "Negative" = col.neg, "NS" = col.ns)

  if (size.by == "none") {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x     = .data$.logp,
      y     = .data$score,
      color = .data$.sig
    )) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x     = .data$.logp,
      y     = .data$score,
      color = .data$.sig,
      size  = .data$.size_var
    )) +
      ggplot2::geom_point(alpha = point.alpha) +
      ggplot2::scale_size_continuous(
        name  = size_lab,
        range = size.range,
        guide = ggplot2::guide_legend(
          override.aes = list(alpha = 1, color = "grey30"),
          order = 2
        )
      )
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = sig_colors,
      name   = "Significance",
      drop   = FALSE,
      guide  = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1),
        order = 1
      )
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = -log10(p.cutoff), linetype = "dashed",
                        color = "grey40", linewidth = 0.5)

  # ── Labels ────────────────────────────────────────────────────────────────
  if (nrow(label_df) > 0) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      message("Install 'ggrepel' for non-overlapping labels.")
      p <- p + ggplot2::geom_text(
        data    = label_df,
        mapping = ggplot2::aes(label = .data$gene),
        color   = label.color,
        size    = label.size,
        hjust   = -0.1, vjust = -0.5
      )
    } else {
      p <- p + ggrepel::geom_text_repel(
        data         = label_df,
        mapping      = ggplot2::aes(label = .data$gene),
        color        = label.color,
        size         = label.size,
        box.padding  = box.padding,
        max.overlaps = max.overlaps,
        segment.color = "grey50",
        show.legend  = FALSE
      )
    }
  }

  # ── Faceting ──────────────────────────────────────────────────────────────
  if (has_group) {
    p <- p + ggplot2::facet_wrap(~ group, ncol = ncol, scales = "free_x")
  }

  # ── Theme ─────────────────────────────────────────────────────────────────
  p <- p +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, size = 14,
                                               face = "bold"),
      axis.title       = ggplot2::element_text(size = 12),
      axis.text        = ggplot2::element_text(size = 10, color = "black"),
      strip.text       = ggplot2::element_text(size = 11, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top",
      legend.title     = ggplot2::element_text(face = "bold", size = 11),
      legend.text      = ggplot2::element_text(size = 10)
    )

  p
}


# ============================================================================
# PlotPropCorrelation
# ============================================================================

#' Cell-Type Proportion vs Gene Expression Correlation
#'
#' For each sample, compute the proportion of a given cell type and the
#' average expression of one or more genes, then draw a scatter-correlation
#' plot.  Internally delegates to \code{\link{PlotScatter}}, so all its
#' features (grouping, marginal plots, regression lines, etc.) are available.
#'
#' @param seu A Seurat object.
#' @param gene Character vector. Gene(s) whose mean expression per sample is
#'   plotted on the y-axis.  When length > 1 the plot is faceted.
#' @param celltype Character vector. Cell type(s) whose proportion per sample
#'   is plotted on the x-axis.  When length > 1 the plot is faceted.
#' @param celltype.col Character. Column in \code{meta.data} containing cell
#'   type annotations.
#' @param sample.col Character. Column in \code{meta.data} containing sample /
#'   patient identifiers.
#' @param expr.in Character vector or \code{NULL}. Cell types in which to
#'   average gene expression.  \code{NULL} (default) = all cells in each
#'   sample.  Set to e.g. \code{"Parathyroid"} to average only within that
#'   cell type.
#' @param group.by Character or \code{NULL}. Metadata column for coloring
#'   points (e.g. \code{"disease"}).  Each sample must have a unique group
#'   label. Default: \code{NULL}.
#' @param assay Character. Seurat assay. Default: \code{DefaultAssay(seu)}.
#' @param layer Character. Data layer. Default: \code{"data"}.
#' @param method Correlation method: \code{"spearman"}, \code{"pearson"}, or
#'   \code{"kendall"}. Default: \code{"spearman"}.
#' @param show.cor Logical. Show correlation statistics. Default: \code{TRUE}.
#' @param show.smooth Logical. Show regression line. Default: \code{TRUE}.
#' @param smooth.method Smoothing method for the trend line. Default:
#'   \code{"lm"}.
#' @param point.size Numeric. Point size. Default: 3.
#' @param point.alpha Numeric. Point transparency. Default: 0.7.
#' @param point.color Character. Point color when \code{group.by = NULL}.
#'   Default: \code{"#984ea3"}.
#' @param cor.size Numeric. Font size for correlation text. Default: 4.
#' @param marginal Character. Marginal plot type passed to
#'   \code{\link{PlotScatter}}: \code{"none"}, \code{"density"},
#'   \code{"histogram"}, etc.  Default: \code{"none"}.
#' @param marginal.size Numeric. Relative size of marginal plots. Default: 5.
#' @param palette Character. Palette name. Default: \code{"Paired"}.
#' @param palcolor Character vector or \code{NULL}. Custom colors.
#' @param title Character or \code{NULL}. Plot title.
#' @param ncol Integer. Facet columns. Default: 3.
#' @param return.data Logical. If \code{TRUE}, return the sample-level
#'   data.frame instead of the plot.  Default: \code{FALSE}.
#'
#' @return A \code{ggplot} / \code{patchwork} object, or a data.frame if
#'   \code{return.data = TRUE}.
#'
#' @examples
#' \dontrun{
#' # Basic: Parathyroid proportion vs PTH expression
#' PlotPropCorrelation(seu, gene = "PTH", celltype = "Parathyroid",
#'                     celltype.col = "celltype", sample.col = "sample_id")
#'
#' # Color by disease group
#' PlotPropCorrelation(seu, gene = "PTH", celltype = "Parathyroid",
#'                     celltype.col = "celltype", sample.col = "sample_id",
#'                     group.by = "disease")
#'
#' # Multiple cell types
#' PlotPropCorrelation(seu, gene = "PTH",
#'                     celltype = c("Parathyroid", "Stromal"),
#'                     celltype.col = "celltype", sample.col = "sample_id")
#'
#' # Multiple genes
#' PlotPropCorrelation(seu, gene = c("PTH", "GCM2"),
#'                     celltype = "Parathyroid",
#'                     celltype.col = "celltype", sample.col = "sample_id",
#'                     marginal = "density")
#'
#' # Return data for custom plotting
#' df <- PlotPropCorrelation(seu, gene = "PTH", celltype = "Parathyroid",
#'                           celltype.col = "celltype",
#'                           sample.col = "sample_id",
#'                           return.data = TRUE)
#' }
#'
#' @seealso \code{\link{PlotScatter}}
#' @export
PlotPropCorrelation <- function(seu,
                                gene,
                                celltype,
                                celltype.col,
                                sample.col,
                                expr.in       = NULL,
                                group.by      = NULL,
                                assay         = NULL,
                                layer         = "data",
                                method        = c("spearman", "pearson",
                                                  "kendall"),
                                show.cor      = TRUE,
                                show.smooth   = TRUE,
                                smooth.method = "lm",
                                point.size    = 3,
                                point.alpha   = 0.7,
                                point.color   = "#984ea3",
                                cor.size      = 4,
                                marginal      = "none",
                                marginal.size = 5,
                                palette       = "Paired",
                                palcolor      = NULL,
                                title         = NULL,
                                ncol          = 3L,
                                return.data   = FALSE) {

  method <- match.arg(method)

  # ── Input validation ──────────────────────────────────────────────────────
  if (!inherits(seu, "Seurat")) {
    stop("'seu' must be a Seurat object.", call. = FALSE)
  }
  meta <- seu@meta.data

  for (col in c(celltype.col, sample.col)) {
    if (!col %in% colnames(meta)) {
      stop("Column '", col, "' not found in meta.data.", call. = FALSE)
    }
  }
  if (!is.null(group.by) && !group.by %in% colnames(meta)) {
    stop("group.by column '", group.by, "' not found in meta.data.",
         call. = FALSE)
  }

  # Resolve assay
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seu)
  expr_mat <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)

  missing_g <- setdiff(gene, rownames(expr_mat))
  if (length(missing_g) > 0) {
    stop("Gene(s) not found in assay '", assay, "': ",
         paste(missing_g, collapse = ", "), call. = FALSE)
  }

  # ── 1. Cell-type proportions per sample ─────────────────────────────────
  ct_table <- table(meta[[sample.col]], meta[[celltype.col]])
  ct_prop  <- prop.table(ct_table, margin = 1)   # row-normalise
  samples  <- rownames(ct_prop)

  missing_ct <- setdiff(celltype, colnames(ct_prop))
  if (length(missing_ct) > 0) {
    stop("Cell type(s) not found in '", celltype.col, "': ",
         paste(missing_ct, collapse = ", "), call. = FALSE)
  }

  # ── 2. Average expression per sample ────────────────────────────────────
  if (!is.null(expr.in)) {
    cells_use <- colnames(seu)[meta[[celltype.col]] %in% expr.in]
    if (length(cells_use) == 0) {
      stop("No cells found for expr.in = ",
           paste(expr.in, collapse = ", "), call. = FALSE)
    }
  } else {
    cells_use <- colnames(seu)
  }

  sample_vec <- meta[cells_use, sample.col]

  avg_list <- lapply(gene, function(g) {
    tapply(as.numeric(expr_mat[g, cells_use]), sample_vec, mean, na.rm = TRUE)
  })
  names(avg_list) <- gene

  # ── 3. Build sample-level data.frame ────────────────────────────────────
  df_list <- list()
  for (ct in celltype) {
    for (g in gene) {
      sub_df <- data.frame(
        sample     = samples,
        proportion = as.numeric(ct_prop[samples, ct]),
        expression = as.numeric(avg_list[[g]][samples]),
        stringsAsFactors = FALSE
      )
      # Facet label
      need_facet <- (length(celltype) > 1) || (length(gene) > 1)
      if (need_facet) {
        if (length(celltype) > 1 && length(gene) > 1) {
          sub_df$.facet <- paste0(ct, " | ", g)
        } else if (length(celltype) > 1) {
          sub_df$.facet <- ct
        } else {
          sub_df$.facet <- g
        }
      }
      df_list[[paste0(ct, "_", g)]] <- sub_df
    }
  }
  plot_df <- do.call(rbind, df_list)
  rownames(plot_df) <- NULL

  # Map sample → group (first occurrence)
  if (!is.null(group.by)) {
    sample_grp <- meta[!duplicated(meta[[sample.col]]),
                       c(sample.col, group.by), drop = FALSE]
    colnames(sample_grp)[1] <- "sample"
    plot_df <- merge(plot_df, sample_grp, by = "sample", all.x = TRUE)
  }

  # Remove samples with NA expression (e.g. expr.in has no cells)
  plot_df <- plot_df[is.finite(plot_df$expression), , drop = FALSE]

  if (nrow(plot_df) == 0) {
    stop("No valid data points. Check sample.col / celltype / gene.",
         call. = FALSE)
  }

  if (return.data) return(plot_df)

  # ── 4. Axis labels ──────────────────────────────────────────────────────
  x_lab <- if (length(celltype) == 1) {
    paste0(celltype, " proportion")
  } else {
    "Cell-type proportion"
  }
  y_lab <- if (length(gene) == 1) {
    paste0(gene, " expression")
  } else {
    "Expression"
  }

  # ── 5. PlotScatter ─────────────────────────────────────────────────────
  split_var <- if (".facet" %in% colnames(plot_df)) ".facet" else NULL

  PlotScatter(
    object        = plot_df,
    var1          = "proportion",
    var2          = "expression",
    group.by      = group.by,
    split.by      = split_var,
    method        = method,
    smooth.method = smooth.method,
    show.cor      = show.cor,
    show.smooth   = show.smooth,
    point.size    = point.size,
    point.alpha   = point.alpha,
    point.color   = point.color,
    cor.size      = cor.size,
    marginal      = marginal,
    marginal.size = marginal.size,
    palette       = palette,
    palcolor      = palcolor,
    title         = title,
    ncol          = ncol
  )
}


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
                        ncol          = 3,
                        palette       = "Paired",
                        palcolor      = NULL,
                        point.color   = "#984ea3",
                        rug.color     = "#7fc97f",
                        title         = NULL,
                        marginal      = c("none", "density", "histogram",
                                          "boxplot", "violin", "densigram"),
                        marginal.size = 5,
                        raster        = NULL,
                        raster.dpi    = 300,
                        ...) {

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
    if (!is.null(group.by)) {
      p <- p + ggplot2::geom_smooth(
        method = smooth.method, se = show.se, na.rm = TRUE,
        linewidth = smooth.size
      )
    } else {
      p <- p + ggplot2::geom_smooth(
        method = smooth.method, se = show.se, na.rm = TRUE,
        linewidth = smooth.size, color = smooth.color
      )
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
      p <- p + ggpubr::stat_cor(
        method  = method,
        digits  = cor.digits,
        size    = cor.size,
        label.x.npc = "left",
        label.y.npc = "top"
      )
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


# ============================================================================
# PlotAnnotation — annotation-bar heatmap for single-cell metadata
# ============================================================================

#' @keywords internal
.auto_anno_colors <- function(meta, columns, palette, palcolor) {
  col_list <- list()
  n_discrete <- sum(!vapply(meta[columns], is.numeric, logical(1)))
  disc_idx <- 0L

  for (col in columns) {
    vals <- meta[[col]]
    if (is.numeric(vals)) {
      rng <- range(vals, na.rm = TRUE)
      col_list[[col]] <- circlize::colorRamp2(rng, c("white", "red"))
    } else {
      disc_idx <- disc_idx + 1L
      lvls <- if (is.factor(vals)) levels(vals) else sort(unique(as.character(vals)))
      if (n_discrete > 1) {
        # rotate palette so each variable gets visually distinct colors
        pal_choices <- c("Paired", "Set1", "Set2", "Set3", "Dark2",
                         "npg", "nejm", "lancet", "jama", "aaas", "d3")
        use_pal <- if (disc_idx == 1) palette else {
          pal_choices[((disc_idx - 1) %% length(pal_choices)) + 1]
        }
      } else {
        use_pal <- palette
      }
      col_list[[col]] <- palette_colors(lvls, palette = use_pal,
                                        palcolor = if (disc_idx == 1) palcolor else NULL)
    }
  }
  col_list
}

#' Annotation Heatmap for Single-Cell Metadata
#'
#' Create a ComplexHeatmap-style annotation heatmap showing cell metadata as
#' colored annotation bars. Optionally includes a gene expression heatmap
#' matrix below the annotations.
#'
#' @param seu A Seurat object.
#' @param columns Character vector of metadata column names to display as
#'   annotation bars (e.g., \code{c("celltype", "sample", "Phase")}).
#' @param sort.by Column name to sort cells by. Default: first element
#'   of \code{columns}.
#' @param features Optional character vector of gene names. If provided,
#'   a scaled expression heatmap is drawn below the annotation bars.
#'   Default: NULL (annotation-only mode).
#' @param anno_point Optional metadata column name to display as a point
#'   annotation at the top (e.g., \code{"nCount_RNA"}).
#' @param downsample Integer; maximum number of cells per group (defined
#'   by \code{sort.by}) to keep. Default: NULL (no downsampling).
#' @param palette Palette name for discrete variables. Default: "Paired".
#' @param palcolor Optional custom color vector (overrides palette for the
#'   first discrete variable).
#' @param use_raster Logical; rasterize the expression heatmap for speed.
#'   Default: TRUE.
#' @param show_column_names Logical; show cell barcode labels. Default: FALSE.
#' @param assay Assay to use for expression data. Default: NULL
#'   (DefaultAssay).
#' @param layer Layer to pull expression from. Default: "data".
#' @param scale_rows Logical; z-score scale each gene across cells.
#'   Default: TRUE.
#' @param ... Additional arguments passed to
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object.
#'
#' @examples
#' \dontrun{
#' # Annotation-only (like clinical correlation heatmap)
#' PlotAnnotation(seu, columns = c("celltype", "sample", "Phase"))
#'
#' # With gene expression heatmap
#' PlotAnnotation(seu, columns = c("celltype", "sample"),
#'                features = c("CD3D", "CD14", "MS4A1"),
#'                sort.by = "celltype")
#'
#' # With point annotation and downsampling
#' PlotAnnotation(seu, columns = c("celltype", "condition"),
#'                anno_point = "nCount_RNA", downsample = 200)
#' }
#'
#' @export
PlotAnnotation <- function(seu,
                           columns,
                           sort.by = NULL,
                           features = NULL,
                           anno_point = NULL,
                           downsample = NULL,
                           palette = "Paired",
                           palcolor = NULL,
                           use_raster = TRUE,
                           show_column_names = FALSE,
                           assay = NULL,
                           layer = "data",
                           scale_rows = TRUE,
                           ...) {

  # --- dependency checks ---
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Package 'ComplexHeatmap' is required. Install with:\n",
         "  BiocManager::install('ComplexHeatmap')")
  if (!inherits(seu, "Seurat"))
    stop("'seu' must be a Seurat object.", call. = FALSE)

  meta <- seu@meta.data

  # validate columns
  missing_cols <- setdiff(columns, colnames(meta))
  if (length(missing_cols) > 0)
    stop("Columns not found in metadata: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  if (!is.null(anno_point) && !anno_point %in% colnames(meta))
    stop("anno_point '", anno_point, "' not found in metadata.", call. = FALSE)

  # --- sort.by default ---
  if (is.null(sort.by)) sort.by <- columns[1]
  if (!sort.by %in% colnames(meta))
    stop("sort.by '", sort.by, "' not found in metadata.", call. = FALSE)

  # --- downsample ---
  if (!is.null(downsample)) {
    grp <- as.character(meta[[sort.by]])
    keep <- unlist(lapply(split(rownames(meta), grp), function(ids) {
      if (length(ids) <= downsample) ids
      else sample(ids, downsample)
    }))
    meta <- meta[keep, , drop = FALSE]
  }

  # --- sort ---
  sort_vals <- meta[[sort.by]]
  if (is.factor(sort_vals)) {
    ord <- order(as.integer(sort_vals))
  } else {
    ord <- order(sort_vals)
  }
  meta <- meta[ord, , drop = FALSE]
  cells <- rownames(meta)

  # --- build annotation colors ---
  col_list <- .auto_anno_colors(meta, columns, palette, palcolor)

  # --- build annotation args ---
  anno_args <- list()
  # point annotation on top
  if (!is.null(anno_point)) {
    anno_args[[anno_point]] <- ComplexHeatmap::anno_points(
      meta[[anno_point]],
      size = grid::unit(0.5, "mm"),
      gp = grid::gpar(col = "grey40")
    )
  }
  # data frame of annotation bars
  anno_args[["df"]] <- meta[, columns, drop = FALSE]
  anno_args[["col"]] <- col_list

  ha <- do.call(ComplexHeatmap::HeatmapAnnotation, anno_args)

  # --- build heatmap ---
  if (is.null(features)) {
    # annotation-only: empty matrix
    ht <- ComplexHeatmap::Heatmap(
      matrix(nrow = 0, ncol = length(cells)),
      top_annotation = ha,
      show_column_names = show_column_names,
      ...
    )
  } else {
    # expression heatmap
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(seu)
    expr_mat <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
    avail <- intersect(features, rownames(expr_mat))
    if (length(avail) == 0)
      stop("None of the requested features found in assay '", assay, "'.",
           call. = FALSE)
    mat <- as.matrix(expr_mat[avail, cells, drop = FALSE])
    if (scale_rows) {
      mat <- t(scale(t(mat)))
      mat[is.nan(mat)] <- 0
    }
    # clamp for better visualization
    cap <- quantile(abs(mat), 0.99, na.rm = TRUE)
    mat[mat > cap] <- cap
    mat[mat < -cap] <- -cap

    hm_col <- circlize::colorRamp2(
      c(-cap, 0, cap), c("blue", "white", "red")
    )
    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = "Expression",
      top_annotation = ha,
      col = hm_col,
      cluster_columns = FALSE,
      cluster_rows = length(avail) > 1,
      show_column_names = show_column_names,
      show_row_names = TRUE,
      use_raster = use_raster,
      ...
    )
  }

  ht
}
