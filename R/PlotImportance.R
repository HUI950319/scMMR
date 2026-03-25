# ============================================================================
# PlotImportance — Importance lollipop / bar plot
# ============================================================================

# --------------------------------------------------------------------------
# Internal helper: clean pathway names for display
# --------------------------------------------------------------------------

#' @keywords internal
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

#' @keywords internal
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
      ## -- Global pathway mode --
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
      ## -- Per-group pathway mode --
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

  ## -- Per-class mode --
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

  ## -- Global mode --
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
