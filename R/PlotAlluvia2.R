# ============================================================================
# PlotAlluvia2 — Alluvial Plot with gaps between strata (ggalluvial-based)
# ============================================================================

#' Alluvial Plot with Gaps Between Cell Types
#'
#' Create an alluvial plot showing cell type composition changes across
#' conditions/groups, with visible gaps between strata (cell types).
#' Based on \pkg{ggalluvial}, identical to \code{\link{PlotAlluvia}} but
#' with spacing between cell type strata for cleaner visual separation.
#' Gaps are implemented by inserting invisible (transparent) strata between
#' real cell types.
#'
#' @param cellmeta A Seurat object or data.frame containing cell metadata.
#' @param by Character string specifying the grouping variable.
#'   Default \code{"group"}.
#' @param fill Character string specifying the cell type variable.
#'   Default \code{"celltype"}.
#' @param palette Character. Colour palette name passed to
#'   \code{\link{palette_colors}}. Default \code{"lancet"}.
#' @param palcolor Character vector. Custom colours overriding
#'   \code{palette}. Default \code{NULL}.
#' @param flow.alpha Numeric (0-1), transparency of alluvial flows.
#'   Default 0.4.
#' @param bar.width Numeric (0-1), width of stratum bars. Default 0.5.
#' @param bar.color Stratum bar border colour. Default \code{"gray50"}.
#' @param gap Numeric. Size of gap between strata as a fraction of total
#'   height. Default 0.005. Set 0 for no gaps (same as PlotAlluvia).
#' @param show.label Logical. Default \code{FALSE}.
#' @param label.size Numeric. Default 3.
#' @param label.color Character. Default \code{"black"}.
#' @param label.box Logical. If \code{TRUE}, boxed labels. Default \code{FALSE}.
#' @param label.fill Background fill for boxed labels. Default \code{"white"}.
#' @param label.alpha Transparency of boxed label background. Default \code{1}.
#' @param show.pct Logical. Default \code{FALSE}.
#' @param label.style Character. Label format: \code{"default"}, \code{"all"},
#'   \code{"n"}, \code{"pct"}, or \code{"name"}.
#' @param label.min.pct Numeric (0-1). Minimum proportion to show label.
#'   Default 0.
#' @param y.type Character. \code{"percent"} (default) or \code{"count"}.
#' @param legend.ncol Integer. Default 1.
#' @param base.size Numeric. Default 15.
#' @param curve.type Character. Flow curve type. Default \code{"sigmoid"}.
#'   Options: \code{"linear"}, \code{"cubic"}, \code{"quintic"},
#'   \code{"sine"}, \code{"arctangent"}, \code{"sigmoid"},
#'   \code{"xspline"}.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' library(ToyData)
#'
#' # --- Basic usage (with default gaps) ---
#' PlotAlluvia2(toy_test)
#'
#' # --- Adjust gap size ---
#' PlotAlluvia2(toy_test, gap = 0.01)   # larger gaps
#' PlotAlluvia2(toy_test, gap = 0)      # no gaps (same as PlotAlluvia)
#'
#' # --- Labels ---
#' PlotAlluvia2(toy_test, label.style = "pct", label.min.pct = 0.05)
#' PlotAlluvia2(toy_test, label.style = "pct", label.min.pct = 0.05,
#'              label.box = TRUE)
#'
#' # --- Custom colours ---
#' PlotAlluvia2(toy_test, palcolor = UtilsR::pal_paraSC)
#' PlotAlluvia2(toy_test, palette = "Paired")
#'
#' # --- Raw counts ---
#' PlotAlluvia2(toy_test, y.type = "count", label.style = "n")
#' }
#'
#' @export
PlotAlluvia2 <- function(cellmeta, by = "group", fill = "celltype",
                         palette = "lancet", palcolor = NULL,
                         flow.alpha = 0.4,
                         bar.width = 0.5, bar.color = "gray50",
                         gap = 0.01,
                         show.label = FALSE, label.size = 3,
                         label.color = "black", label.box = TRUE,
                         label.fill = "white", label.alpha = 1,
                         show.pct = TRUE,
                         label.style = c("all", "default", "n", "pct", "name"),
                         label.min.pct = 0,
                         y.type = c("percent", "count"),
                         legend.ncol = 1, base.size = 15,
                         curve.type = "linear") {

  label.style <- match.arg(label.style)
  y.type <- match.arg(y.type)

  if (!requireNamespace("ggalluvial", quietly = TRUE))
    stop("Package 'ggalluvial' is required. Install with: ",
         "install.packages('ggalluvial')")

  cellmeta <- .extract_cellmeta(cellmeta)
  if (!is.character(by) || !is.character(fill))
    stop("by and fill must be character strings.")
  if (!(by %in% names(cellmeta)) || !(fill %in% names(cellmeta)))
    stop("by and fill must be columns in cellmeta.")

  # ---- Ensure factors ----
  if (!is.factor(cellmeta[[by]])) cellmeta[[by]] <- factor(cellmeta[[by]])
  if (!is.factor(cellmeta[[fill]])) cellmeta[[fill]] <- factor(cellmeta[[fill]])

  # ---- Colours ----
  colors <- palette_colors(levels(cellmeta[[fill]]),
                           palette = palette, palcolor = palcolor)

  # ---- Frequency table ----
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

  # Keep all levels (including 0-count for missing types)
  stat <- stat[stat$Freq > 0 | ave(stat$Freq, stat$Celltype, FUN = sum) > 0, ]

  # ---- Within-group proportions ----
  stat <- stat %>%
    dplyr::group_by(.data$Group) %>%
    dplyr::mutate(
      group_total = sum(.data$Freq),
      pct = ifelse(.data$group_total > 0,
                   .data$Freq / .data$group_total, 0)
    ) %>%
    dplyr::ungroup()

  # ---- Insert gaps ----
  n_types <- nlevels(stat$Celltype)
  total_gap <- gap * (n_types - 1)

  if (y.type == "percent") {
    scale_factor <- (1 - total_gap)
    stat$y_val <- stat$pct * scale_factor
  } else {
    stat$y_val <- stat$Freq
  }

  # Add invisible gap strata between each celltype
  if (gap > 0) {
    gap_rows <- list()
    ct_levels <- levels(stat$Celltype)
    grp_levels <- levels(stat$Group)
    gap_ct_levels <- character(0)

    for (i in seq_len(n_types - 1)) {
      gap_name <- paste0(".gap_", i)
      gap_ct_levels <- c(gap_ct_levels, gap_name)
      for (g in grp_levels) {
        gap_val <- if (y.type == "percent") {
          gap
        } else {
          gap * max(stat$group_total)
        }
        gap_rows <- c(gap_rows, list(data.frame(
          Group = g, Celltype = gap_name, Freq = 0,
          group_total = stat$group_total[stat$Group == g][1],
          pct = 0, y_val = gap_val,
          stringsAsFactors = FALSE
        )))
      }
    }
    gap_df <- do.call(rbind, gap_rows)

    # Interleave: ct1, gap1, ct2, gap2, ..., ctN
    new_levels <- character(0)
    for (i in seq_along(ct_levels)) {
      new_levels <- c(new_levels, ct_levels[i])
      if (i < length(ct_levels)) {
        new_levels <- c(new_levels, gap_ct_levels[i])
      }
    }

    stat <- rbind(stat, gap_df)
    stat$Celltype <- factor(stat$Celltype, levels = new_levels)
    stat$Group <- factor(stat$Group, levels = grp_levels)
  }

  # ---- Labels ----
  stat$label <- ifelse(
    stat$pct >= label.min.pct &
      !grepl("^\\.gap_", stat$Celltype), {
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

  # ---- Build plot ----
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

  # ---- Label layer ----
  # Extract actual stratum positions from ggplot_build, then overlay labels
  use_label <- if (label.style != "default") TRUE else (show.label || show.pct)
  if (use_label) {
    # Build to get computed stratum positions
    built <- ggplot2::ggplot_build(p)
    # geom_stratum is the 2nd layer (after geom_alluvium)
    stratum_data <- built$data[[2]]

    # Map stratum integer back to Celltype factor level
    ct_all_levels <- levels(stat$Celltype)
    stratum_data$Celltype <- ct_all_levels[stratum_data$stratum]

    # Map x integer back to Group
    grp_levels <- levels(stat$Group)
    stratum_data$Group <- grp_levels[as.integer(stratum_data$x)]

    # Compute center y
    stratum_data$y_center <- (stratum_data$ymin + stratum_data$ymax) / 2

    # Merge with labels from stat
    label_lookup <- stat[, c("Group", "Celltype", "label")]
    label_lookup <- label_lookup[!is.na(label_lookup$label), ]
    label_lookup$Group <- as.character(label_lookup$Group)
    label_lookup$Celltype <- as.character(label_lookup$Celltype)

    stratum_data$Group_chr <- as.character(stratum_data$Group)
    stratum_data$Celltype_chr <- as.character(stratum_data$Celltype)

    label_df <- merge(stratum_data, label_lookup,
                      by.x = c("Group_chr", "Celltype_chr"),
                      by.y = c("Group", "Celltype"), all.x = FALSE)

    # Filter gap strata and NA labels
    label_df <- label_df[!grepl("^\\.gap_", label_df$Celltype_chr) &
                           !is.na(label_df$label), ]

    if (nrow(label_df) > 0) {
      if (label.box) {
        p <- p + ggplot2::geom_label(
          data = label_df,
          ggplot2::aes(x = .data$x, y = .data$y_center,
                       label = .data$label),
          inherit.aes = FALSE,
          size = label.size, color = label.color,
          fill = label.fill, alpha = label.alpha
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(x = .data$x, y = .data$y_center,
                       label = .data$label),
          inherit.aes = FALSE,
          size = label.size, color = label.color
        )
      }
    }
  }

  # ---- Colours: hide gap strata ----
  gap_names <- grep("^\\.gap_", levels(stat$Celltype), value = TRUE)
  real_names <- setdiff(levels(stat$Celltype), gap_names)
  gap_cols <- stats::setNames(rep("transparent", length(gap_names)), gap_names)
  real_cols <- stats::setNames(as.character(colors[real_names]), real_names)
  all_cols <- c(real_cols, gap_cols)

  p <- p + ggplot2::scale_fill_manual(
    values = all_cols,
    breaks = real_names
  )

  # ---- Y-axis ----
  if (y.type == "percent") {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::label_percent())
  }

  # ---- Theme ----
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
