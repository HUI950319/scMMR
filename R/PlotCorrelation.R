# ============================================================================
# PlotCorrelation — Volcano-Style Correlation Plot
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
#' @param name.col Character.  Column name to use as the label/name
#'   (e.g. gene or pathway names).  Default: \code{NULL} (auto-detect
#'   from \code{gene}/\code{pathway}/\code{name}/\code{feature},
#'   case-insensitive).
#' @param score.col Character.  Column name to use as the score (y-axis).
#'   Default: \code{NULL} (auto-detect from \code{score}/\code{rho}/
#'   \code{NES}).
#' @param p.col Character.  Column name to use as the p-value (x-axis).
#'   Default: \code{NULL} (auto-detect from \code{padj}/\code{padjust}/
#'   \code{pvalue}/\code{PValue}).  Overrides \code{use.padj} when set.
#' @param use.padj Logical.  Use adjusted p-value (\code{padj} or
#'   \code{padjust}) instead of raw \code{pvalue} for the x-axis and
#'   significance threshold.  Ignored when \code{p.col} is specified.
#'   Default: \code{TRUE}.
#' @param p.cutoff Numeric.  Significance threshold on the chosen p-value
#'   column. Default: 0.05.
#' @param cor.cutoff Numeric.  Minimum absolute correlation to be
#'   considered significant. Default: 0 (any direction).
#' @param label Character or character vector controlling which genes to
#'   label.
#'   \itemize{
#'     \item \code{"top"} (default) — label \code{topn} positive and
#'       \code{topn} negative genes among significant hits.
#'     \item \code{"sig"} — label all significant genes.
#'     \item \code{"all"} — label every point.
#'     \item \code{"none"} — no labels.
#'     \item A character vector of specific gene names.
#'   }
#' @param topn Integer.  Number of top genes to label per direction
#'   (positive and negative) when \code{label = "top"}. Default: 5.
#' @param col.pos Character. Color for significant positive correlations.
#'   Default: \code{"#E64B35"} (red).
#' @param col.neg Character. Color for significant negative correlations.
#'   Default: \code{"#4DBBD5"} (blue).
#' @param col.ns Character. Color for non-significant points.
#'   Default: \code{"grey70"}.
#' @param size.by Character. Variable to map to point size.
#'   Built-in shortcuts: \code{"none"} (default, fixed size),
#'   \code{"cor"} (absolute correlation), \code{"pvalue"}
#'   (\eqn{-\log_{10}(p)}). Or pass any column name in \code{data}.
#' @param point.size Numeric. Fixed point size when \code{size.by = "none"}.
#'   Default: 2.5.
#' @param size.range Numeric vector of length 2. Range of point sizes when
#'   \code{size.by} is not \code{"none"}. Default: \code{c(0.5, 5)}.
#' @param alpha.by Character. Variable to map to point transparency.
#'   Built-in shortcuts: \code{"none"} (default, fixed alpha),
#'   \code{"cor"} (absolute correlation), \code{"pvalue"}
#'   (\eqn{-\log_{10}(p)}). Or pass any column name in \code{data}.
#' @param point.alpha Numeric. Fixed point transparency when
#'   \code{alpha.by = "none"}. Default: 0.7.
#' @param alpha.range Numeric vector of length 2. Range of alpha values
#'   when \code{alpha.by} is not \code{"none"}. Default: \code{c(0.1, 1)}.
#' @param label.size Numeric. Label text size. Default: 3.5.
#' @param label.color Character. Label text color. Default: \code{"black"}.
#' @param box.padding Numeric. Padding around label boxes. Default: 0.5.
#' @param max.overlaps Integer. Maximum overlapping labels. Default: 20.
#' @param title Character. Plot title. Default: \code{NULL} (auto).
#' @param xlab Character. X-axis label. Default: auto.
#' @param ylab Character. Y-axis label. Default: \code{NULL} (auto).
#' @param ncol Integer. Number of facet columns. Default: 3.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' res <- RunCorrelation(seu, target = "PTH")
#' PlotCorrelation(res)                              # top 5 pos + 5 neg
#' PlotCorrelation(res, topn = 10)                   # top 10 pos + 10 neg
#' PlotCorrelation(res, label = c("GCM2", "CASR"))   # specific genes
#'
#' # Custom score and p-value columns
#' PlotCorrelation(res, score.col = "prop_cor", p.col = "prop_padj")
#'
#' # Map size and alpha to custom columns
#' PlotCorrelation(res, size.by = "prop_cor", alpha.by = "de_abs_logFC")
#'
#' # Per-group
#' res <- RunCorrelation(seu, target = "PTH", group.by = "celltype")
#' PlotCorrelation(res, size.by = "cor", alpha.by = "cor")
#'
#' # RunTraceGene output
#' tg <- RunTraceGene(seu, lineages = c("Lineage1", "Lineage2"))
#' PlotCorrelation(tg, size.by = "cor")
#'
#' # RunTraceGSEA output
#' gsea <- RunTraceGSEA(seu, lineages = c("Lineage1"), gene.sets = gmt)
#' PlotCorrelation(gsea, size.by = "pvalue", alpha.by = "pvalue")
#' }
#'
#' @seealso \code{\link{RunCorrelation}}, \code{\link{PlotRankScatter}},
#'   \code{\link{PlotScatter}}
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline
#'   scale_color_manual scale_size_continuous scale_alpha_continuous
#'   facet_wrap labs theme_bw
#'   theme element_text element_blank guide_legend
#' @importFrom rlang .data
#' @export
PlotCorrelation <- function(data,
                            name.col     = NULL,
                            score.col    = NULL,
                            p.col        = NULL,
                            use.padj     = TRUE,
                            p.cutoff     = 0.05,
                            cor.cutoff   = 0,
                            label        = "top",
                            topn         = 5L,
                            col.pos      = "#E64B35",
                            col.neg      = "#4DBBD5",
                            col.ns       = "grey70",
                            size.by      = "none",
                            point.size   = 2.5,
                            size.range   = c(0.5, 5),
                            alpha.by     = "none",
                            point.alpha  = 0.7,
                            alpha.range  = c(0.1, 1),
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

  # Resolve name column: user-specified or auto-detect
  if (!is.null(name.col)) {
    if (!name.col %in% colnames(data)) {
      stop(sprintf("name.col '%s' not found in data. Available: %s",
                   name.col, paste(colnames(data), collapse = ", ")),
           call. = FALSE)
    }
    name_col <- name.col
  } else {
    name_col <- .mc(c("gene", "pathway", "name", "feature"), colnames(data))
    if (is.na(name_col)) {
      stop("No name column auto-detected (gene/pathway/name/feature). ",
           "Set name.col manually.", call. = FALSE)
    }
  }

  # Resolve score column: user-specified or auto-detect
  if (!is.null(score.col)) {
    if (!score.col %in% colnames(data)) {
      stop(sprintf("score.col '%s' not found in data.", score.col), call. = FALSE)
    }
    s_col <- score.col
  } else {
    s_col <- .mc(c("score", "rho", "NES"), colnames(data))
    if (is.na(s_col)) {
      stop("Need a score column ('score', 'rho', or 'NES'), or set score.col.",
           call. = FALSE)
    }
  }

  # Resolve p-value column: user-specified or auto-detect
  if (!is.null(p.col)) {
    if (!p.col %in% colnames(data)) {
      stop(sprintf("p.col '%s' not found in data.", p.col), call. = FALSE)
    }
    pv_col <- p.col
  } else {
    if (isTRUE(use.padj)) {
      pv_col <- .mc(c("padj", "padjust", "AdjPValue"), colnames(data))
    } else {
      pv_col <- NA_character_
    }
    if (is.na(pv_col)) {
      pv_col <- .mc(c("pvalue", "PValue"), colnames(data))
    }
    if (is.na(pv_col)) {
      stop("No p-value column found (pvalue/PValue/padj/padjust/AdjPValue), or set p.col.",
           call. = FALSE)
    }
  }

  # Resolve group column: group > lineage > Lineage
  group_col <- .mc(c("group", "lineage"), colnames(data))
  has_group <- !is.na(group_col)

  # ── Derived columns ───────────────────────────────────────────────────────
  df <- data
  df$gene  <- df[[name_col]]
  df$score <- df[[s_col]]
  if (has_group) df$group <- df[[group_col]]
  df$.p     <- df[[pv_col]]
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
    is_adj <- tolower(pv_col) %in% c("padj", "padjust", "adjpvalue")
    xlab <- bquote(-log[10] ~ italic(.(ifelse(is_adj, "p.adj", "p"))))
  }
  if (is.null(ylab)) {
    is_gsea <- tolower(s_col) == "nes" ||
      ("Method" %in% colnames(data) && any(data[["Method"]] == "GSEA"))
    ylab <- if (is_gsea) "Normalized Enrichment Score" else "Correlation"
  }

  # ── Resolve label subset (top N positive + top N negative) ────────────────
  sig_df <- df[df$.sig != "NS", , drop = FALSE]

  if (length(label) == 1 && label == "top") {
    .pick_top_both <- function(sub, n) {
      pos <- sub[sub$score > 0, , drop = FALSE]
      neg <- sub[sub$score < 0, , drop = FALSE]
      pos <- pos[order(-pos$score), ]
      neg <- neg[order(neg$score), ]
      rbind(head(pos, n), head(neg, n))
    }
    if (has_group) {
      label_df <- do.call(rbind, lapply(
        split(sig_df, sig_df$group), function(sub) .pick_top_both(sub, topn)
      ))
    } else {
      label_df <- .pick_top_both(sig_df, topn)
    }
  } else if (length(label) == 1 && label == "sig") {
    label_df <- sig_df
  } else if (length(label) == 1 && label == "all") {
    label_df <- df
  } else if (length(label) == 1 && label == "none") {
    label_df <- df[0, ]
  } else {
    label_df <- df[df$gene %in% label, , drop = FALSE]
  }

  # ── Helper: resolve aesthetic variable ────────────────────────────────────
  .resolve_aes <- function(by, df, aes_name) {
    if (by == "none") return(list(var = NULL, lab = NULL))
    if (by == "cor") {
      df$.tmp <- abs(df$score)
      return(list(var = ".tmp", lab = "|Correlation|"))
    }
    if (by == "pvalue") {
      df$.tmp <- df$.logp
      return(list(var = ".tmp", lab = bquote(-log[10] ~ italic(p))))
    }
    if (!by %in% colnames(df)) {
      stop(sprintf("Column '%s' not found in data for %s mapping.", by, aes_name),
           call. = FALSE)
    }
    list(var = by, lab = by)
  }

  # ── size.by ──────────────────────────────────────────────────────────────
  size_info <- .resolve_aes(size.by, df, "size.by")
  if (!is.null(size_info$var)) {
    if (size_info$var == ".tmp") {
      if (size.by == "cor") df$.size_var <- abs(df$score)
      else df$.size_var <- df$.logp
      size_info$var <- ".size_var"
    }
    size_lab <- size_info$lab
  }

  # ── alpha.by ─────────────────────────────────────────────────────────────
  alpha_info <- .resolve_aes(alpha.by, df, "alpha.by")
  if (!is.null(alpha_info$var)) {
    if (alpha_info$var == ".tmp") {
      if (alpha.by == "cor") df$.alpha_var <- abs(df$score)
      else df$.alpha_var <- df$.logp
      alpha_info$var <- ".alpha_var"
    }
    alpha_lab <- alpha_info$lab
  }

  # ── Build plot ────────────────────────────────────────────────────────────
  sig_colors <- c("Positive" = col.pos, "Negative" = col.neg, "NS" = col.ns)
  use_size  <- !is.null(size_info$var)
  use_alpha <- !is.null(alpha_info$var)

  if (use_size && use_alpha) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data$.logp, y = .data$score, color = .data$.sig,
      size = .data[[size_info$var]], alpha = .data[[alpha_info$var]]
    )) + ggplot2::geom_point()
  } else if (use_size) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data$.logp, y = .data$score, color = .data$.sig,
      size = .data[[size_info$var]]
    )) + ggplot2::geom_point(alpha = point.alpha)
  } else if (use_alpha) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data$.logp, y = .data$score, color = .data$.sig,
      alpha = .data[[alpha_info$var]]
    )) + ggplot2::geom_point(size = point.size)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data$.logp, y = .data$score, color = .data$.sig
    )) + ggplot2::geom_point(size = point.size, alpha = point.alpha)
  }

  if (use_size) {
    p <- p + ggplot2::scale_size_continuous(
      name = size_lab, range = size.range,
      guide = ggplot2::guide_legend(
        override.aes = list(alpha = 1, color = "grey30"), order = 2))
  }
  if (use_alpha) {
    p <- p + ggplot2::scale_alpha_continuous(
      name = alpha_lab, range = alpha.range,
      guide = ggplot2::guide_legend(
        override.aes = list(size = 3, color = "grey30"), order = 3))
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = sig_colors, name = "Significance", drop = FALSE,
      guide = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1), order = 1)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = -log10(p.cutoff), linetype = "dashed",
                        color = "grey40", linewidth = 0.5)

  # ── Labels ────────────────────────────────────────────────────────────────
  if (nrow(label_df) > 0) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      message("Install 'ggrepel' for non-overlapping labels.")
      p <- p + ggplot2::geom_text(
        data = label_df, mapping = ggplot2::aes(label = .data$gene),
        color = label.color, size = label.size, hjust = -0.1, vjust = -0.5)
    } else {
      p <- p + ggrepel::geom_text_repel(
        data = label_df, mapping = ggplot2::aes(label = .data$gene),
        color = label.color, size = label.size,
        box.padding = box.padding, max.overlaps = max.overlaps,
        segment.color = "grey50", show.legend = FALSE)
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
