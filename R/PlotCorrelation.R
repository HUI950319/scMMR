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
