# ============================================================================
# Dynamic feature visualization along pseudotime
# ============================================================================

# --------------------------------------------------------------------------
# Internal helper: fit a single feature along pseudotime
# --------------------------------------------------------------------------

#' @keywords internal
.fit_dynamic_feature <- function(y, x, fit_method, smooth_k, loess_span,
                                 bspline_knot, family, ...) {
  ord <- order(x)
  x_ord <- x[ord]
  y_ord <- y[ord]
  n <- length(x_ord)
  pvalue <- NA_real_
  r_sq   <- NA_real_

  if (fit_method == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for fit_method = 'gam'.\n",
           "Install with: install.packages('mgcv')")
    }
    k_use <- min(smooth_k, n - 1)
    mod <- mgcv::gam(
      y_val ~ s(x_val, bs = "cs", k = k_use),
      family = family,
      data = data.frame(y_val = y_ord, x_val = x_ord),
      ...
    )
    pre <- stats::predict(mod, type = "link", se.fit = TRUE)
    fit_link <- pre$fit
    se_link <- pre$se.fit
    upr_link <- fit_link + 2 * se_link
    lwr_link <- fit_link - 2 * se_link
    fitted_vals <- mod$family$linkinv(fit_link)
    upr_vals <- mod$family$linkinv(upr_link)
    lwr_vals <- mod$family$linkinv(lwr_link)
    # statistics
    mod_summary <- summary(mod)
    pvalue <- mod_summary$s.table[1, "p-value"]
    r_sq   <- max(0, min(1, mod_summary$r.sq))

  } else if (fit_method == "loess") {
    mod <- stats::loess(y_val ~ x_val,
                        data = data.frame(y_val = y_ord, x_val = x_ord),
                        span = loess_span, ...)
    pre <- stats::predict(mod, se = TRUE)
    fitted_vals <- pre$fit
    se_vals <- pre$se.fit
    upr_vals <- fitted_vals + 2 * se_vals
    lwr_vals <- fitted_vals - 2 * se_vals
    # loess: compute R² and approximate F-test p-value
    SSE <- sum((y_ord - fitted_vals)^2, na.rm = TRUE)
    SST <- sum((y_ord - mean(y_ord, na.rm = TRUE))^2, na.rm = TRUE)
    r_sq <- if (SST > 0) max(0, 1 - SSE / SST) else 0
    enp <- mod$trace.hat  # equivalent number of parameters
    df_model <- enp - 1
    df_resid <- n - enp
    if (df_model > 0 && df_resid > 0 && SSE > 0) {
      f_stat <- ((SST - SSE) / df_model) / (SSE / df_resid)
      pvalue <- stats::pf(f_stat, df1 = df_model, df2 = df_resid,
                           lower.tail = FALSE)
    }

  } else if (fit_method == "bspline") {
    df_use <- bspline_knot + 3
    B <- splines::bs(x_ord, df = df_use, intercept = FALSE)
    B <- B[, apply(B, 2, stats::sd) > 0, drop = FALSE]
    B <- cbind(1, B)
    mod <- stats::lm.fit(B, y_ord)
    fitted_vals <- mod$fitted.values
    res <- mod$residuals
    sigma2 <- sum(res^2) / max(1, n - ncol(B))
    hat_diag <- rowSums(B * (B %*% chol2inv(chol(crossprod(B)))))
    se_vals <- sqrt(sigma2 * hat_diag)
    upr_vals <- fitted_vals + 2 * se_vals
    lwr_vals <- fitted_vals - 2 * se_vals
    # F-test: spline model vs intercept-only
    SSE <- sum(res^2)
    SST <- sum((y_ord - mean(y_ord))^2)
    r_sq <- if (SST > 0) max(0, 1 - SSE / SST) else 0
    df_model <- ncol(B) - 1
    df_resid <- n - ncol(B)
    if (df_model > 0 && df_resid > 0 && SSE > 0) {
      f_stat <- ((SST - SSE) / df_model) / (SSE / df_resid)
      pvalue <- stats::pf(f_stat, df1 = df_model, df2 = df_resid,
                           lower.tail = FALSE)
    }

  } else {
    stop("Unknown fit_method: ", fit_method)
  }

  list(
    curve = data.frame(
      pseudotime = x_ord,
      fitted     = as.numeric(fitted_vals),
      upr        = as.numeric(upr_vals),
      lwr        = as.numeric(lwr_vals),
      stringsAsFactors = FALSE
    ),
    pvalue = pvalue,
    r_sq   = r_sq
  )
}


# ============================================================================
# PlotDynamicFeatures
# ============================================================================

#' Plot Dynamic Features Along Pseudotime
#'
#' Visualise gene expression or metadata features as a function of pseudotime
#' (or any continuous trajectory variable).
#' Supports three fitting methods with user-controllable smoothness:
#' \itemize{
#'   \item \strong{GAM} (\code{mgcv}): control smoothness via \code{smooth_k}
#'         (smaller = smoother).
#'   \item \strong{loess}: control smoothness via \code{loess_span}
#'         (larger = smoother).
#'   \item \strong{B-spline}: control smoothness via \code{bspline_knot}
#'         (more knots = more flexible / less smooth).
#' }
#'
#' @param srt A Seurat object.
#' @param pseudotime Character vector of column name(s) in
#'   \code{srt@@meta.data} containing pseudotime values.
#'   Multiple names enable multi-lineage comparison.
#' @param features Character vector of gene names (looked up in the assay)
#'   and/or column names in \code{srt@@meta.data}.
#' @param group.by Optional character string.
#'   A column in \code{srt@@meta.data} used to colour scatter points.
#' @param assay Character. Seurat assay to use. Default: active assay.
#' @param layer Character. Data layer to use. Default \code{"counts"}.
#' @param fit_method Character, one of \code{"gam"}, \code{"loess"}, or
#'   \code{"bspline"}. Default \code{"gam"}.
#' @param smooth_k Integer. Number of basis dimensions for the GAM smooth
#'   term (\code{mgcv::s(k = ...)}). Smaller values produce smoother
#'   curves. Default 10.
#' @param loess_span Numeric (0, 1]. The \code{span} parameter for
#'   \code{stats::loess()}. Larger values produce smoother curves.
#'   Default 0.75.
#' @param bspline_knot Integer. Number of internal knots for B-spline
#'   fitting. More knots allow a more flexible (less smooth) curve.
#'   Default 3.
#' @param family Character or family object for GAM. If \code{NULL}
#'   (default), \code{"gaussian"} is used.
#' @param exp_method Character, one of \code{"log1p"}, \code{"raw"}, or
#'   \code{"zscore"}. Transformation applied to expression values before
#'   plotting. Default \code{"log1p"}.
#' @param lib_normalize Logical.
#'   Whether to normalise by library size. Default \code{TRUE} when
#'   \code{layer = "counts"}.
#' @param add_point Logical. Show scatter points. Default \code{TRUE}.
#' @param add_line Logical. Show fitted curve. Default \code{TRUE}.
#' @param add_interval Logical. Show 95\% confidence interval ribbon.
#'   Default \code{TRUE}.
#' @param add_rug Logical. Show rug marks along x-axis. Default \code{TRUE}.
#' @param pt.size Numeric. Point size. Default 1.
#' @param line.size Numeric. Line width. Default 1.
#' @param point_palette Character. Palette name for scatter points / rug
#'   coloured by \code{group.by}. See \code{\link{show_palettes}} for
#'   available names. Default \code{"Paired"}.
#' @param point_palcolor Optional character vector of custom colours for
#'   the \code{group.by} variable. Overrides \code{point_palette}.
#' @param line_palette Character. Palette name for fitted lines / ribbons.
#'   Default \code{"Dark2"}.
#' @param line_palcolor Optional character vector of custom colours for
#'   lineages or features. Overrides \code{line_palette}.
#' @param compare_lineages Logical. When multiple \code{pseudotime} columns
#'   are given, overlay them in a single panel per feature. Default \code{TRUE}.
#' @param compare_features Logical. When multiple \code{features} are given,
#'   overlay them in a single panel per lineage. Default \code{FALSE}.
#' @param ncol Integer. Number of columns in the facet layout.
#' @param nrow Integer. Number of rows in the facet layout.
#' @param reverse Logical. Reverse the pseudotime axis. Default \code{FALSE}.
#' @param flip Logical. Flip x and y axes. Default \code{FALSE}.
#' @param seed Integer. Random seed. Default 11.
#' @param ... Additional arguments passed to the fitting function
#'   (\code{mgcv::gam}, \code{stats::loess}, or \code{stats::lm.fit}).
#'
#' @return A \code{ggplot} object (or list of ggplot objects if both
#'   \code{compare_lineages} and \code{compare_features} are \code{FALSE}).
#'
#' @examples
#' \dontrun{
#' # After running Slingshot or any trajectory method:
#' # Default GAM fitting
#' PlotDynamicFeatures(srt, pseudotime = "Lineage1",
#'                     features = c("Gene1", "Gene2"))
#'
#' # Smoother curve (smaller k)
#' PlotDynamicFeatures(srt, pseudotime = "Lineage1",
#'                     features = "Gene1", smooth_k = 5)
#'
#' # More flexible curve (larger k)
#' PlotDynamicFeatures(srt, pseudotime = "Lineage1",
#'                     features = "Gene1", smooth_k = 20)
#'
#' # Use loess with custom span
#' PlotDynamicFeatures(srt, pseudotime = "Lineage1",
#'                     features = "Gene1",
#'                     fit_method = "loess", loess_span = 0.3)
#'
#' # B-spline with 5 knots
#' PlotDynamicFeatures(srt, pseudotime = "Lineage1",
#'                     features = "Gene1",
#'                     fit_method = "bspline", bspline_knot = 5)
#'
#' # Compare two lineages
#' PlotDynamicFeatures(srt,
#'                     pseudotime = c("Lineage1", "Lineage2"),
#'                     features = c("Gene1", "Gene2"),
#'                     compare_lineages = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon geom_rug
#'   scale_color_manual scale_fill_manual facet_wrap facet_grid labs theme
#'   theme_classic element_text scale_x_continuous scale_y_continuous
#'   expansion coord_flip guides guide_legend
#' @importFrom rlang .data
#' @export
PlotDynamicFeatures <- function(
    srt,
    pseudotime,
    features,
    group.by      = NULL,
    assay         = NULL,
    layer         = "counts",
    fit_method    = c("gam", "loess", "bspline"),
    smooth_k      = 10,
    loess_span    = 0.75,
    bspline_knot  = 3,
    family        = NULL,
    exp_method    = c("log1p", "raw", "zscore"),
    lib_normalize = (layer == "counts"),
    add_point     = TRUE,
    add_line      = TRUE,
    add_interval  = TRUE,
    add_rug       = TRUE,
    pt.size       = 1,
    line.size     = 1,
    point_palette  = "Paired",
    point_palcolor = NULL,
    line_palette   = "Dark2",
    line_palcolor  = NULL,
    compare_lineages = TRUE,
    compare_features = FALSE,
    ncol          = NULL,
    nrow          = NULL,
    reverse       = FALSE,
    flip          = FALSE,
    seed          = 11,
    ...) {

  set.seed(seed)
  fit_method <- match.arg(fit_method)
  exp_method <- match.arg(exp_method)

  if (is.null(family)) family <- "gaussian"

  # --- validate pseudotime columns ---
  missing_pt <- pseudotime[!pseudotime %in% colnames(srt@meta.data)]
  if (length(missing_pt) > 0) {
    stop("Pseudotime column(s) not found in meta.data: ",
         paste(missing_pt, collapse = ", "))
  }

  # --- resolve assay ---
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(srt)
  }

  # --- separate gene vs meta features ---
  gene_features <- features[features %in% rownames(srt[[assay]])]
  meta_features <- features[features %in% colnames(srt@meta.data)]
  features <- unique(c(gene_features, meta_features))
  if (length(features) == 0) {
    stop("No valid features found in assay or meta.data.")
  }

  # --- extract expression matrix (genes x cells) ---
  if (length(gene_features) > 0) {
    expr_mat <- as.matrix(
      SeuratObject::GetAssayData(srt, assay = assay, layer = layer)[
        gene_features, , drop = FALSE
      ]
    )
  } else {
    expr_mat <- matrix(nrow = 0, ncol = ncol(srt))
    colnames(expr_mat) <- colnames(srt)
  }

  # --- library size normalisation ---
  if (isTRUE(lib_normalize) && length(gene_features) > 0) {
    lib_sizes <- Matrix::colSums(
      SeuratObject::GetAssayData(srt, assay = assay, layer = "counts")
    )
    med_lib <- stats::median(lib_sizes)
    size_factors <- lib_sizes / med_lib
    expr_mat <- sweep(expr_mat, 2, size_factors, "/")
  }

  # --- expression transformation ---
  if (exp_method == "log1p" && length(gene_features) > 0) {
    expr_mat <- log1p(expr_mat)
  } else if (exp_method == "zscore" && length(gene_features) > 0) {
    row_means <- rowMeans(expr_mat)
    row_sds   <- apply(expr_mat, 1, stats::sd)
    row_sds[row_sds == 0] <- 1
    expr_mat <- sweep(sweep(expr_mat, 1, row_means, "-"), 1, row_sds, "/")
  }

  # --- build per-cell data for each (lineage, feature) ---
  df_list  <- list()
  stat_list <- list()

  for (lin in pseudotime) {
    # cells with finite pseudotime on this lineage
    pt_vals <- srt@meta.data[[lin]]
    valid_cells <- colnames(srt)[is.finite(pt_vals)]
    if (length(valid_cells) < 10) {
      warning("Lineage '", lin, "' has < 10 valid cells; skipping.")
      next
    }
    pt_valid <- pt_vals[match(valid_cells, colnames(srt))]

    for (feat in features) {
      # get expression / metadata values
      if (feat %in% gene_features) {
        y_vals <- expr_mat[feat, valid_cells]
      } else {
        y_vals <- as.numeric(srt@meta.data[valid_cells, feat])
        if (exp_method == "log1p") y_vals <- log1p(y_vals)
        if (exp_method == "zscore") {
          m <- mean(y_vals, na.rm = TRUE)
          s <- stats::sd(y_vals, na.rm = TRUE)
          if (is.na(s) || s == 0) s <- 1
          y_vals <- (y_vals - m) / s
        }
      }

      # remove NA
      keep <- is.finite(y_vals) & is.finite(pt_valid)
      if (sum(keep) < 10) next
      x_use <- pt_valid[keep]
      y_use <- y_vals[keep]
      cells_use <- valid_cells[keep]

      # fit
      fit_result <- tryCatch(
        .fit_dynamic_feature(
          y = y_use, x = x_use,
          fit_method = fit_method,
          smooth_k = smooth_k,
          loess_span = loess_span,
          bspline_knot = bspline_knot,
          family = family, ...
        ),
        error = function(e) {
          warning("Fitting failed for ", feat, " on ", lin, ": ", e$message)
          NULL
        }
      )
      if (is.null(fit_result)) next

      # raw points data
      raw_df <- data.frame(
        Cell       = cells_use,
        Pseudotime = x_use,
        Expression = y_use,
        Feature    = feat,
        Lineage    = lin,
        stringsAsFactors = FALSE
      )

      # fitted curve data
      curve_df <- fit_result$curve
      curve_df$Feature <- feat
      curve_df$Lineage <- lin

      # collect statistics
      key <- paste(lin, feat, sep = "|")
      stat_list[[key]] <- data.frame(
        Feature = feat, Lineage = lin,
        pvalue = fit_result$pvalue, R2 = fit_result$r_sq,
        stringsAsFactors = FALSE
      )

      df_list[[key]] <- list(raw = raw_df, fit = curve_df)
    }
  }

  if (length(df_list) == 0) {
    stop("No valid feature-lineage combinations could be fitted.")
  }

  raw_all  <- do.call(rbind, lapply(df_list, `[[`, "raw"))
  fit_all  <- do.call(rbind, lapply(df_list, `[[`, "fit"))
  stat_all <- do.call(rbind, stat_list)
  rownames(raw_all) <- NULL
  rownames(fit_all) <- NULL
  rownames(stat_all) <- NULL

  # --- add group info ---
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(srt@meta.data)) {
      stop("group.by '", group.by, "' not found in meta.data.")
    }
    raw_all$Group <- srt@meta.data[raw_all$Cell, group.by]
    if (!is.factor(raw_all$Group)) {
      raw_all$Group <- factor(raw_all$Group)
    }
  }

  raw_all$Feature <- factor(raw_all$Feature, levels = features)
  raw_all$Lineage <- factor(raw_all$Lineage, levels = pseudotime)
  fit_all$Feature <- factor(fit_all$Feature, levels = features)
  fit_all$Lineage <- factor(fit_all$Lineage, levels = pseudotime)

  # --- determine y-axis label ---
  y_label <- switch(exp_method,
    log1p  = paste0("log1p(", layer, ")"),
    raw    = layer,
    zscore = paste0("zscore(", layer, ")")
  )

  # --- decide faceting and colour-by logic ---
  n_lin  <- length(unique(fit_all$Lineage))
  n_feat <- length(unique(fit_all$Feature))

  # determine what the line color represents
  if (isTRUE(compare_lineages) && n_lin > 1) {
    line_color_by <- "Lineage"
  } else if (isTRUE(compare_features) && n_feat > 1) {
    line_color_by <- "Feature"
  } else {
    line_color_by <- NULL
  }

  # determine facet formula
  if (isTRUE(compare_lineages) && isTRUE(compare_features)) {
    facet_formula <- NULL
  } else if (isTRUE(compare_lineages) && !isTRUE(compare_features)) {
    facet_formula <- if (n_feat > 1) ~ Feature else NULL
  } else if (!isTRUE(compare_lineages) && isTRUE(compare_features)) {
    facet_formula <- if (n_lin > 1) ~ Lineage else NULL
  } else {
    # neither compare: facet by both
    if (n_lin > 1 && n_feat > 1) {
      facet_formula <- Lineage ~ Feature
    } else if (n_feat > 1) {
      facet_formula <- ~ Feature
    } else if (n_lin > 1) {
      facet_formula <- ~ Lineage
    } else {
      facet_formula <- NULL
    }
  }

  # --- build plot ---
  p <- ggplot2::ggplot()

  # --- resolve color palettes using palette_colors() ---
  # point/rug colors (group.by)
  if (!is.null(group.by)) {
    point_cols <- palette_colors(
      raw_all$Group,
      palette  = point_palette,
      palcolor = point_palcolor
    )
  }

  # line/ribbon colors (lineage or feature)
  if (!is.null(line_color_by)) {
    line_cols <- palette_colors(
      fit_all[[line_color_by]],
      palette  = line_palette,
      palcolor = line_palcolor
    )
  }
  # default single color for line/ribbon when no multi-comparison
  default_line_col <- palette_colors(
    "default",
    palette  = line_palette,
    palcolor = line_palcolor
  )[1]

  # scatter points + rug (share the same Group color scale)
  if (!is.null(group.by)) {
    # points colored by group
    if (isTRUE(add_point)) {
      p <- p + ggplot2::geom_point(
        data = raw_all,
        ggplot2::aes(
          x     = .data$Pseudotime,
          y     = .data$Expression,
          color = .data$Group
        ),
        size  = pt.size,
        alpha = 0.6
      )
    }
    # rug colored by group (same color scale as points)
    if (isTRUE(add_rug)) {
      p <- p + ggplot2::geom_rug(
        data = raw_all,
        ggplot2::aes(
          x     = .data$Pseudotime,
          color = .data$Group
        ),
        alpha  = 0.5,
        length = grid::unit(0.04, "npc"),
        show.legend = FALSE
      )
    }
    # apply palette_colors result for points + rug
    p <- p + ggplot2::scale_color_manual(
      values = point_cols,
      name   = group.by
    )
    # start new color scale for fitted lines
    p <- p + ggnewscale::new_scale_color()
  } else {
    # no group: grey points
    if (isTRUE(add_point)) {
      p <- p + ggplot2::geom_point(
        data = raw_all,
        ggplot2::aes(
          x = .data$Pseudotime,
          y = .data$Expression
        ),
        size  = pt.size,
        alpha = 0.4,
        color = "grey60"
      )
    }
    # rug without group color
    if (isTRUE(add_rug)) {
      p <- p + ggplot2::geom_rug(
        data = raw_all,
        ggplot2::aes(x = .data$Pseudotime),
        alpha  = 0.3,
        length = grid::unit(0.04, "npc"),
        show.legend = FALSE
      )
    }
  }

  # confidence interval ribbon
  if (isTRUE(add_interval)) {
    if (!is.null(line_color_by)) {
      p <- p + ggplot2::geom_ribbon(
        data = fit_all,
        ggplot2::aes(
          x     = .data$pseudotime,
          ymin  = .data$lwr,
          ymax  = .data$upr,
          fill  = .data[[line_color_by]],
          group = interaction(.data$Lineage, .data$Feature)
        ),
        alpha = 0.2
      )
    } else {
      p <- p + ggplot2::geom_ribbon(
        data = fit_all,
        ggplot2::aes(
          x    = .data$pseudotime,
          ymin = .data$lwr,
          ymax = .data$upr
        ),
        fill  = default_line_col,
        alpha = 0.2
      )
    }
  }

  # fitted line
  if (isTRUE(add_line)) {
    if (!is.null(line_color_by)) {
      p <- p + ggplot2::geom_line(
        data = fit_all,
        ggplot2::aes(
          x     = .data$pseudotime,
          y     = .data$fitted,
          color = .data[[line_color_by]],
          group = interaction(.data$Lineage, .data$Feature)
        ),
        linewidth = line.size
      )
    } else {
      p <- p + ggplot2::geom_line(
        data = fit_all,
        ggplot2::aes(
          x = .data$pseudotime,
          y = .data$fitted
        ),
        color     = default_line_col,
        linewidth = line.size
      )
    }
  }

  # apply line/fill palette for multi-lineage / multi-feature comparisons
  if (!is.null(line_color_by)) {
    p <- p + ggplot2::scale_color_manual(
      values = line_cols,
      name   = line_color_by
    )
    p <- p + ggplot2::scale_fill_manual(
      values = line_cols,
      name   = line_color_by
    )
  }

  # x-axis transformation
  x_trans <- ifelse(isTRUE(reverse), "reverse", "identity")
  p <- p +
    ggplot2::scale_x_continuous(trans = x_trans,
                                expand = ggplot2::expansion(c(0.02, 0.02))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0.05, 0.05)))

  # faceting
  if (!is.null(facet_formula)) {
    if (inherits(facet_formula, "formula") &&
        length(facet_formula) == 3) {
      # two-sided formula: Lineage ~ Feature
      p <- p + ggplot2::facet_grid(facet_formula, scales = "free_y")
    } else {
      p <- p + ggplot2::facet_wrap(facet_formula, scales = "free_y",
                                   ncol = ncol, nrow = nrow)
    }
  }

  # --- annotate p-value and R² ---
  # build label for each facet panel
  stat_all$label <- vapply(seq_len(nrow(stat_all)), function(i) {
    pv <- stat_all$pvalue[i]
    r2 <- stat_all$R2[i]
    pv_txt <- if (is.na(pv)) {
      "P = NA"
    } else if (pv < 2.2e-16) {
      "P < 2.2e-16"
    } else if (pv < 0.001) {
      sprintf("P = %.2e", pv)
    } else {
      sprintf("P = %.4f", pv)
    }
    r2_txt <- if (is.na(r2)) "R\u00b2 = NA" else sprintf("R\u00b2 = %.3f", r2)
    paste0(pv_txt, "\n", r2_txt)
  }, character(1))

  stat_all$Feature <- factor(stat_all$Feature, levels = features)
  stat_all$Lineage <- factor(stat_all$Lineage, levels = pseudotime)

  # compute per-panel label positions (top-left corner)
  label_df <- merge(
    stat_all,
    do.call(rbind, lapply(split(raw_all, list(raw_all$Feature, raw_all$Lineage),
                                drop = TRUE), function(d) {
      data.frame(
        Feature = d$Feature[1], Lineage = d$Lineage[1],
        x_pos = min(d$Pseudotime, na.rm = TRUE),
        y_pos = max(d$Expression, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })),
    by = c("Feature", "Lineage")
  )

  # when compare_lineages: multiple lineages share panel, show combined label
  if (isTRUE(compare_lineages) && length(pseudotime) > 1) {
    label_df$label <- vapply(seq_len(nrow(label_df)), function(i) {
      paste0(label_df$Lineage[i], ": ", label_df$label[i])
    }, character(1))
    # stack labels for same feature
    label_df <- do.call(rbind, lapply(
      split(label_df, label_df$Feature), function(d) {
        d$label_combined <- paste(d$label, collapse = "\n")
        d[1, , drop = FALSE]
      }
    ))
    label_df$label <- label_df$label_combined
    label_df$x_pos <- ave(raw_all$Pseudotime, raw_all$Feature,
                          FUN = min)[match(label_df$Feature, raw_all$Feature)]
    label_df$y_pos <- ave(raw_all$Expression, raw_all$Feature,
                          FUN = function(x) max(x, na.rm = TRUE))[
      match(label_df$Feature, raw_all$Feature)
    ]
  }

  p <- p + ggplot2::geom_text(
    data = label_df,
    ggplot2::aes(x = .data$x_pos, y = .data$y_pos, label = .data$label),
    hjust = 0, vjust = 1, size = 3.2, color = "black",
    inherit.aes = FALSE
  )

  # labels and theme
  p <- p +
    ggplot2::labs(x = "Pseudotime", y = y_label) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  if (isTRUE(flip)) {
    p <- p + ggplot2::coord_flip()
  }

  p
}
