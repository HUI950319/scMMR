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
#' @param stat_method Character. Method for computing the annotation
#'   statistics displayed on each panel. \code{"fit"} (default) shows
#'   the p-value and R-squared from the fitted curve. \code{"spearman"}
#'   shows the Spearman correlation coefficient (rho) and FDR-adjusted
#'   p-value, matching \code{\link{RunTraceGene}} output.
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
#' @param raster Logical or \code{NULL}. Whether to rasterise scatter
#'   points via \code{\link{rasterise_layer}}. When \code{NULL}
#'   (default), rasterisation is enabled automatically if the number of
#'   cells exceeds 100 000. Set \code{TRUE} / \code{FALSE} to override.
#' @param raster.dpi Numeric. DPI for rasterised points. Default 300.
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
    stat_method   = c("fit", "spearman"),
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
    raster        = NULL,
    raster.dpi    = 300,
    seed          = 11,
    ...) {

  set.seed(seed)
  fit_method <- match.arg(fit_method)
  exp_method <- match.arg(exp_method)
  stat_method <- match.arg(stat_method)

  # ---- Auto raster ----
  n_cells <- base::ncol(srt)
  raster <- raster %||% (n_cells > 1e5)
  if (!is.numeric(raster.dpi) || length(raster.dpi) != 1) {
    stop("'raster.dpi' must be a numeric scalar.", call. = FALSE)
  }

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
    expr_mat <- matrix(nrow = 0, ncol = base::ncol(srt))
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
      if (stat_method == "spearman") {
        ct <- stats::cor.test(x_use, y_use, method = "spearman", exact = FALSE)
        stat_list[[key]] <- data.frame(
          Feature = feat, Lineage = lin,
          pvalue = ct$p.value, R2 = fit_result$r_sq,
          rho = unname(ct$estimate),
          stringsAsFactors = FALSE)
      } else {
        stat_list[[key]] <- data.frame(
          Feature = feat, Lineage = lin,
          pvalue = fit_result$pvalue, R2 = fit_result$r_sq,
          rho = NA_real_,
          stringsAsFactors = FALSE)
      }

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

  # adjust p-values when using spearman correlation
  if (stat_method == "spearman") {
    stat_all$pvalue <- stats::p.adjust(stat_all$pvalue, method = "fdr")
  }

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

  # --- resolve color palettes ---
  if (!is.null(group.by)) {
    point_cols <- palette_colors(
      raw_all$Group, palette = point_palette, palcolor = point_palcolor
    )
  }
  if (!is.null(line_color_by)) {
    line_cols <- palette_colors(
      fit_all[[line_color_by]], palette = line_palette, palcolor = line_palcolor
    )
  }
  default_line_col <- palette_colors(
    "default", palette = line_palette, palcolor = line_palcolor
  )[1]

  # --- stat_cor style annotation (ggpubr-compatible) ---
  .format_cor_label <- function(pv, r2, rho, stat_method) {
    # Format like ggpubr::stat_cor: italic R, italic p
    if (stat_method == "spearman") {
      r_txt <- if (is.na(rho)) "NA" else sprintf("%.3f", rho)
      p_txt <- if (is.na(pv)) "NA"
               else if (pv < 2.2e-16) "2.2e-16"
               else if (pv < 0.001)   sprintf("%.2e", pv)
               else sprintf("%.3g", pv)
      # Use rho symbol for spearman
      paste0("italic(rho) == ", r_txt, "*\",\"~~italic(p) ",
             if (!is.na(pv) && pv < 2.2e-16) "< " else "== ", p_txt)
    } else {
      r2_txt <- if (is.na(r2)) "NA" else sprintf("%.3f", r2)
      p_txt  <- if (is.na(pv)) "NA"
                else if (pv < 2.2e-16) "2.2e-16"
                else if (pv < 0.001)   sprintf("%.2e", pv)
                else sprintf("%.3g", pv)
      paste0("italic(R)^2 == ", r2_txt, "*\",\"~~italic(p) ",
             if (!is.na(pv) && pv < 2.2e-16) "< " else "== ", p_txt)
    }
  }

  stat_all$Feature <- factor(stat_all$Feature, levels = features)
  stat_all$Lineage <- factor(stat_all$Lineage, levels = pseudotime)

  # --- build per-feature plots and combine with patchwork ---
  plot_list <- list()

  for (feat in features) {
    raw_sub  <- raw_all[raw_all$Feature == feat, , drop = FALSE]
    fit_sub  <- fit_all[fit_all$Feature == feat, , drop = FALSE]
    stat_sub <- stat_all[stat_all$Feature == feat, , drop = FALSE]

    if (nrow(raw_sub) == 0) next

    # Start ggplot with raw_sub as default data (enables p[[i]]$data access)
    pi <- ggplot2::ggplot(data = raw_sub)

    # --- Points ---
    if (!is.null(group.by) && isTRUE(add_point)) {
      .pt_l <- ggplot2::geom_point(
        data = raw_sub,
        ggplot2::aes(x = .data$Pseudotime, y = .data$Expression,
                     color = .data$Group),
        size = pt.size, alpha = 0.6
      )
      pi <- pi + if (isTRUE(raster)) rasterise_layer(.pt_l, dpi = raster.dpi) else .pt_l
    } else if (isTRUE(add_point)) {
      .pt_l <- ggplot2::geom_point(
        data = raw_sub,
        ggplot2::aes(x = .data$Pseudotime, y = .data$Expression),
        size = pt.size, alpha = 0.4, color = "grey60"
      )
      pi <- pi + if (isTRUE(raster)) rasterise_layer(.pt_l, dpi = raster.dpi) else .pt_l
    }

    # --- Rug ---
    if (isTRUE(add_rug)) {
      if (!is.null(group.by)) {
        pi <- pi + ggplot2::geom_rug(
          data = raw_sub,
          ggplot2::aes(x = .data$Pseudotime, color = .data$Group),
          alpha = 0.5, length = grid::unit(0.04, "npc"), show.legend = FALSE
        )
      } else {
        pi <- pi + ggplot2::geom_rug(
          data = raw_sub,
          ggplot2::aes(x = .data$Pseudotime),
          alpha = 0.3, length = grid::unit(0.04, "npc"), show.legend = FALSE
        )
      }
    }

    # --- Point color scale + new scale for lines ---
    if (!is.null(group.by)) {
      pi <- pi +
        ggplot2::scale_color_manual(values = point_cols, name = group.by) +
        ggnewscale::new_scale_color()
    }

    # --- Ribbon ---
    if (isTRUE(add_interval)) {
      if (!is.null(line_color_by)) {
        pi <- pi + ggplot2::geom_ribbon(
          data = fit_sub,
          ggplot2::aes(x = .data$pseudotime, ymin = .data$lwr, ymax = .data$upr,
                       fill = .data[[line_color_by]],
                       group = interaction(.data$Lineage, .data$Feature)),
          alpha = 0.2
        )
      } else {
        pi <- pi + ggplot2::geom_ribbon(
          data = fit_sub,
          ggplot2::aes(x = .data$pseudotime, ymin = .data$lwr, ymax = .data$upr),
          fill = default_line_col, alpha = 0.2
        )
      }
    }

    # --- Line ---
    if (isTRUE(add_line)) {
      if (!is.null(line_color_by)) {
        pi <- pi + ggplot2::geom_line(
          data = fit_sub,
          ggplot2::aes(x = .data$pseudotime, y = .data$fitted,
                       color = .data[[line_color_by]],
                       group = interaction(.data$Lineage, .data$Feature)),
          linewidth = line.size
        )
      } else {
        pi <- pi + ggplot2::geom_line(
          data = fit_sub,
          ggplot2::aes(x = .data$pseudotime, y = .data$fitted),
          color = default_line_col, linewidth = line.size
        )
      }
    }

    # --- Line/fill scale ---
    if (!is.null(line_color_by)) {
      pi <- pi +
        ggplot2::scale_color_manual(values = line_cols, name = line_color_by) +
        ggplot2::scale_fill_manual(values = line_cols, name = line_color_by)
    }

    # --- Stat annotation (ggpubr stat_cor style) ---
    if (nrow(stat_sub) > 0) {
      if (isTRUE(compare_lineages) && length(pseudotime) > 1) {
        # Multiple lineages: stack labels
        lab_lines <- vapply(seq_len(nrow(stat_sub)), function(j) {
          paste0(stat_sub$Lineage[j], ": ",
                 .format_cor_label(stat_sub$pvalue[j], stat_sub$R2[j],
                                   stat_sub$rho[j], stat_method))
        }, character(1))
        lab_expr <- paste(lab_lines, collapse = "\n")
        pi <- pi + ggplot2::annotate(
          "text",
          x = min(raw_sub$Pseudotime, na.rm = TRUE),
          y = max(raw_sub$Expression, na.rm = TRUE),
          label = lab_expr, hjust = 0, vjust = 1, size = 3.5,
          color = "black"
        )
      } else {
        lab_expr <- .format_cor_label(stat_sub$pvalue[1], stat_sub$R2[1],
                                      stat_sub$rho[1], stat_method)
        pi <- pi + ggplot2::annotate(
          "text",
          x = min(raw_sub$Pseudotime, na.rm = TRUE),
          y = max(raw_sub$Expression, na.rm = TRUE),
          label = parse(text = lab_expr), hjust = 0, vjust = 1, size = 3.5,
          color = "black"
        )
      }
    }

    # --- Axes ---
    x_trans <- ifelse(isTRUE(reverse), "reverse", "identity")
    pi <- pi +
      ggplot2::scale_x_continuous(trans = x_trans,
                                  expand = ggplot2::expansion(c(0.02, 0.02))) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0.05, 0.05)))

    # --- Theme (matches PlotScatter) ---
    pi <- pi +
      ggplot2::labs(x = "Pseudotime", y = y_label, title = feat) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title       = ggplot2::element_text(size = 12),
        axis.text        = ggplot2::element_text(size = 10, color = "black"),
        strip.text       = ggplot2::element_text(size = 10),
        panel.grid.minor = ggplot2::element_blank(),
        legend.title     = ggplot2::element_text(face = "bold", size = 11),
        legend.text      = ggplot2::element_text(size = 10)
      )

    if (isTRUE(flip)) pi <- pi + ggplot2::coord_flip()

    plot_list[[feat]] <- pi
  }

  # --- Combine with patchwork ---
  if (length(plot_list) == 1) return(plot_list[[1]])

  pw <- patchwork::wrap_plots(plot_list, ncol = ncol, nrow = nrow)
  pw
}
