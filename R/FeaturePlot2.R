# ============================================================================
# FeaturePlot2.R -- Feature expression on 2-D reduction plot
# ============================================================================

# ---------------------------------------------------------------------------
# Internal colour helpers
# ---------------------------------------------------------------------------

#' Lighten a colour toward white
#' @param color Character vector of colour(s).
#' @param factor Numeric 0–1. 0 = white, 1 = original colour.
#' @noRd
.fp2_lighten <- function(color, factor) {
  rgb_mat <- grDevices::col2rgb(color) / 255
  rgb_new <- rgb_mat + (1 - rgb_mat) * (1 - factor)
  grDevices::rgb(rgb_new[1, ], rgb_new[2, ], rgb_new[3, ])
}

#' Blend two hex colours using one photographic mode
#' @noRd
.fp2_blend_two <- function(c1, c2, mode) {
  m1 <- as.numeric(grDevices::col2rgb(c1, alpha = FALSE)) / 255
  m2 <- as.numeric(grDevices::col2rgb(c2, alpha = FALSE)) / 255
  out <- switch(mode,
    blend    = (m1 + m2) / 2,
    average  = (m1 + m2) / 2,
    screen   = 1 - (1 - m1) * (1 - m2),
    multiply = m1 * m2
  )
  out <- pmax(0, pmin(1, out))
  grDevices::rgb(out[1], out[2], out[3])
}

#' Blend a vector of per-cell colours across features
#' @param colors_mat Character matrix (n_cells × n_features); NA = background.
#' @param mode One of "blend", "average", "screen", "multiply".
#' @noRd
.fp2_blendcolors <- function(colors_mat, mode = "blend") {
  apply(colors_mat, 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) == 0) return(NA_character_)
    if (length(row) == 1) return(row)
    result <- row[1]
    for (i in seq(2, length(row))) {
      result <- .fp2_blend_two(result, row[i], mode)
    }
    result
  })
}

#' Map numeric feature values to hex colours via a two-colour gradient
#' @param values Numeric vector (may contain NAs).
#' @param lo_color Hex colour for low end.
#' @param hi_color Hex colour for high end.
#' @param lo_val,hi_val Numeric range for clipping.
#' @param bg_color Colour returned for NA values.
#' @noRd
.fp2_map_to_color <- function(values, lo_color, hi_color,
                              lo_val, hi_val, bg_color = "grey80") {
  result <- rep(NA_character_, length(values))
  valid  <- !is.na(values)
  if (!any(valid)) return(rep(bg_color, length(values)))

  v <- pmax(lo_val, pmin(hi_val, values[valid]))
  if (hi_val == lo_val) {
    frac <- rep(1, sum(valid))
  } else {
    frac <- (v - lo_val) / (hi_val - lo_val)
  }
  lo_rgb <- grDevices::col2rgb(lo_color) / 255
  hi_rgb <- grDevices::col2rgb(hi_color) / 255
  blended_r <- lo_rgb[1] + frac * (hi_rgb[1] - lo_rgb[1])
  blended_g <- lo_rgb[2] + frac * (hi_rgb[2] - lo_rgb[2])
  blended_b <- lo_rgb[3] + frac * (hi_rgb[3] - lo_rgb[3])
  result[valid] <- grDevices::rgb(blended_r, blended_g, blended_b)
  result[!valid] <- bg_color
  result
}


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

#' @title Feature Expression on 2-D Reduction Plot
#'
#' @description
#' Visualise feature values (gene expression, metadata scores, reduction
#' embeddings) on a 2-D reduction plot.  Accepts Seurat objects or plain
#' data.frames.  Key options include:
#' \itemize{
#'   \item **Multiple features / split panels** with unified colour scale
#'     (\code{keep_scale}).
#'   \item **Co-expression** geometric mean across gene features
#'     (\code{calculate_coexp = TRUE}).
#'   \item **Compare mode** – blend all features per cell into a single colour
#'     (\code{compare_features = TRUE}).
#' }
#'
#' @md
#' @param seu A Seurat object or a \code{data.frame}.  For data.frames the
#'   columns specified by \code{dims} are used as x/y coordinates; all feature
#'   columns must also be present in the data.frame.
#' @param features Character vector (or named character vector) of features to
#'   plot.  Features can be gene names (Seurat assay rows), metadata columns,
#'   reduction embedding columns, or data.frame columns.  Named features use
#'   the name as subtitle.  A list is also accepted and is flattened.
#' @param reduction Character. Reduction to use for coordinates.  \code{NULL}
#'   auto-detects UMAP/tSNE.  Ignored for data.frame input.
#' @param dims Integer vector of length 2.  Dimension indices within the
#'   reduction (Seurat) or column indices (data.frame).  Default \code{c(1,2)}.
#' @param split.by Column name to split into facets.  \code{NULL} = no split.
#' @param cells Character vector of cell/row names to include.  \code{NULL} =
#'   all.
#' @param layer Seurat assay layer to use.  Default \code{"data"}.
#' @param assay Seurat assay name.  \code{NULL} = default assay.
#' @param show_stat Logical.  Show positive-cell statistics in subtitle.
#'   \code{NULL} = auto (\code{FALSE} for blank theme).
#' @param palette Colour palette name for the gradient.  Default
#'   \code{"Spectral"} (standard) or \code{"Set1"} (compare mode).
#' @param palcolor Custom colour vector.  \code{NULL} = use palette.
#' @param pt.size Point size.  \code{NULL} = auto.
#' @param pt.alpha Point transparency.  Default \code{1}.
#' @param bg_cutoff Numeric.  Feature values ≤ cutoff are treated as
#'   background (shown in \code{bg_color}).  Default \code{0}.
#' @param bg_color Background colour for low/NA values.  Default
#'   \code{"grey80"}.
#' @param keep_scale How to align colour scales across panels:
#'   \itemize{
#'     \item \code{NULL} – each panel scaled independently to its own range.
#'     \item \code{"feature"} (default) – panels for the same feature share a
#'       scale.
#'     \item \code{"all"} – all panels share a single global scale.
#'   }
#' @param lower_quantile,upper_quantile Quantile bounds for colour scale.
#'   Default \code{0} / \code{0.99}.
#' @param lower_cutoff,upper_cutoff Explicit value bounds (override quantile).
#'   \code{NULL} = use quantiles.
#' @param add_density Logical.  Overlay 2-D density contour.  Default
#'   \code{FALSE}.
#' @param density_color Contour colour.  Default \code{"grey80"}.
#' @param density_filled Logical.  Filled density raster.  Default
#'   \code{FALSE}.
#' @param density_filled_palette Palette for filled density.  Default
#'   \code{"Greys"}.
#' @param density_filled_palcolor Custom colours for filled density.
#'   \code{NULL}.
#' @param cells.highlight Cells to highlight (names or \code{TRUE} for all
#'   expressed cells).  \code{NULL} = none.
#' @param cols.highlight Border colour of highlighted cells.  Default
#'   \code{"black"}.
#' @param sizes.highlight Size of highlighted points.  Default \code{1}.
#' @param alpha.highlight Alpha of highlighted points.  Default \code{1}.
#' @param stroke.highlight Border width around highlighted points.  Default
#'   \code{0.5}.
#' @param calculate_coexp Logical.  Replace gene features with their geometric
#'   mean co-expression value.  Default \code{FALSE}.
#' @param compare_features Logical.  Blend all features per cell into one
#'   colour.  Default \code{FALSE}.
#' @param color_blend_mode Blend mode when \code{compare_features = TRUE}.
#'   One of \code{"blend"}, \code{"average"}, \code{"screen"},
#'   \code{"multiply"}.
#' @param label Logical.  Label the high-expression region.  Default
#'   \code{FALSE}.
#' @param label.size Label text size.  Default \code{4}.
#' @param label.fg Label foreground colour.  Default \code{"white"}.
#' @param label.bg Label background colour.  Default \code{"black"}.
#' @param label.bg.r Label background ratio.  Default \code{0.1}.
#' @param label_insitu Logical.  Use feature names as labels (insitu) instead
#'   of numbers.  Default \code{FALSE}.
#' @param label_repel Logical.  Use \code{ggrepel} for labels.  Default
#'   \code{FALSE}.
#' @param label_repulsion Repulsion force.  Default \code{20}.
#' @param label_point_size Center point size for repelled labels.  Default
#'   \code{1}.
#' @param label_point_color Center point colour.  Default \code{"black"}.
#' @param label_segment_color Segment colour.  Default \code{"black"}.
#' @param lineages Character vector of metadata column names containing
#'   pseudotime values (e.g. from \code{RunSlingshot}).  \code{NULL} = none.
#' @param lineages_trim Quantile range for LOESS fit. Default \code{c(0.01,0.99)}.
#' @param lineages_span LOESS span. Default \code{0.75}.
#' @param lineages_palette Palette for lineage colours. Default \code{"Dark2"}.
#' @param lineages_palcolor Custom colour vector. Default \code{NULL}.
#' @param lineages_arrow Arrow for path end. Default
#'   \code{grid::arrow(length = grid::unit(0.1, "inches"))}.
#' @param lineages_linewidth Foreground linewidth. Default \code{1}.
#' @param lineages_line_bg Outline colour. Default \code{"white"}.
#' @param lineages_line_bg_stroke Outline extra width. Default \code{0.5}.
#' @param lineages_whiskers Draw cell-to-curve segments. Default \code{FALSE}.
#' @param lineages_whiskers_linewidth Whisker linewidth. Default \code{0.5}.
#' @param lineages_whiskers_alpha Whisker alpha. Default \code{0.5}.
#' @param raster Logical.  Rasterise points.  \code{NULL} = auto (> 100k).
#' @param raster_method Character. Rasterisation backend:
#'   \code{"rasterise"} (default, faithful colours via \code{rasterise_layer})
#'   or \code{"scattermore"} (faster C-level rendering, may shift colours).
#' @param raster.dpi Numeric scalar or 2-length vector.  Default
#'   \code{c(512,512)}.
#' @param aspect.ratio Aspect ratio.  Default \code{1}.
#' @param title Plot title.  \code{NULL} = none.
#' @param subtitle Subtitle(s).  \code{NULL} = auto from feature names.
#'   Single value or vector of \code{length(features)}.
#' @param xlab,ylab Axis labels.  \code{NULL} = auto.
#' @param legend.position Legend position.  Default \code{"right"}.
#' @param legend.direction Legend direction.  Default \code{"vertical"}.
#' @param legend.title Legend title.  \code{NULL} = feature name.
#' @param theme_use Theme specification (function, character, or
#'   \code{theme} object).  \code{NULL} → \code{UtilsR::theme_blank}.
#' @param theme_args Extra arguments for \code{theme_use}.  Default
#'   \code{list()}.
#' @param combine Logical.  Combine panels with patchwork.  Default
#'   \code{TRUE}.
#' @param nrow,ncol Layout dimensions.  \code{NULL} = auto.
#' @param byrow Logical.  Fill by row.  Default \code{TRUE}.
#' @param force Logical.  Force plot when > 50 features.  Default
#'   \code{FALSE}.
#' @param seed Random seed.  Default \code{11}.
#'
#' @return A ggplot or patchwork object (or list when
#'   \code{combine = FALSE}).
#'
#' @examples
#' \dontrun{
#' # Seurat input
#' FeaturePlot2(seu, features = "CD3D")
#' FeaturePlot2(seu, features = c("CD3D", "CD8A"), keep_scale = "all")
#' FeaturePlot2(seu, features = c("CD3D", "CD8A"),
#'              compare_features = TRUE, color_blend_mode = "blend")
#'
#' # data.frame input
#' df <- data.frame(UMAP1 = rnorm(500), UMAP2 = rnorm(500),
#'                  GeneA = abs(rnorm(500)), GeneB = abs(rnorm(500)))
#' FeaturePlot2(df, features = "GeneA", dims = c(1, 2))
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn
#'   scale_color_identity scale_x_continuous scale_y_continuous
#'   labs theme guide_colorbar facet_grid stat_density_2d
#'   scale_fill_gradientn geom_density_2d
#' @export
FeaturePlot2 <- function(
    seu,
    features,
    reduction       = NULL,
    dims            = c(1, 2),
    split.by        = NULL,
    cells           = NULL,
    layer           = "data",
    assay           = NULL,
    show_stat       = NULL,
    palette         = NULL,
    palcolor        = NULL,
    pt.size         = NULL,
    pt.alpha        = 1,
    bg_cutoff       = 0,
    bg_color        = "grey80",
    keep_scale      = "feature",
    lower_quantile  = 0,
    upper_quantile  = 0.99,
    lower_cutoff    = NULL,
    upper_cutoff    = NULL,
    add_density     = FALSE,
    density_color   = "grey80",
    density_filled  = FALSE,
    density_filled_palette  = "Greys",
    density_filled_palcolor = NULL,
    cells.highlight = NULL,
    cols.highlight  = "black",
    sizes.highlight = 1,
    alpha.highlight = 1,
    stroke.highlight = 0.5,
    calculate_coexp = FALSE,
    compare_features = FALSE,
    color_blend_mode = c("blend", "average", "screen", "multiply"),
    label           = FALSE,
    label.size      = 4,
    label.fg        = "white",
    label.bg        = "black",
    label.bg.r      = 0.1,
    label_insitu    = FALSE,
    label_repel     = FALSE,
    label_repulsion = 20,
    label_point_size = 1,
    label_point_color = "black",
    label_segment_color = "black",
    lineages                    = NULL,
    lineages_trim               = c(0.01, 0.99),
    lineages_span               = 0.75,
    lineages_palette            = "Dark2",
    lineages_palcolor           = NULL,
    lineages_arrow              = grid::arrow(length = grid::unit(0.1, "inches")),
    lineages_linewidth          = 1,
    lineages_line_bg            = "white",
    lineages_line_bg_stroke     = 0.5,
    lineages_whiskers           = FALSE,
    lineages_whiskers_linewidth = 0.5,
    lineages_whiskers_alpha     = 0.5,
    raster          = NULL,
    raster_method   = c("scattermore", "rasterise"),
    raster.dpi      = c(512, 512),
    aspect.ratio    = 1,
    title           = NULL,
    subtitle        = NULL,
    xlab            = NULL,
    ylab            = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    legend.title    = NULL,
    theme_use       = NULL,
    theme_args      = list(),
    combine         = TRUE,
    nrow            = NULL,
    ncol            = NULL,
    byrow           = TRUE,
    force           = FALSE,
    seed            = 11) {

  set.seed(seed)
  color_blend_mode <- match.arg(color_blend_mode)
  if (!is.null(keep_scale)) {
    keep_scale <- match.arg(keep_scale, c("feature", "all"))
  }

  # ---- Normalise features argument ----
  feature_input <- features
  if (is.list(features)) {
    if (is.null(names(features))) {
      features <- unlist(features)
    } else {
      features <- stats::setNames(
        unlist(features),
        nm = rep(names(features), vapply(features, length, 1L))
      )
    }
  }
  if (!is.character(features)) stop("'features' must be a character vector.", call. = FALSE)

  # ---- Resolve palette default ----
  if (is.null(palette)) {
    palette <- if (isTRUE(compare_features)) "Set1" else "Spectral"
  }

  # =========================================================================
  # 1. Extract coordinates and metadata
  # =========================================================================
  if (inherits(seu, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Package 'SeuratObject' required for Seurat input.", call. = FALSE)
    }
    assay <- assay %||% SeuratObject::DefaultAssay(seu)

    # Resolve reduction
    if (is.null(reduction)) {
      red_names <- names(seu@reductions)
      umap_idx  <- grep("umap|tsne", tolower(red_names))
      reduction <- if (length(umap_idx) > 0) red_names[umap_idx[1]] else tail(red_names, 1)
    }
    if (!reduction %in% names(seu@reductions)) {
      stop("Reduction '", reduction, "' not found. Available: ",
           paste(names(seu@reductions), collapse = ", "), call. = FALSE)
    }
    reduction_key <- seu@reductions[[reduction]]@key
    dat_dim       <- seu@reductions[[reduction]]@cell.embeddings
    colnames(dat_dim) <- paste0(reduction_key, seq_len(ncol(dat_dim)))
    rownames(dat_dim) <- rownames(dat_dim) %||% colnames(seu)
    dat_meta <- seu@meta.data

    # Feature sources
    emb_cols <- unlist(lapply(seu@reductions,
                               function(r) colnames(r@cell.embeddings)))
    valid_features <- c(rownames(seu@assays[[assay]]),
                        colnames(dat_meta), emb_cols)

  } else if (is.data.frame(seu)) {
    reduction_key <- "Dim"
    reduction     <- "custom"
    if (any(dims > ncol(seu))) {
      stop("'dims' indices out of range for the data.frame (ncol = ",
           ncol(seu), ").", call. = FALSE)
    }
    dat_dim_mat <- as.matrix(seu[, dims, drop = FALSE])
    colnames(dat_dim_mat) <- paste0(reduction_key, seq_len(2))
    rownames(dat_dim_mat) <- rownames(seu)
    dat_dim  <- dat_dim_mat
    dat_meta <- seu
    assay    <- NULL
    emb_cols <- character(0)
    valid_features <- colnames(seu)
  } else {
    stop("'seu' must be a Seurat object or a data.frame.", call. = FALSE)
  }

  dim_x <- paste0(reduction_key, dims[1])
  dim_y <- paste0(reduction_key, dims[2])

  # ---- Validate / filter features ----
  features_unique  <- unique(features)
  features_invalid <- features_unique[!features_unique %in% valid_features]
  if (length(features_invalid) > 0) {
    warning("Features not found and dropped: ",
            paste(features_invalid, collapse = ", "), call. = FALSE)
    keep <- features %in% valid_features
    features      <- features[keep]
    feature_input <- feature_input[keep]
  }
  if (length(features) == 0) stop("No valid features remaining.", call. = FALSE)

  # ---- Force check ----
  if (length(features) > 50 && isFALSE(force)) {
    warning("More than 50 features requested. Use force = TRUE to proceed.",
            call. = FALSE)
    return(invisible(NULL))
  }

  # =========================================================================
  # 2. Extract feature data
  # =========================================================================
  if (inherits(seu, "Seurat")) {
    features_gene  <- features[features %in% rownames(seu@assays[[assay]])]
    features_meta  <- features[features %in% colnames(dat_meta)]
    features_emb   <- features[features %in% emb_cols]

    # co-expression: geometric mean across gene features
    if (isTRUE(calculate_coexp) && length(features_gene) > 0) {
      if (length(features_meta) > 0) {
        warning("Metadata features ignored when calculate_coexp = TRUE.",
                call. = FALSE)
      }
      mat_gene <- as.matrix(
        SeuratObject::GetAssayData(seu, assay = assay, layer = layer)[
          features_gene, , drop = FALSE]
      )
      coexp <- apply(mat_gene, 2, function(x) exp(mean(log(x + 1e-9))))
      dat_meta[["CoExp"]] <- coexp
      features      <- c("CoExp")
      features_meta <- "CoExp"
      features_gene <- features_emb <- character(0)
    }

    fetch_vars <- unique(c(features_gene, features_meta, features_emb))
    dat_exp <- as.matrix(
      SeuratObject::FetchData(seu, vars = fetch_vars, layer = layer)
    )
  } else {
    dat_exp <- as.matrix(seu[, features, drop = FALSE])
  }

  # Keep only valid rows (intersection with dat_dim)
  shared <- intersect(rownames(dat_dim), rownames(dat_exp))
  dat_dim <- dat_dim[shared, , drop = FALSE]
  dat_exp <- dat_exp[shared, features, drop = FALSE]
  dat_meta <- dat_meta[shared, , drop = FALSE]

  # Apply background cutoff
  dat_exp[dat_exp <= bg_cutoff] <- NA

  # =========================================================================
  # 3. Split column
  # =========================================================================
  if (is.null(split.by)) {
    split.by <- ".split_all"
    dat_meta[[split.by]] <- factor("")
  } else {
    if (!split.by %in% colnames(dat_meta)) {
      stop("'split.by' column '", split.by, "' not found.", call. = FALSE)
    }
    if (!is.factor(dat_meta[[split.by]])) {
      dat_meta[[split.by]] <- factor(dat_meta[[split.by]],
                                     levels = unique(dat_meta[[split.by]]))
    }
  }
  split_use <- split.by

  # =========================================================================
  # 4. Combine into working data frame
  # =========================================================================
  dat_use <- cbind(
    as.data.frame(dat_dim),
    dat_meta[, split_use, drop = FALSE],
    dat_exp
  )
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }

  # =========================================================================
  # 4b. Lineage layers (computed once from full coordinates)
  # =========================================================================
  lineage_layers <- NULL
  if (!is.null(lineages)) {
    if (!requireNamespace("ggnewscale", quietly = TRUE)) {
      stop("Package 'ggnewscale' is required for lineages.", call. = FALSE)
    }
    missing_lin <- setdiff(lineages, colnames(dat_meta))
    if (length(missing_lin) > 0) {
      stop("Lineage column(s) not found in metadata: ",
           paste(missing_lin, collapse = ", "), call. = FALSE)
    }
    lin_dat <- as.data.frame(dat_dim)
    for (l in lineages) lin_dat[[l]] <- dat_meta[rownames(lin_dat), l]
    if (!is.null(cells)) {
      lin_dat <- lin_dat[intersect(rownames(lin_dat), cells), , drop = FALSE]
    }
    lineage_layers <- c(
      list(ggnewscale::new_scale_color()),
      .compute_lineage_layers(
        dat                = lin_dat,
        x_col              = dim_x,
        y_col              = dim_y,
        lineages           = lineages,
        trim               = lineages_trim,
        span               = lineages_span,
        palette            = lineages_palette,
        palcolor           = lineages_palcolor,
        lineages_arrow     = lineages_arrow,
        linewidth          = lineages_linewidth,
        line_bg            = lineages_line_bg,
        line_bg_stroke     = lineages_line_bg_stroke,
        whiskers           = lineages_whiskers,
        whiskers_linewidth = lineages_whiskers_linewidth,
        whiskers_alpha     = lineages_whiskers_alpha
      )
    )
  }

  # =========================================================================
  # 5. Housekeeping
  # =========================================================================
  if (is.null(pt.size)) pt.size <- min(3000 / nrow(dat_use), 0.5)

  raster        <- raster %||% (nrow(dat_use) > 1e5)
  raster_method <- match.arg(raster_method)
  if (!is.null(raster.dpi) && length(raster.dpi) == 2) raster.dpi <- max(raster.dpi)
  raster.dpi <- raster.dpi %||% 300

  # ---- Theme ----
  if (is.null(theme_use)) {
    if (requireNamespace("UtilsR", quietly = TRUE)) {
      theme_fn <- UtilsR::theme_blank; is_blank <- TRUE
    } else {
      theme_fn <- ggplot2::theme_bw;   is_blank <- FALSE
    }
  } else if (is.character(theme_use)) {
    theme_fn <- tryCatch(match.fun(theme_use), error = function(e) {
      if (requireNamespace("UtilsR", quietly = TRUE))
        tryCatch(utils::getFromNamespace(theme_use, "UtilsR"),
                 error = function(e2) ggplot2::theme_bw)
      else ggplot2::theme_bw
    })
    is_blank <- identical(theme_use, "theme_blank")
  } else if (is.function(theme_use)) {
    theme_fn <- theme_use
    is_blank  <- isTRUE(attr(theme_use, "blank")) ||
      (requireNamespace("UtilsR", quietly = TRUE) &&
         identical(theme_fn, UtilsR::theme_blank))
  } else if (inherits(theme_use, "theme")) {
    theme_fn <- function(...) theme_use; is_blank <- FALSE
  } else {
    stop("'theme_use' must be NULL, character, function, or theme object.", call. = FALSE)
  }
  if (is.null(show_stat)) show_stat <- !is_blank

  # ---- Axis labels ----
  xlab <- xlab %||% dim_x
  ylab <- ylab %||% dim_y
  if (isTRUE(is_blank)) {
    theme_args[["xlab"]] <- xlab
    theme_args[["ylab"]] <- ylab
  }

  # ---- Subtitle lookup ----
  # Build named vector: feature_value -> subtitle_text
  if (is.null(subtitle)) {
    if (!is.null(names(feature_input))) {
      sub_lookup <- stats::setNames(names(feature_input), nm = feature_input)
    } else {
      sub_lookup <- stats::setNames(rep(NA_character_, length(features)),
                                    nm = features)
    }
  } else {
    if (length(subtitle) == 1) {
      sub_lookup <- stats::setNames(rep(subtitle, length(features)), nm = features)
    } else if (length(subtitle) == length(features)) {
      sub_lookup <- stats::setNames(subtitle, nm = features)
    } else {
      stop("'subtitle' length must be 1 or length(features).", call. = FALSE)
    }
  }

  # ---- Axis limits (from full embedding) ----
  xlim_range <- range(dat_dim[, 1], na.rm = TRUE)
  ylim_range <- range(dat_dim[, 2], na.rm = TRUE)

  # ---- Point geom helper ----
  # raster_method = "rasterise" : faithful colours (default)
  # raster_method = "scattermore": faster but may shift colours
  .pt <- function(data, mapping, size_val, alpha_val, color = NULL, ...) {
    if (isTRUE(raster) && raster_method == "scattermore") {
      if (!requireNamespace("scattermore", quietly = TRUE))
        stop("Package 'scattermore' is required for raster_method = 'scattermore'.",
             call. = FALSE)
      args <- list(data = data, mapping = mapping,
                   pointsize = ceiling(size_val), alpha = alpha_val,
                   pixels = c(raster.dpi, raster.dpi), ...)
      if (!is.null(color)) args$color <- color
      return(do.call(scattermore::geom_scattermore, args))
    }
    args <- list(data = data, mapping = mapping, size = size_val,
                 alpha = alpha_val, ...)
    if (!is.null(color)) args$color <- color
    layer <- do.call(ggplot2::geom_point, args)
    if (isTRUE(raster)) rasterise_layer(layer, dpi = raster.dpi) else layer
  }

  # ---- Density layer helper ----
  .density_layer <- function(dat_sub) {
    if (!isTRUE(add_density)) return(NULL)
    if (isTRUE(density_filled)) {
      if (!requireNamespace("ggnewscale", quietly = TRUE)) {
        warning("'ggnewscale' required for density_filled; skipping.", call. = FALSE)
        return(NULL)
      }
      filled_cols <- palette_colors(palette  = density_filled_palette,
                                    palcolor = density_filled_palcolor)
      list(
        ggplot2::stat_density_2d(
          data = dat_sub,
          ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                       fill = ggplot2::after_stat(density)),
          geom = "raster", contour = FALSE,
          inherit.aes = FALSE, show.legend = FALSE),
        ggplot2::scale_fill_gradientn(name = "Density", colours = filled_cols),
        ggnewscale::new_scale_fill()
      )
    } else {
      ggplot2::geom_density_2d(
        data = dat_sub,
        ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
        color = density_color, inherit.aes = FALSE, show.legend = FALSE)
    }
  }

  # ---- Label layer helper ----
  .label_layer <- function(p, dat_sub, feature_name, label_colors = NULL) {
    if (!isTRUE(label)) return(p)
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      warning("'ggrepel' required for label = TRUE; skipping.", call. = FALSE)
      return(p)
    }
    vals <- dat_sub[["value"]]
    if (all(is.na(vals))) return(p)
    q_lo <- stats::quantile(vals[is.finite(vals)], 0.95, na.rm = TRUE)
    q_hi <- stats::quantile(vals[is.finite(vals)], 0.99, na.rm = TRUE)
    lab_df <- dat_sub[is.finite(vals) & vals >= q_lo & vals <= q_hi, , drop = FALSE]
    if (nrow(lab_df) == 0) return(p)
    lab_summary <- data.frame(
      label = feature_name,
      x     = stats::median(lab_df[["x"]]),
      y     = stats::median(lab_df[["y"]])
    )
    if (!isTRUE(label_insitu)) lab_summary[["label"]] <- "★"

    if (isTRUE(label_repel)) {
      p <- p +
        ggplot2::geom_point(
          data = lab_summary,
          ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
          color = label_point_color, size = label_point_size,
          inherit.aes = FALSE) +
        ggrepel::geom_text_repel(
          data = lab_summary,
          ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                       label = .data[["label"]]),
          fontface = "bold", min.segment.length = 0,
          segment.color = label_segment_color,
          point.size = label_point_size, max.overlaps = 100,
          force = label_repulsion,
          color = label.fg, bg.color = label.bg, bg.r = label.bg.r,
          size = label.size, inherit.aes = FALSE, show.legend = FALSE)
    } else {
      p <- p +
        ggrepel::geom_text_repel(
          data = lab_summary,
          ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                       label = .data[["label"]]),
          fontface = "bold", min.segment.length = 0,
          segment.color = label_segment_color,
          point.size = NA, max.overlaps = 100, force = 0,
          color = label.fg, bg.color = label.bg, bg.r = label.bg.r,
          size = label.size, inherit.aes = FALSE, show.legend = FALSE)
    }
    p
  }

  # ---- Gradient colours (for standard path) ----
  colors_gradient <- palette_colors(palette = palette, palcolor = palcolor)

  # =========================================================================
  # 6A. COMPARE path
  # =========================================================================
  if (isTRUE(compare_features) && length(features) > 1) {

    # Discrete colour per feature (saturated endpoint)
    feat_colors <- palette_colors(features, palette = palette, palcolor = palcolor)

    # Split data by split column
    dat_all <- cbind(dat_use, x = dat_use[[dim_x]], y = dat_use[[dim_y]])
    dat_all[["x"]] <- dat_all[[dim_x]]
    dat_all[["y"]] <- dat_all[[dim_y]]
    dat_split <- split(dat_all, dat_all[[split_use]])

    plist <- lapply(
      stats::setNames(levels(dat_use[[split_use]]),
                      levels(dat_use[[split_use]])),
      function(s) {
        dat <- dat_split[[s]]
        if (is.null(dat) || nrow(dat) == 0) dat <- dat_all

        # Per-feature, per-cell colour mapping
        colors_mat <- do.call(cbind, lapply(seq_along(features), function(i) {
          f        <- features[i]
          hi_col   <- feat_colors[i]
          lo_col   <- .fp2_lighten(hi_col, 0.1)
          vals     <- dat[[f]]
          lo_val   <- lower_cutoff %||%
            stats::quantile(vals[is.finite(vals)], lower_quantile, na.rm = TRUE)
          hi_val   <- upper_cutoff %||%
            stats::quantile(vals[is.finite(vals)], upper_quantile, na.rm = TRUE)
          if (is.na(lo_val) || is.na(hi_val)) return(rep(NA_character_, nrow(dat)))
          .fp2_map_to_color(vals, lo_col, hi_col, lo_val, hi_val, bg_color = NA)
        }))
        colnames(colors_mat) <- features

        # Blend per cell
        dat[["color_blend"]] <- .fp2_blendcolors(colors_mat, mode = color_blend_mode)
        all_na <- rowSums(!is.na(colors_mat)) == 0
        dat[["color_blend"]][all_na] <- bg_color

        # Compute brightness for ordering (bg-colour cells first)
        is_bg  <- dat[["color_blend"]] == bg_color | is.na(dat[["color_blend"]])
        dat[["color_blend"]][is.na(dat[["color_blend"]])] <- bg_color
        dat <- dat[order(is_bg, decreasing = TRUE), , drop = FALSE]

        # highlight
        cells.highlight_use <- cells.highlight
        if (isTRUE(cells.highlight_use)) {
          cells.highlight_use <- rownames(dat)[!is_bg[order(is_bg, decreasing = TRUE)]]
        }

        subtitle_use <- paste(sub_lookup[features], collapse = "|")
        if (all(is.na(sub_lookup[features]))) subtitle_use <- NULL

        p <- ggplot2::ggplot(dat) +
          .density_layer(dat[!is_bg, , drop = FALSE]) +
          .pt(dat, ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                color = .data[["color_blend"]]),
              size_val = pt.size, alpha_val = pt.alpha) +
          ggplot2::scale_color_identity() +
          ggplot2::scale_x_continuous(limits = xlim_range) +
          ggplot2::scale_y_continuous(limits = ylim_range) +
          ggplot2::labs(title = title, subtitle = subtitle_use,
                        x = xlab, y = ylab) +
          do.call(theme_fn, theme_args) +
          ggplot2::theme(aspect.ratio = aspect.ratio,
                         legend.position = "none")

        # Highlight layer
        if (!is.null(cells.highlight_use) && length(cells.highlight_use) > 0) {
          cell_df <- dat[rownames(dat) %in% cells.highlight_use, , drop = FALSE]
          if (nrow(cell_df) > 0) {
            p <- p +
              .pt(cell_df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                  size_val = sizes.highlight + stroke.highlight,
                  alpha_val = alpha.highlight, color = cols.highlight) +
              .pt(cell_df,
                  ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                               color = .data[["color_blend"]]),
                  size_val = sizes.highlight, alpha_val = alpha.highlight) +
              ggplot2::scale_color_identity()
          }
        }

        # == Lineage overlay (compare path) ==
        if (!is.null(lineage_layers)) {
          p <- p + lineage_layers
        }

        # Label layer (one label per feature)
        for (f in features) {
          dat[["value"]] <- dat[[f]]
          p <- .label_layer(p, dat, f)
        }

        # ---- Build per-feature legend plots then combine ----
        legend_plots <- lapply(seq_along(features), function(i) {
          f      <- features[i]
          hi_col <- feat_colors[i]
          lo_col <- .fp2_lighten(hi_col, 0.1)
          vals   <- dat_all[[f]]
          lo_val <- lower_cutoff %||%
            stats::quantile(vals[is.finite(vals)], lower_quantile, na.rm = TRUE)
          hi_val <- upper_cutoff %||%
            stats::quantile(vals[is.finite(vals)], upper_quantile, na.rm = TRUE)
          if (is.na(lo_val) || is.na(hi_val)) {
            lo_val <- 0; hi_val <- 1
          }
          leg_title <- legend.title %||% f
          dummy_df  <- data.frame(x = 0, y = 0, v = seq(lo_val, hi_val, length.out = 10))
          ggplot2::ggplot(dummy_df, ggplot2::aes(x = .data[["x"]],
                                                  y = .data[["y"]],
                                                  color = .data[["v"]])) +
            ggplot2::geom_point(size = 0) +
            ggplot2::scale_color_gradientn(
              name     = leg_title,
              colours  = c(lo_col, hi_col),
              limits   = c(lo_val, hi_val),
              na.value = bg_color,
              guide    = ggplot2::guide_colorbar(
                frame.colour = "black", ticks.colour = "black",
                title.hjust = 0, order = i)
            ) +
            ggplot2::theme_void() +
            ggplot2::theme(
              legend.position = "right",
              legend.direction = legend.direction,
              legend.title = ggplot2::element_text(size = 9),
              legend.text  = ggplot2::element_text(size = 8)
            )
        })

        if (!requireNamespace("patchwork", quietly = TRUE)) {
          return(p)
        }
        if (legend.position == "none") return(p)
        legend_col <- patchwork::wrap_plots(legend_plots, ncol = 1) &
          patchwork::plot_layout(guides = "keep")
        p + patchwork::inset_element(legend_col,
                                     left = 1, bottom = 0, right = 1.3, top = 1,
                                     align_to = "plot")
      }
    )
    names(plist) <- paste0(
      levels(dat_use[[split_use]]), ":", paste(features, collapse = "|")
    )

  } else {
    # =========================================================================
    # 6B. STANDARD path (per feature × per split)
    # =========================================================================
    comb <- expand.grid(
      split   = levels(dat_use[[split_use]]),
      feature = features,
      stringsAsFactors = FALSE
    )
    rownames(comb) <- paste0(comb[["split"]], ":", comb[["feature"]])

    dat_split <- split(dat_use, dat_use[[split_use]])

    plist <- lapply(
      stats::setNames(rownames(comb), rownames(comb)),
      function(i) {
        f <- comb[i, "feature"]
        s <- comb[i, "split"]
        dat <- dat_split[[s]]
        if (is.null(dat) || nrow(dat) == 0) dat <- dat_use

        # Subset to needed columns + set x/y/value
        dat[["x"]]     <- dat[[dim_x]]
        dat[["y"]]     <- dat[[dim_y]]
        dat[["value"]] <- dat[[f]]

        # Sort: NA first (background), expressed cells sorted ascending (high on top)
        dat <- dat[order(dat[["value"]], method = "radix",
                         decreasing = FALSE, na.last = FALSE), , drop = FALSE]

        # Highlight
        cells.highlight_use <- cells.highlight
        if (isTRUE(cells.highlight_use)) {
          cells.highlight_use <- rownames(dat)[!is.na(dat[["value"]])]
        }

        # Color scale range
        all_vals <- dat_use[[f]]  # full data (before split) for keep_scale
        if (all(is.na(dat[["value"]]))) {
          lo_val <- 0; hi_val <- 1
        } else if (is.null(keep_scale)) {
          # Per-split range
          lo_val <- lower_cutoff %||%
            stats::quantile(dat[["value"]][is.finite(dat[["value"]])],
                            lower_quantile, na.rm = TRUE)
          hi_val <- upper_cutoff %||%
            stats::quantile(dat[["value"]][is.finite(dat[["value"]])],
                            upper_quantile, na.rm = TRUE) + 1e-9
        } else if (keep_scale == "feature") {
          lo_val <- lower_cutoff %||%
            stats::quantile(all_vals[is.finite(all_vals)],
                            lower_quantile, na.rm = TRUE)
          hi_val <- upper_cutoff %||%
            stats::quantile(all_vals[is.finite(all_vals)],
                            upper_quantile, na.rm = TRUE) + 1e-9
        } else { # "all"
          all_matrix <- as.vector(as.matrix(dat_use[, features, drop = FALSE]))
          lo_val <- lower_cutoff %||%
            stats::quantile(all_matrix[is.finite(all_matrix)],
                            lower_quantile, na.rm = TRUE)
          hi_val <- upper_cutoff %||%
            stats::quantile(all_matrix[is.finite(all_matrix)],
                            upper_quantile, na.rm = TRUE) + 1e-9
        }
        # Clamp
        dat[["value"]] <- pmax(lo_val, pmin(hi_val, dat[["value"]]))

        # Subtitle
        sub_name <- sub_lookup[[f]]
        if (is.na(sub_name)) sub_name <- NULL
        if (isTRUE(show_stat) && split_use != ".split_all") {
          n_pos   <- sum(dat[["value"]] > 0, na.rm = TRUE)
          pct_pos <- round(n_pos / nrow(dat) * 100, 2)
          subtitle_use <- sub_name %||%
            paste0(s, " nPos:", n_pos, ", ", pct_pos, "%")
        } else {
          subtitle_use <- sub_name
        }

        legend_title_use <- legend.title %||% f

        # == Base plot ==
        p <- ggplot2::ggplot(dat) +
          .density_layer(dat[!is.na(dat[["value"]]), , drop = FALSE]) +
          ggplot2::scale_x_continuous(limits = xlim_range) +
          ggplot2::scale_y_continuous(limits = ylim_range) +
          ggplot2::labs(title = title, subtitle = subtitle_use,
                        x = xlab, y = ylab) +
          do.call(theme_fn, theme_args) +
          ggplot2::theme(aspect.ratio = aspect.ratio,
                         legend.position = legend.position,
                         legend.direction = legend.direction)

        # == Point layer ==
        p <- p +
          .pt(dat,
              ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                           color = .data[["value"]]),
              size_val = pt.size, alpha_val = pt.alpha)

        # == Highlight ==
        if (!is.null(cells.highlight_use) && length(cells.highlight_use) > 0) {
          cell_df <- dat[rownames(dat) %in% cells.highlight_use, , drop = FALSE]
          if (nrow(cell_df) > 0) {
            p <- p +
              .pt(cell_df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
                  size_val = sizes.highlight + stroke.highlight,
                  alpha_val = alpha.highlight, color = cols.highlight) +
              .pt(cell_df,
                  ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                               color = .data[["value"]]),
                  size_val = sizes.highlight, alpha_val = alpha.highlight)
          }
        }

        # == Color scale ==
        if (all(is.na(dat[["value"]]))) {
          p <- p + ggplot2::scale_color_gradient(
            name     = legend_title_use,
            na.value = bg_color
          )
        } else {
          p <- p + ggplot2::scale_color_gradientn(
            name     = legend_title_use,
            colours  = colors_gradient,
            values   = scales::rescale(seq(lo_val, hi_val, length.out = length(colors_gradient))),
            limits   = c(lo_val, hi_val),
            na.value = bg_color,
            guide    = ggplot2::guide_colorbar(
              frame.colour = "black", ticks.colour = "black",
              title.hjust = 0, order = 1
            )
          )
        }

        # == Lineage overlay ==
        if (!is.null(lineage_layers)) {
          p <- p + lineage_layers
        }

        # == Facet ==
        if (split_use != ".split_all") {
          p <- p + ggplot2::facet_grid(
            stats::formula(paste0(split_use, " ~ features"))
          )
        } else {
          p <- p + ggplot2::facet_grid(. ~ features)
        }

        # Stash features column for facet label
        p$data[["features"]] <- f

        # == Label ==
        p <- .label_layer(p, dat, f)

        p
      }
    )
  }

  # =========================================================================
  # 7. Combine
  # =========================================================================
  if (isTRUE(combine)) {
    if (length(plist) == 1) return(plist[[1]])
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' required for combining plots.", call. = FALSE)
    }
    return(patchwork::wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol,
                                  byrow = byrow))
  } else {
    return(plist)
  }
}
