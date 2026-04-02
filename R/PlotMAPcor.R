# ============================================================================
# PlotMAPcor.R -- Cell type correlation heatmap (simplified from scop)
# ============================================================================

#' @title Cell Type Correlation Heatmap
#'
#' @description
#' Visualise the similarity between cell types in query and reference Seurat
#' objects using a ComplexHeatmap-based heatmap.  Cell types are defined by
#' collapsing (averaging) cells within each group, then computing pairwise
#' similarity (cosine, Pearson, or Spearman).
#'
#' When \code{ref = NULL}, the function computes self-correlation within the
#' query object.
#'
#' @md
#' @param query A Seurat object (query dataset).
#' @param ref A Seurat object (reference dataset). \code{NULL} = self-vs-self.
#' @param query_group Column name in query metadata for cell type grouping.
#' @param ref_group Column name in ref metadata for cell type grouping.
#' @param query_assay Assay name for query. \code{NULL} = DefaultAssay.
#' @param ref_assay Assay name for ref. \code{NULL} = DefaultAssay.
#' @param features Character vector of features. \code{NULL} = auto HVF.
#' @param nfeatures Max number of HVF features. Default 2000.
#' @param method Similarity method: \code{"cosine"}, \code{"pearson"}, or
#'   \code{"spearman"}. Default \code{"cosine"}.
#' @param query_annotation Character vector of categorical metadata columns
#'   in query to show as annotation tracks. Default \code{NULL}.
#' @param ref_annotation Character vector of categorical metadata columns
#'   in ref to show as annotation tracks. Default \code{NULL}.
#' @param heatmap_palette Palette for heatmap body. Default \code{"RdBu"}.
#' @param heatmap_palcolor Custom colours for heatmap body. Default \code{NULL}.
#' @param query_palette Palette for query group block. Default \code{"Paired"}.
#' @param query_palcolor Custom colours for query group. Default \code{NULL}.
#' @param ref_palette Palette for ref group block. Default \code{"Set1"}.
#' @param ref_palcolor Custom colours for ref group. Default \code{NULL}.
#' @param query_annotation_palette Palette(s) for query annotation tracks.
#'   Default \code{"Set2"}.
#' @param query_annotation_palcolor Custom colour(s) for query annotations.
#' @param ref_annotation_palette Palette(s) for ref annotation tracks.
#'   Default \code{"Set3"}.
#' @param ref_annotation_palcolor Custom colour(s) for ref annotations.
#' @param border Logical. Draw cell borders. Default \code{TRUE}.
#' @param limits Numeric vector of length 2 for colour scale limits.
#'   \code{NULL} = auto.
#' @param cluster_rows Logical. Cluster rows. Default \code{FALSE}.
#' @param cluster_columns Logical. Cluster columns. Default \code{FALSE}.
#' @param show_row_names Logical. Default \code{TRUE}.
#' @param show_column_names Logical. Default \code{TRUE}.
#' @param row_names_side Side for row names. Default \code{"left"}.
#' @param column_names_side Side for column names. Default \code{"top"}.
#' @param nlabel Integer. Top-N similarity values to label per row/column.
#'   Default \code{0}.
#' @param label_cutoff Numeric. Min similarity to show label. Default \code{0}.
#' @param label_by Label selection dimension: \code{"row"}, \code{"column"},
#'   or \code{"both"}. Default \code{"row"}.
#' @param label_size Numeric. Label font size. Default \code{10}.
#' @param title Character. Heatmap title. Default \code{NULL} = auto.
#' @param width Numeric. Heatmap body width in inches. \code{NULL} = auto.
#' @param height Numeric. Heatmap body height in inches. \code{NULL} = auto.
#' @param seed Integer. Random seed. Default \code{11}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{plot}{A patchwork/ggplot object.}
#'   \item{simil_matrix}{Similarity matrix (query groups x ref groups).}
#'   \item{features}{Features used for computation.}
#' }
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # ---- simulate test data ----
#' set.seed(42)
#' ng <- 200; nc <- 100
#' counts <- matrix(0, nrow = ng, ncol = nc)
#' rownames(counts) <- paste0("Gene", 1:ng)
#' colnames(counts) <- paste0("Cell", 1:nc)
#' ct <- rep(c("A", "B", "C", "D", "E"), each = 20)
#' for (i in 1:nc) {
#'   base <- rpois(ng, lambda = 3)
#'   gi <- match(ct[i], c("A", "B", "C", "D", "E"))
#'   ms <- (gi - 1) * 20 + 1
#'   base[ms:(ms + 19)] <- base[ms:(ms + 19)] + rpois(20, 10)
#'   counts[, i] <- base
#' }
#' srt <- CreateSeuratObject(counts = counts)
#' srt <- NormalizeData(srt, verbose = FALSE)
#' srt <- AddMetaData(srt,
#'   factor(ct, levels = c("A", "B", "C", "D", "E")), "celltype")
#' srt <- AddMetaData(srt,
#'   factor(sample(c("G1", "S", "G2M"), nc, replace = TRUE)), "Phase")
#' srt <- AddMetaData(srt,
#'   factor(sample(c("B1", "B2"), nc, replace = TRUE)), "batch")
#'
#' # ---- 1. Self-correlation (cosine, default) ----
#' res1 <- PlotMAPcor(query = srt, query_group = "celltype")
#' res1$plot
#'
#' # ---- 2. Switch similarity method ----
#' res2 <- PlotMAPcor(query = srt, query_group = "celltype",
#'                     method = "pearson")
#'
#' # ---- 3. Add categorical annotation tracks ----
#' res3 <- PlotMAPcor(query = srt, query_group = "celltype",
#'                     query_annotation = c("Phase", "batch"))
#' res3$plot
#'
#' # ---- 4. Show top-N similarity labels ----
#' res4 <- PlotMAPcor(query = srt, query_group = "celltype",
#'                     nlabel = 2, label_cutoff = 0.9)
#'
#' # ---- 5. Hierarchical clustering ----
#' res5 <- PlotMAPcor(query = srt, query_group = "celltype",
#'                     cluster_rows = TRUE, cluster_columns = TRUE)
#'
#' # ---- 6. Custom colours and fixed limits ----
#' res6 <- PlotMAPcor(query = srt, query_group = "celltype",
#'                     heatmap_palette = "Spectral",
#'                     query_palette = "npg",
#'                     limits = c(0, 1))
#'
#' # ---- 7. Custom feature set ----
#' res7 <- PlotMAPcor(query = srt, query_group = "celltype",
#'                     features = paste0("Gene", 1:50))
#'
#' # ---- 8. Query vs Reference ----
#' counts2 <- matrix(0, nrow = ng, ncol = 60)
#' rownames(counts2) <- paste0("Gene", 1:ng)
#' colnames(counts2) <- paste0("R", 1:60)
#' ct2 <- rep(c("CT1", "CT2", "CT3"), each = 20)
#' for (i in 1:60) {
#'   base <- rpois(ng, lambda = 3)
#'   gi <- match(ct2[i], c("CT1", "CT2", "CT3"))
#'   ms <- (gi - 1) * 30 + 1
#'   base[ms:(ms + 29)] <- base[ms:(ms + 29)] + rpois(30, 8)
#'   counts2[, i] <- base
#' }
#' srt2 <- CreateSeuratObject(counts = counts2)
#' srt2 <- NormalizeData(srt2, verbose = FALSE)
#' srt2 <- AddMetaData(srt2, factor(ct2), "celltype")
#' srt2 <- AddMetaData(srt2,
#'   factor(sample(c("10X", "SS2"), 60, replace = TRUE)), "tech")
#'
#' res8 <- PlotMAPcor(
#'   query = srt, ref = srt2,
#'   query_group = "celltype", ref_group = "celltype",
#'   query_annotation = "Phase",
#'   ref_annotation = "tech",
#'   method = "spearman",
#'   nlabel = 1
#' )
#' res8$plot
#'
#' # ---- 9. Save output ----
#' ggplot2::ggsave("heatmap.pdf", res8$plot, width = 9, height = 7)
#'
#' # ---- 10. Access similarity matrix ----
#' res1$simil_matrix   # query_groups x ref_groups
#' res1$features       # features used
#' }
#'
#' @export
PlotMAPcor <- function(
    query,
    ref                       = NULL,
    query_group               = NULL,
    ref_group                 = NULL,
    query_assay               = NULL,
    ref_assay                 = NULL,
    # --- feature selection ---
    features                  = NULL,
    nfeatures                 = 2000,
    # --- similarity ---
    method                    = c("cosine", "pearson", "spearman"),
    # --- annotation ---
    query_annotation          = NULL,
    ref_annotation            = NULL,
    # --- palette ---
    heatmap_palette           = "RdBu",
    heatmap_palcolor          = NULL,
    query_palette             = "Paired",
    query_palcolor            = NULL,
    ref_palette               = "Set1",
    ref_palcolor              = NULL,
    query_annotation_palette  = "Set2",
    query_annotation_palcolor = NULL,
    ref_annotation_palette    = "Set3",
    ref_annotation_palcolor   = NULL,
    # --- heatmap options ---
    border                    = TRUE,
    limits                    = NULL,
    cluster_rows              = FALSE,
    cluster_columns           = FALSE,
    show_row_names            = TRUE,
    show_column_names         = TRUE,
    row_names_side            = "left",
    column_names_side         = "top",
    nlabel                    = 0,
    label_cutoff              = 0,
    label_by                  = "row",
    label_size                = 10,
    # --- output ---
    title                     = NULL,
    width                     = NULL,
    height                    = NULL,
    seed                      = 11) {

  set.seed(seed)

  # --- check dependencies ---
  for (pkg in c("ComplexHeatmap", "circlize", "patchwork", "Seurat", "SeuratObject")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required. Install it first."))
    }
  }

  method <- match.arg(method)

  # ========================================================================
  # 1. Input handling
  # ========================================================================
  self_cor <- FALSE
  if (is.null(ref)) {
    ref <- query
    ref_group <- query_group
    ref_assay <- query_assay
    ref_palette <- query_palette
    ref_palcolor <- query_palcolor
    ref_annotation <- query_annotation
    ref_annotation_palette <- query_annotation_palette
    ref_annotation_palcolor <- query_annotation_palcolor
    self_cor <- TRUE
  }

  query_assay <- query_assay %||% SeuratObject::DefaultAssay(query)
  ref_assay   <- ref_assay   %||% SeuratObject::DefaultAssay(ref)

  if (is.null(query_group)) stop("'query_group' is required.")
  if (is.null(ref_group))   stop("'ref_group' is required.")

  # ensure factor with original level order
  if (!is.factor(query[[query_group, drop = TRUE]])) {
    query@meta.data[[query_group]] <- factor(
      query[[query_group, drop = TRUE]],
      levels = unique(query[[query_group, drop = TRUE]])
    )
  }
  if (!is.factor(ref[[ref_group, drop = TRUE]])) {
    ref@meta.data[[ref_group]] <- factor(
      ref[[ref_group, drop = TRUE]],
      levels = unique(ref[[ref_group, drop = TRUE]])
    )
  }

  # ========================================================================
  # 2. Feature selection (HVF intersection)
  # ========================================================================
  if (is.null(features)) {
    query <- Seurat::FindVariableFeatures(
      query, assay = query_assay, nfeatures = nfeatures, verbose = FALSE
    )
    ref <- Seurat::FindVariableFeatures(
      ref, assay = ref_assay, nfeatures = nfeatures, verbose = FALSE
    )
    hvf_q <- SeuratObject::VariableFeatures(query, assay = query_assay)
    hvf_r <- SeuratObject::VariableFeatures(ref, assay = ref_assay)
    features <- union(hvf_q, hvf_r)
  }
  features <- intersect(
    features,
    intersect(rownames(query[[query_assay]]), rownames(ref[[ref_assay]]))
  )
  if (length(features) < 5) {
    stop("Too few shared features (", length(features), "). Check assay names.")
  }

  # ========================================================================
  # 3. AverageExpression + similarity
  # ========================================================================
  query_avg <- Seurat::AverageExpression(
    query, features = features, layer = "data",
    assays = query_assay, group.by = query_group, verbose = FALSE
  )[[1]]
  query_avg <- t(as.matrix(log1p(query_avg)))

  ref_avg <- Seurat::AverageExpression(
    ref, features = features, layer = "data",
    assays = ref_assay, group.by = ref_group, verbose = FALSE
  )[[1]]
  ref_avg <- t(as.matrix(log1p(ref_avg)))

  # compute similarity
  if (method == "cosine") {
    norm_q <- sqrt(rowSums(query_avg^2))
    norm_r <- sqrt(rowSums(ref_avg^2))
    norm_q[norm_q == 0] <- 1
    norm_r[norm_r == 0] <- 1
    simil_matrix <- (query_avg %*% t(ref_avg)) /
      (outer(norm_q, norm_r, "*"))
    simil_name <- "Cosine similarity"
  } else {
    # pearson or spearman
    simil_matrix <- stats::cor(
      t(query_avg), t(ref_avg), method = method
    )
    simil_name <- paste0(
      toupper(substring(method, 1, 1)),
      substring(method, 2),
      " correlation"
    )
  }

  simil_matrix[is.na(simil_matrix)] <- 0
  exp_name <- title %||% simil_name

  query_levels <- levels(query[[query_group, drop = TRUE]])
  ref_levels   <- levels(ref[[ref_group, drop = TRUE]])
  simil_matrix <- simil_matrix[
    intersect(query_levels, rownames(simil_matrix)),
    intersect(ref_levels, colnames(simil_matrix)),
    drop = FALSE
  ]

  # ========================================================================
  # 4. Colour mapping
  # ========================================================================
  if (is.null(limits)) {
    color_fun <- circlize::colorRamp2(
      seq(min(simil_matrix, na.rm = TRUE),
          max(simil_matrix, na.rm = TRUE), length.out = 100),
      palette_colors(palette = heatmap_palette, palcolor = heatmap_palcolor)
    )
  } else {
    color_fun <- circlize::colorRamp2(
      seq(limits[1], limits[2], length.out = 100),
      palette_colors(palette = heatmap_palette, palcolor = heatmap_palcolor)
    )
  }

  # ========================================================================
  # 5. Legend list
  # ========================================================================
  lgd <- list()
  lgd[["ht"]] <- ComplexHeatmap::Legend(
    title = exp_name, col_fun = color_fun, border = TRUE
  )

  # ========================================================================
  # 6. Query annotation (left side): group block + categorical annotations
  # ========================================================================
  ha_left <- NULL

  # --- query group block ---
  query_group_cols <- palette_colors(
    query_levels, palette = query_palette, palcolor = query_palcolor
  )
  block_fun_query <- function(index, levels) {
    lv <- levels[1]
    grid::grid.rect(
      gp = grid::gpar(
        fill = query_group_cols[lv],
        col  = if (border) "black" else NA
      )
    )
  }
  query_group_label <- paste0("Query:", query_group)
  anno_query_list <- list()
  anno_query_list[[query_group_label]] <- ComplexHeatmap::anno_block(
    align_to = split(
      seq_along(query_levels), query_levels
    ),
    panel_fun = block_fun_query,
    which = "row", show_name = FALSE
  )
  qblock_args <- anno_query_list
  qblock_args[["which"]] <- "row"
  qblock_args[["show_annotation_name"]] <- TRUE
  qblock_args[["annotation_name_side"]] <- "bottom"
  qblock_args[["border"]] <- TRUE
  ha_left <- do.call(ComplexHeatmap::HeatmapAnnotation, qblock_args)
  lgd[[query_group_label]] <- ComplexHeatmap::Legend(
    title = query_group_label,
    labels = query_levels,
    legend_gp = grid::gpar(fill = query_group_cols[query_levels]),
    border = TRUE
  )

  # --- query categorical annotation tracks ---
  if (!is.null(query_annotation)) {
    if (length(query_annotation_palette) == 1) {
      query_annotation_palette <- rep(
        query_annotation_palette, length(query_annotation)
      )
    }
    ha_anno_list <- list()
    for (i in seq_along(query_annotation)) {
      col_name <- query_annotation[i]
      pal <- query_annotation_palette[i]
      palcol <- if (!is.null(query_annotation_palcolor) &&
                    length(query_annotation_palcolor) >= i) {
        query_annotation_palcolor[[i]]
      } else {
        NULL
      }
      # compute per-group majority vote for categorical annotation
      anno_vals <- .collapse_annotation(
        query, query_group, col_name, query_levels
      )
      anno_cols <- palette_colors(
        levels(anno_vals), palette = pal, palcolor = palcol
      )
      ha_anno_list[[col_name]] <- ComplexHeatmap::anno_simple(
        x = as.character(anno_vals[rownames(simil_matrix)]),
        col = anno_cols,
        which = "row", na_col = "transparent", border = TRUE
      )
      lgd[[paste0("Query:", col_name)]] <- ComplexHeatmap::Legend(
        title = paste0("Query:", col_name),
        labels = levels(anno_vals),
        legend_gp = grid::gpar(fill = anno_cols[levels(anno_vals)]),
        border = TRUE
      )
    }
    anno_args <- ha_anno_list
    anno_args[["which"]] <- "row"
    anno_args[["show_annotation_name"]] <- TRUE
    anno_args[["annotation_name_side"]] <- "bottom"
    anno_args[["width"]] <- grid::unit(10, "mm")
    anno_args[["border"]] <- TRUE
    ha_query_anno <- do.call(ComplexHeatmap::HeatmapAnnotation, anno_args)
    ha_left <- c(ha_left, ha_query_anno)
  }

  # ========================================================================
  # 7. Reference annotation (top): group block + categorical annotations
  # ========================================================================
  ha_top <- NULL

  ref_group_cols <- palette_colors(
    ref_levels, palette = ref_palette, palcolor = ref_palcolor
  )
  block_fun_ref <- function(index, levels) {
    lv <- levels[1]
    grid::grid.rect(
      gp = grid::gpar(
        fill = ref_group_cols[lv],
        col  = if (border) "black" else NA
      )
    )
  }
  ref_group_label <- paste0("Ref:", ref_group)
  anno_ref_list <- list()
  anno_ref_list[[ref_group_label]] <- ComplexHeatmap::anno_block(
    align_to = split(
      seq_along(ref_levels), ref_levels
    ),
    panel_fun = block_fun_ref,
    which = "column", show_name = FALSE
  )
  rblock_args <- anno_ref_list
  rblock_args[["which"]] <- "column"
  rblock_args[["show_annotation_name"]] <- TRUE
  rblock_args[["annotation_name_side"]] <- "left"
  rblock_args[["border"]] <- TRUE
  ha_top <- do.call(ComplexHeatmap::HeatmapAnnotation, rblock_args)
  if (!self_cor) {
    lgd[[ref_group_label]] <- ComplexHeatmap::Legend(
      title = ref_group_label,
      labels = ref_levels,
      legend_gp = grid::gpar(fill = ref_group_cols[ref_levels]),
      border = TRUE
    )
  }

  # --- ref categorical annotation tracks ---
  if (!is.null(ref_annotation)) {
    if (length(ref_annotation_palette) == 1) {
      ref_annotation_palette <- rep(
        ref_annotation_palette, length(ref_annotation)
      )
    }
    ha_ref_anno_list <- list()
    for (i in seq_along(ref_annotation)) {
      col_name <- ref_annotation[i]
      pal <- ref_annotation_palette[i]
      palcol <- if (!is.null(ref_annotation_palcolor) &&
                    length(ref_annotation_palcolor) >= i) {
        ref_annotation_palcolor[[i]]
      } else {
        NULL
      }
      anno_vals <- .collapse_annotation(
        ref, ref_group, col_name, ref_levels
      )
      anno_cols <- palette_colors(
        levels(anno_vals), palette = pal, palcolor = palcol
      )
      ha_ref_anno_list[[col_name]] <- ComplexHeatmap::anno_simple(
        x = as.character(anno_vals[colnames(simil_matrix)]),
        col = anno_cols,
        which = "column", na_col = "transparent", border = TRUE
      )
      if (!self_cor) {
        lgd[[paste0("Ref:", col_name)]] <- ComplexHeatmap::Legend(
          title = paste0("Ref:", col_name),
          labels = levels(anno_vals),
          legend_gp = grid::gpar(fill = anno_cols[levels(anno_vals)]),
          border = TRUE
        )
      }
    }
    ranno_args <- ha_ref_anno_list
    ranno_args[["which"]] <- "column"
    ranno_args[["show_annotation_name"]] <- TRUE
    ranno_args[["annotation_name_side"]] <- "left"
    ranno_args[["height"]] <- grid::unit(10, "mm")
    ranno_args[["border"]] <- TRUE
    ha_ref_anno <- do.call(ComplexHeatmap::HeatmapAnnotation, ranno_args)
    ha_top <- c(ha_top, ha_ref_anno)
  }

  # ========================================================================
  # 8. layer_fun: label top-N similarity values
  # ========================================================================
  layer_fun <- NULL
  if (nlabel > 0) {
    layer_fun <- function(j, i, x, y, w, h, fill) {
      mat <- simil_matrix
      value <- ComplexHeatmap::pindex(mat, i, j)
      ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)

      inds <- NULL
      if (label_by %in% c("row", "both")) {
        for (row in seq_len(nrow(ind_mat))) {
          ind <- ind_mat[row, ]
          cutoff <- max(
            c(sort(value[ind], decreasing = TRUE)[nlabel]),
            na.rm = TRUE
          )
          ind <- ind[value[ind] >= cutoff & value[ind] >= label_cutoff]
          inds <- c(inds, ind)
        }
      }
      if (label_by %in% c("column", "both")) {
        for (col in seq_len(ncol(ind_mat))) {
          ind <- ind_mat[, col]
          cutoff <- max(
            c(sort(value[ind], decreasing = TRUE)[nlabel]),
            na.rm = TRUE
          )
          ind <- ind[value[ind] >= cutoff & value[ind] >= label_cutoff]
          inds <- c(inds, ind)
        }
      }
      if (label_by == "both") inds <- inds[duplicated(inds)]
      if (length(inds) > 0) {
        # white outline + black text
        theta <- seq(pi / 8, 2 * pi, length.out = 16)
        for (th in theta) {
          grid::grid.text(
            round(value[inds], 2),
            x = x[inds] + grid::unit(cos(th) * label_size / 30, "mm"),
            y = y[inds] + grid::unit(sin(th) * label_size / 30, "mm"),
            gp = grid::gpar(fontsize = label_size, col = "white")
          )
        }
        grid::grid.text(
          round(value[inds], 2),
          x[inds], y[inds],
          gp = grid::gpar(fontsize = label_size, col = "black")
        )
      }
    }
  }

  # ========================================================================
  # 9. Heatmap assembly
  # ========================================================================
  ht_args <- list(
    name             = exp_name,
    matrix           = simil_matrix,
    col              = color_fun,
    row_title        = paste0("Query:", query_group),
    row_title_side   = "left",
    row_title_rot    = 90,
    column_title     = paste0("Ref:", ref_group),
    column_title_side = "top",
    column_title_rot = 0,
    cluster_rows     = cluster_rows,
    cluster_columns  = cluster_columns,
    show_row_names   = show_row_names,
    show_column_names = show_column_names,
    row_names_side   = row_names_side,
    column_names_side = column_names_side,
    row_names_rot    = 0,
    column_names_rot = 90,
    top_annotation   = ha_top,
    left_annotation  = ha_left,
    show_heatmap_legend = FALSE,
    border           = border
  )
  if (!is.null(layer_fun)) ht_args[["layer_fun"]] <- layer_fun
  if (!is.null(width))  ht_args[["width"]]  <- grid::unit(width, "inch")
  if (!is.null(height)) ht_args[["height"]] <- grid::unit(height, "inch")

  ht <- do.call(ComplexHeatmap::Heatmap, ht_args)

  # ========================================================================
  # 10. Render and return
  # ========================================================================
  # estimate render size
  nr <- nrow(simil_matrix)
  nc <- ncol(simil_matrix)
  ht_w <- width  %||% max(nc * 0.6 + 2.5, 6)
  ht_h <- height %||% max(nr * 0.6 + 2.5, 6)

  g_tree <- grid::grid.grabExpr(
    ComplexHeatmap::draw(ht, annotation_legend_list = lgd),
    width  = grid::unit(ht_w, "inch"),
    height = grid::unit(ht_h, "inch"),
    wrap = TRUE, wrap.grobs = TRUE
  )
  p <- patchwork::wrap_plots(g_tree)

  list(
    plot          = p,
    simil_matrix  = simil_matrix,
    features      = features
  )
}


# --------------------------------------------------------------------------
# Internal: collapse categorical annotation by majority vote per group
# --------------------------------------------------------------------------
#' @keywords internal
.collapse_annotation <- function(srt, group_col, anno_col, group_levels) {
  meta <- srt@meta.data
  if (!anno_col %in% colnames(meta)) {
    stop("Annotation column '", anno_col, "' not found in metadata.")
  }
  vals <- meta[[anno_col]]
  grp  <- meta[[group_col]]
  if (!is.factor(vals)) vals <- factor(vals, levels = unique(vals))
  # majority vote: most frequent category in each group
  result <- vapply(group_levels, function(g) {
    sub_vals <- vals[grp == g]
    if (length(sub_vals) == 0) return(NA_character_)
    tb <- sort(table(sub_vals), decreasing = TRUE)
    names(tb)[1]
  }, character(1))
  result <- factor(result, levels = levels(vals))
  names(result) <- group_levels
  result
}
