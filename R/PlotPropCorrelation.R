# ============================================================================
# PlotPropCorrelation — Cell-Type Proportion vs Gene Expression Correlation
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
