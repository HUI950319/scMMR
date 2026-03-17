# ============================================================================
# CytoTRACE2 wrapper — cellular potency prediction
# ============================================================================


#' Run CytoTRACE 2 Cellular Potency Prediction
#'
#' Predict single-cell developmental potential (stemness / differentiation)
#' using CytoTRACE 2.  Wraps \code{CytoTRACE2::cytotrace2()} with a
#' cleaner interface, automatic species handling, and Seurat v5
#' compatibility.
#'
#' @param object A Seurat object, matrix (genes x cells), or data.frame.
#' @param species Character: \code{"human"} or \code{"mouse"}.
#'   Default \code{"human"}.
#' @param assay Character. Seurat assay to use (default: active assay).
#'   Ignored for matrix input.
#' @param layer Character: \code{"counts"} (default) or \code{"data"}.
#'   Ignored for matrix input.
#' @param batch_size Integer or \code{NULL}. Number of cells per batch
#'   for KNN smoothing.  Default 10000 (recommended for >10K cells).
#'   Set \code{NULL} to disable subsampling.
#' @param smooth_batch_size Integer or \code{NULL}. Cells to subsample
#'   within each batch for diffusion smoothing.  Default 1000.
#' @param cores Integer. Number of parallel cores.
#'   Default \code{NULL} (auto-detect, use half).
#' @param seed Integer. Random seed for reproducibility (default 14).
#' @param ... Additional arguments passed to
#'   \code{CytoTRACE2::cytotrace2()}.
#'
#' @return
#' For Seurat input: the Seurat object with metadata columns added:
#' \describe{
#'   \item{CytoTRACE2_Score}{Predicted potency score (0-1, higher = more potent)}
#'   \item{CytoTRACE2_Potency}{Potency category (Differentiated, Unipotent,
#'     Oligopotent, Multipotent, Pluripotent, Totipotent)}
#'   \item{CytoTRACE2_Relative}{Relative order normalised to 0-1}
#'   \item{preKNN_CytoTRACE2_Score}{Score before KNN smoothing}
#'   \item{preKNN_CytoTRACE2_Potency}{Category before KNN smoothing}
#' }
#'
#' For matrix/data.frame input: a data.frame with cell IDs as row names
#' and the same columns as above.
#'
#' @details
#' CytoTRACE 2 predicts cellular potency in three steps:
#' \enumerate{
#'   \item Preprocess input data
#'   \item Predict potency (discrete categories + continuous score 0-1)
#'   \item Smooth predictions via diffusion + KNN rescaling
#' }
#'
#' If CytoTRACE2 is not installed, the function will prompt you with
#' installation instructions.
#'
#' @examples
#' \dontrun{
#' # Seurat object
#' seu <- RunCytoTRACE2(seu, species = "human")
#' FeaturePlot(seu, features = "CytoTRACE2_Score")
#' VlnPlot(seu, features = "CytoTRACE2_Score", group.by = "cell_type")
#'
#' # Raw matrix (genes x cells)
#' result <- RunCytoTRACE2(counts_matrix, species = "mouse")
#' }
#'
#' @export
RunCytoTRACE2 <- function(object, ...) {
  UseMethod("RunCytoTRACE2")
}


#' @rdname RunCytoTRACE2
#' @method RunCytoTRACE2 Seurat
#' @export
RunCytoTRACE2.Seurat <- function(object,
                                 species          = c("human", "mouse"),
                                 assay            = NULL,
                                 layer            = c("counts", "data"),
                                 batch_size       = 10000,
                                 smooth_batch_size = 1000,
                                 cores            = NULL,
                                 seed             = 14,
                                 ...) {

  .check_cytotrace2()

  species <- match.arg(species)
  layer   <- match.arg(layer)
  assay   <- assay %||% SeuratObject::DefaultAssay(object)

  message("[CytoTRACE2] Running on Seurat object (",
          ncol(object), " cells, species = ", species, ") ...")

  ct2_args <- list(
    input                 = object,
    species               = species,
    is_seurat             = TRUE,
    slot_type             = layer,
    batch_size            = batch_size,
    smooth_batch_size     = smooth_batch_size,
    parallelize_models    = TRUE,
    parallelize_smoothing = TRUE,
    ncores                = cores,
    seed                  = seed
  )
  ct2_args <- c(ct2_args, list(...))

  result <- do.call(CytoTRACE2::cytotrace2, ct2_args)

  message("[CytoTRACE2] Done. Potency scores added to metadata.")
  result
}


#' @rdname RunCytoTRACE2
#' @method RunCytoTRACE2 default
#' @export
RunCytoTRACE2.default <- function(object,
                                  species          = c("human", "mouse"),
                                  batch_size       = 10000,
                                  smooth_batch_size = 1000,
                                  cores            = NULL,
                                  seed             = 14,
                                  ...) {

  .check_cytotrace2()

  species <- match.arg(species)

  # Convert Matrix to data.frame (CytoTRACE2 requirement)
  if (inherits(object, "Matrix")) {
    object <- as.matrix(object)
  }
  if (is.matrix(object)) {
    object <- as.data.frame(object)
  }
  if (!is.data.frame(object))
    stop("'object' must be a Seurat object, matrix, or data.frame ",
         "(genes as rows, cells as columns).")

  if (is.null(rownames(object)))
    stop("'object' must have row names (gene names).")
  if (is.null(colnames(object)))
    stop("'object' must have column names (cell IDs).")

  message("[CytoTRACE2] Running on matrix (",
          nrow(object), " genes x ", ncol(object),
          " cells, species = ", species, ") ...")

  ct2_args <- list(
    input                 = object,
    species               = species,
    is_seurat             = FALSE,
    batch_size            = batch_size,
    smooth_batch_size     = smooth_batch_size,
    parallelize_models    = TRUE,
    parallelize_smoothing = TRUE,
    ncores                = cores,
    seed                  = seed
  )
  ct2_args <- c(ct2_args, list(...))

  result <- do.call(CytoTRACE2::cytotrace2, ct2_args)

  message("[CytoTRACE2] Done.")
  result
}


# --------------------------------------------------------------------------
# Internal: check CytoTRACE2 availability
# --------------------------------------------------------------------------

.check_cytotrace2 <- function() {
  if (!requireNamespace("CytoTRACE2", quietly = TRUE))
    stop("Package 'CytoTRACE2' is required but not installed.\n",
         "Install with:\n",
         '  devtools::install_github("digitalcytometry/cytotrace2",\n',
         '    subdir = "cytotrace2_r")\n',
         "or:\n",
         '  pak::pak("digitalcytometry/cytotrace2/cytotrace2_r")',
         call. = FALSE)
}


# ============================================================================
# CytoTRACE2 Visualization
# ============================================================================

#' CytoTRACE 2 Potency Plot
#'
#' Violin + box plot of CytoTRACE2 potency scores grouped by cell type,
#' with optional statistical comparisons.
#'
#' @param object A Seurat object with CytoTRACE2 results (from
#'   \code{RunCytoTRACE2}).
#' @param group.by Character. Metadata column for grouping (e.g.
#'   \code{"cell_type"}).  Default uses active idents.
#' @param score Character. Which score to plot:
#'   \code{"CytoTRACE2_Score"} (default) or \code{"CytoTRACE2_Relative"}.
#' @param order Logical. If \code{TRUE} (default), reorder groups by
#'   median score (most differentiated to most potent).
#' @param colors Character vector of colours, or \code{NULL} for default.
#' @param base_size Numeric. Base font size (default 12).
#' @param point.size Numeric. Jitter point size (default 0.3).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' seu <- RunCytoTRACE2(seu, species = "human")
#' PlotCytoTRACE2(seu, group.by = "cell_type")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_jitter
#'   scale_fill_manual labs coord_flip theme_bw theme element_text
#'   element_blank position_jitter
#' @export
PlotCytoTRACE2 <- function(object,
                            group.by   = NULL,
                            score      = "CytoTRACE2_Score",
                            order      = TRUE,
                            colors     = NULL,
                            base_size  = 12,
                            point.size = 0.3) {

  if (!score %in% colnames(object@meta.data))
    stop("Score column '", score, "' not found in metadata. ",
         "Run RunCytoTRACE2() first.")

  # Build plot data
  if (is.null(group.by)) {
    df <- data.frame(
      Group = as.character(Seurat::Idents(object)),
      Score = object@meta.data[[score]],
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Group = as.character(object@meta.data[[group.by]]),
      Score = object@meta.data[[score]],
      stringsAsFactors = FALSE
    )
  }

  # Order by median score
  if (order) {
    med_order <- stats::aggregate(Score ~ Group, data = df, FUN = median)
    med_order <- med_order[order(med_order$Score), ]
    df$Group  <- factor(df$Group, levels = med_order$Group)
  }

  n_groups <- length(unique(df$Group))

  # Default colour palette
  if (is.null(colors)) {
    if (n_groups <= 8) {
      colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
    } else {
      colors <- colorspace::qualitative_hcl(n_groups, palette = "Dark 3")
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Group,
                                         y = .data$Score,
                                         fill = .data$Group)) +
    ggplot2::geom_violin(alpha = 0.7, scale = "width", linewidth = 0.3) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA,
                          fill = "white", alpha = 0.6) +
    ggplot2::geom_jitter(size = point.size, alpha = 0.3, width = 0.1) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL,
                  y = gsub("_", " ", score),
                  title = "CytoTRACE2 Cellular Potency") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text     = ggplot2::element_text(color = "black"),
      plot.title    = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}
