# ══════════════════════════════════════════════════════════════════════════════
# DNN_deconv_predict.R — Predict cell type proportions from bulk RNA-seq
# ══════════════════════════════════════════════════════════════════════════════

#' Predict Cell Type Proportions from Bulk RNA-seq
#'
#' Uses a trained deconvolution DNN model (from \code{\link{DNN_deconv_train}})
#' to estimate cell type proportions from bulk RNA-seq expression data.
#'
#' The function handles gene alignment automatically: genes present in both
#' the bulk data and the model are used; missing genes are zero-filled.
#' Expression is normalized to log1p(CPM) to match the training procedure.
#'
#' When \code{adaptive = TRUE}, the model runs a TAPE-style adaptive stage
#' that alternately fine-tunes the decoder and encoder on the target bulk
#' data, producing tissue-adapted cell-type-specific gene expression
#' profiles (GEP / sigmatrix).
#'
#' @param bulk_expr  Bulk expression data. Accepts multiple formats:
#'   \itemize{
#'     \item A \strong{matrix} or \strong{data.frame} (genes x samples).
#'       Rownames = gene names, colnames = sample names.
#'     \item A \strong{CSV file path} (genes x samples, first column =
#'       gene names or rownames).
#'     \item An \strong{h5ad file path}.
#'     \item A \strong{Seurat object} (uses RNA assay counts).
#'   }
#' @param model_path  Path to the trained deconvolution model (\code{.pt}
#'   file created by \code{\link{DNN_deconv_train}}).
#' @param device  Compute device: \code{"auto"} (default), \code{"cpu"},
#'   or \code{"cuda"}.
#' @param adaptive  Logical; if \code{TRUE}, run TAPE-style adaptive
#'   stage to refine proportions and extract cell-type-specific GEP
#'   (default \code{FALSE}).
#' @param adaptive_mode  Adaptive mode: \code{"overall"} (default) adapts
#'   on all samples jointly (one shared GEP), or \code{"high-resolution"}
#'   adapts per-sample (individual GEPs, slower).
#' @param adaptive_max_iter  Number of alternating decoder/encoder rounds
#'   (default 3).
#' @param adaptive_steps  Gradient steps per phase per round (default 300).
#' @param adaptive_lr  Learning rate for adaptive optimization
#'   (default 1e-4).
#'
#' @return When \code{adaptive = FALSE}: a \code{data.frame} with
#'   \code{sample_id} and one column per cell type (backward compatible).
#'
#'   When \code{adaptive = TRUE}: a named \code{list} with:
#'   \describe{
#'     \item{proportions}{data.frame of cell type proportions.}
#'     \item{sigmatrix}{Cell-type-specific GEP.
#'       For \code{"overall"} mode: a data.frame (cell_types x genes).
#'       For \code{"high-resolution"} mode: a 3D array
#'       (samples x cell_types x genes).}
#'   }
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#' use_scMMR_python(condaenv = "scMMR")
#'
#' # Standard prediction (no adaptive)
#' props <- DNN_deconv_predict(
#'   bulk_expr  = bulk_matrix,
#'   model_path = "models/deconv_model.pt"
#' )
#'
#' # Adaptive prediction with GEP extraction
#' result <- DNN_deconv_predict(
#'   bulk_expr      = bulk_matrix,
#'   model_path     = "models/deconv_model.pt",
#'   adaptive       = TRUE,
#'   adaptive_mode  = "overall"
#' )
#' result$proportions   # cell type proportions
#' result$sigmatrix     # cell-type-specific GEP (K x genes)
#' }
#'
#' @seealso \code{\link{DNN_deconv_train}} for training the model.
#'
#' @export
DNN_deconv_predict <- function(bulk_expr,
                                model_path,
                                device              = "auto",
                                adaptive            = FALSE,
                                adaptive_mode       = "overall",
                                adaptive_max_iter   = 3L,
                                adaptive_steps      = 300L,
                                adaptive_lr         = 1e-4) {

  .ensure_python()

  # ── 0. Set device ───────────────────────────────────────────────────────────
  .scMMR_env$deconv_set_device(device)

  # ── 1. Load model ──────────────────────────────────────────────────────────
  model_path <- normalizePath(model_path, mustWork = TRUE)
  message("\n== Loading deconvolution model ==")
  model <- .scMMR_env$deconv_load_model(model_path)
  cell_types <- reticulate::py_to_r(model$cell_type_names)
  n_genes    <- reticulate::py_to_r(model$input_size)
  message(sprintf("  %d cell types: %s",
                  length(cell_types), paste(cell_types, collapse = ", ")))
  message(sprintf("  %d genes in model", n_genes))

  # ── 2. Prepare bulk expression ─────────────────────────────────────────────
  message("\n== Preparing bulk expression ==")
  bulk_data <- .prepare_bulk_input(bulk_expr)
  message(sprintf("  Input: %d genes x %d samples",
                  length(bulk_data$gene_names), length(bulk_data$sample_names)))

  X_bulk <- .scMMR_env$deconv_prepare_bulk(
    bulk_matrix  = reticulate::r_to_py(bulk_data$matrix),
    model        = model,
    gene_names   = bulk_data$gene_names,
    sample_names = bulk_data$sample_names
  )

  # ── 3. Predict ─────────────────────────────────────────────────────────────
  if (isTRUE(adaptive)) {
    message(sprintf("\n== Adaptive prediction (mode=%s, max_iter=%d, steps=%d) ==",
                    adaptive_mode, adaptive_max_iter, adaptive_steps))
  } else {
    message("\n== Predicting proportions ==")
  }

  pred <- .scMMR_env$deconv_predict(
    model    = model,
    X_bulk   = X_bulk,
    adaptive = adaptive,
    mode     = adaptive_mode,
    max_iter = as.integer(adaptive_max_iter),
    step     = as.integer(adaptive_steps),
    lr       = adaptive_lr
  )
  pred_r <- reticulate::py_to_r(pred)

  # ── 4. Build output ────────────────────────────────────────────────────────
  proportions <- as.data.frame(pred_r$proportions, stringsAsFactors = FALSE)
  colnames(proportions) <- cell_types
  proportions <- cbind(
    sample_id = bulk_data$sample_names,
    proportions,
    stringsAsFactors = FALSE
  )

  message(sprintf("\n  Predicted: %d samples x %d cell types",
                  nrow(proportions), length(cell_types)))

  if (!isTRUE(adaptive)) {
    # ── Backward compatible: return data.frame directly ──
    message("== Done ==\n")
    return(proportions)
  }

  # ── Adaptive mode: return list with proportions + sigmatrix ──
  result <- list(proportions = proportions)

  if (!is.null(pred_r$sigmatrix)) {
    gep <- pred_r$sigmatrix
    var_genes <- pred_r$var_genes

    if (adaptive_mode == "overall") {
      # K × n_genes matrix
      gep <- as.data.frame(gep, stringsAsFactors = FALSE)
      if (!is.null(var_genes) && length(var_genes) == ncol(gep)) {
        colnames(gep) <- var_genes
      }
      rownames(gep) <- cell_types
      result$sigmatrix <- gep
    } else if (adaptive_mode == "high-resolution") {
      # N × K × n_genes 3D array — keep as-is
      # dimnames: samples × cell_types × genes
      if (is.array(gep) && length(dim(gep)) == 3) {
        dimnames(gep) <- list(
          bulk_data$sample_names,
          cell_types,
          if (!is.null(var_genes)) var_genes else NULL
        )
      }
      result$sigmatrix <- gep
    }
  }

  message("== Done ==\n")
  result
}
