# ══════════════════════════════════════════════════════════════════════════════
# DNN_predict.R — Predict cell types using a trained multi-task model
# ══════════════════════════════════════════════════════════════════════════════

#' Predict Cell Types Using a Trained Multi-Task DNN Model
#'
#' Loads a trained multi-task model, runs predictions on query data,
#' and optionally computes gene importance via Integrated Gradients.
#' When a GMT file is provided, per-cell pathway importance scores are
#' also computed by aggregating gene-level IG attributions to pathways.
#'
#' @param query       Path to an h5ad file or a Seurat object containing the
#'                    query single-cell data.
#' @param model_path  Path to the saved model file (\code{.pt}).
#' @param save_path   Path to save results as HDF5 (\code{.h5}).
#'                    Set to \code{NULL} to skip saving.
#' @param true_label_col  Column in \code{obs} with ground-truth labels
#'                        for accuracy evaluation. Set to \code{NULL} if
#'                        no labels available.
#' @param explain     Logical. If \code{TRUE}, compute gene importance
#'                    via Integrated Gradients (default \code{FALSE}).
#' @param top_k_global  Number of top genes for global importance
#'                      (default 15).
#' @param top_k_class   Number of top genes per cell type (default 10).
#' @param n_cells_explain  Number of cells for attribution baseline
#'                         (speed tradeoff, default 50).
#' @param pathway_gmt  Path to a GMT file for pathway-level scoring.
#'   Requires \code{explain = TRUE}. When provided, computes per-cell
#'   pathway importance via \code{|IG attributions| × GMT mask / pathway_size}.
#'   Built-in GMT files are available in the package
#'   (\code{system.file("extdata/gmt", package = "scMMR")}):
#'   \code{reactome.gmt}, \code{GO_bp.gmt}, \code{TF.gmt},
#'   \code{immune.gmt} (human), plus \code{m_reactome.gmt},
#'   \code{m_GO_bp.gmt}, \code{m_TF.gmt} (mouse).
#'   Set to \code{NULL} to skip pathway scoring.
#' @param pathway_min_genes  Minimum number of overlapping genes between
#'   a pathway and the model's variable genes for the pathway to be
#'   retained (default 5).
#' @param pathway_n_steps  Number of Integrated Gradients interpolation
#'   steps for pathway scoring (default 20). Lower values are faster
#'   but less precise. For gene importance \code{n_steps=50} is used;
#'   pathway scoring tolerates fewer steps since attributions are
#'   aggregated across many genes.
#' @param return_embedding Logical. If \code{TRUE}, include the 512-dim
#'                    shared embedding in the output (default \code{FALSE}).
#' @param device      Compute device: \code{"auto"} (default),
#'                    \code{"cpu"}, or \code{"cuda"}.
#'
#' @return A named list with:
#'   \describe{
#'     \item{predictions}{data.frame with columns: cell_id, cell_type_pred,
#'           confidence, is_ood, umap_1_pred, umap_2_pred.
#'           \code{cell_type_pred} always contains the model's best-guess
#'           label (never \code{"Unknown"}). Use \code{is_ood} or
#'           \code{confidence} for downstream filtering.}
#'     \item{shared_embedding}{Matrix of shape (n_cells, 512) — shared
#'           representation from the model backbone. Only included when
#'           \code{return_embedding = TRUE}.}
#'     \item{imp_global}{data.frame (gene, importance) if \code{explain=TRUE},
#'           otherwise \code{NULL}. Computed via Integrated Gradients.}
#'     \item{imp_per_class}{data.frame (cell_type, n_cells, rank, gene,
#'           importance) if \code{explain=TRUE}, otherwise \code{NULL}.}
#'     \item{pathway_scores}{Matrix of shape (n_cells, n_pathways) with
#'           per-cell pathway importance scores. Rownames are cell IDs;
#'           column names are pathway names. Only included when
#'           \code{pathway_gmt} is provided and \code{explain = TRUE}.
#'           Scores represent each pathway's discriminative contribution
#'           (via Integrated Gradients), not expression enrichment.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#' use_scMMR_python(condaenv = "scanpy-env")
#'
#' # Basic prediction
#' pred <- DNN_predict(
#'   query      = "query.h5ad",
#'   model_path = "models/my_model.pt",
#'   explain    = TRUE,
#'   device     = "cpu"
#' )
#' head(pred$predictions)
#' pred$imp_global
#'
#' # With pathway scoring (GMT file)
#' gmt <- system.file("extdata/gmt", "reactome.gmt", package = "scMMR")
#' pred2 <- DNN_predict(
#'   query       = "query.h5ad",
#'   model_path  = "models/my_model.pt",
#'   explain     = TRUE,
#'   pathway_gmt = gmt,
#'   device      = "cpu"
#' )
#' dim(pred2$pathway_scores)  # n_cells x n_pathways
#' }
#'
#' @export
DNN_predict <- function(query,
                        model_path,
                        save_path          = NULL,
                        true_label_col     = "cell_type",
                        explain            = FALSE,
                        top_k_global       = 15L,
                        top_k_class        = 10L,
                        n_cells_explain    = 50L,
                        pathway_gmt        = NULL,
                        pathway_min_genes  = 5L,
                        pathway_n_steps    = 20L,
                        return_embedding   = FALSE,
                        device             = "auto") {

  .ensure_python()

  # ── 0. Set device ───────────────────────────────────────────────────────────
  .scMMR_env$mt_set_device(device)

  # ── 1. Load model ──────────────────────────────────────────────────────────
  model_path <- normalizePath(model_path, mustWork = TRUE)
  message("\n== Loading model ==")
  message("  Path: ", model_path)
  model   <- .scMMR_env$mt_load_model(model_path)
  n_cls   <- reticulate::py_to_r(model$num_classes)
  n_genes <- reticulate::py_to_r(model$input_size)
  message(sprintf("  %d classes | %d genes", n_cls, n_genes))

  # ── 2. Load query data ─────────────────────────────────────────────────────
  message("\n== Loading query data ==")
  query_adata <- .resolve_input(query, embedding_key = "umap")

  # ── 3. Prepare feature matrix ──────────────────────────────────────────────
  message("\n== Preparing feature matrix ==")
  X <- .scMMR_env$mt_prepare_query(query_adata, model)

  # ── 4. Run predictions ─────────────────────────────────────────────────────
  message("\n== Running predictions ==")
  pred_out  <- .scMMR_env$mt_predict(model, query_adata, X)
  cell_list <- reticulate::py_to_r(pred_out[["cell_level"]])
  cell_df   <- as.data.frame(lapply(cell_list, as.vector),
                             stringsAsFactors = FALSE)
  cell_df$is_ood <- as.logical(cell_df$is_ood)

  # Shared embedding (only compute if requested)
  shared_emb <- NULL
  if (return_embedding) {
    shared_emb <- reticulate::py_to_r(pred_out[["shared_embedding"]])
    rownames(shared_emb) <- cell_df$cell_id
  }

  # ── Summary ────────────────────────────────────────────────────────────────
  message("\n== Summary ==")
  n_ood <- sum(cell_df$is_ood)
  n_tot <- nrow(cell_df)
  conf  <- cell_df$confidence
  message(sprintf("  Cells          : %d", n_tot))
  message(sprintf("  Low-confidence : %d/%d (%.1f%%) [is_ood flag]",
                  n_ood, n_tot, 100 * n_ood / n_tot))
  message(sprintf("  Confidence     : [%.4f, %.4f]  mean=%.4f",
                  min(conf), max(conf), mean(conf)))

  # Accuracy (if ground truth available)
  if (!is.null(true_label_col)) {
    obs_meta <- as.data.frame(reticulate::py_to_r(query_adata$obs),
                              stringsAsFactors = FALSE)
    if (true_label_col %in% names(obs_meta)) {
      true_labels <- as.character(obs_meta[[true_label_col]])
      acc <- mean(cell_df$cell_type_pred == true_labels)
      message(sprintf("  Accuracy       : %.2f%%", acc * 100))
    }
  }

  # ── 5. Gene importance (optional) ──────────────────────────────────────────
  imp_global    <- NULL
  imp_per_class <- NULL

  if (explain) {
    message("\n== Computing gene importance ==")

    # Global (classification)
    message("  Global importance (classification) ...")
    global_result <- .scMMR_env$mt_explain_global(
      model, X, task = "cls",
      top_k   = as.integer(top_k_global),
      n_cells = as.integer(n_cells_explain)
    )
    global_r <- reticulate::py_to_r(global_result)
    imp_global <- data.frame(
      gene       = as.character(global_r[["gene"]]),
      importance = as.numeric(global_r[["importance"]]),
      stringsAsFactors = FALSE
    )
    message("  Top genes: ", paste(head(imp_global$gene, 5), collapse = ", "))

    # Per-class
    message("  Per-class importance ...")
    # Use true labels if available, else predicted
    if (!is.null(true_label_col)) {
      obs_meta <- as.data.frame(reticulate::py_to_r(query_adata$obs),
                                stringsAsFactors = FALSE)
      if (true_label_col %in% names(obs_meta)) {
        use_labels <- as.character(obs_meta[[true_label_col]])
      } else {
        use_labels <- cell_df$cell_type_pred
      }
    } else {
      use_labels <- cell_df$cell_type_pred
    }

    pc_result <- .scMMR_env$mt_explain_per_class(
      model, X,
      labels = use_labels,
      top_k  = as.integer(top_k_class)
    )
    pc_r <- reticulate::py_to_r(pc_result)
    imp_per_class <- data.frame(
      cell_type  = as.character(pc_r[["cell_type"]]),
      n_cells    = as.integer(pc_r[["n_cells"]]),
      rank       = as.integer(pc_r[["rank"]]),
      gene       = as.character(pc_r[["gene"]]),
      importance = as.numeric(pc_r[["importance"]]),
      stringsAsFactors = FALSE
    )
  }

  # ── 6. Pathway scoring (optional: GMT + IG) ────────────────────────────────
  pathway_scores <- NULL

  if (!is.null(pathway_gmt)) {
    if (!explain) {
      warning("pathway_gmt is provided but explain=FALSE. ",
              "Setting explain=TRUE for pathway scoring.")
      explain <- TRUE
    }

    gmt_path <- normalizePath(pathway_gmt, mustWork = TRUE)
    message("\n== Computing pathway scores ==")
    message("  GMT: ", gmt_path)

    pw_result <- .scMMR_env$mt_pathway_score(
      model, X, gmt_path,
      task       = "cls",
      min_genes  = as.integer(pathway_min_genes),
      n_steps    = as.integer(pathway_n_steps),
      batch_size = 256L
    )
    pw_r <- reticulate::py_to_r(pw_result)

    pathway_scores <- pw_r[["scores"]]
    pw_names <- as.character(pw_r[["pathway_names"]])
    pw_sizes <- as.integer(pw_r[["pathway_sizes"]])

    rownames(pathway_scores) <- cell_df$cell_id
    colnames(pathway_scores) <- pw_names

    message(sprintf("  Result: %d cells × %d pathways",
                    nrow(pathway_scores), ncol(pathway_scores)))
    message(sprintf("  Pathway gene counts: [%d, %d]  median=%d",
                    min(pw_sizes), max(pw_sizes),
                    as.integer(stats::median(pw_sizes))))
  }

  # ── 7. Save results (optional) ─────────────────────────────────────────────
  if (!is.null(save_path)) {
    message("\n== Saving results ==")
    cell_dict <- reticulate::r_to_py(as.list(cell_df))
    if (!is.null(imp_per_class)) {
      pclass_dict <- reticulate::r_to_py(as.list(imp_per_class))
    } else {
      # Empty dict if no importance computed
      pclass_dict <- reticulate::r_to_py(list(
        cell_type  = character(0),
        n_cells    = integer(0),
        rank       = integer(0),
        gene       = character(0),
        importance = numeric(0)
      ))
    }
    .scMMR_env$mt_save_h5(cell_dict, pclass_dict, save_path)
    message("  Saved: ", save_path)
  }

  message("\n== Done ==")

  # Build output list — only include non-NULL components
  out <- list(predictions = cell_df)
  if (return_embedding)          out$shared_embedding <- shared_emb
  if (!is.null(imp_global))      out$imp_global       <- imp_global
  if (!is.null(imp_per_class))   out$imp_per_class    <- imp_per_class
  if (!is.null(pathway_scores))  out$pathway_scores   <- pathway_scores
  invisible(out)
}
