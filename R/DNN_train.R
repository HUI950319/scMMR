# ══════════════════════════════════════════════════════════════════════════════
# DNN_train.R — Train a multi-task DNN model
# ══════════════════════════════════════════════════════════════════════════════

#' Train a Multi-Task DNN Model
#'
#' One-step function that loads reference data, selects highly variable genes,
#' creates a multi-task neural network model, trains it, configures the OOD
#' confidence threshold, and saves the trained model. The model jointly learns
#' cell type classification
#' and embedding (e.g. UMAP) regression.
#'
#' @param input       Path to an h5ad file or a Seurat object containing the
#'                    reference single-cell data. Must include cell type labels
#'                    in \code{obs} and embedding coordinates in \code{obsm}.
#' @param label_column  Column name in \code{obs} for cell type labels
#'                      (e.g. \code{"cell_type"}).
#' @param embedding_key Key in \code{obsm} for embedding coordinates
#'                      (e.g. \code{"umap"}, \code{"X_umap"}).
#' @param save_path   Path to save the trained model (\code{.pt} file).
#'                    Parent directories are created automatically.
#' @param n_top_genes Number of highly variable genes to select (default 6000).
#' @param batch_key   Column in \code{obs} for batch-aware HVG selection.
#'                    Set to \code{NULL} for no batch correction.
#' @param ood_threshold OOD (out-of-distribution) confidence threshold
#'                      (default 0.5). Cells with confidence below this are
#'                      flagged as \code{is_ood = TRUE} in predictions.
#'                      Cell type labels always retain the model's best-guess;
#'                      users can filter by \code{is_ood} or \code{confidence}
#'                      downstream.
#' @param hidden_size   Width of the shared backbone (default 512).
#' @param num_epochs    Maximum training epochs (default 50).
#' @param learning_rate Learning rate (default 0.001).
#' @param batch_size    Training batch size (default 256).
#' @param device        Compute device: \code{"auto"} (default, uses GPU if
#'                      available), \code{"cpu"} (force CPU), or \code{"cuda"}
#'                      (require GPU).
#' @param head_hidden_size   Width of task-specific heads (default 256).
#' @param n_resnet_blocks    Number of ResNet blocks (default 4).
#' @param input_dropout_rate Dropout rate for input layer (default 0.25).
#' @param resnet_dropout_rate Dropout rate for ResNet blocks (default 0.1).
#' @param head_dropout_rate   Dropout rate for task heads (default 0.1).
#' @param test_size           Fraction of data for validation (default 0.2).
#' @param early_stopping_patience Epochs to wait before early stopping
#'                                (default 10).
#' @param early_stopping_monitor  Metric to monitor: \code{"cls"} or
#'                                \code{"reg"} (default \code{"cls"}).
#' @param use_gradnorm      Use GradNorm for multi-task balancing (default TRUE).
#' @param gradnorm_alpha    GradNorm alpha parameter (default 0.15).
#' @param label_smoothing   Label smoothing factor (default 0.1).
#' @param use_lr_scheduler  Use cosine annealing LR scheduler (default TRUE).
#' @param use_uncertainty   Use uncertainty-aware regression (default FALSE).
#' @param random_state      Random seed (default 42).
#'
#' @return A named list with:
#'   \describe{
#'     \item{model}{Trained Python MultiTaskModel object.}
#'     \item{var_genes}{Character vector of selected variable genes.}
#'     \item{history}{data.frame with per-epoch training metrics.}
#'     \item{best_val_acc}{Best validation accuracy (\%).}
#'     \item{best_val_rmse}{Best validation RMSE.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#' use_scMMR_python(condaenv = "scanpy-env")
#'
#' result <- DNN_train(
#'   input         = "reference.h5ad",
#'   label_column  = "cell_type",
#'   embedding_key = "umap",
#'   save_path     = "models/my_model.pt",
#'   device        = "auto"
#' )
#' print(result$best_val_acc)
#' print(result$history)
#' }
#'
#' @export
DNN_train <- function(input,
                      label_column,
                      embedding_key,
                      save_path              = "model.pt",
                      n_top_genes            = 6000L,
                      batch_key              = NULL,
                      ood_threshold          = 0.5,
                      hidden_size            = 512L,
                      num_epochs             = 50L,
                      learning_rate          = 0.001,
                      batch_size             = 256L,
                      device                 = "auto",
                      head_hidden_size       = 256L,
                      n_resnet_blocks        = 4L,
                      input_dropout_rate     = 0.25,
                      resnet_dropout_rate    = 0.1,
                      head_dropout_rate      = 0.1,
                      test_size              = 0.2,
                      early_stopping_patience = 10L,
                      early_stopping_monitor = "cls",
                      use_gradnorm           = TRUE,
                      gradnorm_alpha         = 0.15,
                      label_smoothing        = 0.1,
                      use_lr_scheduler       = TRUE,
                      use_uncertainty        = FALSE,
                      random_state           = 42L) {

  .ensure_python()

  # ── 0. Set device ───────────────────────────────────────────────────────────
  .scMMR_env$mt_set_device(device)

  # ── 1. Load reference data ─────────────────────────────────────────────────
  message("\n== Loading reference data ==")
  adata <- .resolve_input(input, embedding_key = embedding_key)

  # ── 2. Select highly variable genes ────────────────────────────────────────
  message("\n== Selecting highly variable genes ==")
  var_genes <- .select_hvg(adata, n_top = n_top_genes, batch_key = batch_key)
  message(sprintf("  %d HVGs selected", length(var_genes)))

  # ── 3. Auto-detect dimensions ──────────────────────────────────────────────
  input_size  <- length(var_genes)
  obs_meta    <- as.data.frame(reticulate::py_to_r(adata$obs),
                               stringsAsFactors = FALSE)
  num_classes <- length(unique(obs_meta[[label_column]]))
  emb_mat     <- reticulate::py_to_r(adata$obsm[[embedding_key]])
  embedding_dim <- ncol(emb_mat)

  message(sprintf("\n  Cell types : %d", num_classes))
  message(sprintf("  Embedding  : %dD (%s)", embedding_dim, embedding_key))
  message("  Class distribution:")
  ct_table <- sort(table(obs_meta[[label_column]]), decreasing = TRUE)
  for (i in seq_along(ct_table)) {
    message(sprintf("    %-30s %d", names(ct_table)[i], ct_table[i]))
  }

  # ── 4. Create model ────────────────────────────────────────────────────────
  message("\n== Creating MultiTaskModel ==")
  model <- .scMMR_env$mt_create_model(
    input_size              = as.integer(input_size),
    num_classes             = as.integer(num_classes),
    embedding_dim           = as.integer(embedding_dim),
    hidden_size             = as.integer(hidden_size),
    head_hidden_size        = as.integer(head_hidden_size),
    n_resnet_blocks         = as.integer(n_resnet_blocks),
    input_dropout_rate      = input_dropout_rate,
    resnet_dropout_rate     = resnet_dropout_rate,
    head_dropout_rate       = head_dropout_rate,
    batch_size              = as.integer(batch_size),
    test_size               = test_size,
    num_epochs              = as.integer(num_epochs),
    early_stopping_patience = as.integer(early_stopping_patience),
    early_stopping_monitor  = early_stopping_monitor,
    learning_rate           = learning_rate,
    use_gradnorm            = use_gradnorm,
    gradnorm_alpha          = gradnorm_alpha,
    label_smoothing         = label_smoothing,
    use_lr_scheduler        = use_lr_scheduler,
    use_uncertainty         = use_uncertainty,
    random_state            = as.integer(random_state)
  )
  message(sprintf("  Model: %d inputs -> %d classes + %dD embedding",
                  input_size, num_classes, embedding_dim))

  # ── 5. Train ───────────────────────────────────────────────────────────────
  message("\n== Training ==")
  result <- .scMMR_env$mt_train(model, adata,
                                 label_column  = label_column,
                                 embedding_key = embedding_key,
                                 var_genes     = var_genes)
  result_r <- reticulate::py_to_r(result)
  message(sprintf("\n  Training complete (%d epochs)", result_r$n_epochs))
  if (!is.null(result_r$best_val_acc)) {
    message(sprintf("    Best val accuracy : %.2f%%", result_r$best_val_acc))
  }
  if (!is.null(result_r$best_val_rmse)) {
    message(sprintf("    Best val RMSE     : %.4f", result_r$best_val_rmse))
  }

  # ── 6. Set OOD confidence threshold (advisory is_ood flag) ─────────────────
  .scMMR_env$mt_set_ood(model, ood_threshold)
  message(sprintf("  OOD confidence threshold: %.2f (is_ood flag only)", ood_threshold))

  # ── 7. Save model ──────────────────────────────────────────────────────────
  message("\n== Saving model ==")
  .scMMR_env$mt_save_model(model, save_path)
  message("  Saved: ", save_path)

  # ── 8. Training history ────────────────────────────────────────────────────
  hist_dict <- .scMMR_env$mt_get_history(model)
  hist_r    <- reticulate::py_to_r(hist_dict)
  history   <- if (length(hist_r) > 0) {
    as.data.frame(hist_r, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  invisible(list(
    model         = model,
    var_genes     = var_genes,
    history       = history,
    best_val_acc  = result_r$best_val_acc,
    best_val_rmse = result_r$best_val_rmse
  ))
}
