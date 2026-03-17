# ══════════════════════════════════════════════════════════════════════════════
# DNN_deconv_train.R — Train a deconvolution DNN model
# ══════════════════════════════════════════════════════════════════════════════

#' Train a Deconvolution DNN Model
#'
#' Trains a deep neural network to estimate cell type proportions from bulk
#' RNA-seq expression data. The model is trained on simulated pseudo-bulk
#' samples generated from a single-cell RNA-seq reference dataset.
#'
#' The model architecture reuses the same ResNet backbone as
#' \code{\link{DNN_train}} but with key differences:
#' \itemize{
#'   \item \strong{Input}: continuous log-CPM expression (not binary)
#'   \item \strong{Output}: cell type proportions (softmax, sum = 1)
#'   \item \strong{Training data}: simulated pseudo-bulk from scRNA-seq
#' }
#'
#' @param reference  Path to an h5ad file or a Seurat object containing the
#'   scRNA-seq reference. Must include raw counts and cell type labels.
#' @param label_column  Column name in \code{obs}/\code{meta.data} for cell
#'   type labels (e.g. \code{"cell_type"}).
#' @param save_path  Path to save the trained model (\code{.pt} file).
#'   Parent directories are created automatically.
#' @param n_pseudobulk  Number of pseudo-bulk samples to generate for
#'   training (default 5000). More samples generally improve accuracy.
#' @param n_cells_per_sample  Number of cells to mix per pseudo-bulk sample
#'   (default 500).
#' @param n_top_genes  Number of highly variable genes to select (default 6000).
#' @param batch_key  Column in \code{obs} for batch-aware HVG selection.
#'   Set to \code{NULL} for no batch correction.
#' @param hidden_size  Width of the shared backbone (default 512).
#' @param num_epochs  Maximum training epochs (default 50).
#' @param learning_rate  Learning rate (default 0.001).
#' @param batch_size  Training batch size (default 256).
#' @param loss_function  Loss function: \code{"mse"} (default) or
#'   \code{"kl"} (KL divergence).
#' @param device  Compute device: \code{"auto"} (default, uses GPU if
#'   available), \code{"cpu"}, or \code{"cuda"}.
#' @param head_hidden_size  Width of the proportion head (default 256).
#' @param n_resnet_blocks  Number of ResNet blocks (default 4).
#' @param input_dropout_rate  Dropout rate for input layer (default 0.25).
#' @param resnet_dropout_rate  Dropout rate for ResNet blocks (default 0.1).
#' @param head_dropout_rate  Dropout rate for proportion head (default 0.1).
#' @param test_size  Fraction of pseudo-bulk for validation (default 0.1).
#' @param early_stopping_patience  Epochs to wait before early stopping
#'   (default 10).
#' @param use_lr_scheduler  Use cosine annealing LR scheduler (default TRUE).
#' @param random_state  Random seed (default 42).
#' @param recon_weight  Weight for reconstruction loss relative to proportion
#'   loss (default 0.1). Higher values push the decoder to reconstruct input
#'   expression more faithfully. Set to 0 to disable reconstruction loss.
#'
#' @return A named list with:
#'   \describe{
#'     \item{model}{Trained Python DeconvModel object.}
#'     \item{var_genes}{Character vector of selected variable genes.}
#'     \item{cell_types}{Character vector of cell type names (output order).}
#'     \item{history}{data.frame with per-epoch training metrics
#'       (train_loss, val_loss, train_corr, val_corr).}
#'     \item{best_val_loss}{Best validation loss achieved.}
#'     \item{best_val_corr}{Best validation Pearson correlation.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#' use_scMMR_python(condaenv = "scMMR")
#'
#' # Train deconvolution model from scRNA-seq reference
#' result <- DNN_deconv_train(
#'   reference    = "reference.h5ad",
#'   label_column = "cell_type",
#'   save_path    = "models/deconv_model.pt",
#'   n_pseudobulk = 5000L,
#'   device       = "auto"
#' )
#' print(result$cell_types)
#' print(result$best_val_corr)
#' }
#'
#' @seealso \code{\link{DNN_deconv_predict}} for predicting proportions,
#'   \code{\link{DNN_train}} for cell type classification.
#'
#' @export
DNN_deconv_train <- function(reference,
                              label_column,
                              save_path               = "deconv_model.pt",
                              n_pseudobulk            = 5000L,
                              n_cells_per_sample      = 500L,
                              n_top_genes             = 6000L,
                              batch_key               = NULL,
                              hidden_size             = 512L,
                              num_epochs              = 50L,
                              learning_rate           = 0.001,
                              batch_size              = 256L,
                              loss_function           = "mse",
                              device                  = "auto",
                              head_hidden_size        = 256L,
                              n_resnet_blocks         = 4L,
                              input_dropout_rate      = 0.25,
                              resnet_dropout_rate     = 0.1,
                              head_dropout_rate       = 0.1,
                              test_size               = 0.1,
                              early_stopping_patience = 10L,
                              use_lr_scheduler        = TRUE,
                              random_state            = 42L,
                              recon_weight            = 0.1) {

  .ensure_python()

  # ── 0. Set device ───────────────────────────────────────────────────────────
  .scMMR_env$deconv_set_device(device)

  # ── 1. Load reference data ─────────────────────────────────────────────────
  message("\n== Loading reference data ==")
  adata <- .resolve_input(reference)

  # ── 2. Select highly variable genes ────────────────────────────────────────
  message("\n== Selecting highly variable genes ==")
  var_genes <- .select_hvg(adata, n_top = n_top_genes, batch_key = batch_key)
  message(sprintf("  %d HVGs selected", length(var_genes)))

  # ── 3. Show cell type distribution ─────────────────────────────────────────
  obs_meta <- as.data.frame(reticulate::py_to_r(adata$obs),
                            stringsAsFactors = FALSE)
  ct_table <- sort(table(obs_meta[[label_column]]), decreasing = TRUE)
  n_types  <- length(ct_table)
  message(sprintf("\n  Cell types : %d", n_types))
  message("  Class distribution:")
  for (i in seq_along(ct_table)) {
    message(sprintf("    %-30s %d", names(ct_table)[i], ct_table[i]))
  }

  # ── 4. Generate pseudo-bulk ─────────────────────────────────────────────────
  message("\n== Generating pseudo-bulk training data ==")
  pb_result <- .scMMR_env$deconv_generate_pseudobulk(
    adata              = adata,
    label_column       = label_column,
    n_samples          = as.integer(n_pseudobulk),
    n_cells_per_sample = as.integer(n_cells_per_sample),
    var_genes          = var_genes,
    random_state       = as.integer(random_state)
  )
  pb_r <- reticulate::py_to_r(pb_result)
  cell_types <- as.character(pb_r$cell_type_names)
  var_genes_used <- as.character(pb_r$var_genes)
  message(sprintf("  Generated %d pseudo-bulk samples", n_pseudobulk))
  message(sprintf("  Cell types (%d): %s", length(cell_types),
                  paste(cell_types, collapse = ", ")))

  # ── 5. Create model ─────────────────────────────────────────────────────────
  message("\n== Creating DeconvModel ==")
  input_size     <- length(var_genes_used)
  num_cell_types <- length(cell_types)
  model <- .scMMR_env$deconv_create_model(
    input_size              = as.integer(input_size),
    num_cell_types          = as.integer(num_cell_types),
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
    learning_rate           = learning_rate,
    loss_function           = loss_function,
    use_lr_scheduler        = use_lr_scheduler,
    random_state            = as.integer(random_state),
    recon_weight            = recon_weight
  )

  # Store var_genes on the model
  model$var_genes <- var_genes_used

  # ── 6. Train ────────────────────────────────────────────────────────────────
  message("\n== Training ==")
  result <- .scMMR_env$deconv_train(
    model           = model,
    X_pb            = pb_r$X_pb,
    y_pb            = pb_r$y_pb,
    cell_type_names = cell_types
  )
  result_r <- reticulate::py_to_r(result)
  message(sprintf("\n  Training complete (%d epochs)", result_r$n_epochs))
  message(sprintf("    Best val loss : %.5f", result_r$best_val_loss))
  message(sprintf("    Best val corr : %.4f", result_r$best_val_corr))

  # ── 7. Save model ──────────────────────────────────────────────────────────
  message("\n== Saving model ==")
  save_dir <- dirname(normalizePath(save_path, mustWork = FALSE))
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  .scMMR_env$deconv_save_model(model, save_path)
  message("  Saved: ", save_path)

  # ── 8. Training history ────────────────────────────────────────────────────
  hist_dict <- .scMMR_env$deconv_get_history(model)
  hist_r    <- reticulate::py_to_r(hist_dict)
  history   <- if (length(hist_r) > 0 && length(hist_r$train_loss) > 0) {
    as.data.frame(hist_r, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  invisible(list(
    model         = model,
    var_genes     = var_genes_used,
    cell_types    = cell_types,
    history       = history,
    best_val_loss = result_r$best_val_loss,
    best_val_corr = result_r$best_val_corr
  ))
}
