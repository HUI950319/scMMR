# ══════════════════════════════════════════════════════════════════════════════
# test_scMMR.R — Quick test for scMMR package
# ══════════════════════════════════════════════════════════════════════════════

library(scMMR)

# Set Python environment
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# ── Test 1: DNN_train ────────────────────────────────────────────────────────
cat("\n========== TEST 1: DNN_train ==========\n")

result <- DNN_train(
  input         = "/home/oyh/project/para/python_pacakge/data/toy_train.h5ad",
  label_column  = "cell_type",
  embedding_key = "umap",
  save_path     = "/tmp/scMMR_test_model.pt",
  n_top_genes   = 6000L,
  batch_key     = "sample",
  num_epochs    = 50L,
  device        = "auto"
)

cat("\n--- Training Results ---\n")
cat(sprintf("Best val accuracy: %.2f%%\n", result$best_val_acc))
cat(sprintf("Best val RMSE: %.4f\n", result$best_val_rmse))
cat(sprintf("Epochs: %d\n", nrow(result$history)))
cat(sprintf("HVGs: %d\n", length(result$var_genes)))
cat("\nHistory (last 5 rows):\n")
print(tail(result$history[, c("epoch", "train_acc", "val_acc",
                               "train_rmse", "val_rmse")], 5))

# ── Test 2: DNN_predict ──────────────────────────────────────────────────────
cat("\n========== TEST 2: DNN_predict ==========\n")

pred <- DNN_predict(
  query      = "/home/oyh/project/para/python_pacakge/data/toy_test.h5ad",
  model_path = "/tmp/scMMR_test_model.pt",
  explain    = FALSE,
  device     = "auto"
)

cat("\n--- Prediction Results ---\n")
cat(sprintf("Cells predicted: %d\n", nrow(pred$predictions)))
cat(sprintf("Shared embedding dim: %s\n",
            paste(dim(pred$shared_embedding), collapse = " x ")))
cat("\nPredictions (first 6 rows):\n")
print(head(pred$predictions[, c("cell_id", "cell_type_pred",
                                 "confidence", "is_ood")]))

cat("\n========== ALL TESTS PASSED ==========\n")
