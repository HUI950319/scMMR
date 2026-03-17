# ============================================================================
# scMMR Demo: Bulk RNA-seq Deconvolution
# ============================================================================
#
# This script demonstrates the scMMR deconvolution workflow:
#   1. Load scRNA-seq reference data (Baron pancreas)
#   2. Train a deconvolution DNN model (DNN_deconv_train)
#   3. Generate test pseudo-bulk with known proportions
#   4. Predict cell type proportions (DNN_deconv_predict)
#   5. Evaluate: correlation, RMSE, visualization
#
# Requirements:
#   - scMMR installed
#   - Python environment configured (conda: scanpy-env)
#   - scRNAseq Bioconductor package (for Baron pancreas data)
#   - ggplot2, tidyr, pheatmap (for visualization)
#
# Expected run time: ~3-5 minutes (CPU), ~1-2 minutes (GPU)
# ============================================================================

library(scMMR)
library(ggplot2)

# ── 0. Setup Python environment ─────────────────────────────────────────────
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# ── 1. Prepare reference scRNA-seq data ─────────────────────────────────────
# Using Baron human pancreas dataset from Bioconductor
if (!requireNamespace("scRNAseq", quietly = TRUE)) {
  stop("Install scRNAseq: BiocManager::install('scRNAseq')")
}

message("\n========================================")
message("  Step 1: Loading reference data")
message("========================================")

sce <- scRNAseq::BaronPancreasData("human")

# Filter low-quality cells (keep cell types with >= 50 cells)
ct_counts <- table(SummarizedExperiment::colData(sce)[["label"]])
keep_types <- names(ct_counts[ct_counts >= 50])
keep_idx <- which(SummarizedExperiment::colData(sce)[["label"]] %in% keep_types)
sce <- sce[, keep_idx]

cat(sprintf("  Reference: %d cells x %d genes\n", ncol(sce), nrow(sce)))
cat(sprintf("  Cell types (%d): %s\n",
            length(keep_types), paste(keep_types, collapse = ", ")))
cat("\n  Cell type distribution:\n")
ct_table <- sort(table(SummarizedExperiment::colData(sce)[["label"]]),
                 decreasing = TRUE)
for (i in seq_along(ct_table)) {
  cat(sprintf("    %-20s %d\n", names(ct_table)[i], ct_table[i]))
}

# Save as h5ad for scMMR (Seurat also works - see Alternative below)
tmp_h5ad <- tempfile(fileext = ".h5ad")

# Convert SCE → AnnData via reticulate
counts_mat <- as.matrix(SingleCellExperiment::counts(sce))
obs_df <- as.data.frame(SummarizedExperiment::colData(sce))
anndata <- reticulate::import("anndata")
adata <- anndata$AnnData(
  X = reticulate::r_to_py(t(counts_mat)),   # cells x genes
  obs = reticulate::r_to_py(obs_df)
)
adata$var_names <- rownames(counts_mat)
adata$obs_names <- colnames(counts_mat)
adata$write_h5ad(tmp_h5ad)
cat(sprintf("\n  Saved temp h5ad: %s\n", tmp_h5ad))

# ── Alternative: Seurat object as input ──────────────────────────────────────
# If you have a Seurat object, pass it directly:
#   result <- DNN_deconv_train(reference = seurat_obj, ...)

# ============================================================================
# Step 2: Train deconvolution model
# ============================================================================
message("\n========================================")
message("  Step 2: Training deconvolution model")
message("========================================")

model_path <- file.path(tempdir(), "deconv_baron.pt")

result <- DNN_deconv_train(
  reference          = tmp_h5ad,
  label_column       = "label",
  save_path          = model_path,

  # Pseudo-bulk settings
  n_pseudobulk       = 3000L,        # samples for training
  n_cells_per_sample = 500L,         # cells mixed per pseudo-bulk

  # Gene selection
  n_top_genes        = 3000L,        # HVGs

  # Model architecture (defaults work well)
  hidden_size        = 512L,
  n_resnet_blocks    = 4L,

  # Training
  num_epochs         = 30L,
  learning_rate      = 0.001,
  batch_size         = 256L,
  loss_function      = "mse",        # "mse" or "kl"

  # Device
  device             = "auto"        # auto-detect GPU
)

# ── Training summary ────────────────────────────────────────────────────────
cat("\n== Training Results ==\n")
cat(sprintf("  Model saved: %s\n", model_path))
cat(sprintf("  HVGs used: %d\n", length(result[["var_genes"]])))
cat(sprintf("  Cell types: %s\n",
            paste(result[["cell_types"]], collapse = ", ")))
cat(sprintf("  Best val loss: %.5f\n", result[["best_val_loss"]]))
cat(sprintf("  Best val corr: %.4f\n", result[["best_val_corr"]]))

# ── Plot training curves ───────────────────────────────────────────────────
history <- result[["history"]]
if (nrow(history) > 0) {
  history$epoch <- seq_len(nrow(history))

  # Loss curve
  p_loss <- ggplot(history, aes(x = epoch)) +
    geom_line(aes(y = train_loss, color = "Train"), linewidth = 0.8) +
    geom_line(aes(y = val_loss, color = "Validation"), linewidth = 0.8) +
    scale_color_manual(values = c("Train" = "#2196F3", "Validation" = "#F44336")) +
    labs(title = "Deconvolution Training Loss",
         x = "Epoch", y = "Loss (MSE)", color = NULL) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8))
  print(p_loss)

  # Correlation curve
  p_corr <- ggplot(history, aes(x = epoch)) +
    geom_line(aes(y = train_corr, color = "Train"), linewidth = 0.8) +
    geom_line(aes(y = val_corr, color = "Validation"), linewidth = 0.8) +
    scale_color_manual(values = c("Train" = "#2196F3", "Validation" = "#F44336")) +
    labs(title = "Deconvolution Training Correlation",
         x = "Epoch", y = "Pearson r", color = NULL) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.2))
  print(p_corr)
}


# ============================================================================
# Step 3: Generate test pseudo-bulk with KNOWN proportions
# ============================================================================
message("\n========================================")
message("  Step 3: Generating test pseudo-bulk")
message("========================================")

# Use Python helper to generate test samples (separate from training data)
.ensure_python <- scMMR:::.ensure_python
.ensure_python()
env <- scMMR:::.scMMR_env

# Load reference again for generating test data
adata_test <- env$mt_load_query(tmp_h5ad)

# Generate 20 test pseudo-bulk with different seed
test_pb <- env$deconv_generate_pseudobulk(
  adata              = adata_test,
  label_column       = "label",
  n_samples          = 20L,
  n_cells_per_sample = 500L,
  var_genes          = result[["var_genes"]],
  random_state       = 999L        # different seed from training!
)

test_pb_r <- reticulate::py_to_r(test_pb)

cat(sprintf("  Test samples: %d\n", nrow(test_pb_r$X_pb)))
cat(sprintf("  True proportions matrix: %d x %d\n",
            nrow(test_pb_r$y_pb), ncol(test_pb_r$y_pb)))

# Build a bulk expression matrix (genes x samples) for DNN_deconv_predict
# Note: X_pb is already log1p(CPM). DNN_deconv_predict will re-normalize
# its input to log1p(CPM), so we need to pass CPM values (expm1 reverses log1p)
test_matrix <- t(expm1(test_pb_r$X_pb))   # CPM, genes x samples
rownames(test_matrix) <- test_pb_r$var_genes
colnames(test_matrix) <- paste0("sample_", seq_len(ncol(test_matrix)))

# Save as CSV for testing CSV input path
test_csv <- file.path(tempdir(), "test_bulk.csv")
write.csv(test_matrix, test_csv)
cat(sprintf("  Saved test CSV: %s\n", test_csv))


# ============================================================================
# Step 4: Predict cell type proportions
# ============================================================================
message("\n========================================")
message("  Step 4: Predicting proportions")
message("========================================")

# ── 4a. Predict from matrix ────────────────────────────────────────────────
props_mat <- DNN_deconv_predict(
  bulk_expr  = test_matrix,
  model_path = model_path,
  device     = "auto"
)

cat("\n  Predicted proportions (first 5 samples):\n")
print(head(props_mat, 5))

# ── 4b. Predict from CSV file ─────────────────────────────────────────────
props_csv <- DNN_deconv_predict(
  bulk_expr  = test_csv,
  model_path = model_path,
  device     = "auto"
)

cat("\n  CSV input also works — row sums:\n")
cat(sprintf("    %s\n", paste(round(rowSums(props_csv[, -1]), 4), collapse = ", ")))


# ============================================================================
# Step 5: Evaluation
# ============================================================================
message("\n========================================")
message("  Step 5: Evaluation")
message("========================================")

cell_types <- result[["cell_types"]]

# ── 5a. Check row sums = 1 ────────────────────────────────────────────────
row_sums <- rowSums(props_mat[, -1])
cat(sprintf("\n  Row sum range: [%.6f, %.6f] (should be ~1.0)\n",
            min(row_sums), max(row_sums)))
stopifnot(all(abs(row_sums - 1.0) < 1e-4))
cat("  ✓ All row sums ≈ 1.0\n")

# ── 5b. Per-cell-type Pearson correlation ─────────────────────────────────
true_props <- test_pb_r$y_pb
pred_props <- as.matrix(props_mat[, -1])    # remove sample_id column

cat("\n  Per-cell-type Pearson correlation:\n")
corr_vec <- numeric(length(cell_types))
rmse_vec <- numeric(length(cell_types))
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  r <- cor(true_props[, i], pred_props[, i])
  rmse <- sqrt(mean((true_props[, i] - pred_props[, i])^2))
  corr_vec[i] <- r
  rmse_vec[i] <- rmse
  cat(sprintf("    %-20s  r = %.4f  RMSE = %.4f\n", ct, r, rmse))
}

# Overall
overall_r <- cor(as.vector(true_props), as.vector(pred_props))
overall_rmse <- sqrt(mean((true_props - pred_props)^2))
cat(sprintf("\n  Overall:  r = %.4f,  RMSE = %.4f\n", overall_r, overall_rmse))

# Summary table
eval_df <- data.frame(
  cell_type = cell_types,
  pearson_r = round(corr_vec, 4),
  RMSE      = round(rmse_vec, 4),
  stringsAsFactors = FALSE
)
cat("\n  Summary:\n")
print(eval_df)


# ============================================================================
# Step 6: Visualization
# ============================================================================
message("\n========================================")
message("  Step 6: Visualization")
message("========================================")

# ── 6a. Scatter plot: predicted vs true (all cell types) ─────────────────
scatter_df <- data.frame(
  true      = as.vector(true_props),
  predicted = as.vector(pred_props),
  cell_type = rep(cell_types, each = nrow(true_props))
)

p_scatter <- ggplot(scatter_df, aes(x = true, y = predicted, color = cell_type)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = sprintf("Deconvolution: Predicted vs True (r = %.3f)", overall_r),
       x = "True Proportion", y = "Predicted Proportion", color = "Cell Type") +
  scale_x_continuous(limits = c(0, max(scatter_df$true, scatter_df$predicted) * 1.05)) +
  scale_y_continuous(limits = c(0, max(scatter_df$true, scatter_df$predicted) * 1.05)) +
  theme_classic() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14))
print(p_scatter)

# ── 6b. Per-cell-type scatter (faceted) ──────────────────────────────────
p_facet <- ggplot(scatter_df, aes(x = true, y = predicted)) +
  geom_point(aes(color = cell_type), alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ cell_type, scales = "free") +
  labs(title = "Per-Cell-Type Deconvolution Accuracy",
       x = "True Proportion", y = "Predicted Proportion") +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#f0f0f0"),
        plot.title = element_text(hjust = 0.5))
print(p_facet)

# ── 6c. Heatmap: true vs predicted side by side ─────────────────────────
# Combine into long format for heatmap
heatmap_true <- reshape2::melt(true_props)
heatmap_pred <- reshape2::melt(pred_props)
heatmap_true$type <- "True"
heatmap_pred$type <- "Predicted"
heatmap_true$cell_type <- rep(cell_types, each = nrow(true_props))
heatmap_pred$cell_type <- rep(cell_types, each = nrow(pred_props))
heatmap_true$sample <- rep(paste0("S", seq_len(nrow(true_props))),
                           times = length(cell_types))
heatmap_pred$sample <- rep(paste0("S", seq_len(nrow(pred_props))),
                           times = length(cell_types))
heatmap_df <- rbind(heatmap_true[, c("sample", "cell_type", "value", "type")],
                    heatmap_pred[, c("sample", "cell_type", "value", "type")])
heatmap_df$type <- factor(heatmap_df$type, levels = c("True", "Predicted"))

p_heat <- ggplot(heatmap_df, aes(x = sample, y = cell_type, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ type) +
  scale_fill_gradient2(low = "white", mid = "#F4A460", high = "#8B0000",
                       midpoint = 0.2, name = "Proportion") +
  labs(title = "Cell Type Proportions: True vs Predicted",
       x = "Sample", y = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(hjust = 0.5, size = 14),
        strip.text = element_text(size = 12, face = "bold"))
print(p_heat)

# ── 6d. Bar chart: per-cell-type accuracy metrics ───────────────────────
p_bar <- ggplot(eval_df, aes(x = reorder(cell_type, -pearson_r), y = pearson_r)) +
  geom_col(fill = "#4CAF50", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", pearson_r)),
            vjust = -0.3, size = 3) +
  labs(title = "Per-Cell-Type Deconvolution Accuracy",
       x = "Cell Type", y = "Pearson r") +
  ylim(0, 1.05) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14))
print(p_bar)

# ── 6e. Stacked bar: predicted proportions per sample ────────────────────
prop_long <- tidyr::pivot_longer(
  props_mat,
  cols      = -sample_id,
  names_to  = "cell_type",
  values_to = "proportion"
)

p_stack <- ggplot(prop_long, aes(x = sample_id, y = proportion, fill = cell_type)) +
  geom_col(width = 0.8) +
  labs(title = "Predicted Cell Type Proportions per Sample",
       x = "Sample", y = "Proportion", fill = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(hjust = 0.5, size = 14))
print(p_stack)


# ============================================================================
# Step 7: Compare with different loss function (optional)
# ============================================================================
# Uncomment to compare MSE vs KL divergence loss:
#
# result_kl <- DNN_deconv_train(
#   reference     = tmp_h5ad,
#   label_column  = "label",
#   save_path     = file.path(tempdir(), "deconv_baron_kl.pt"),
#   n_pseudobulk  = 3000L,
#   n_top_genes   = 3000L,
#   num_epochs    = 30L,
#   loss_function = "kl",
#   device        = "auto"
# )
# cat(sprintf("KL: val_corr = %.4f vs MSE: val_corr = %.4f\n",
#             result_kl[["best_val_corr"]], result[["best_val_corr"]]))

# ============================================================================
# Step 8: Clean up
# ============================================================================
message("\n========================================")
message("  Cleanup")
message("========================================")

file.remove(tmp_h5ad)
file.remove(test_csv)
cat("  Temp files removed.\n")
cat(sprintf("  Model kept at: %s\n", model_path))

message("\n========================================")
message("  Demo complete!")
message(sprintf("  Overall: Pearson r = %.4f, RMSE = %.4f", overall_r, overall_rmse))
message("========================================\n")

