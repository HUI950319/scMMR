# ============================================================================
# Demo: Comprehensive Classification Benchmark
# ============================================================================
#
# Compares scMMR DNN against 8 ML models + 3 single-cell methods:
#
#   Deep Learning:
#     - scMMR DNN (pre-trained multi-task ResNet)
#
#   ML models (mlr3verse, 5-fold stratified CV, PCA 50-dim):
#     - Random Forest, XGBoost, SVM, Elastic Net
#     - LDA, Naive Bayes, KNN, Single-layer NNet
#
#   Single-cell methods (5-fold CV, same splits):
#     - SingleR        (correlation-based)
#     - Seurat LT      (CCA anchor-based label transfer)
#     - CellTypist      (logistic regression, Python)
#
# Usage:
#   source("demo_benchmark_classification.R")
# ============================================================================

library(scMMR)
library(Seurat)
library(qs)
library(ggplot2)
library(mlr3verse)
library(data.table)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

cat("=== Demo: Comprehensive Classification Benchmark ===\n\n")


# ── Helper ───────────────────────────────────────────────────────────────────

compute_metrics <- function(true_labels, pred_labels, ct_levels) {
  acc <- mean(true_labels == pred_labels, na.rm = TRUE)
  per_class <- sapply(ct_levels, function(ct) {
    idx <- true_labels == ct
    if (sum(idx) == 0) return(NA)
    mean(pred_labels[idx] == ct, na.rm = TRUE)
  })
  list(acc = acc, bacc = mean(per_class, na.rm = TRUE),
       ce = 1 - acc, per_class = per_class)
}


# ── 1. Load data ─────────────────────────────────────────────────────────────

test_path  <- system.file("extdata", "toy_test.qs", package = "scMMR")
model_path <- system.file("extdata", "model.pt",    package = "scMMR")

toy_test <- qread(test_path)
cat("Data:", ncol(toy_test), "cells x", nrow(toy_test), "genes\n")
cat("Cell types:", length(unique(toy_test@meta.data[["cell_type"]])), "\n")
print(sort(table(toy_test@meta.data[["cell_type"]]), decreasing = TRUE))
cat("\n")


# ── 2. DNN prediction ───────────────────────────────────────────────────────

# Save true labels by cell barcode BEFORE any processing
true_label_map <- setNames(
  as.character(toy_test@meta.data[["cell_type"]]),
  colnames(toy_test)
)

cat("=== [1/7] DNN prediction (pre-trained) ===\n")
pred_dnn <- DNN_predict(
  query = toy_test, model_path = model_path,
  true_label_col = "cell_type", device = "cpu"
)
# Save DNN predictions keyed by cell barcode
dnn_pred_map <- setNames(
  as.character(pred_dnn$predictions$cell_type_pred),
  colnames(toy_test)
)


# ── 3. Prepare features ─────────────────────────────────────────────────────

cat("=== [2/7] Preparing PCA features ===\n")
toy_test <- NormalizeData(toy_test, verbose = FALSE)
toy_test <- FindVariableFeatures(toy_test, nfeatures = 2000, verbose = FALSE)
toy_test <- ScaleData(toy_test, verbose = FALSE)
toy_test <- RunPCA(toy_test, npcs = 50, verbose = FALSE)
pca_mat <- Embeddings(toy_test, "pca")
cat("PCA:", nrow(pca_mat), "x", ncol(pca_mat), "\n")

# Reconstruct true labels and DNN predictions in POST-PCA cell order
dnn_true <- true_label_map[colnames(toy_test)]
dnn_pred <- dnn_pred_map[colnames(toy_test)]
cat("DNN accuracy:", round(mean(dnn_true == dnn_pred), 4), "\n\n")


# ── 4. mlr3verse benchmark (8 ML models) ────────────────────────────────────

cat("=== [3/7] mlr3verse benchmark (8 ML models, 5-fold CV) ===\n")

labels <- factor(dnn_true)
task_df <- as.data.frame(pca_mat)
task_df$cell_type <- labels

task <- TaskClassif$new(
  id = "cell_type", backend = task_df, target = "cell_type"
)
task$col_roles$stratum <- "cell_type"

learners <- list(
  lrn("classif.ranger",     id = "RF",         num.trees = 500),
  lrn("classif.xgboost",    id = "XGBoost",    nrounds = 100, max_depth = 6, verbose = 0),
  lrn("classif.svm",        id = "SVM",        type = "C-classification", kernel = "radial"),
  lrn("classif.glmnet",     id = "ElasticNet",  alpha = 0.5),
  lrn("classif.lda",        id = "LDA"),
  lrn("classif.naive_bayes", id = "NaiveBayes"),
  lrn("classif.kknn",       id = "KNN",        k = 15),
  lrn("classif.nnet",       id = "NNet",       size = 20, MaxNWts = 50000, maxit = 200, trace = FALSE)
)

cv5 <- rsmp("cv", folds = 5)
design <- benchmark_grid(task, learners, cv5)
cat("Running 8 learners x 5 folds...\n")
suppressWarnings(suppressMessages(bmr <- benchmark(design, store_models = FALSE)))

measures <- msrs(c("classif.acc", "classif.bacc", "classif.ce"))
agg <- bmr$aggregate(measures)
cat("mlr3 complete.\n\n")

# Extract CV fold assignments
rs_dt <- as.data.table(bmr$resamplings)
rs <- rs_dt$resampling[[1]]

# Extract combined predictions for per-class metrics
scores_dt <- as.data.table(bmr$score(msr("classif.acc")))
ml_combined_preds <- list()
for (lid in unique(scores_dt$learner_id)) {
  fold_sc <- scores_dt[learner_id == lid]
  preds <- character(ncol(toy_test))
  for (j in seq_len(nrow(fold_sc))) {
    p <- fold_sc$prediction[[j]]
    ti <- fold_sc$iteration[j]
    idx <- fold_sc$resampling[[j]]$test_set(ti)
    preds[idx] <- as.character(p$response)
  }
  ml_combined_preds[[lid]] <- preds
}


# ── 5. Single-cell methods (5-fold CV, same folds) ──────────────────────────

norm_data <- GetAssayData(toy_test, layer = "data")

# Use named vectors keyed by cell barcode to guarantee correct pairing
cell_names_all <- colnames(toy_test)
n_cells <- ncol(toy_test)

singler_preds    <- setNames(character(n_cells), cell_names_all)
seurat_lt_preds  <- setNames(character(n_cells), cell_names_all)
celltypist_preds <- setNames(character(n_cells), cell_names_all)
# Collect true labels during CV (from each query subset, guaranteeing pairing)
sc_true_labels   <- setNames(character(n_cells), cell_names_all)

# Check availability
singler_ok    <- requireNamespace("SingleR", quietly = TRUE) &&
                 requireNamespace("SingleCellExperiment", quietly = TRUE)
celltypist_ok <- tryCatch({
  reticulate::import("celltypist"); TRUE
}, error = function(e) FALSE)

# ── 5a. SingleR ─────────────────────────────────────────────────────────────

if (singler_ok) {
  cat("=== [4/7] SingleR (5-fold CV) ===\n")
  library(SingleR)
  library(SingleCellExperiment)

  for (fold in 1:5) {
    cat(sprintf("  fold %d/5...\n", fold))
    ti <- rs$test_set(fold); tri <- rs$train_set(fold)

    ref_sce <- SingleCellExperiment(assays = list(logcounts = norm_data[, tri]))
    ref_sce$label <- as.character(dnn_true[tri])
    query_sce <- SingleCellExperiment(assays = list(logcounts = norm_data[, ti]))

    sr <- SingleR(test = query_sce, ref = ref_sce,
                  labels = ref_sce$label, de.method = "wilcox")
    query_cells <- colnames(toy_test)[ti]
    singler_preds[query_cells] <- sr$labels
    sc_true_labels[query_cells] <- as.character(dnn_true[ti])
  }
  cat("  done.\n\n")
}

# ── 5b. Seurat Label Transfer ───────────────────────────────────────────────

cat("=== [5/7] Seurat Label Transfer (5-fold CV) ===\n")
for (fold in 1:5) {
  cat(sprintf("  fold %d/5...\n", fold))
  ti <- rs$test_set(fold); tri <- rs$train_set(fold)

  ref_s   <- toy_test[, tri]
  query_s <- toy_test[, ti]

  anchors <- FindTransferAnchors(
    reference = ref_s, query = query_s,
    dims = 1:30, verbose = FALSE
  )
  ref_labels <- as.character(ref_s@meta.data[["cell_type"]])
  transferred <- TransferData(
    anchorset = anchors, refdata = ref_labels,
    dims = 1:30, verbose = FALSE
  )
  preds_fold <- as.character(transferred$predicted.id)
  true_fold  <- as.character(query_s@meta.data[["cell_type"]])
  # Use cell barcodes from the query subset
  query_cells <- colnames(query_s)
  seurat_lt_preds[query_cells] <- preds_fold
  sc_true_labels[query_cells] <- true_fold
}
cat("  done.\n\n")

# ── 5c. CellTypist ─────────────────────────────────────────────────────────

if (celltypist_ok) {
  cat("=== [6/7] CellTypist (5-fold CV) ===\n")
  ct <- reticulate::import("celltypist")
  sc <- reticulate::import("scanpy")
  np <- reticulate::import("numpy")
  pd <- reticulate::import("pandas")

  for (fold in 1:5) {
    cat(sprintf("  fold %d/5...\n", fold))
    ti <- rs$test_set(fold); tri <- rs$train_set(fold)

    query_cells <- cell_names_all[ti]
    ref_true <- as.character(dnn_true[tri])
    query_true <- as.character(dnn_true[ti])

    # Build AnnData for reference
    ref_mat <- as.matrix(t(norm_data[, tri]))  # cells x genes
    ref_adata <- sc$AnnData(
      X = ref_mat,
      obs = pd$DataFrame(list(cell_type = ref_true))
    )
    ref_adata$var_names <- rownames(norm_data)

    # Build AnnData for query
    query_mat <- as.matrix(t(norm_data[, ti]))
    query_adata <- sc$AnnData(
      X = query_mat,
      obs = pd$DataFrame(list(cell_type = query_true))
    )
    query_adata$var_names <- rownames(norm_data)

    tryCatch({
      # Train custom model (use_SGD=FALSE for better convergence on small data)
      model <- ct$train(
        ref_adata,
        labels = "cell_type",
        use_SGD = FALSE,
        n_jobs = 4L,
        max_iter = 200L
      )

      # Predict
      preds <- ct$annotate(
        query_adata,
        model = model,
        majority_voting = FALSE
      )

      # Extract predictions: reticulate::py_to_r converts to data.frame
      pred_df <- reticulate::py_to_r(preds$predicted_labels)
      celltypist_preds[query_cells] <- as.character(pred_df[["predicted_labels"]])
    }, error = function(e) {
      cat("    CellTypist error:", e$message, "\n")
    })
  }
  cat("  done.\n\n")
}


# ── 6. Compile all results ──────────────────────────────────────────────────

cat("=== [7/7] Compiling results ===\n\n")

ct_levels <- sort(unique(dnn_true))

# DNN metrics
m_dnn <- compute_metrics(dnn_true, dnn_pred, ct_levels)

# ML metrics
ml_metrics <- lapply(ml_combined_preds, function(p) {
  compute_metrics(dnn_true, p, ct_levels)
})

# Single-cell metrics (use sc_true_labels collected during CV for correct pairing)
sc_methods <- list()
if (singler_ok) {
  sc_methods[["SingleR"]] <- compute_metrics(sc_true_labels, singler_preds, ct_levels)
}
sc_methods[["Seurat_LT"]] <- compute_metrics(sc_true_labels, seurat_lt_preds, ct_levels)
if (celltypist_ok && any(celltypist_preds != "")) {
  sc_methods[["CellTypist"]] <- compute_metrics(sc_true_labels, celltypist_preds, ct_levels)
}


# ── 7. Results table ─────────────────────────────────────────────────────────


cat("================================================================\n")
cat("     Comprehensive Classification Benchmark Results\n")
cat("================================================================\n\n")

results <- data.frame(
  Method = "DNN (scMMR)", Accuracy = m_dnn$acc,
  Balanced_Acc = m_dnn$bacc, Error_Rate = m_dnn$ce,
  Type = "Deep Learning", stringsAsFactors = FALSE
)

for (lid in names(ml_metrics)) {
  results <- rbind(results, data.frame(
    Method = lid, Accuracy = ml_metrics[[lid]]$acc,
    Balanced_Acc = ml_metrics[[lid]]$bacc, Error_Rate = ml_metrics[[lid]]$ce,
    Type = "ML (CV)"
  ))
}

for (nm in names(sc_methods)) {
  results <- rbind(results, data.frame(
    Method = nm, Accuracy = sc_methods[[nm]]$acc,
    Balanced_Acc = sc_methods[[nm]]$bacc, Error_Rate = sc_methods[[nm]]$ce,
    Type = "Single-cell"
  ))
}

results <- results[order(-results$Accuracy), ]
results$Rank <- seq_len(nrow(results))

cat(sprintf("%-20s %-15s %10s %14s %6s\n",
            "Method", "Type", "Accuracy", "Balanced_Acc", "Rank"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("%-20s %-15s %10.4f %14.4f %6d\n",
              results$Method[i], results$Type[i],
              results$Accuracy[i], results$Balanced_Acc[i],
              results$Rank[i]))
}

cat(sprintf("\nTotal: %d methods compared\n", nrow(results)))
cat("  DNN: pre-trained on separate training data (all genes)\n")
cat("  ML:  5-fold stratified CV on PCA 50-dim\n")
cat("  SC:  5-fold CV, training fold as reference\n\n")


# ── 8. Per-class recall ─────────────────────────────────────────────────────

cat("--- Per-class recall (top methods) ---\n\n")

all_per_class <- list(DNN = m_dnn$per_class)
for (lid in names(ml_metrics)) all_per_class[[lid]] <- ml_metrics[[lid]]$per_class
for (nm in names(sc_methods))  all_per_class[[nm]]  <- sc_methods[[nm]]$per_class

# Show top 4 ML + all single-cell
top_ml <- names(sort(sapply(ml_metrics, function(m) m$acc), decreasing = TRUE))[1:4]
display <- c("DNN", top_ml, names(sc_methods))

cat(sprintf("%-22s %4s", "Cell Type", "N"))
for (m in display) cat(sprintf(" %10s", m))
cat("\n")
cat(strrep("-", 26 + 10 * length(display)), "\n")
for (ct in ct_levels) {
  cat(sprintf("%-22s %4d", ct, sum(dnn_true == ct)))
  for (m in display) cat(sprintf(" %10.3f", all_per_class[[m]][ct]))
  cat("\n")
}


# ── 9. Confusion matrix (DNN) ───────────────────────────────────────────────

cat("\n--- DNN Confusion Matrix ---\n")
print(table(True = dnn_true, Predicted = dnn_pred))


# ── 10. Visualization ───────────────────────────────────────────────────────

cat("\n--- Generating plots ---\n")

# Color scheme
type_colors <- c("Deep Learning" = "#E74C3C", "ML (CV)" = "#3498DB",
                 "Single-cell" = "#2ECC71")

# 10a. Accuracy bar chart
plot_df <- results
plot_df$Method <- factor(plot_df$Method,
                         levels = plot_df$Method[order(plot_df$Accuracy)])

p_acc <- ggplot(plot_df, aes(x = Method, y = Accuracy, fill = Type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Accuracy)), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = type_colors) +
  scale_y_continuous(limits = c(0, 1.1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Cell Type Classification Benchmark",
       subtitle = sprintf("%d methods | %d types | %d cells",
                          nrow(results), length(ct_levels), length(dnn_true)),
       x = NULL, y = "Accuracy") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
print(p_acc)

# 10b. Balanced accuracy bar chart
p_bacc <- ggplot(plot_df, aes(x = Method, y = Balanced_Acc, fill = Type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Balanced_Acc)), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = type_colors) +
  scale_y_continuous(limits = c(0, 1.1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Balanced Accuracy (mean per-class recall)",
       x = NULL, y = "Balanced Accuracy") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
print(p_bacc)

# 10c. Per-class heatmap
heat_list <- lapply(display, function(m) {
  data.frame(cell_type = ct_levels, method = m,
             recall = all_per_class[[m]], stringsAsFactors = FALSE)
})
heat_df <- do.call(rbind, heat_list)
heat_df$method <- factor(heat_df$method, levels = display)
heat_df$text_col <- ifelse(heat_df$recall > 0.6, "white", "black")

p_heat <- ggplot(heat_df, aes(x = method, y = cell_type, fill = recall)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", recall)),
            color = heat_df$text_col, size = 2.5) +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#fee0d2", "#fc9272", "#de2d26", "#a50f15"),
    limits = c(0, 1), name = "Recall"
  ) +
  labs(title = "Per-class Recall Heatmap", x = "Method", y = "Cell Type") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
print(p_heat)

# 10d. Acc vs BAcc scatter
p_scatter <- ggplot(results, aes(x = Accuracy, y = Balanced_Acc,
                                  color = Type, label = Method)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 20) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = type_colors) +
  labs(title = "Accuracy vs Balanced Accuracy",
       x = "Accuracy", y = "Balanced Accuracy") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
print(p_scatter)


# ── 11. Summary ──────────────────────────────────────────────────────────────

cat("\n================================================================\n")
cat("                         Summary\n")
cat("================================================================\n\n")

dnn_rank <- results$Rank[results$Method == "DNN (scMMR)"]
cat(sprintf("Best method: %s (Acc=%.4f)\n", results$Method[1], results$Accuracy[1]))
cat(sprintf("DNN rank: %d / %d\n\n", dnn_rank, nrow(results)))

ml_res <- results[results$Type == "ML (CV)", ]
sc_res <- results[results$Type == "Single-cell", ]
cat("Category winners:\n")
cat(sprintf("  Best ML:          %-15s Acc=%.4f\n", ml_res$Method[1], ml_res$Accuracy[1]))
cat(sprintf("  Best Single-cell: %-15s Acc=%.4f\n", sc_res$Method[1], sc_res$Accuracy[1]))
cat(sprintf("\nDNN vs best ML:     %+.1f%%\n", (m_dnn$acc - ml_res$Accuracy[1]) * 100))
cat(sprintf("DNN vs best SC:     %+.1f%%\n", (m_dnn$acc - sc_res$Accuracy[1]) * 100))

cat("\nDNN unique advantages:\n")
cat("  1. 512-dim shared embedding for downstream analysis\n")
cat("  2. Perturbation ranking (RankPerturbation)\n")
cat("  3. Differential abundance (RankPercent)\n")
cat("  4. Multi-task learning: classification + embedding regression\n")

cat("\n=== Benchmark Complete ===\n")
