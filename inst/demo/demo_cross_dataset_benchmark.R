# ============================================================================
# Demo: Cross-Dataset Pancreas Benchmark
# ============================================================================
#
# Evaluates scMMR DNN against 11 other methods on 4 independent pancreas
# datasets (Baron, Muraro, Segerstolpe, Xin). Each dataset is benchmarked
# independently using 5-fold stratified CV.
#
# Methods (12 total):
#   Deep Learning:  DNN (scMMR)
#   ML (mlr3):      RF, XGBoost, SVM, ElasticNet, LDA, NaiveBayes, KNN, NNet
#   Single-cell:    SingleR, Seurat_LT, CellTypist
#
# Output: 4 × 12 accuracy heatmap
#
# Usage:
#   source("demo_cross_dataset_benchmark.R")
# ============================================================================

library(scMMR)
library(Seurat)
library(qs)
library(ggplot2)
library(mlr3verse)
library(data.table)
library(scRNAseq)
library(SingleR)
library(SingleCellExperiment)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# ── Quick test mode: set to TRUE to subsample each dataset for fast debugging
QUICK_TEST <- FALSE      # TRUE → 500 cells/dataset, 3-fold, num_epochs=10
QUICK_N    <- 500L
QUICK_FOLDS <- 3L
QUICK_EPOCHS <- 10L

cat("=== Demo: Cross-Dataset Pancreas Benchmark ===\n")
if (QUICK_TEST) cat("  *** QUICK TEST MODE (subsample) ***\n")
cat("\n")


# ── Helper functions ─────────────────────────────────────────────────────────

compute_accuracy <- function(true_labels, pred_labels) {
  valid <- !is.na(true_labels) & !is.na(pred_labels)
  if (sum(valid) == 0) return(NA_real_)
  mean(true_labels[valid] == pred_labels[valid])
}

#' Convert SCE to Seurat with harmonized cell type labels
sce_to_seurat <- function(sce, dataset_name) {
  # Get raw counts (fallback to first available assay, e.g. rpkm for Xin)
  if ("counts" %in% assayNames(sce)) {
    cts <- counts(sce)
  } else {
    avail <- assayNames(sce)[1]
    cat(sprintf("    Note: no 'counts' assay, using '%s' instead\n", avail))
    cts <- assay(sce, avail)
  }

  # Harmonize gene names to symbols
  if (dataset_name %in% c("Muraro", "Xin")) {
    # Muraro: rownames are "GENE__chrN", Xin: rownames are Entrez IDs
    # Both need rowData$symbol for proper gene names
    syms <- as.character(rowData(sce)$symbol)
    # Seurat converts underscores to dashes, so do this first
    syms <- gsub("_", "-", syms)
    valid <- !is.na(syms) & syms != "" & !duplicated(syms)
    cts <- cts[valid, ]
    rownames(cts) <- syms[valid]
  }

  # Get cell type labels
  label_map <- list(
    Baron = list(
      col = "label",
      map = c(
        "acinar" = "acinar", "alpha" = "alpha", "beta" = "beta",
        "delta" = "delta", "ductal" = "ductal", "endothelial" = "endothelial",
        "epsilon" = "epsilon", "gamma" = "gamma",
        "activated_stellate" = "stellate", "quiescent_stellate" = "stellate",
        "macrophage" = "macrophage", "mast" = "mast",
        "schwann" = "schwann", "t_cell" = "t_cell"
      )
    ),
    Muraro = list(
      col = "label",
      map = c(
        "acinar" = "acinar", "alpha" = "alpha", "beta" = "beta",
        "delta" = "delta", "duct" = "ductal", "endothelial" = "endothelial",
        "epsilon" = "epsilon", "pp" = "gamma", "mesenchymal" = "mesenchyme"
      )
    ),
    Segerstolpe = list(
      col = "cell type",
      map = c(
        "acinar cell" = "acinar", "alpha cell" = "alpha", "beta cell" = "beta",
        "delta cell" = "delta", "ductal cell" = "ductal",
        "endothelial cell" = "endothelial", "epsilon cell" = "epsilon",
        "gamma cell" = "gamma", "PSC cell" = "stellate", "mast cell" = "mast"
      )
    ),
    Xin = list(
      col = "cell.type",
      map = c(
        "alpha" = "alpha", "beta" = "beta",
        "delta" = "delta", "PP" = "gamma"
      )
    )
  )

  info <- label_map[[dataset_name]]
  raw_labels <- as.character(colData(sce)[[info$col]])
  mapped <- info$map[raw_labels]

  # Keep only cells with valid mapped labels
  keep <- !is.na(mapped)
  cts <- cts[, keep]
  mapped <- mapped[keep]

  # Seurat converts underscores to dashes in feature names, which can

  # create duplicates.  Pre-emptively deduplicate after the same conversion.
  rn <- gsub("_", "-", rownames(cts))
  dup <- duplicated(rn)
  if (any(dup)) {
    cts <- cts[!dup, ]
    rn  <- rn[!dup]
  }
  rownames(cts) <- rn
  seu <- CreateSeuratObject(counts = cts, min.cells = 0, min.features = 0)
  seu <- AddMetaData(seu, metadata = as.character(mapped), col.name = "cell_type")

  cat(sprintf("  %s: %d cells, %d genes, %d types\n",
              dataset_name, ncol(seu), nrow(seu), length(unique(seu@meta.data[["cell_type"]]))))
  cat("    ", paste(names(table(seu@meta.data[["cell_type"]])), collapse = ", "), "\n")
  seu
}


# ── 1. Load and prepare datasets ────────────────────────────────────────────

cat("=== [1/4] Loading pancreas datasets ===\n\n")

datasets <- list()

cat("Loading Baron...\n")
datasets[["Baron"]] <- sce_to_seurat(BaronPancreasData("human"), "Baron")

cat("Loading Muraro...\n")
datasets[["Muraro"]] <- sce_to_seurat(MuraroPancreasData(), "Muraro")

cat("Loading Segerstolpe...\n")
datasets[["Segerstolpe"]] <- sce_to_seurat(SegerstolpePancreasData(), "Segerstolpe")

cat("Loading Xin...\n")
datasets[["Xin"]] <- sce_to_seurat(XinPancreasData(), "Xin")

# Subsample for quick testing
if (QUICK_TEST) {
  cat("\n--- Subsampling to", QUICK_N, "cells per dataset ---\n")
  for (nm in names(datasets)) {
    seu_tmp <- datasets[[nm]]
    n <- ncol(seu_tmp)
    if (n > QUICK_N) {
      set.seed(42)
      # Stratified subsample
      ct_vec <- seu_tmp@meta.data[["cell_type"]]
      keep_idx <- unlist(lapply(split(seq_len(n), ct_vec), function(idx) {
        n_keep <- max(2, round(length(idx) / n * QUICK_N))
        sample(idx, min(n_keep, length(idx)))
      }))
      datasets[[nm]] <- seu_tmp[, keep_idx]
      cat(sprintf("  %s: %d → %d cells\n", nm, n, length(keep_idx)))
    }
  }
}

cat("\nDatasets loaded.\n\n")

n_folds <- if (QUICK_TEST) QUICK_FOLDS else 5L
n_epochs <- if (QUICK_TEST) QUICK_EPOCHS else 50L

# ── 2. Run benchmark on each dataset ────────────────────────────────────────

cat(sprintf("=== [2/4] Running %d-fold CV benchmark on each dataset ===\n\n", n_folds))

method_names <- c("DNN", "RF", "XGBoost", "SVM", "ElasticNet",
                   "LDA", "NaiveBayes", "KNN", "NNet",
                   "SingleR", "Seurat_LT", "CellTypist")
dataset_names <- names(datasets)
results_matrix <- matrix(NA_real_, nrow = length(dataset_names),
                          ncol = length(method_names),
                          dimnames = list(dataset_names, method_names))

# Check CellTypist availability
celltypist_ok <- tryCatch({
  reticulate::import("celltypist"); TRUE
}, error = function(e) FALSE)

for (ds_name in dataset_names) {

  cat(sprintf("\n════════════════════════════════════════════════════════\n"))
  cat(sprintf("  Dataset: %s (%d cells, %d types)\n",
              ds_name, ncol(datasets[[ds_name]]),
              length(unique(datasets[[ds_name]]@meta.data[["cell_type"]]))))
  cat(sprintf("════════════════════════════════════════════════════════\n\n"))

  seu <- datasets[[ds_name]]

  # ── Pre-process ──────────────────────────────────────────────────────────
  # Save true labels by barcode
  true_label_map <- setNames(
    as.character(seu@meta.data[["cell_type"]]),
    colnames(seu)
  )

  # Remove classes with too few cells for CV (need at least n_folds * 2)
  ct_counts <- table(true_label_map)
  min_needed <- n_folds * 2
  rare_types <- names(ct_counts[ct_counts < min_needed])
  if (length(rare_types) > 0) {
    keep_cells <- colnames(seu)[!true_label_map %in% rare_types]
    cat(sprintf("  Removing %d rare types with < %d cells: %s\n",
                length(rare_types), min_needed, paste(rare_types, collapse = ", ")))
    seu <- seu[, keep_cells]
    true_label_map <- true_label_map[keep_cells]
  }

  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  pca_mat <- Embeddings(seu, "pca")

  # Reconstruct labels in post-PCA order
  true_labels <- true_label_map[colnames(seu)]
  ct_levels <- sort(unique(true_labels))

  # ── DNN (CV) ────────────────────────────────────────────────────────────
  cat(sprintf("  [DNN] Training + predicting (%d-fold CV)...\n", n_folds))
  tryCatch({
    # Create stratified folds
    labels_f <- factor(true_labels)
    task_tmp <- TaskClassif$new(
      id = "tmp", backend = data.frame(y = labels_f, x = 1), target = "y")
    task_tmp$col_roles$stratum <- "y"
    cv5 <- rsmp("cv", folds = n_folds)
    cv5$instantiate(task_tmp)

    cell_names_all <- colnames(seu)
    dnn_preds <- setNames(character(length(cell_names_all)), cell_names_all)

    for (fold in seq_len(n_folds)) {
      cat(sprintf("    fold %d/%d...\n", fold, n_folds))
      tri <- cv5$train_set(fold)
      ti  <- cv5$test_set(fold)

      train_seu <- seu[, tri]
      test_seu  <- seu[, ti]

      model_path <- tempfile(fileext = ".pt")
      suppressMessages({
        DNN_train(
          input         = train_seu,
          label_column  = "cell_type",
          embedding_key = "umap",
          save_path     = model_path,
          n_top_genes   = 6000L,
          num_epochs    = n_epochs,
          batch_size    = 256L,
          device        = "auto"
        )
        pred <- DNN_predict(
          query          = test_seu,
          model_path     = model_path,
          true_label_col = "cell_type",
          device         = "auto"
        )
      })
      query_cells <- colnames(test_seu)
      dnn_preds[query_cells] <- as.character(pred$predictions$cell_type_pred)
      file.remove(model_path)
    }
    results_matrix[ds_name, "DNN"] <- compute_accuracy(true_labels, dnn_preds)
    cat(sprintf("    DNN accuracy: %.4f\n", results_matrix[ds_name, "DNN"]))
  }, error = function(e) {
    cat(sprintf("    DNN error: %s\n", e$message))
  })

  # ── ML models (mlr3verse, CV) ──────────────────────────────────────────
  cat(sprintf("  [ML] Running 8 learners (%d-fold CV)...\n", n_folds))
  tryCatch({
    task_df <- as.data.frame(pca_mat)
    task_df$cell_type <- factor(true_labels)
    task <- TaskClassif$new(id = "ct", backend = task_df, target = "cell_type")
    task$col_roles$stratum <- "cell_type"

    learners <- list(
      lrn("classif.ranger",      id = "RF",         num.trees = 500),
      lrn("classif.xgboost",     id = "XGBoost",    nrounds = 100, max_depth = 6, verbose = 0),
      lrn("classif.svm",         id = "SVM",        type = "C-classification", kernel = "radial"),
      lrn("classif.glmnet",      id = "ElasticNet",  alpha = 0.5),
      lrn("classif.lda",         id = "LDA"),
      lrn("classif.naive_bayes", id = "NaiveBayes"),
      lrn("classif.kknn",        id = "KNN",        k = 15),
      lrn("classif.nnet",        id = "NNet",       size = 20, MaxNWts = 50000, maxit = 200, trace = FALSE)
    )

    cv5_ml <- rsmp("cv", folds = n_folds)
    design <- benchmark_grid(task, learners, cv5_ml)
    suppressWarnings(suppressMessages(
      bmr <- benchmark(design, store_models = FALSE)
    ))
    agg <- bmr$aggregate(msr("classif.acc"))

    for (j in seq_len(nrow(agg))) {
      lid <- agg$learner_id[j]
      results_matrix[ds_name, lid] <- agg$classif.acc[j]
    }
    cat(sprintf("    ML done. Best: %s (%.4f)\n",
                agg$learner_id[which.max(agg$classif.acc)],
                max(agg$classif.acc)))

    # Extract fold assignments for single-cell methods
    rs_dt <- as.data.table(bmr$resamplings)
    rs <- rs_dt$resampling[[1]]
  }, error = function(e) {
    cat(sprintf("    ML error: %s\n", e$message))
    # Create fallback folds
    labels_f <- factor(true_labels)
    task_fb <- TaskClassif$new(
      id = "fb", backend = data.frame(y = labels_f, x = 1), target = "y")
    task_fb$col_roles$stratum <- "y"
    rs <<- rsmp("cv", folds = n_folds)
    rs$instantiate(task_fb)
  })

  # ── Single-cell methods (same folds) ────────────────────────────────────
  norm_data <- GetAssayData(seu, layer = "data")
  cell_names_all <- colnames(seu)
  n_cells <- ncol(seu)

  singler_preds    <- setNames(character(n_cells), cell_names_all)
  seurat_lt_preds  <- setNames(character(n_cells), cell_names_all)
  celltypist_preds <- setNames(character(n_cells), cell_names_all)
  sc_true_labels   <- setNames(character(n_cells), cell_names_all)

  # SingleR
  cat(sprintf("  [SingleR] %d-fold CV...\n", n_folds))
  tryCatch({
    for (fold in seq_len(n_folds)) {
      ti <- rs$test_set(fold); tri <- rs$train_set(fold)
      ref_sce <- SingleCellExperiment(assays = list(logcounts = norm_data[, tri]))
      ref_sce$label <- as.character(true_labels[tri])
      query_sce <- SingleCellExperiment(assays = list(logcounts = norm_data[, ti]))
      sr <- SingleR(test = query_sce, ref = ref_sce,
                    labels = ref_sce$label, de.method = "wilcox")
      query_cells <- cell_names_all[ti]
      singler_preds[query_cells] <- sr$labels
      sc_true_labels[query_cells] <- as.character(true_labels[ti])
    }
    results_matrix[ds_name, "SingleR"] <- compute_accuracy(sc_true_labels, singler_preds)
    cat(sprintf("    SingleR: %.4f\n", results_matrix[ds_name, "SingleR"]))
  }, error = function(e) {
    cat(sprintf("    SingleR error: %s\n", e$message))
  })

  # Seurat Label Transfer
  cat(sprintf("  [Seurat LT] %d-fold CV...\n", n_folds))
  tryCatch({
    for (fold in seq_len(n_folds)) {
      ti <- rs$test_set(fold); tri <- rs$train_set(fold)
      ref_s   <- seu[, tri]
      query_s <- seu[, ti]
      anchors <- FindTransferAnchors(
        reference = ref_s, query = query_s, dims = 1:30, verbose = FALSE)
      ref_labels <- as.character(ref_s@meta.data[["cell_type"]])
      transferred <- TransferData(
        anchorset = anchors, refdata = ref_labels, dims = 1:30, verbose = FALSE)
      query_cells <- colnames(query_s)
      seurat_lt_preds[query_cells] <- as.character(transferred$predicted.id)
      sc_true_labels[query_cells] <- as.character(query_s@meta.data[["cell_type"]])
    }
    results_matrix[ds_name, "Seurat_LT"] <- compute_accuracy(sc_true_labels, seurat_lt_preds)
    cat(sprintf("    Seurat LT: %.4f\n", results_matrix[ds_name, "Seurat_LT"]))
  }, error = function(e) {
    cat(sprintf("    Seurat LT error: %s\n", e$message))
  })

  # CellTypist
  if (celltypist_ok) {
    cat(sprintf("  [CellTypist] %d-fold CV...\n", n_folds))
    tryCatch({
      ct <- reticulate::import("celltypist")
      sc_py <- reticulate::import("scanpy")
      pd <- reticulate::import("pandas")

      for (fold in seq_len(n_folds)) {
        ti <- rs$test_set(fold); tri <- rs$train_set(fold)
        query_cells <- cell_names_all[ti]

        ref_mat <- as.matrix(t(norm_data[, tri]))
        ref_adata <- sc_py$AnnData(
          X = ref_mat,
          obs = pd$DataFrame(list(cell_type = as.character(true_labels[tri])))
        )
        ref_adata$var_names <- rownames(norm_data)

        query_mat <- as.matrix(t(norm_data[, ti]))
        query_adata <- sc_py$AnnData(X = query_mat)
        query_adata$var_names <- rownames(norm_data)

        model <- ct$train(ref_adata, labels = "cell_type",
                          use_SGD = FALSE, n_jobs = 4L, max_iter = 200L)
        preds <- ct$annotate(query_adata, model = model, majority_voting = FALSE)
        pred_df <- reticulate::py_to_r(preds$predicted_labels)
        celltypist_preds[query_cells] <- as.character(pred_df[["predicted_labels"]])
      }
      results_matrix[ds_name, "CellTypist"] <- compute_accuracy(
        sc_true_labels, celltypist_preds)
      cat(sprintf("    CellTypist: %.4f\n", results_matrix[ds_name, "CellTypist"]))
    }, error = function(e) {
      cat(sprintf("    CellTypist error: %s\n", e$message))
    })
  }

  cat(sprintf("\n  %s complete.\n", ds_name))
}


# ── 3. Results ──────────────────────────────────────────────────────────────

cat("\n\n=== [3/4] Results ===\n\n")
cat("================================================================\n")
cat("    Cross-Dataset Pancreas Classification Benchmark\n")
cat(sprintf("    (%d-fold stratified CV, Accuracy)\n", n_folds))
cat("================================================================\n\n")

# Print table
cat(sprintf("%-15s", "Dataset"))
for (m in method_names) cat(sprintf(" %10s", m))
cat("\n")
cat(strrep("-", 15 + 11 * length(method_names)), "\n")
for (ds in dataset_names) {
  cat(sprintf("%-15s", ds))
  for (m in method_names) {
    val <- results_matrix[ds, m]
    cat(sprintf(" %10s", ifelse(is.na(val), "  NA", sprintf("%.4f", val))))
  }
  cat("\n")
}
cat(strrep("-", 15 + 11 * length(method_names)), "\n")
cat(sprintf("%-15s", "Mean"))
for (m in method_names) {
  cat(sprintf(" %10.4f", mean(results_matrix[, m], na.rm = TRUE)))
}
cat("\n\n")

# Ranking by mean accuracy
mean_acc <- colMeans(results_matrix, na.rm = TRUE)
ranked <- sort(mean_acc, decreasing = TRUE)
cat("--- Ranking (by mean accuracy) ---\n")
for (i in seq_along(ranked)) {
  cat(sprintf("  %2d. %-12s %.4f\n", i, names(ranked)[i], ranked[i]))
}


# ── 4. Heatmap visualization ────────────────────────────────────────────────

cat("\n\n=== [4/4] Generating heatmap ===\n")

# Convert to long format
heat_df <- as.data.frame(as.table(results_matrix))
colnames(heat_df) <- c("Dataset", "Method", "Accuracy")
heat_df$Accuracy <- as.numeric(heat_df$Accuracy)

# Order methods by mean accuracy (best on right)
method_order <- names(sort(mean_acc))
heat_df$Method <- factor(heat_df$Method, levels = method_order)
heat_df$Dataset <- factor(heat_df$Dataset, levels = rev(dataset_names))

# Text color
heat_df$text_col <- ifelse(
  is.na(heat_df$Accuracy), "black",
  ifelse(heat_df$Accuracy > 0.85, "white", "black")
)

# Dynamic color scale
min_val <- floor(min(results_matrix, na.rm = TRUE) * 20) / 20
max_val <- 1.0

p_heat <- ggplot(heat_df, aes(x = Method, y = Dataset, fill = Accuracy)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(
    aes(label = ifelse(is.na(Accuracy), "NA", sprintf("%.3f", Accuracy))),
    color = heat_df$text_col, size = 3.8, fontface = "bold"
  ) +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
               "#FDDBC7", "#EF8A62", "#B2182B"),
    limits = c(min_val, max_val),
    na.value = "grey80",
    name = "Accuracy"
  ) +
  labs(
    title = "Cross-Dataset Pancreas Classification Benchmark",
    subtitle = "5-fold stratified CV on 4 independent datasets",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.text.y   = element_text(face = "bold", size = 12),
    panel.grid     = element_blank(),
    legend.position = "right"
  )
print(p_heat)

# Save results
out_dir <- file.path(dirname(dirname(dirname(getwd()))), "output_01")
if (!dir.exists(out_dir)) out_dir <- tempdir()

results_path <- file.path(out_dir, "pancreas_benchmark_results.qs")
qsave(list(matrix = results_matrix, heatmap = p_heat), results_path)
cat(sprintf("\nResults saved to: %s\n", results_path))

heatmap_path <- file.path(out_dir, "pancreas_benchmark_heatmap.pdf")
ggsave(heatmap_path, p_heat, width = 14, height = 5)
cat(sprintf("Heatmap saved to: %s\n", heatmap_path))

cat("\n=== Cross-Dataset Benchmark Complete ===\n")
