# ============================================================================
# Demo: EvaluateEmbedding — 评估 DNN 512维嵌入的信息量与 PCA 一致性
# ============================================================================
#
# 使用 toy_test 数据集演示 EvaluateEmbedding() 和 PlotEmbeddingEval()。
#
# Requirements:
#   - scMMR installed (with model.pt and toy_test.qs in inst/extdata/)
#   - Python environment configured (scanpy-env)
#   - Packages: qs, Seurat, FNN, cluster
#
# Usage:
#   source("demo_evaluate_embedding.R")
# ============================================================================

library(scMMR)
library(Seurat)
library(qs)
library(ggplot2)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

cat("=== Demo: EvaluateEmbedding ===\n\n")


# ── 1. 加载数据与模型 ─────────────────────────────────────────────────────────

test_path  <- system.file("extdata", "toy_test.qs", package = "scMMR")
model_path <- system.file("extdata", "model.pt",    package = "scMMR")

toy_test <- qread(test_path)
cat("Query data:", ncol(toy_test), "cells x", nrow(toy_test), "genes\n")
cat("Cell types:", paste(unique(toy_test$cell_type), collapse = ", "), "\n")
cat("Groups:", paste(unique(toy_test$group), collapse = ", "), "\n\n")


# ── 2. DNN 预测 (return_embedding = TRUE) ────────────────────────────────────

pred <- DNN_predict(
  query            = toy_test,
  model_path       = model_path,
  true_label_col   = "cell_type",
  return_embedding = TRUE,
  device           = "cpu"
)

cat("Shared embedding:", nrow(pred$shared_embedding), "cells x",
    ncol(pred$shared_embedding), "dims\n\n")

# 合并预测结果到 Seurat 对象
q1 <- AddMetaData(toy_test, pred$predictions)


# ── 3. 确保有 PCA ────────────────────────────────────────────────────────────

# 如果 toy_test 还没有 PCA，先计算
if (is.null(Reductions(q1, "pca"))) {
  cat("Computing PCA on Seurat object ...\n")
  q1 <- NormalizeData(q1, verbose = FALSE)
  q1 <- FindVariableFeatures(q1, nfeatures = 2000, verbose = FALSE)
  q1 <- ScaleData(q1, verbose = FALSE)
  q1 <- RunPCA(q1, npcs = 50, verbose = FALSE)
}
cat("PCA dims:", ncol(Embeddings(q1, "pca")), "\n\n")


# ── 4. 评估嵌入质量 (传 Seurat 对象) ─────────────────────────────────────────

cat("=== Test 1: EvaluateEmbedding with Seurat object ===\n\n")

eval_res <- EvaluateEmbedding(
  embedding     = pred$shared_embedding,
  seurat_obj    = q1,
  cell_type_col = "cell_type",
  n_pcs         = 50L,
  k_values      = c(10L, 20L, 50L, 100L),
  n_pairs       = 5000L,
  seed          = 42L,
  verbose       = TRUE
)

cat("\n--- Results ---\n")
cat(eval_res$summary, "\n\n")

# 内在维度
cat("Effective dim (knee):", eval_res$pca_of_embedding$effective_dim, "\n")
cat("Cumulative variance (first 10 PCs):\n")
print(round(eval_res$pca_of_embedding$cum_var[1:10], 3))

# KNN overlap
cat("\nKNN overlap:\n")
print(eval_res$consistency$knn_overlap)

# 距离相关
cat("\nDistance correlation:", eval_res$consistency$dist_cor$correlation, "\n")
cat("RV coefficient:", eval_res$consistency$rv_coef, "\n")

# Silhouette
if (!is.null(eval_res$consistency$silhouette)) {
  cat("\nSilhouette:\n")
  cat("  DNN mean:", eval_res$consistency$silhouette$emb_mean, "\n")
  cat("  PCA mean:", eval_res$consistency$silhouette$ref_mean, "\n")
  print(eval_res$consistency$silhouette$per_type)
}


# ── 5. 可视化 ────────────────────────────────────────────────────────────────

cat("\n=== Test 2: PlotEmbeddingEval ===\n")

# 5a. 获取所有图
plots <- PlotEmbeddingEval(eval_res, which = "all")
cat("Available plots:", paste(names(plots), collapse = ", "), "\n")

# 5b. 逐个绘制
print(PlotEmbeddingEval(eval_res, which = "elbow"))
print(PlotEmbeddingEval(eval_res, which = "cumvar"))
print(PlotEmbeddingEval(eval_res, which = "knn"))
print(PlotEmbeddingEval(eval_res, which = "distcor"))

if (!is.null(eval_res$consistency$silhouette)) {
  print(PlotEmbeddingEval(eval_res, which = "silhouette"))
}

# 5c. 使用 patchwork 组合 (可选)
# library(patchwork)
# (plots$elbow + plots$cumvar) / (plots$knn + plots$distcor)
# ggsave("embedding_eval.pdf", width = 14, height = 10)


# ── 6. 仅评估内在维度 (无 PCA reference) ─────────────────────────────────────

cat("\n=== Test 3: Intrinsic dimensionality only ===\n\n")

eval_intrinsic <- EvaluateEmbedding(
  embedding = pred$shared_embedding,
  n_pcs     = 50L,
  verbose   = TRUE
)

cat("Summary:\n", eval_intrinsic$summary, "\n")
cat("Consistency is NULL:", is.null(eval_intrinsic$consistency), "\n")

plots2 <- PlotEmbeddingEval(eval_intrinsic, which = "all")
cat("Available plots (no ref):", paste(names(plots2), collapse = ", "), "\n")


# ── 7. 用 pca_ref 直接传矩阵 + 标签向量 ─────────────────────────────────────

cat("\n=== Test 4: Using pca_ref + label vector ===\n\n")

pca_mat <- Embeddings(q1, "pca")
labels  <- setNames(as.character(q1$cell_type), colnames(q1))

eval_ref <- EvaluateEmbedding(
  embedding     = pred$shared_embedding,
  pca_ref       = pca_mat,
  cell_type_col = labels,
  k_values      = c(20L, 50L),
  verbose       = TRUE
)
cat("Summary:\n", eval_ref$summary, "\n")


# ── 8. 下游任务对比: DNN 嵌入 vs PCA 空间的扰动检测能力 ─────────────────────

cat("\n=== Test 5: RankPerturbation — DNN vs PCA ===\n\n")

# 取两个条件进行对比
conditions <- c("PT", "PH")

# 8a. 在 DNN 512维嵌入空间上做扰动排名
cat("--- RankPerturbation on DNN embedding (512d) ---\n")
res_dnn <- RankPerturbation(
  embedding     = pred$shared_embedding,
  cell_meta     = q1@meta.data,
  cell_type_col = "cell_type",
  condition_col = "group",
  conditions    = conditions,
  method        = "wasserstein",
  n_permutations = 500L,
  seed          = 42L,
  verbose       = TRUE
)

# 8b. 在 PCA 50维空间上做扰动排名
cat("\n--- RankPerturbation on PCA embedding (50d) ---\n")
pca_emb <- Embeddings(q1, "pca")
res_pca <- RankPerturbation(
  embedding     = pca_emb,
  cell_meta     = q1@meta.data,
  cell_type_col = "cell_type",
  condition_col = "group",
  conditions    = conditions,
  method        = "wasserstein",
  n_permutations = 500L,
  n_pcs         = NULL,      # 不再降维，直接用 50 PCs
  seed          = 42L,
  verbose       = TRUE
)

# 8c. 对比结果
cat("\n========================================\n")
cat("  DNN vs PCA Perturbation Comparison\n")
cat(sprintf("  Conditions: %s vs %s\n", conditions[1], conditions[2]))
cat("========================================\n\n")

# 显著性对比
n_sig_dnn <- sum(res_dnn$results$p_adj < 0.05, na.rm = TRUE)
n_sig_pca <- sum(res_pca$results$p_adj < 0.05, na.rm = TRUE)
cat(sprintf("Significant cell types (p_adj < 0.05):\n"))
cat(sprintf("  DNN: %d / %d\n", n_sig_dnn, nrow(res_dnn$results)))
cat(sprintf("  PCA: %d / %d\n", n_sig_pca, nrow(res_pca$results)))

# 效应量对比
cat(sprintf("\nMean perturbation score:\n"))
cat(sprintf("  DNN: %.4f\n", mean(res_dnn$results$score)))
cat(sprintf("  PCA: %.4f\n", mean(res_pca$results$score)))

# 排名一致性
common_ct <- intersect(res_dnn$results$cell_type, res_pca$results$cell_type)
if (length(common_ct) >= 3) {
  rank_dnn <- setNames(res_dnn$results$rank, res_dnn$results$cell_type)
  rank_pca <- setNames(res_pca$results$rank, res_pca$results$cell_type)
  rho <- cor(rank_dnn[common_ct], rank_pca[common_ct], method = "spearman")
  cat(sprintf("\nRank correlation (Spearman): %.3f\n", rho))
  cat("  (> 0.8 = highly consistent, < 0.5 = divergent)\n")
}

# 并排展示排名
cat("\n--- Side-by-side ranking ---\n")
merged <- merge(
  res_dnn$results[, c("cell_type", "score", "p_adj", "rank")],
  res_pca$results[, c("cell_type", "score", "p_adj", "rank")],
  by = "cell_type", suffixes = c("_DNN", "_PCA")
)
merged <- merged[order(merged$rank_DNN), ]
cat(sprintf("%-25s %8s %8s %8s %8s\n",
            "Cell Type", "DNN_rank", "PCA_rank", "DNN_padj", "PCA_padj"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(merged))) {
  cat(sprintf("%-25s %8d %8d %8.3f %8.3f\n",
              merged$cell_type[i],
              merged$rank_DNN[i], merged$rank_PCA[i],
              merged$p_adj_DNN[i], merged$p_adj_PCA[i]))
}

# 8d. 可视化对比
p_dnn <- PlotPerturbation(res_dnn) +
  ggplot2::labs(title = "DNN Embedding (512d)")
p_pca <- PlotPerturbation(res_pca) +
  ggplot2::labs(title = "PCA Embedding (50d)")
print(p_dnn)
print(p_pca)

# 组合展示 (patchwork)
# library(patchwork)
# p_dnn + p_pca + plot_layout(ncol = 2)
# ggsave("perturbation_DNN_vs_PCA.pdf", width = 16, height = 6)

cat("\n=== Interpretation ===\n")
cat("If DNN detects MORE significant cell types with LARGER effect sizes,\n")
cat("it suggests the DNN embedding provides a more sensitive space for\n")
cat("perturbation analysis — beyond just better cell type separation.\n")


cat("\n=== Demo Complete ===\n")
