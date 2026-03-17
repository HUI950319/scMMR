# ============================================================================
# scMMR Demo: Prediction & Visualization
# ============================================================================
#
# This script demonstrates the full scMMR workflow:
#   1. Load test data and trained model
#   2. Run DNN_predict() for cell type prediction
#   3. Merge predictions into Seurat object
#   4. Visualize with PlotAlluvia, PlotGroupPreference, PlotMAP
#
# Requirements:
#   - scMMR installed
#   - Python environment configured (use_scMMR_python)
#   - qs package for loading test data
# ============================================================================

library(scMMR)
library(Seurat)
library(tidyverse)
library(qs)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# ── 0. Setup Python environment ─────────────────────────────────────────────
# Uncomment and modify to match your setup:
# use_scMMR_python(condaenv = "scanpy-env")

# ── 1. Load test data, model & reference UMAP ────────────────────────────────
test_path  <- system.file("extdata", "toy_test.qs",  package = "scMMR")
model_path <- system.file("extdata", "model.pt",     package = "scMMR")
ref_path   <- system.file("extdata", "ref_umap.qs",  package = "scMMR")

toy_test <- qread(test_path)
ref_data <- qread(ref_path)   # data.frame: umap_1, umap_2, cell_type

cat("Query data:", ncol(toy_test), "cells x", nrow(toy_test), "genes\n")
cat("Reference UMAP:", nrow(ref_data), "cells\n")
# cat("Groups:", paste(levels(factor(toy_test$group)), collapse = ", "), "\n")

# ── 2. Run prediction (with gene importance + pathway scoring) ───────────────
gmt_path <- system.file("extdata/gmt", "reactome.gmt", package = "scMMR")

pred <- DNN_predict(
  query          = toy_test,
  model_path     = model_path,
  true_label_col = "celltype",
  explain        = TRUE,
  top_k_class = 50L,
  pathway_gmt    = gmt_path,       # optional: GMT for pathway scoring
  device         = "gpu")

head(pred$predictions)
# columns: cell_id, cell_type_pred, confidence, is_ood, umap_1_pred, umap_2_pred
# cell_type_pred = always best-guess label; use is_ood / confidence to filter

# Pathway scores matrix: n_cells × n_pathways
cat("Pathway scores:", nrow(pred$pathway_scores), "cells ×",
    ncol(pred$pathway_scores), "pathways\n")
cat("Top 5 pathways:\n")
cat(head(colnames(pred$pathway_scores), 5), sep = "\n")

# ── 3. Merge predictions into Seurat ────────────────────────────────────────
q1 <- AddMetaData(toy_test, pred$predictions)

# Check merged columns
cat("\nNew meta columns:\n")
cat(setdiff(colnames(q1@meta.data), colnames(toy_test@meta.data)), sep = ", ")
cat("\n")


colors = ggsci::pal_d3("category20")(length(unique(q1$cell_type_pred)))
names(colors) = unique(q1$cell_type_pred)
# ── 4. PlotAlluvia: Cell composition across groups ──────────────────────────
# Shows predicted cell type proportions changing across PH / PT / SH
p1 <- PlotAlluvia(
  cellmeta = q1@meta.data,
  by       = "group",
  fill     = "cell_type_pred",
  colors   = colors)

print(p1)
# ggsave("alluvia_plot.pdf", p1, width = 8, height = 6)

# ── 5. PlotGroupPreference: O/E ratio heatmap ──────────────────────────────
# Highlights which cell types are enriched (O/E > 1) in each group
p2 <- PlotGroupPreference(
  cellmeta      = q1@meta.data,
  group.by      = "group",
  preference.on = "cell_type_pred",
  palette = "Blues",
  column_names_rot = 45
)
print(p2)
# pdf("group_preference.pdf", width = 6, height = 8); print(p2); dev.off()

# ── 6. PlotMAP: UMAP projection ────────────────────────────────────────────
# Overlay predicted query coordinates onto reference UMAP (170k cells)

# 6a. Basic projection (all query cells together)
p3<- PlotMAP(
  ref        = ref_data,
  query_meta = q1@meta.data,
  color_by   = "celltype")

print(p3)

# 6b. Faceted by group
p4 <- PlotMAP(
  ref        = ref_data,
  query_meta = q1@meta.data,
  color_by   = "celltype",
  facet_by   = "group")

print(p4)
# ggsave("umap_projection.pdf", p4, width = 14, height = 5)









# ── 7. PlotImportance: Gene importance ────────────────────────────────────────
# 7a. Global importance (top 15 genes)
p5 <- PlotImportance(pred$imp_global, top_k = 15)
print(p5)

# 7b. Per-class importance (top 10 genes per cell type)
p6 <- PlotImportance(pred$imp_per_class, top_k = 10, ncol = 4)
print(p6)
# ggsave("importance_per_class.pdf", p6, width = 16, height = 12)

# 7c. Bar style with different palette
p7 <- PlotImportance(pred$imp_global, display = "bar", palette = "Blues 3")
print(p7)

# ── 7d. Pathway importance: global (top 20 pathways) ─────────────────────
p8 <- PlotImportance(pred$pathway_scores, top_k = 20)
print(p8)
# ggsave("pathway_global.pdf", p8, width = 10, height = 8)

# ── 7e. Pathway importance: per cell type ────────────────────────────────
p9 <- PlotImportance(
  pred$pathway_scores,
  top_k    = 10,
  group_by = q1$cell_type_pred,
  ncol     = 3
)
print(p9)
# ggsave("pathway_per_celltype.pdf", p9, width = 16, height = 14)

# ── 7f. Pathway importance: per sample group ─────────────────────────────
p10 <- PlotImportance(
  pred$pathway_scores,
  top_k    = 10,
  group_by = q1$group,
  ncol     = 3,
  palette  = "Purples 3"
)
print(p10)
# ggsave("pathway_per_group.pdf", p10, width = 14, height = 6)

# ── 8. Comparison: true vs predicted labels ─────────────────────────────────
# Confusion table
cat("\n== Prediction accuracy ==\n")
acc <- mean(q1$cell_type_pred == q1$cell_type)
cat(sprintf("Overall accuracy: %.2f%%\n", acc * 100))

cat("\n== Confusion (first 10 types) ==\n")
ct <- table(True = q1$cell_type, Predicted = q1$cell_type_pred)
print(ct[1:min(10, nrow(ct)), 1:min(10, ncol(ct))])

# Low-confidence summary (is_ood flag)
cat(sprintf("\nLow-confidence cells: %d / %d (%.1f%%) [is_ood flag]\n",
            sum(q1$is_ood), ncol(q1), 100 * mean(q1$is_ood)))
cat(sprintf("Confidence: mean=%.3f, median=%.3f\n",
            mean(q1$confidence), median(q1$confidence)))

# ── 9. Pathway scoring analysis ──────────────────────────────────────────────
# pred$pathway_scores is a matrix [n_cells × n_pathways]
# Each value = mean absolute IG attribution aggregated to that pathway
# This captures the model's discriminative use of each pathway (not expression)

pw_scores <- pred$pathway_scores

# 9a. Top pathways by mean score across all cells
mean_pw <- colMeans(pw_scores)
top_pw  <- sort(mean_pw, decreasing = TRUE)[1:min(20, length(mean_pw))]
cat("\n== Top 20 pathways (mean score across all cells) ==\n")
for (i in seq_along(top_pw)) {
  cat(sprintf("  %2d. %-50s %.4f\n", i, names(top_pw)[i], top_pw[i]))
}

# 9b. Per-group pathway comparison (e.g., detect group-specific pathway activity)
# Merge pathway scores into metadata for group comparisons
pw_df <- as.data.frame(pw_scores)
pw_df$group     <- q1$group[match(rownames(pw_scores), colnames(q1))]
pw_df$cell_type <- q1$cell_type_pred[match(rownames(pw_scores), colnames(q1))]

# # Example: Compare a specific pathway across groups (within one cell type)
# if ("Macrophages" %in% unique(pw_df$cell_type) && ncol(pw_scores) > 0) {
#   target_ct <- "Macrophages"
#   target_pw <- names(top_pw)[1]  # top pathway
# 
#   sub_df <- pw_df[pw_df$cell_type == target_ct, ]
#   cat(sprintf("\n== Pathway '%s' in %s across groups ==\n", target_pw, target_ct))
#   for (g in sort(unique(sub_df$group))) {
#     vals <- sub_df[sub_df$group == g, target_pw]
#     if (length(vals) > 0) {
#       cat(sprintf("  %s: mean=%.4f (n=%d)\n", g, mean(vals), length(vals)))
#     }
#   }
# }

# 9c. Use different GMT file (e.g., GO biological process)
# gmt_go <- system.file("extdata/gmt", "GO_bp.gmt", package = "scMMR")
# pred_go <- DNN_predict(
#   query = toy_test, model_path = model_path,
#   explain = TRUE, pathway_gmt = gmt_go, device = "cpu"
# )
# dim(pred_go$pathway_scores)

# Available GMT files:
# Human: reactome.gmt, GO_bp.gmt, TF.gmt, immune.gmt
# Mouse: m_reactome.gmt, m_GO_bp.gmt, m_TF.gmt

# ── 10. RankPerturbation: Which cell types shift most? ─────────────────────
# Requires return_embedding = TRUE in DNN_predict()
pred_emb <- DNN_predict(
  query          = toy_test,
  model_path     = model_path,
  true_label_col = "cell_type",
  return_embedding = TRUE,
  device         = "gpu"
)
q2 <- AddMetaData(toy_test, pred_emb$predictions)

# 10a. Wasserstein distance (default): PH vs SH
perturb_ws <- RankPerturbation(
  embedding     = pred_emb$shared_embedding,
  cell_meta     = q2@meta.data,
  cell_type_col = "cell_type_pred",
  condition_col = "group",
  conditions    = c("PT", "PH"),
  method        = "wasserstein"
)
cat("\n== Perturbation ranking (Wasserstein, PH vs SH) ==\n")
print(perturb_ws$results)


p11 <- PlotPerturbation(perturb_ws)
print(p11)
# ggsave("perturbation_wasserstein.pdf", p11, width = 8, height = 6)

# 10b. MMD distance: PH vs PT
perturb_mmd <- RankPerturbation(
  embedding     = pred_emb$shared_embedding,
  cell_meta     = q2@meta.data,
  conditions    = c("PH", "PT"),
  method        = "mmd"
)
p12 <- PlotPerturbation(perturb_mmd, display = "bar", palette = "Blues 3")
print(p12)

# 10c. Energy distance (fast, parameter-free)
perturb_energy <- RankPerturbation(
  embedding     = pred_emb$shared_embedding,
  cell_meta     = q2@meta.data,
  conditions    = c("PT", "PH"),
  method        = "energy"
)
p12b <- PlotPerturbation(perturb_energy, palette = "Greens 3")
print(p12b)
# ggsave("perturbation_energy.pdf", p12b, width = 8, height = 6)

# 10d. Classifier AUC (Augur-like, now uses LDA by default)
perturb_auc <- RankPerturbation(
  embedding     = pred_emb$shared_embedding,
  cell_meta     = q2@meta.data,
  conditions    = c("PH", "SH"),
  method        = "classifier",
  n_permutations = 200   # faster for demo
)
p13 <- PlotPerturbation(perturb_auc)
print(p13)

# 10e. With balanced downsampling and parallel computation
# perturb_bal <- RankPerturbation(
#   embedding     = pred_emb$shared_embedding,
#   cell_meta     = q2@meta.data,
#   conditions    = c("PH", "SH"),
#   method        = "wasserstein",
#   balance_cells = TRUE,
#   n_cores       = 4
# )
# PlotPerturbation(perturb_bal)

# ── 11. RankPercent: Differential abundance testing ────────────────────────
# miloR-style KNN neighborhood DA testing on DNN embedding

# 11a. Fisher's exact test
da <- RankPercent(
  embedding     = pred_emb$shared_embedding,
  cell_meta     = q2@meta.data,
  cell_type_col = "cell_type_pred",
  condition_col = "group",
  conditions    = c( "PT", "PH"),
  k             = 30,
  prop_sample   = 0.2
  )


cat("\n== Differential abundance summary ==\n")
print(da$cell_type_summary)

# Add per-cell DA score to Seurat
q2 <- AddMetaData(q2, da$cell_da_scores, col.name = "da_score")

# # 11b. Beeswarm plot (miloR-style)
# p14 <- PlotPercent(da,fdr_threshold = 0)
# print(p14)
# # ggsave("da_beeswarm.pdf", p14, width = 8, height = 6)

# 11c. Stricter FDR threshold
p15 <- PlotPercent(da, fdr_threshold = 0.1, show_boxplot = T)
print(p15)



