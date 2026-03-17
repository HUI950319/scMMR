library(scMMR)
library(qs)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# Load data
test_path  <- system.file("extdata", "toy_test.qs", package = "scMMR")
model_path <- system.file("extdata", "model.pt", package = "scMMR")
gmt_path   <- system.file("extdata/gmt", "reactome.gmt", package = "scMMR")

toy_test <- qread(test_path)
cat("Data:", ncol(toy_test), "cells\n\n")

# ── Test 1: CPU ──────────────────────────────────────────────────────────────
cat("========== CPU Test ==========\n")
t_cpu <- system.time({
  pred_cpu <- DNN_predict(
    query          = toy_test,
    model_path     = model_path,
    true_label_col = "cell_type",
    explain        = TRUE,
    pathway_gmt    = gmt_path,
    pathway_n_steps = 20L,
    device         = "cpu"
  )
})
cat(sprintf("\nCPU total time: %.1f seconds\n", t_cpu["elapsed"]))
cat(sprintf("CPU pathway_scores: %d x %d\n\n",
            nrow(pred_cpu$pathway_scores), ncol(pred_cpu$pathway_scores)))

# ── Test 2: GPU ──────────────────────────────────────────────────────────────
cat("========== GPU Test ==========\n")
t_gpu <- system.time({
  pred_gpu <- DNN_predict(
    query          = toy_test,
    model_path     = model_path,
    true_label_col = "cell_type",
    explain        = TRUE,
    pathway_gmt    = gmt_path,
    pathway_n_steps = 20L,
    device         = "auto"
  )
})
cat(sprintf("\nGPU total time: %.1f seconds\n", t_gpu["elapsed"]))
cat(sprintf("GPU pathway_scores: %d x %d\n\n",
            nrow(pred_gpu$pathway_scores), ncol(pred_gpu$pathway_scores)))

# ── Summary ──────────────────────────────────────────────────────────────────
cat("========== Speed Comparison ==========\n")
cat(sprintf("CPU: %.1f seconds\n", t_cpu["elapsed"]))
cat(sprintf("GPU: %.1f seconds\n", t_gpu["elapsed"]))
cat(sprintf("Speedup: %.1fx\n", t_cpu["elapsed"] / t_gpu["elapsed"]))

# Verify results match
cat("\nResults verification:\n")
cat(sprintf("  Accuracy match: %s\n",
            identical(pred_cpu$predictions$cell_type_pred,
                      pred_gpu$predictions$cell_type_pred)))
cat(sprintf("  Pathway dim match: %s\n",
            identical(dim(pred_cpu$pathway_scores),
                      dim(pred_gpu$pathway_scores))))

# Top 5 pathways comparison
cat("\nTop 5 pathways (GPU):\n")
ms <- sort(colMeans(pred_gpu$pathway_scores), decreasing = TRUE)
for (i in 1:5) {
  cat(sprintf("  %d. %-45s %.6f\n", i, names(ms)[i], ms[i]))
}
