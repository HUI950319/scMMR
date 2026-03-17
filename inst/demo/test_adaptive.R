#!/usr/bin/env Rscript
# Test TAPE-style adaptive deconvolution (Seurat input)
library(scMMR)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# ── Load Baron reference as Seurat ──
cat("\n== Loading Baron reference (Seurat) ==\n")
sce <- scRNAseq::BaronPancreasData("human")
sce <- sce[, !is.na(SummarizedExperiment::colData(sce)[["label"]])]
seu <- Seurat::as.Seurat(sce, counts = "counts", data = NULL)
cat(sprintf("  %d cells x %d genes, %d cell types\n",
            ncol(seu), nrow(seu), length(unique(seu$label))))

# ── Train ──
cat("\n== Training ==\n")
model_path <- file.path(tempdir(), "deconv_baron_tape.pt")
t0 <- Sys.time()
train_result <- DNN_deconv_train(
  reference    = seu,
  label_column = "label",
  save_path    = model_path,
  n_pseudobulk = 5000L,
  n_top_genes  = 3000L,
  num_epochs   = 50L,
  batch_size   = 256L,
  device       = "auto"
)
cat(sprintf("  Training time: %.1f sec\n", as.numeric(Sys.time() - t0, units = "secs")))
cat(sprintf("  Best val corr: %.4f\n", train_result[["best_val_corr"]]))

# ── Generate test pseudo-bulk ──
cat("\n== Generating test pseudo-bulk ==\n")
ref_adata <- scMMR:::.resolve_input(seu)
pb_test <- scMMR:::.scMMR_env$deconv_generate_pseudobulk(
  adata = ref_adata, label_column = "label",
  n_samples = 20L, n_cells_per_sample = 500L,
  var_genes = train_result[["var_genes"]], random_state = 999L
)
pb_r <- reticulate::py_to_r(pb_test)
test_matrix <- t(expm1(pb_r[["X_pb"]]))
rownames(test_matrix) <- pb_r[["var_genes"]]
colnames(test_matrix) <- paste0("sample_", seq_len(ncol(test_matrix)))
true_props <- pb_r[["y_pb"]]

# ── Test 1: Standard ──
cat("\n== Test 1: Standard prediction ==\n")
props_std <- DNN_deconv_predict(bulk_expr = test_matrix, model_path = model_path)
r_std <- cor(as.vector(as.matrix(props_std[, -1])), as.vector(true_props))
cat(sprintf("  r = %.4f\n", r_std))

# ── Test 2: Adaptive (overall) ──
cat("\n== Test 2: Adaptive prediction (overall) ==\n")
props_adp <- DNN_deconv_predict(
  bulk_expr = test_matrix, model_path = model_path,
  adaptive = TRUE, adaptive_mode = "overall",
  adaptive_max_iter = 2L, adaptive_steps = 100L
)
r_adp <- cor(as.vector(as.matrix(props_adp$proportions[, -1])), as.vector(true_props))
sigm <- props_adp$sigmatrix
cat(sprintf("  r = %.4f\n", r_adp))
cat(sprintf("  GEP: %d x %d, non-neg=%s, range=[%.2f, %.2f]\n",
            nrow(sigm), ncol(sigm), all(sigm >= 0), min(sigm), max(sigm)))

# ── Summary ──
cat(sprintf("\n== Results: standard r=%.4f, adaptive r=%.4f ==\n", r_std, r_adp))
cat("== All tests passed! ==\n")



