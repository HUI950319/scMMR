#!/usr/bin/env Rscript
# Benchmark DNN_deconv_train() timing
library(scMMR)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")

# ── Stage 1: Load reference ──
cat("\n== Stage 1: Loading Baron reference ==\n")
t1 <- Sys.time()
sce <- scRNAseq::BaronPancreasData("human")
sce <- sce[, !is.na(SummarizedExperiment::colData(sce)[["label"]])]
ref_path <- tempfile(fileext = ".h5ad")
sceasy::convertFormat(sce, from = "sce", to = "anndata", outFile = ref_path)
t2 <- Sys.time()
cat(sprintf("  Stage 1 (data load + convert): %.1f sec\n", as.numeric(difftime(t2, t1, units = "secs"))))

# ── Stage 2: Full training pipeline ──
cat("\n== Stage 2: DNN_deconv_train() ==\n")
t3 <- Sys.time()
train_result <- DNN_deconv_train(
  reference    = ref_path,
  label_column = "label",
  save_path    = file.path(tempdir(), "deconv_baron_tape.pt"),
  n_pseudobulk = 5000L,
  n_top_genes  = 3000L,
  num_epochs   = 50L,
  batch_size   = 256L,
  device       = "auto"
)
t4 <- Sys.time()
cat(sprintf("  Stage 2 (DNN_deconv_train total): %.1f sec\n", as.numeric(difftime(t4, t3, units = "secs"))))

# ── Summary ──
cat("\n== Timing Summary ==\n")
cat(sprintf("  Data load + convert:  %.1f sec\n", as.numeric(difftime(t2, t1, units = "secs"))))
cat(sprintf("  DNN_deconv_train():   %.1f sec\n", as.numeric(difftime(t4, t3, units = "secs"))))
cat(sprintf("  Total:                %.1f sec\n", as.numeric(difftime(t4, t1, units = "secs"))))
cat(sprintf("  Best val corr: %.4f\n", train_result[["best_val_corr"]]))
cat("\n== Done ==\n")
