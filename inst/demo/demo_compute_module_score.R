# ============================================================================
# scMMR Demo: ComputeModuleScore
# ============================================================================
#
# This script demonstrates the ComputeModuleScore function:
#   1. Load toy Seurat test data
#   2. Test three gene.sets input formats: named list, data.frame, GMT file
#   3. Test three scoring methods: AUCell, Seurat, UCell
#   4. Test two storage modes: assay vs metadata
#
# Requirements:
#   - scMMR installed
#   - Seurat, AUCell, UCell packages
#   - qs package for loading test data
# ============================================================================

library(scMMR)
library(Seurat)
library(qs)

# ── 1. Load toy test data ────────────────────────────────────────────────────
test_path <- system.file("extdata", "toy_test.qs", package = "scMMR")
seu <- qread(test_path)
cat("Loaded Seurat object:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Ensure RNA assay has counts and normalized data
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)

# ── 2. Prepare gene sets in three formats ────────────────────────────────────

# Get genes actually present in the data for building test gene sets
all_genes <- rownames(seu)

## 2a. Named list format
gene_list <- list(
  GeneSet_A = head(all_genes, 30),
  GeneSet_B = all_genes[31:60],
  GeneSet_C = all_genes[61:90]
)
cat("\n[Format: named list]", length(gene_list), "gene sets\n")

## 2b. data.frame format (term-to-gene)
gene_df <- data.frame(
  term = rep(names(gene_list), sapply(gene_list, length)),
  gene = unlist(gene_list, use.names = FALSE)
)
cat("[Format: data.frame]", nrow(gene_df), "rows,", length(unique(gene_df$term)), "terms\n")

## 2c. GMT file format
gmt_path <- system.file("extdata", "gmt", "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")
cat("[Format: GMT file]", gmt_path, "\n")


# ── 3. Test gene.sets input formats ─────────────────────────────────────────
cat("\n========== Test gene.sets input formats ==========\n")

# 3a. Named list
cat("\n--- 3a. Named list input ---\n")
seu_test <- ComputeModuleScore(seu, gene.sets = gene_list, method = "AUCell",cores = 10,
                               min.size = 20, store = "assay")
cat("AUCell assay features:", rownames(seu_test[["AUCell"]]), "\n")

# 3b. data.frame
cat("\n--- 3b. data.frame input ---\n")
seu_test <- ComputeModuleScore(seu, gene.sets = gene_df, method = "AUCell",
                               min.size = 20, store = "assay")
cat("AUCell assay features:", rownames(seu_test[["AUCell"]]), "\n")

# 3c. GMT file
cat("\n--- 3c. GMT file input ---\n")
gmt_sets <- parse_gene_sets(gmt_path)
cat("Parsed", length(gmt_sets), "gene sets from GMT\n")
cat("First 3 sets:", paste(head(names(gmt_sets), 3), collapse = ", "), "\n")
# Use a smaller min.size since toy data has fewer genes
seu_test <- ComputeModuleScore(seu, gene.sets = gmt_path, method = "AUCell",
                               min.size = 5, store = "metadata")
cat("Scored gene sets:", nrow(seu_test[["AUCell"]]), "\n")


# ── 4. Test three scoring methods ────────────────────────────────────────────
cat("\n========== Test scoring methods ==========\n")

# 4a. AUCell (uses counts, batch processing)
cat("\n--- 4a. AUCell method ---\n")
seu_test <- ComputeModuleScore(seu, gene.sets = gene_list, method = "AUCell",
                               min.size = 20, batch.size = 100,
                               store = "assay", assay.name = "AUCell")
cat("AUCell scores range:",
    round(range(GetAssayData(seu_test, assay = "AUCell")), 4), "\n")

# 4b. Seurat method (uses normalized data, AddModuleScore-like)
cat("\n--- 4b. Seurat method ---\n")
seu_test <- ComputeModuleScore(seu_test, gene.sets = gene_list, method = "Seurat",
                               min.size = 20, nbin = 10, ctrl = 50,
                               store = "assay", assay.name = "SeuratScore")
cat("Seurat scores range:",
    round(range(GetAssayData(seu_test, assay = "SeuratScore")), 4), "\n")

# 4c. UCell method
cat("\n--- 4c. UCell method ---\n")
seu_test <- ComputeModuleScore(seu_test, gene.sets = gene_list, method = "UCell",
                               min.size = 20,
                               store = "assay", assay.name = "UCell")
cat("UCell scores range:",
    round(range(GetAssayData(seu_test, assay = "UCell")), 4), "\n")

cat("\nAll assays:", Assays(seu_test), "\n")


# ── 5. Test storage modes ────────────────────────────────────────────────────
cat("\n========== Test storage modes ==========\n")

# 5a. Store as assay (default)
cat("\n--- 5a. store = 'assay' ---\n")
seu_assay <- ComputeModuleScore(seu, gene.sets = gene_list, method = "AUCell",
                                min.size = 20, store = "assay", assay.name = "MyScore")
cat("Assays:", Assays(seu_assay), "\n")
cat("Score assay features:", rownames(seu_assay[["MyScore"]]), "\n")

# 5b. Store as metadata
cat("\n--- 5b. store = 'metadata' ---\n")
seu_meta <- ComputeModuleScore(seu, gene.sets = gene_list, method = "AUCell",
                               min.size = 20, store = "metadata", prefix = "AUC")
score_cols <- grep("^AUC_", colnames(seu_meta@meta.data), value = TRUE)
cat("New metadata columns:", score_cols, "\n")
cat("Score values (first cell):\n")
print(seu_meta@meta.data[1, score_cols, drop = FALSE])

# 5c. Store as metadata with Seurat method
cat("\n--- 5c. Seurat method -> metadata ---\n")
seu_meta2 <- ComputeModuleScore(seu, gene.sets = gene_list, method = "Seurat",
                                min.size = 20, store = "metadata", prefix = "Mod")
score_cols2 <- grep("^Mod_", colnames(seu_meta2@meta.data), value = TRUE)
cat("New metadata columns:", score_cols2, "\n")


# ── 6. Test on raw matrix (default method) ───────────────────────────────────
cat("\n========== Test on raw matrix ==========\n")
mat <- GetAssayData(seu, layer = "counts")
score_mat <- ComputeModuleScore(mat, gene.sets = gene_list, method = "AUCell",
                                min.size = 20, batch.size = 100)
cat("Score matrix dims:", dim(score_mat), "(gene_sets x cells)\n")
cat("Row names:", rownames(score_mat), "\n")
cat("Score range:", round(range(score_mat), 4), "\n")


# ── 7. Summary ───────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("All tests completed successfully!\n")
cat("========================================\n")
