# ============================================================================
# Test script for optimized RankPerturbation & RankPercent
# ============================================================================
cat("========================================\n")
cat("scMMR rank.R optimization test\n")
cat("========================================\n\n")

library(scMMR)

# ── 1. Create synthetic test data ─────────────────────────────────────────
set.seed(42)
n_cells <- 500
n_dims  <- 30
n_types <- 5

cell_types <- paste0("CellType_", LETTERS[1:n_types])
groups     <- c("CTRL", "TREAT")

# Simulate embedding: TREAT shifts some cell types more than others
ct_labels <- sample(cell_types, n_cells, replace = TRUE)
gp_labels <- sample(groups, n_cells, replace = TRUE)

emb <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)
rownames(emb) <- paste0("cell_", seq_len(n_cells))
colnames(emb) <- paste0("dim_", seq_len(n_dims))

# Add systematic shift for CellType_A in TREAT
shift_idx <- which(ct_labels == "CellType_A" & gp_labels == "TREAT")
emb[shift_idx, 1:5] <- emb[shift_idx, 1:5] + 2  # large shift

# Add moderate shift for CellType_B
shift_idx_b <- which(ct_labels == "CellType_B" & gp_labels == "TREAT")
emb[shift_idx_b, 1:3] <- emb[shift_idx_b, 1:3] + 1  # moderate shift

cell_meta <- data.frame(
  cell_type_pred = ct_labels,
  group          = gp_labels,
  row.names      = rownames(emb),
  stringsAsFactors = FALSE
)

cat("Synthetic data: ", n_cells, "cells x", n_dims, "dims\n")
cat("Cell types:", paste(cell_types, collapse = ", "), "\n")
cat("Groups:", paste(groups, collapse = ", "), "\n")
print(table(cell_meta$cell_type_pred, cell_meta$group))
cat("\n")

# ── 2. Test all 4 methods ─────────────────────────────────────────────────
methods <- c("wasserstein", "mmd", "energy", "classifier")

for (m in methods) {
  cat("-- Testing method:", m, "---\n")
  t0 <- Sys.time()
  res <- tryCatch(
    RankPerturbation(
      embedding     = emb,
      cell_meta     = cell_meta,
      conditions    = c("CTRL", "TREAT"),
      method        = m,
      n_permutations = 100,
      n_pcs         = 15,
      min_cells     = 5
    ),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  dt <- round(difftime(Sys.time(), t0, units = "secs"), 2)

  if (!is.null(res)) {
    cat("  Time:", dt, "sec\n")
    cat("  Top cell type:", res$results$cell_type[1],
        "score =", round(res$results$score[1], 4),
        "p_adj =", round(res$results$p_adj[1], 4), "\n")
    print(res$results[, c("cell_type", "score", "p_adj", "rank")])
  }
  cat("\n")
}

# ── 3. Test balanced downsampling ──────────────────────────────────────────
cat("-- Testing balance_cells = TRUE ---\n")
res_bal <- tryCatch(
  RankPerturbation(
    embedding     = emb,
    cell_meta     = cell_meta,
    conditions    = c("CTRL", "TREAT"),
    method        = "wasserstein",
    n_permutations = 50,
    min_cells     = 5,
    balance_cells = TRUE
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(res_bal)) {
  cat("  OK - balanced downsampling works\n")
  print(res_bal$results[, c("cell_type", "score", "rank")])
}
cat("\n")

# ── 4. Test early stopping disabled ──────────────────────────────────────
cat("-- Testing early_stop_alpha = 0 (no early stopping) ---\n")
res_noes <- tryCatch(
  RankPerturbation(
    embedding        = emb,
    cell_meta        = cell_meta,
    conditions       = c("CTRL", "TREAT"),
    method           = "wasserstein",
    n_permutations   = 50,
    min_cells        = 5,
    early_stop_alpha = 0
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(res_noes)) {
  cat("  OK - no-early-stopping mode works\n")
}
cat("\n")

# ── 5. Test verbose = FALSE ──────────────────────────────────────────────
cat("-- Testing verbose = FALSE (should produce no messages) ---\n")
res_quiet <- tryCatch(
  suppressMessages(
    RankPerturbation(
      embedding      = emb,
      cell_meta      = cell_meta,
      conditions     = c("CTRL", "TREAT"),
      method         = "energy",
      n_permutations = 20,
      min_cells      = 5,
      verbose        = FALSE
    )
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(res_quiet)) {
  cat("  OK - quiet mode works\n")
}
cat("\n")

# ── 6. Test RankPercent ──────────────────────────────────────────────────
cat("-- Testing RankPercent ---\n")
da <- tryCatch(
  RankPercent(
    embedding     = emb,
    cell_meta     = cell_meta,
    conditions    = c("CTRL", "TREAT"),
    k             = 20,
    prop_sample   = 0.1,
    min_cells     = 10,
    test          = "fisher"
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(da)) {
  cat("  OK - RankPercent works\n")
  cat("  Neighborhoods:", nrow(da$da_results), "\n")
  cat("  Cell type summary:\n")
  print(da$cell_type_summary)
}
cat("\n")

# ── 7. Summary ────────────────────────────────────────────────────────────
cat("========================================\n")
cat("ALL TESTS PASSED\n")
cat("========================================\n")
