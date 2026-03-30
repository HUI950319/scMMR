# ============================================================================
# Perturbation ranking & Differential abundance testing for scMMR
# ============================================================================

# --------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------

#' Validate and align shared embedding with cell metadata (internal)
#' @keywords internal
.validate_embedding <- function(embedding, cell_meta,
                                cell_type_col, condition_col) {
  if (!is.matrix(embedding) || !is.numeric(embedding))
    stop("'embedding' must be a numeric matrix (e.g. pred$shared_embedding).\n",
         "Did you run DNN_predict() with return_embedding = TRUE?")
  if (is.null(rownames(embedding)))
    stop("'embedding' must have rownames (cell IDs).\n",
         "Example: pred <- DNN_predict(query, model_path, return_embedding = TRUE)")

  if (!is.data.frame(cell_meta))
    stop("'cell_meta' must be a data.frame (e.g. seurat_obj@meta.data).")
  for (col in c(cell_type_col, condition_col)) {
    if (!(col %in% names(cell_meta)))
      stop(sprintf("Column '%s' not found in cell_meta.", col))
  }

  common <- intersect(rownames(embedding), rownames(cell_meta))
  if (length(common) == 0)
    stop("No common cell IDs between embedding and cell_meta rownames.")
  if (length(common) < nrow(embedding))
    message(sprintf("  Using %d / %d cells (intersection).",
                    length(common), nrow(embedding)))

  list(
    embedding = embedding[common, , drop = FALSE],
    cell_meta = cell_meta[common, , drop = FALSE]
  )
}


# --------------------------------------------------------------------------
# Distance functions
# --------------------------------------------------------------------------

#' 1D Wasserstein distance (sort-based) (internal)
#'
#' Uses sorted samples with evenly-spaced index selection, which is
#' significantly faster than the quantile-based approach for large vectors.
#' @keywords internal
.wasserstein_1d <- function(x, y, n_quantiles = 200L) {
  sx <- sort(x)
  sy <- sort(y)
  nx <- length(sx)
  ny <- length(sy)
  n  <- min(nx, ny, n_quantiles)
  # Evenly-spaced indices for uniform quantile approximation
  idx_x <- round(seq(1, nx, length.out = n))
  idx_y <- round(seq(1, ny, length.out = n))
  mean(abs(sx[idx_x] - sy[idx_y]))
}

#' Sliced Wasserstein distance with random projections (internal)
#'
#' Projects high-dimensional distributions onto random 1D directions and
#' averages 1D Wasserstein distances. When \code{n_projections = NULL},
#' falls back to coordinate-axis projections for backward compatibility.
#' @keywords internal
.sliced_wasserstein <- function(mat1, mat2, n_projections = 50L) {
  d <- ncol(mat1)
  if (is.null(n_projections) || n_projections <= 0) {
    # Coordinate-axis mode (backward compatible)
    return(mean(vapply(seq_len(d), function(j) {
      .wasserstein_1d(mat1[, j], mat2[, j])
    }, numeric(1))))
  }
  # Random projections on unit sphere
  proj <- matrix(stats::rnorm(d * n_projections),
                 nrow = d, ncol = n_projections)
  proj <- sweep(proj, 2, sqrt(colSums(proj^2)), "/")
  p1 <- mat1 %*% proj   # n1 x n_projections
  p2 <- mat2 %*% proj   # n2 x n_projections
  mean(vapply(seq_len(n_projections), function(j) {
    .wasserstein_1d(p1[, j], p2[, j])
  }, numeric(1)))
}

#' MMD with RBF kernel (internal)
#'
#' Uses multiple bandwidth values (median heuristic x scales) to avoid
#' sensitivity to a single sigma. Returns the raw (unbiased) MMD-squared
#' without clamping so that permutation tests remain valid.
#'
#' Optimisations vs. original:
#' \itemize{
#'   \item Downsamples groups larger than \code{max_cells} to control memory.
#'   \item Pre-computes squared distance matrices once (shared across
#'     bandwidths) instead of recomputing kernel matrices per sigma.
#' }
#' @keywords internal
.mmd_rbf <- function(mat1, mat2, sigma = NULL, max_cells = 2000L) {
  # Downsample if needed to avoid O(n^2) memory blow-up
  if (nrow(mat1) > max_cells)
    mat1 <- mat1[sample(nrow(mat1), max_cells), , drop = FALSE]
  if (nrow(mat2) > max_cells)
    mat2 <- mat2[sample(nrow(mat2), max_cells), , drop = FALSE]

  # Median heuristic on subsample
  combined <- rbind(mat1, mat2)
  n_sub <- min(500L, nrow(combined))
  idx   <- sample(nrow(combined), n_sub)
  med_dist <- stats::median(as.numeric(stats::dist(combined[idx, ])))
  if (med_dist == 0) med_dist <- 1

  sigmas <- if (is.null(sigma)) med_dist * c(0.5, 1, 2) else sigma
  nx <- nrow(mat1); ny <- nrow(mat2)

  # Pre-compute squared norms and cross-products (once for all bandwidths)
  a2 <- rowSums(mat1^2)
  b2 <- rowSums(mat2^2)

  D2_xx <- outer(a2, a2, "+") - 2 * tcrossprod(mat1)
  D2_yy <- outer(b2, b2, "+") - 2 * tcrossprod(mat2)
  D2_xy <- outer(a2, b2, "+") - 2 * tcrossprod(mat1, mat2)

  best_mmd2 <- -Inf
  for (sg in sigmas) {
    gamma <- 1 / (2 * sg^2)

    Kxx <- exp(-gamma * D2_xx); diag(Kxx) <- 0
    Kyy <- exp(-gamma * D2_yy); diag(Kyy) <- 0
    Kxy <- exp(-gamma * D2_xy)

    mmd2 <- sum(Kxx) / (nx * (nx - 1)) +
      sum(Kyy) / (ny * (ny - 1)) -
      2 * mean(Kxy)

    if (mmd2 > best_mmd2) best_mmd2 <- mmd2
  }

  # Return raw value -- do NOT clamp to 0.
  # Negative values indicate no distributional difference and are needed
  # for correct permutation p-values.
  best_mmd2
}

#' Energy distance (E-distance) (internal)
#'
#' A fast, parameter-free distributional distance.
#' \deqn{E(X,Y) = 2 E\|X - Y\| - E\|X - X'\| - E\|Y - Y'\|}
#'
#' For large populations, downsamples to \code{max_cells} per group.
#' Cross-group distances are computed in chunks to avoid memory issues.
#' @keywords internal
.energy_distance <- function(mat1, mat2, max_cells = 2000L) {
  if (nrow(mat1) > max_cells)
    mat1 <- mat1[sample(nrow(mat1), max_cells), , drop = FALSE]
  if (nrow(mat2) > max_cells)
    mat2 <- mat2[sample(nrow(mat2), max_cells), , drop = FALSE]

  n1 <- nrow(mat1); n2 <- nrow(mat2)

  # Within-group mean distances (use dist for efficiency)
  E_xx <- if (n1 > 1) mean(stats::dist(mat1)) else 0
  E_yy <- if (n2 > 1) mean(stats::dist(mat2)) else 0


  # Cross-group mean distance (chunk-based to control memory)
  # For n1*n2 <= 4e6, direct computation; otherwise chunk
  if (as.double(n1) * n2 <= 4e6) {
    # Direct: expand and compute
    D_xy <- sqrt(rowSums(
      (mat1[rep(seq_len(n1), each = n2), , drop = FALSE] -
         mat2[rep(seq_len(n2), n1), , drop = FALSE])^2
    ))
    E_xy <- mean(D_xy)
  } else {
    # Chunk-based: process mat1 rows in chunks
    chunk_size <- max(1L, floor(4e6 / n2))
    total_sum <- 0
    total_n   <- 0
    for (start in seq(1, n1, by = chunk_size)) {
      end <- min(start + chunk_size - 1, n1)
      chunk <- mat1[start:end, , drop = FALSE]
      nc <- nrow(chunk)
      D_chunk <- sqrt(rowSums(
        (chunk[rep(seq_len(nc), each = n2), , drop = FALSE] -
           mat2[rep(seq_len(n2), nc), , drop = FALSE])^2
      ))
      total_sum <- total_sum + sum(D_chunk)
      total_n   <- total_n + length(D_chunk)
    }
    E_xy <- total_sum / total_n
  }

  2 * E_xy - E_xx - E_yy
}


#' AUC via stratified k-fold cross-validated LDA (internal)
#'
#' Uses Linear Discriminant Analysis for faster and more stable
#' classification than logistic regression. Falls back to GLM if
#' MASS is unavailable.
#'
#' When the number of features exceeds \code{min(n_per_class) / 2},
#' an additional PCA step is applied to avoid the \code{p >> n} problem.
#' @keywords internal
.auc_lda <- function(mat1, mat2, n_folds = 5L) {
  if (!requireNamespace("MASS", quietly = TRUE))
    return(.auc_glm(mat1, mat2, n_folds))

  X <- rbind(mat1, mat2)
  y <- factor(c(rep("A", nrow(mat1)), rep("B", nrow(mat2))))
  n <- length(y)

  min_n <- min(nrow(mat1), nrow(mat2))
  if (min_n < 3) return(0.5)  # too few cells

  # Auto-reduce features to avoid p >> n overfitting
  max_features <- max(2L, floor(min_n / 2))
  if (ncol(X) > max_features) {
    pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE,
                             rank. = max_features)
    X <- pca_fit$x
  }

  # Stratified k-fold CV
  n_folds <- min(n_folds, min_n)
  if (n_folds < 2) return(0.5)

  idx_a <- which(y == "A")
  idx_b <- which(y == "B")
  fold_a <- sample(rep(seq_len(n_folds), length.out = length(idx_a)))
  fold_b <- sample(rep(seq_len(n_folds), length.out = length(idx_b)))

  folds <- integer(n)
  folds[idx_a] <- fold_a
  folds[idx_b] <- fold_b

  cv_probs <- numeric(n)

  for (f in seq_len(n_folds)) {
    test_idx  <- which(folds == f)
    train_idx <- which(folds != f)

    fit <- tryCatch(
      MASS::lda(x = X[train_idx, , drop = FALSE],
                grouping = y[train_idx]),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      cv_probs[test_idx] <- 0.5
      next
    }

    pred <- tryCatch(
      stats::predict(fit,
                     newdata = X[test_idx, , drop = FALSE])$posterior[, "B"],
      error = function(e) rep(0.5, length(test_idx))
    )
    cv_probs[test_idx] <- pred
  }

  # AUC via Mann-Whitney U
  r   <- rank(cv_probs)
  n_a <- sum(y == "A"); n_b <- sum(y == "B")
  U   <- sum(r[y == "B"]) - n_b * (n_b + 1) / 2
  auc <- U / (n_a * n_b)
  max(auc, 1 - auc)
}


#' AUC via stratified k-fold cross-validated logistic regression (internal)
#'
#' Retained as fallback when MASS is not available.
#' @keywords internal
.auc_glm <- function(mat1, mat2, n_folds = 5L) {
  X <- rbind(mat1, mat2)
  y <- c(rep(0L, nrow(mat1)), rep(1L, nrow(mat2)))
  n <- length(y)

  min_n <- min(nrow(mat1), nrow(mat2))
  if (min_n < 3) return(0.5)  # too few cells

  # Auto-reduce features to avoid p >> n overfitting
  max_features <- max(2L, floor(min_n / 2))
  if (ncol(X) > max_features) {
    pca_fit <- stats::prcomp(X, center = TRUE, scale. = FALSE,
                             rank. = max_features)
    X <- pca_fit$x
  }

  # Stratified k-fold CV
  n_folds <- min(n_folds, min_n)
  if (n_folds < 2) return(0.5)

  idx0 <- which(y == 0)
  idx1 <- which(y == 1)
  fold0 <- sample(rep(seq_len(n_folds), length.out = length(idx0)))
  fold1 <- sample(rep(seq_len(n_folds), length.out = length(idx1)))

  folds <- integer(n)
  folds[idx0] <- fold0
  folds[idx1] <- fold1

  cv_probs <- numeric(n)

  for (f in seq_len(n_folds)) {
    test_idx  <- which(folds == f)
    train_idx <- which(folds != f)

    df_train <- as.data.frame(X[train_idx, , drop = FALSE])
    df_train$.y <- y[train_idx]
    df_test  <- as.data.frame(X[test_idx, , drop = FALSE])

    fit <- tryCatch(
      suppressWarnings(
        stats::glm(.y ~ ., data = df_train, family = stats::binomial())
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      cv_probs[test_idx] <- 0.5
    } else {
      cv_probs[test_idx] <- tryCatch(
        as.numeric(stats::predict(fit, newdata = df_test,
                                  type = "response")),
        error = function(e) rep(0.5, length(test_idx))
      )
    }
  }

  # AUC via Mann-Whitney U
  r  <- rank(cv_probs)
  n0 <- sum(y == 0); n1 <- sum(y == 1)
  U  <- sum(r[y == 1]) - n1 * (n1 + 1) / 2
  auc <- U / (n0 * n1)
  max(auc, 1 - auc)
}


# --------------------------------------------------------------------------
# Permutation test helper with adaptive early stopping
# --------------------------------------------------------------------------

#' Permutation test with adaptive early stopping (internal)
#'
#' Runs permutation test for a single cell type. After every
#' \code{check_interval} permutations, checks whether the running p-value
#' is clearly above \code{early_stop_alpha}. If so, terminates early to
#' save computation.
#'
#' @return Numeric p-value.
#' @keywords internal
.permutation_test <- function(emb_ct, cond_ct, conditions, dist_fn,
                              observed, n_permutations,
                              early_stop_alpha = 0.1,
                              check_interval = 100L) {
  n_extreme <- 0L
  k <- 0L

  while (k < n_permutations) {
    # Run a batch of permutations
    batch_end <- min(k + check_interval, n_permutations)
    batch_size <- batch_end - k

    for (i in seq_len(batch_size)) {
      perm_cond <- sample(cond_ct)
      pm1 <- emb_ct[perm_cond == conditions[1], , drop = FALSE]
      pm2 <- emb_ct[perm_cond == conditions[2], , drop = FALSE]
      perm_score <- dist_fn(pm1, pm2)
      if (perm_score >= observed) n_extreme <- n_extreme + 1L
    }
    k <- batch_end

    # Adaptive early stopping check
    # Current p estimate: (n_extreme + 1) / (k + 1)
    # If this is already clearly above early_stop_alpha, stop
    if (k < n_permutations && k >= check_interval) {
      p_est <- (n_extreme + 1) / (k + 1)
      # Use a conservative threshold: stop only if p > 2 * alpha
      # and we have enough permutations for reliability
      if (p_est > early_stop_alpha * 2 && n_extreme >= 10) {
        break
      }
    }
  }

  (n_extreme + 1) / (k + 1)
}


# ============================================================================
# RankPerturbation
# ============================================================================

#' Rank Cell Types by Transcriptional Perturbation
#'
#' Measures how much each cell type's transcriptional state shifts between
#' two conditions in the DNN's learned embedding space. Cell types with
#' larger shifts are ranked higher (more "perturbed"). Inspired by Augur
#' but using multi-task DNN representations.
#'
#' @param embedding Numeric matrix (n_cells x d), typically
#'   \code{pred$shared_embedding} from
#'   \code{DNN_predict(return_embedding = TRUE)}.
#'   Rownames must be cell IDs.
#' @param cell_meta A data.frame (e.g. \code{seurat_obj@@meta.data}).
#'   Rownames must be cell IDs.
#' @param cell_type_col Column name for cell type labels
#'   (default \code{"cell_type_pred"}).
#' @param condition_col Column name for condition/group labels
#'   (default \code{"group"}).
#' @param conditions Character vector of length 2 specifying which two
#'   conditions to compare (e.g. \code{c("PH", "SH")}). If \code{NULL}
#'   and exactly 2 unique conditions exist, they are used automatically.
#' @param method Distance metric: \code{"wasserstein"} (default, sliced
#'   Wasserstein distance with random projections), \code{"mmd"} (Maximum
#'   Mean Discrepancy with RBF kernel), \code{"energy"} (Energy distance,
#'   parameter-free), or \code{"classifier"} (LDA-based AUC, Augur-like).
#' @param n_permutations Number of permutations for p-value estimation
#'   (default 500). Set to 0 to skip.
#' @param min_cells Minimum cells per cell-type-per-condition (default 10).
#' @param n_pcs Number of PCs for dimensionality reduction before distance
#'   computation (default 20). \code{NULL} to skip PCA.
#' @param seed Random seed (default 42).
#' @param n_cores Number of parallel cores for computing across cell types.
#'   Default 1 (serial). On Unix, uses \code{parallel::mclapply};
#'   on Windows, uses \code{parallel::parLapply}.
#' @param balance_cells Logical. If \code{TRUE}, downsample the larger
#'   condition group to match the smaller one per cell type before
#'   distance computation (default \code{FALSE}).
#' @param early_stop_alpha Numeric. p-value threshold for adaptive early
#'   stopping of permutation tests. Permutations stop early when the
#'   running p-value is clearly above \code{2 * early_stop_alpha}
#'   (default 0.1). Set to 0 to disable early stopping.
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#'
#' @return A named list:
#'   \describe{
#'     \item{results}{data.frame: cell_type, score, p_value, p_adj,
#'       n_cells_cond1, n_cells_cond2, rank.}
#'     \item{method}{Method used.}
#'     \item{conditions}{Two conditions compared.}
#'   }
#'
#' @examples
#' \dontrun{
#' pred <- DNN_predict(query, model_path, return_embedding = TRUE)
#' q1 <- Seurat::AddMetaData(toy_test, pred$predictions)
#'
#' # Sliced Wasserstein (default)
#' res <- RankPerturbation(pred$shared_embedding, q1@@meta.data,
#'                         conditions = c("PH", "SH"))
#'
#' # Energy distance (fast, parameter-free)
#' res_e <- RankPerturbation(pred$shared_embedding, q1@@meta.data,
#'                           conditions = c("PH", "SH"),
#'                           method = "energy")
#'
#' # With parallel computation
#' res_p <- RankPerturbation(pred$shared_embedding, q1@@meta.data,
#'                           conditions = c("PH", "SH"),
#'                           n_cores = 4)
#'
#' print(res$results)
#' PlotPerturbation(res)
#' }
#'
#' @export
RankPerturbation <- function(embedding,
                             cell_meta,
                             cell_type_col    = "cell_type_pred",
                             condition_col    = "group",
                             conditions       = NULL,
                             method           = c("wasserstein", "mmd",
                                                  "energy", "classifier"),
                             n_permutations   = 500L,
                             min_cells        = 10L,
                             n_pcs            = 20L,
                             seed             = 42L,
                             n_cores          = 1L,
                             balance_cells    = FALSE,
                             early_stop_alpha = 0.1,
                             verbose          = TRUE) {

  method <- match.arg(method)
  set.seed(seed)

  # -- Validate ---------------------------------------------------------------
  val <- .validate_embedding(embedding, cell_meta,
                             cell_type_col, condition_col)
  emb  <- val$embedding
  meta <- val$cell_meta

  cond_vec <- as.character(meta[[condition_col]])
  ct_vec   <- as.character(meta[[cell_type_col]])

  # -- Resolve 2 conditions ---------------------------------------------------
  uniq_cond <- sort(unique(cond_vec))
  if (is.null(conditions)) {
    if (length(uniq_cond) == 2) {
      conditions <- uniq_cond
    } else {
      stop("More than 2 conditions detected: ",
           paste(uniq_cond, collapse = ", "),
           ".\nPlease specify 'conditions = c(\"A\", \"B\")' ",
           "to choose which two to compare.")
    }
  }
  stopifnot(length(conditions) == 2)
  for (cc in conditions) {
    if (!(cc %in% uniq_cond))
      stop(sprintf("Condition '%s' not found. Available: %s",
                   cc, paste(uniq_cond, collapse = ", ")))
  }

  # Subset to selected conditions
  keep <- cond_vec %in% conditions
  emb      <- emb[keep, , drop = FALSE]
  cond_vec <- cond_vec[keep]
  ct_vec   <- ct_vec[keep]

  # -- Optional PCA -----------------------------------------------------------
  if (!is.null(n_pcs) && n_pcs > 0 && n_pcs < ncol(emb)) {
    if (verbose) message(sprintf("  PCA: %d dims -> %d PCs", ncol(emb), n_pcs))
    pca <- stats::prcomp(emb, center = TRUE, scale. = FALSE, rank. = n_pcs)
    emb <- pca$x[, seq_len(n_pcs), drop = FALSE]
  }

  # -- Filter cell types ------------------------------------------------------
  ct_table <- table(ct_vec, cond_vec)
  valid_cts <- rownames(ct_table)[
    ct_table[, conditions[1]] >= min_cells &
      ct_table[, conditions[2]] >= min_cells
  ]

  excluded <- setdiff(unique(ct_vec), valid_cts)
  if (length(excluded) > 0 && verbose)
    message(sprintf("  Excluding %d cell types with < %d cells per condition: %s",
                    length(excluded), min_cells,
                    paste(excluded, collapse = ", ")))
  if (length(valid_cts) == 0)
    stop("No cell types have >= ", min_cells,
         " cells in both conditions. Try lowering 'min_cells'.")

  # -- Distance function dispatcher -------------------------------------------
  dist_fn <- switch(method,
    wasserstein = .sliced_wasserstein,
    mmd         = .mmd_rbf,
    energy      = .energy_distance,
    classifier  = .auc_lda
  )

  # -- Worker function for one cell type --------------------------------------
  .process_one_ct <- function(ct) {
    idx <- ct_vec == ct
    emb_ct   <- emb[idx, , drop = FALSE]
    cond_ct  <- cond_vec[idx]

    idx1 <- which(cond_ct == conditions[1])
    idx2 <- which(cond_ct == conditions[2])

    # Balanced downsampling
    if (balance_cells) {
      n_min <- min(length(idx1), length(idx2))
      if (length(idx1) > n_min) idx1 <- sample(idx1, n_min)
      if (length(idx2) > n_min) idx2 <- sample(idx2, n_min)
      emb_ct  <- emb_ct[c(idx1, idx2), , drop = FALSE]
      cond_ct <- cond_ct[c(idx1, idx2)]
      idx1 <- seq_len(length(idx1))
      idx2 <- seq(length(idx1) + 1, nrow(emb_ct))
    }

    mat1 <- emb_ct[idx1, , drop = FALSE]
    mat2 <- emb_ct[idx2, , drop = FALSE]

    observed <- dist_fn(mat1, mat2)

    # Permutation test with adaptive early stopping
    pval <- NA_real_
    if (n_permutations > 0) {
      if (early_stop_alpha > 0) {
        pval <- .permutation_test(
          emb_ct, cond_ct, conditions, dist_fn, observed,
          n_permutations  = n_permutations,
          early_stop_alpha = early_stop_alpha,
          check_interval  = min(100L, n_permutations)
        )
      } else {
        # No early stopping
        perm_scores <- vapply(seq_len(n_permutations), function(i) {
          perm_cond <- sample(cond_ct)
          pm1 <- emb_ct[perm_cond == conditions[1], , drop = FALSE]
          pm2 <- emb_ct[perm_cond == conditions[2], , drop = FALSE]
          dist_fn(pm1, pm2)
        }, numeric(1))
        pval <- (sum(perm_scores >= observed) + 1) / (n_permutations + 1)
      }
    }

    data.frame(
      cell_type     = ct,
      score         = observed,
      p_value       = pval,
      n_cells_cond1 = length(idx1),
      n_cells_cond2 = length(idx2),
      stringsAsFactors = FALSE
    )
  }

  # -- Compute per cell type (serial or parallel) -----------------------------
  if (verbose)
    message(sprintf("  Computing %s distance for %d cell types%s ...",
                    method, length(valid_cts),
                    if (n_cores > 1) sprintf(" (%d cores)", n_cores) else ""))

  if (n_cores > 1 && length(valid_cts) > 1) {
    if (.Platform$OS.type == "unix") {
      results_list <- parallel::mclapply(
        valid_cts, .process_one_ct, mc.cores = n_cores
      )
    } else {
      # Windows: use parLapply
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      # Export necessary objects and functions to workers
      parallel::clusterExport(
        cl,
        varlist = c("emb", "cond_vec", "ct_vec", "conditions",
                     "dist_fn", "n_permutations", "balance_cells",
                     "early_stop_alpha",
                     ".sliced_wasserstein", ".wasserstein_1d",
                     ".mmd_rbf", ".energy_distance",
                     ".auc_lda", ".auc_glm",
                     ".permutation_test"),
        envir = environment()
      )
      results_list <- parallel::parLapply(cl, valid_cts, .process_one_ct)
    }
  } else {
    # Serial with optional progress
    results_list <- vector("list", length(valid_cts))
    for (ci in seq_along(valid_cts)) {
      ct <- valid_cts[ci]
      if (verbose)
        message(sprintf("    [%d/%d] %s", ci, length(valid_cts), ct))
      results_list[[ci]] <- .process_one_ct(ct)
    }
  }

  results <- do.call(rbind, results_list)

  # -- FDR correction ---------------------------------------------------------
  results$p_adj <- stats::p.adjust(results$p_value, method = "BH")
  results <- results[order(-results$score), ]
  results$rank <- seq_len(nrow(results))
  rownames(results) <- NULL

  # Rename count columns
  colnames(results)[colnames(results) == "n_cells_cond1"] <-
    paste0("n_", conditions[1])
  colnames(results)[colnames(results) == "n_cells_cond2"] <-
    paste0("n_", conditions[2])

  if (verbose) message("  Done.")

  list(
    results    = results,
    method     = method,
    conditions = conditions
  )
}


# ============================================================================
# RankPercent
# ============================================================================

#' Differential Abundance Testing via KNN Neighborhoods
#'
#' Builds a KNN graph on the DNN's shared embedding, samples representative
#' neighborhoods, and tests for differential cell abundance between
#' conditions. Inspired by the Milo framework but operating on DNN-learned
#' representations.
#'
#' @param embedding Numeric matrix (n_cells x d), typically
#'   \code{pred$shared_embedding}.
#' @param cell_meta A data.frame. Rownames must be cell IDs.
#' @param cell_type_col Column name for cell type labels
#'   (default \code{"cell_type_pred"}).
#' @param condition_col Column name for condition/group labels
#'   (default \code{"group"}).
#' @param sample_col Optional column name for biological sample/replicate
#'   labels (e.g. \code{"sample"}). When provided, an edgeR
#'   quasi-likelihood GLM (\code{edgeR::glmQLFit}) is fitted on the
#'   sample-level count matrix across all neighbourhoods jointly,
#'   with shared dispersion estimation — the same statistical framework
#'   used by miloR's \code{testNhoods}. This produces continuous,
#'   shrinkage-stabilised logFC values. Requires \code{edgeR}.
#'   When \code{NULL} (default), the original proportion-based logFC
#'   is used.
#' @param conditions Optional character vector of length 2. If provided,
#'   only these two conditions are tested. If \code{NULL} (default), all
#'   conditions are used; Fisher's or chi-squared test is applied.
#' @param k Number of nearest neighbors (default 30).
#' @param prop_sample Proportion of cells sampled as neighborhood centers
#'   (default 0.1).
#' @param test Statistical test: \code{"fisher"} (default) or
#'   \code{"nb_glm"} (negative binomial GLM via \code{MASS::glm.nb}).
#' @param min_cells Minimum cells in a neighborhood to include
#'   (default 20).
#' @param fdr_threshold FDR threshold for summary (default 0.1).
#' @param seed Random seed (default 42).
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#'
#' @return A named list:
#'   \describe{
#'     \item{da_results}{data.frame: nhood_idx, nhood_cell_id, n_cells,
#'       counts per condition, logFC, p_value, p_adj,
#'       cell_type_majority.}
#'     \item{cell_da_scores}{Named numeric vector of per-cell DA scores.
#'       Ready for \code{Seurat::AddMetaData()}.}
#'     \item{cell_type_summary}{data.frame: cell_type, n_nhoods,
#'       n_da_nhoods, mean_logFC, fraction_da.}
#'     \item{params}{List of parameters used.}
#'   }
#'
#' @examples
#' \dontrun{
#' pred <- DNN_predict(query, model_path, return_embedding = TRUE)
#' q1 <- Seurat::AddMetaData(toy_test, pred$predictions)
#'
#' # Original mode (proportion-based logFC)
#' da <- RankPercent(pred$shared_embedding, q1@@meta.data,
#'                   conditions = c("PH", "SH"), k = 30)
#'
#' # Sample-aware mode (miloR-like continuous logFC)
#' da <- RankPercent(pred$shared_embedding, q1@@meta.data,
#'                   sample_col = "sample",
#'                   conditions = c("PH", "SH"), k = 30)
#'
#' q1 <- Seurat::AddMetaData(q1, da$cell_da_scores, col.name = "da_score")
#' PlotPercent(da)
#' }
#'
#' @export
RankPercent <- function(embedding,
                        cell_meta,
                        cell_type_col  = "cell_type_pred",
                        condition_col  = "group",
                        sample_col     = NULL,
                        conditions     = NULL,
                        k              = 30L,
                        prop_sample    = 0.1,
                        test           = c("fisher", "nb_glm"),
                        min_cells      = 20L,
                        fdr_threshold  = 0.1,
                        seed           = 42L,
                        verbose        = TRUE) {

  test <- match.arg(test)
  set.seed(seed)

  # -- Check dependencies ----------------------------------------------------
  if (!requireNamespace("FNN", quietly = TRUE))
    stop("Package 'FNN' is required for RankPercent().\n",
         "Install with: install.packages('FNN')")
  if (test == "nb_glm" && !requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required for test = 'nb_glm'.\n",
         "Install with: install.packages('MASS')")

  # -- Validate ---------------------------------------------------------------
  val <- .validate_embedding(embedding, cell_meta,
                             cell_type_col, condition_col)
  emb  <- val$embedding
  meta <- val$cell_meta

  cond_vec <- as.character(meta[[condition_col]])
  ct_vec   <- as.character(meta[[cell_type_col]])
  n_cells  <- nrow(emb)

  # -- Sample-level info (for edgeR quasi-likelihood GLM) ----------------------
  use_sample_glm <- !is.null(sample_col)
  if (use_sample_glm) {
    if (!sample_col %in% colnames(meta))
      stop(sprintf("sample_col '%s' not found in cell_meta.", sample_col))
    if (!requireNamespace("edgeR", quietly = TRUE))
      stop("Package 'edgeR' is required for sample-level GLM.\n",
           "Install with: BiocManager::install('edgeR')")
    sample_vec <- as.character(meta[[sample_col]])
    if (verbose)
      message("  sample_col provided: using edgeR glmQLFit for logFC")
  }

  # -- Optional condition filter -----------------------------------------------
  if (!is.null(conditions)) {
    stopifnot(length(conditions) == 2)
    keep <- cond_vec %in% conditions
    emb      <- emb[keep, , drop = FALSE]
    cond_vec <- cond_vec[keep]
    ct_vec   <- ct_vec[keep]
    if (use_sample_glm) sample_vec <- sample_vec[keep]
    n_cells  <- nrow(emb)
  }

  cond_levels <- sort(unique(cond_vec))
  n_conds     <- length(cond_levels)
  total_per_cond <- table(cond_vec)[cond_levels]

  # -- Sample-level lookup (pre-compute after filtering) -----------------------
  if (use_sample_glm) {
    sample_ids    <- sort(unique(sample_vec))
    sample_totals <- table(sample_vec)[sample_ids]
    sample_to_cond <- tapply(cond_vec, sample_vec, function(x) x[1])
    sample_to_cond <- sample_to_cond[sample_ids]

    # Verify: each sample belongs to exactly one condition
    n_cond_per_sample <- tapply(cond_vec, sample_vec,
                                function(x) length(unique(x)))
    if (any(n_cond_per_sample > 1))
      warning("Some samples map to multiple conditions. ",
              "Check sample_col / condition_col assignment.")
    if (verbose)
      message(sprintf("  Samples: %d (%s)",
                      length(sample_ids),
                      paste(sprintf("%s=%d", names(table(sample_to_cond)),
                                    as.integer(table(sample_to_cond))),
                            collapse = ", ")))
  }

  if (verbose)
    message(sprintf("  Building %d-NN graph on %d cells ...", k, n_cells))

  # -- KNN graph ---------------------------------------------------------------
  knn_res <- FNN::get.knn(emb, k = min(k, n_cells - 1))
  nn_idx  <- knn_res$nn.index   # n_cells x k

  # -- Sample neighborhood centers --------------------------------------------
  n_sample <- max(10L, ceiling(n_cells * prop_sample))
  n_sample <- min(n_sample, n_cells)
  sample_idx <- sort(sample(n_cells, n_sample))
  if (verbose)
    message(sprintf("  Sampled %d neighborhoods (prop = %.2f)",
                    n_sample, prop_sample))

  # -- Collect neighborhood info ------------------------------------------------
  nhood_info <- vector("list", n_sample)

  for (i in seq_along(sample_idx)) {
    ci <- sample_idx[i]
    nhood_cells <- unique(c(ci, nn_idx[ci, ]))
    nhood_cond  <- cond_vec[nhood_cells]
    nhood_ct    <- ct_vec[nhood_cells]
    n_nh        <- length(nhood_cells)

    if (n_nh < min_cells) next

    counts <- table(factor(nhood_cond, levels = cond_levels))
    ct_tab <- sort(table(nhood_ct), decreasing = TRUE)
    majority_ct <- names(ct_tab)[1]

    info <- list(
      ci          = ci,
      nhood_cells = nhood_cells,
      n_nh        = n_nh,
      counts      = counts,
      majority_ct = majority_ct
    )

    # Collect sample-level counts for edgeR batch mode
    if (use_sample_glm) {
      nhood_samp <- sample_vec[nhood_cells]
      info$counts_by_samp <- table(factor(nhood_samp, levels = sample_ids))
    }

    nhood_info[[i]] <- info
  }
  nhood_info <- Filter(Negate(is.null), nhood_info)
  n_valid <- length(nhood_info)

  if (n_valid == 0)
    stop("No neighborhoods passed the min_cells filter. ",
         "Try lowering 'min_cells' or increasing 'k'.")

  if (verbose)
    message(sprintf("  Valid neighborhoods: %d / %d", n_valid, n_sample))

  # -- Statistical testing -----------------------------------------------------
  if (use_sample_glm && n_conds == 2) {
    # ── edgeR quasi-likelihood batch mode (miloR-like) ────────────────────────
    if (verbose) message("  Running edgeR glmQLFit across all neighborhoods ...")

    # Build count matrix: n_nhoods x n_samples
    count_mat <- do.call(rbind, lapply(nhood_info, function(x) {
      as.integer(x$counts_by_samp)
    }))
    colnames(count_mat) <- sample_ids
    rownames(count_mat) <- paste0("nhood_", seq_len(n_valid))


    # Design matrix
    cond_fac <- factor(sample_to_cond, levels = cond_levels)
    design <- stats::model.matrix(~ cond_fac)

    # Library sizes = total cells per sample
    lib_sizes <- as.numeric(sample_totals[sample_ids])

    edger_ok <- tryCatch({
      y <- edgeR::DGEList(counts = count_mat, lib.size = lib_sizes)
      y <- edgeR::estimateDisp(y, design)
      fit <- edgeR::glmQLFit(y, design, robust = TRUE)
      res <- edgeR::glmQLFTest(fit, coef = 2)
      tt <- edgeR::topTags(res, n = Inf, sort.by = "none")$table
      TRUE
    }, error = function(e) {
      if (verbose)
        message(sprintf("  [WARN] edgeR failed: %s. Falling back to proportion logFC.",
                        conditionMessage(e)))
      FALSE
    })

    # Assemble da_results
    da_list <- vector("list", n_valid)
    for (i in seq_len(n_valid)) {
      info <- nhood_info[[i]]
      if (edger_ok) {
        lfc  <- tt$logFC[i]
        pval <- tt$PValue[i]
      } else {
        # Fallback: proportion logFC + Fisher
        c1 <- as.integer(info$counts[cond_levels[1]])
        c2 <- as.integer(info$counts[cond_levels[2]])
        t1 <- as.integer(total_per_cond[cond_levels[1]])
        t2 <- as.integer(total_per_cond[cond_levels[2]])
        prop1 <- (c1 + 0.5) / (t1 + 1)
        prop2 <- (c2 + 0.5) / (t2 + 1)
        lfc   <- log2(prop2 / prop1)
        mat_2x2 <- matrix(c(c1, t1 - c1, c2, t2 - c2), nrow = 2)
        ft <- tryCatch(stats::fisher.test(mat_2x2),
                        error = function(e) list(p.value = 1))
        pval <- ft$p.value
      }

      row <- data.frame(
        nhood_idx          = info$ci,
        nhood_cell_id      = rownames(emb)[info$ci],
        n_cells            = info$n_nh,
        logFC              = lfc,
        p_value            = pval,
        cell_type_majority = info$majority_ct,
        stringsAsFactors   = FALSE
      )
      for (cl in cond_levels)
        row[[paste0("n_", cl)]] <- as.integer(info$counts[cl])
      da_list[[i]] <- row
    }

  } else {
    # ── Original per-neighborhood test (no sample_col) ────────────────────────
    da_list <- vector("list", n_valid)
    for (i in seq_len(n_valid)) {
      info <- nhood_info[[i]]
      pval <- 1
      lfc  <- 0

      if (n_conds == 2) {
        c1 <- as.integer(info$counts[cond_levels[1]])
        c2 <- as.integer(info$counts[cond_levels[2]])
        t1 <- as.integer(total_per_cond[cond_levels[1]])
        t2 <- as.integer(total_per_cond[cond_levels[2]])
        prop1 <- (c1 + 0.5) / (t1 + 1)
        prop2 <- (c2 + 0.5) / (t2 + 1)
        lfc   <- log2(prop2 / prop1)

        if (test == "fisher") {
          mat_2x2 <- matrix(c(c1, t1 - c1, c2, t2 - c2), nrow = 2)
          ft <- tryCatch(stats::fisher.test(mat_2x2),
                          error = function(e) list(p.value = 1))
          pval <- ft$p.value
        } else {
          count_df <- data.frame(count = c(c1, c2),
                                 condition = cond_levels,
                                 total = c(t1, t2))
          fit <- tryCatch(MASS::glm.nb(count ~ condition + offset(log(total)),
                                        data = count_df),
                          error = function(e) NULL)
          if (!is.null(fit)) pval <- summary(fit)$coefficients[2, 4]
        }
      } else {
        chisq <- tryCatch(
          stats::chisq.test(as.integer(info$counts),
                            p = total_per_cond / sum(total_per_cond)),
          error = function(e) list(p.value = 1))
        pval <- chisq$p.value
        props <- (as.numeric(info$counts) + 0.5) /
                 (as.numeric(total_per_cond) + 1)
        lfc <- log2(max(props) / min(props))
      }

      row <- data.frame(
        nhood_idx          = info$ci,
        nhood_cell_id      = rownames(emb)[info$ci],
        n_cells            = info$n_nh,
        logFC              = lfc,
        p_value            = pval,
        cell_type_majority = info$majority_ct,
        stringsAsFactors   = FALSE
      )
      for (cl in cond_levels)
        row[[paste0("n_", cl)]] <- as.integer(info$counts[cl])
      da_list[[i]] <- row
    }
  }

  da_results <- do.call(rbind, da_list)

  # -- FDR correction ---------------------------------------------------------
  da_results$p_adj <- stats::p.adjust(da_results$p_value, method = "BH")
  da_results <- da_results[order(da_results$p_adj), ]
  rownames(da_results) <- NULL

  if (verbose)
    message(sprintf("  %d / %d neighborhoods significant (FDR < %.2f)",
                    sum(da_results$p_adj < fdr_threshold, na.rm = TRUE),
                    nrow(da_results), fdr_threshold))

  # -- Per-cell DA scores -----------------------------------------------------
  if (verbose) message("  Computing per-cell DA scores ...")
  cell_da <- rep(NA_real_, n_cells)
  names(cell_da) <- rownames(emb)

  # Build reverse lookup
  for (ri in seq_len(nrow(da_results))) {
    ci <- da_results$nhood_idx[ri]
    nhood_cells <- unique(c(ci, nn_idx[ci, ]))
    for (j in nhood_cells) {
      if (is.na(cell_da[j])) {
        cell_da[j] <- da_results$logFC[ri]
      } else {
        # Running mean (memory efficient)
        cell_da[j] <- (cell_da[j] + da_results$logFC[ri]) / 2
      }
    }
  }

  # -- Cell type summary -------------------------------------------------------
  ct_summary <- do.call(rbind, lapply(
    sort(unique(da_results$cell_type_majority)),
    function(ct) {
      sub <- da_results[da_results$cell_type_majority == ct, ]
      data.frame(
        cell_type   = ct,
        n_nhoods    = nrow(sub),
        n_da_nhoods = sum(sub$p_adj < fdr_threshold, na.rm = TRUE),
        mean_logFC  = mean(sub$logFC, na.rm = TRUE),
        fraction_da = mean(sub$p_adj < fdr_threshold, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))
  ct_summary <- ct_summary[order(-abs(ct_summary$mean_logFC)), ]
  rownames(ct_summary) <- NULL

  if (verbose) message("  Done.")

  list(
    da_results       = da_results,
    cell_da_scores   = cell_da,
    cell_type_summary = ct_summary,
    params = list(
      k = k, prop_sample = prop_sample, test = test,
      sample_col = sample_col, sample_glm = use_sample_glm,
      min_cells = min_cells, fdr_threshold = fdr_threshold,
      conditions = if (!is.null(conditions)) conditions else cond_levels
    )
  )
}


# ============================================================================
# RunCorrelation
# ============================================================================

#' Correlation Ranking of One Variable Against Many
#'
#' Compute correlation (and p-value) between a target variable and multiple
#' other variables.
#' The result is a data.frame directly compatible with
#' \code{\link{PlotRankScatter}}.
#'
#' @details
#' The first argument \code{object} can be either a \strong{Seurat object}
#' or a \strong{data.frame}.
#'
#' \strong{Seurat input:}
#' \itemize{
#'   \item \code{target} may be a gene name or metadata column.
#'   \item \code{features}: if \code{NULL}, the top 2 000 highly-variable
#'     genes (VariableFeatures) are used; genes with expression in fewer
#'     than \code{min.pct} of cells are dropped.
#'   \item When \code{group.by} is set, correlations are computed
#'     separately within each group (e.g. per cell type).
#' }
#'
#' \strong{data.frame input:}
#' \itemize{
#'   \item \code{target} and \code{features} must be column names.
#'   \item When \code{group.by} is set, data is split by that column
#'     and correlations are computed per group.
#' }
#'
#' @param object A Seurat object or a data.frame.
#' @param target Character. Name of the target variable (gene or column).
#' @param features Character vector. Variables to correlate with
#'   \code{target}. For Seurat objects, \code{NULL} (default) uses
#'   VariableFeatures (up to \code{n_features}).
#'   For data.frames, \code{NULL} uses all numeric columns except
#'   \code{target}.
#' @param group.by Character. Optional column for per-group correlations.
#'   When set, the output contains a \code{group} column and can be
#'   plotted as a faceted \code{PlotRankScatter}. Default: \code{NULL}.
#' @param method Correlation method: \code{"spearman"} (default),
#'   \code{"pearson"}, or \code{"kendall"}.
#' @param assay Character. Seurat assay to use (Seurat only).
#'   Default: \code{DefaultAssay(object)}.
#' @param layer Character. Data layer to extract (Seurat only).
#'   Default: \code{"data"}.
#' @param n_features Integer. Maximum number of variable features to use
#'   when \code{features = NULL} (Seurat only). Default: 2000.
#' @param min.pct Numeric (0–1). Minimum fraction of cells expressing a
#'   gene for it to be included (Seurat only). Default: 0.05.
#' @param min.cells Integer. Minimum number of cells per group for
#'   correlation to be computed. Default: 20.
#' @param cl Integer. Number of parallel workers for per-group
#'   computation (via \code{pbapply::pblapply}). Only used when
#'   \code{group.by} is set and \code{cl > 1}. Default: 1 (serial).
#' @param padjust.method Character. p-value adjustment method passed to
#'   \code{p.adjust()}. Default: \code{"fdr"}.
#' @param verbose Logical. Print progress messages. Default: \code{TRUE}.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{gene}{Feature name.}
#'     \item{score}{Correlation coefficient (Spearman rho, Pearson r, or
#'       Kendall tau).}
#'     \item{pvalue}{Raw p-value from \code{cor.test()}.}
#'     \item{padj}{Adjusted p-value.}
#'     \item{abs_score}{Absolute value of \code{score}, for ranking.}
#'     \item{n_cells}{Number of cells / observations used.}
#'     \item{group}{(Only when \code{group.by} is set) Group label.}
#'   }
#'   Rows are sorted by \code{abs_score} (descending).  The data.frame is
#'   directly compatible with \code{\link{PlotRankScatter}}: the
#'   \code{gene} column maps to the name column and \code{score} maps to
#'   the score column.
#'
#' @examples
#' \dontrun{
#' # --- Seurat: correlate PTH with all HVG ---
#' res <- RunCorrelation(seu, target = "PTH")
#' PlotRankScatter(res)
#'
#' # --- Seurat: per cell-type ---
#' res <- RunCorrelation(seu, target = "PTH", group.by = "celltype")
#' PlotRankScatter(res)
#'
#' # --- Seurat: specific genes ---
#' res <- RunCorrelation(seu, target = "PTH",
#'                features = c("GCM2", "CASR", "VDR", "CYP27B1"))
#'
#' # --- data.frame ---
#' df <- data.frame(x = rnorm(100), a = rnorm(100),
#'                  b = rnorm(100), c = rnorm(100))
#' res <- RunCorrelation(df, target = "x")
#' PlotRankScatter(res)
#' }
#'
#' @export
RunCorrelation <- function(object,
                    target,
                    features       = NULL,
                    group.by       = NULL,
                    method         = c("spearman", "pearson", "kendall"),
                    assay          = NULL,
                    layer          = "data",
                    n_features     = 2000L,
                    min.pct        = 0.05,
                    min.cells      = 20L,
                    cl             = 1L,
                    padjust.method = "fdr",
                    verbose        = TRUE) {

  method <- match.arg(method)

  is_seurat <- inherits(object, "Seurat")
  is_df     <- is.data.frame(object)

  if (!is_seurat && !is_df) {
    stop("'object' must be a Seurat object or a data.frame.", call. = FALSE)
  }

  # ── Build data matrix ──────────────────────────────────────────────────────
  if (is_df) {
    # --- data.frame path ---
    if (!target %in% colnames(object)) {
      stop("'target' column '", target, "' not found.", call. = FALSE)
    }
    if (!is.null(group.by) && !group.by %in% colnames(object)) {
      stop("'group.by' column '", group.by, "' not found.", call. = FALSE)
    }
    if (is.null(features)) {
      # All numeric columns except target (and group.by)
      num_cols <- names(object)[vapply(object, is.numeric, logical(1))]
      features <- setdiff(num_cols, c(target, group.by))
    }
    missing_f <- setdiff(features, colnames(object))
    if (length(missing_f) > 0) {
      warning("Features not found: ",
              paste(head(missing_f, 5), collapse = ", "),
              if (length(missing_f) > 5) "...", call. = FALSE)
      features <- intersect(features, colnames(object))
    }
    if (length(features) == 0) {
      stop("No valid features to compute correlation.", call. = FALSE)
    }
    dat <- object

  } else {
    # --- Seurat path ---
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
    meta_cols <- colnames(object@meta.data)
    all_genes <- rownames(object[[assay]])

    # Validate target
    target_in_gene <- target %in% all_genes
    target_in_meta <- target %in% meta_cols
    if (!target_in_gene && !target_in_meta) {
      stop("'target' ('", target, "') not found in assay or meta.data.",
           call. = FALSE)
    }

    # Resolve features
    if (is.null(features)) {
      features <- tryCatch(
        SeuratObject::VariableFeatures(object, assay = assay),
        error = function(e) NULL
      )
      if (is.null(features) || length(features) == 0) {
        features <- all_genes
      }
      features <- head(features, n_features)
      features <- setdiff(features, target)
      if (verbose) message("  Using ", length(features), " features.")
    } else {
      valid_f <- features %in% all_genes | features %in% meta_cols
      if (any(!valid_f)) {
        warning("Features not found: ",
                paste(head(features[!valid_f], 5), collapse = ", "),
                if (sum(!valid_f) > 5) "...", call. = FALSE)
      }
      features <- features[valid_f]
      features <- setdiff(features, target)
    }
    if (length(features) == 0) {
      stop("No valid features to compute correlation.", call. = FALSE)
    }

    # group.by validation
    if (!is.null(group.by) && !group.by %in% meta_cols) {
      stop("'group.by' column '", group.by,
           "' not found in meta.data.", call. = FALSE)
    }

    # Separate gene vs metadata variables
    all_vars <- unique(c(target, features, group.by))
    gene_vars <- intersect(all_vars, all_genes)
    meta_vars <- setdiff(all_vars, gene_vars)
    if (!is.null(group.by) && group.by %in% gene_vars) {
      gene_vars <- setdiff(gene_vars, group.by)
      meta_vars <- union(meta_vars, group.by)
    }

    cells <- colnames(object)

    # Extract gene expression via matrix slice (fast, no per-gene loop)
    if (length(gene_vars) > 0) {
      if (verbose) message("  Extracting ", length(gene_vars),
                           " gene expressions ...")
      expr_mat <- SeuratObject::GetAssayData(object, assay = assay,
                                             layer = layer)
      gene_df <- as.data.frame(
        as.matrix(Matrix::t(expr_mat[gene_vars, cells, drop = FALSE]))
      )
    } else {
      gene_df <- data.frame(row.names = cells)
    }

    # Add metadata columns
    if (length(meta_vars) > 0) {
      meta_df <- object@meta.data[cells, meta_vars, drop = FALSE]
      dat <- cbind(gene_df, meta_df)
    } else {
      dat <- gene_df
    }

    # Filter genes by min.pct (vectorized on sparse matrix)
    gene_features <- intersect(features, gene_vars)
    if (length(gene_features) > 0 && min.pct > 0) {
      n_total <- length(cells)
      pct <- Matrix::colSums(
        expr_mat[gene_features, cells, drop = FALSE] > 0
      ) / n_total
      drop <- names(pct)[pct < min.pct]
      if (length(drop) > 0 && verbose) {
        message("  Filtered ", length(drop),
                " genes below min.pct = ", min.pct)
      }
      features <- setdiff(features, drop)
    }

    if (length(features) == 0) {
      stop("No features passed the min.pct filter.", call. = FALSE)
    }
  }

  # ── Vectorized correlation computation ────────────────────────────────────
  # For Pearson / Spearman: t = r * sqrt(n-2) / sqrt(1 - r^2), p = 2*pt(...)
  # For Kendall: normal approx z = tau / sqrt(2*(2n+5) / (9*n*(n-1)))
  .cor_pvalue <- function(r, n, method) {
    if (method %in% c("pearson", "spearman")) {
      t_stat <- r * sqrt(n - 2) / sqrt(pmax(1 - r^2, .Machine$double.eps))
      2 * stats::pt(-abs(t_stat), df = n - 2)
    } else {
      # Kendall normal approximation (valid for n >= 10)
      sigma <- sqrt(2 * (2 * n + 5) / (9 * n * (n - 1)))
      z <- r / pmax(sigma, .Machine$double.eps)
      2 * stats::pnorm(-abs(z))
    }
  }

  .cor_one_group <- function(df_sub, grp_label = NULL) {
    target_vec <- df_sub[[target]]
    if (!is.numeric(target_vec)) {
      stop("'target' ('", target, "') must be numeric.", call. = FALSE)
    }

    # Filter to numeric features only
    feat_use <- features[vapply(features, function(f) {
      is.numeric(df_sub[[f]])
    }, logical(1))]
    if (length(feat_use) == 0) return(NULL)

    # Build matrix: target + features
    # Use complete cases for target, pairwise for features
    ok_target <- is.finite(target_vec)
    if (sum(ok_target) < min.cells) return(NULL)

    target_clean <- target_vec[ok_target]
    mat <- as.matrix(df_sub[ok_target, feat_use, drop = FALSE])
    n_obs <- nrow(mat)

    # Rank transform for Spearman
    if (method == "spearman") {
      target_clean <- rank(target_clean)
      mat <- apply(mat, 2, rank)
    }

    # Vectorized cor() — one call for all features
    rho <- as.numeric(stats::cor(target_clean, mat, use = "pairwise.complete.obs"))
    names(rho) <- feat_use

    # Per-feature n (pairwise complete observations)
    n_per <- colSums(is.finite(mat))

    # Filter: enough observations & valid correlation
    keep <- !is.na(rho) & n_per >= max(3L, min.cells)
    if (sum(keep) == 0) return(NULL)

    rho_k   <- rho[keep]
    n_k     <- n_per[keep]
    feat_k  <- feat_use[keep]

    pval <- .cor_pvalue(rho_k, n_k,
                        ifelse(method == "spearman", "pearson", method))

    res <- data.frame(
      gene    = feat_k,
      score   = rho_k,
      pvalue  = pval,
      n_cells = as.integer(n_k),
      stringsAsFactors = FALSE
    )
    if (!is.null(grp_label)) res$group <- grp_label
    res
  }

  # ── With / without grouping ────────────────────────────────────────────────
  if (is.null(group.by)) {
    if (verbose) message("  Computing correlation for ",
                         length(features), " features ...")
    result <- .cor_one_group(dat)
  } else {
    groups <- unique(dat[[group.by]])
    groups <- groups[!is.na(groups)]
    if (verbose) message("  Computing correlation for ",
                         length(features), " features x ",
                         length(groups), " groups ...")

    # Split data once (avoid repeated subsetting)
    split_dat <- split(dat, dat[[group.by]])
    split_dat <- split_dat[groups]

    .worker <- function(g) {
      .cor_one_group(split_dat[[g]], grp_label = as.character(g))
    }

    if (cl > 1L && length(groups) > 1L &&
        requireNamespace("pbapply", quietly = TRUE)) {
      # Parallel via pbapply
      if (verbose) message("  Using ", cl, " parallel workers ...")
      res_list <- pbapply::pblapply(groups, .worker, cl = cl)
    } else {
      # Serial
      res_list <- vector("list", length(groups))
      for (gi in seq_along(groups)) {
        g <- groups[gi]
        if (verbose) message("    [", gi, "/", length(groups), "] ", g,
                             " (n = ", nrow(split_dat[[g]]), ")")
        res_list[[gi]] <- .worker(g)
      }
    }
    result <- do.call(rbind, res_list)
  }

  if (is.null(result) || nrow(result) == 0) {
    warning("No valid correlation results.", call. = FALSE)
    return(data.frame(gene = character(), score = numeric(),
                      pvalue = numeric(), padj = numeric(),
                      abs_score = numeric(), n_cells = integer()))
  }

  # ── p-value adjustment & sorting ───────────────────────────────────────────
  result$padj      <- stats::p.adjust(result$pvalue, method = padjust.method)
  result$abs_score <- abs(result$score)
  result <- result[order(-result$abs_score), ]
  rownames(result) <- NULL

  if (verbose) {
    n_sig <- sum(result$padj < 0.05, na.rm = TRUE)
    message("  Done. ", nrow(result), " pairs; ",
            n_sig, " significant (padj < 0.05).")
  }

  # Store target info as attribute for downstream use
  attr(result, "target") <- target
  attr(result, "method") <- method

  result
}
