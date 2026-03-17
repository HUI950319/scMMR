# ============================================================================
# Embedding quality evaluation — DNN embedding vs PCA comparison
# ============================================================================


# --------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------

#' Find the knee / elbow point using second-derivative method (internal)
#' @keywords internal
.find_knee <- function(var_explained) {
  n <- length(var_explained)
  if (n < 3) return(1L)
  d2 <- diff(diff(var_explained))
  # The knee is where the curvature changes most
  knee <- which.max(abs(d2)) + 1L
  knee
}


#' KNN Jaccard overlap between two spaces (internal)
#'
#' For each cell, compute Jaccard index of its k-nearest neighbors
#' in embedding vs reference space.
#' @keywords internal
.knn_jaccard <- function(emb, ref, k) {
  knn_emb <- FNN::get.knn(emb, k = k)$nn.index
  knn_ref <- FNN::get.knn(ref, k = k)$nn.index
  n <- nrow(emb)
  jaccards <- vapply(seq_len(n), function(i) {
    s1 <- knn_emb[i, ]
    s2 <- knn_ref[i, ]
    length(intersect(s1, s2)) / length(union(s1, s2))
  }, numeric(1))
  c(mean = mean(jaccards), sd = stats::sd(jaccards))
}


#' Distance correlation between two spaces (internal)
#'
#' Randomly sample cell pairs, compute pairwise Euclidean distance
#' in both spaces, and report Spearman rank correlation.
#' @keywords internal
.dist_correlation <- function(emb, ref, n_pairs = 5000L, seed = 42L) {
  set.seed(seed)
  n <- nrow(emb)
  i1 <- sample(n, n_pairs, replace = TRUE)
  i2 <- sample(n, n_pairs, replace = TRUE)
  # Avoid self-pairs
  same <- i1 == i2
  i2[same] <- ((i2[same]) %% n) + 1L

  d_emb <- sqrt(rowSums((emb[i1, , drop = FALSE] - emb[i2, , drop = FALSE])^2))
  d_ref <- sqrt(rowSums((ref[i1, , drop = FALSE] - ref[i2, , drop = FALSE])^2))

  ct <- stats::cor.test(d_emb, d_ref, method = "spearman")
  list(
    correlation = as.numeric(ct$estimate),
    p_value     = ct$p.value,
    method      = "spearman",
    d_emb       = d_emb,
    d_ref       = d_ref
  )
}


#' RV coefficient via trace trick (internal)
#'
#' Multivariate generalization of R-squared. Uses O(n*p*q) instead
#' of O(n^2) by computing through cross-product matrices.
#' @keywords internal
.rv_coefficient_fast <- function(X, Y) {
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  XtY <- crossprod(X, Y)   # p x q
  XtX <- crossprod(X)       # p x p
  YtY <- crossprod(Y)       # q x q
  num   <- sum(XtY^2)
  denom <- sqrt(sum(XtX^2) * sum(YtY^2))
  if (denom == 0) return(0)
  num / denom
}


#' Silhouette comparison between two spaces (internal)
#'
#' Compute mean silhouette score per cell type in both embedding
#' and reference spaces.  Falls back to a simple implementation
#' if the 'cluster' package is unavailable.
#' @keywords internal
.silhouette_comparison <- function(emb, ref, labels, max_cells = 5000L) {
  unique_labels <- unique(labels)
  if (length(unique_labels) < 2) return(NULL)

  # Subsample for performance
  if (length(labels) > max_cells) {
    # Stratified subsample
    idx <- unlist(lapply(unique_labels, function(lab) {
      ii <- which(labels == lab)
      n_take <- max(10L, round(max_cells * length(ii) / length(labels)))
      if (length(ii) <= n_take) ii else sample(ii, n_take)
    }))
    emb    <- emb[idx, , drop = FALSE]
    ref    <- ref[idx, , drop = FALSE]
    labels <- labels[idx]
  }

  label_int <- as.integer(factor(labels))

  if (requireNamespace("cluster", quietly = TRUE)) {
    sil_emb <- cluster::silhouette(label_int, stats::dist(emb))
    sil_ref <- cluster::silhouette(label_int, stats::dist(ref))

    # Per cell type
    lab_levels <- levels(factor(labels))
    per_type <- data.frame(
      cell_type = lab_levels,
      sil_emb   = vapply(seq_along(lab_levels), function(i) {
        mean(sil_emb[label_int == i, 3])
      }, numeric(1)),
      sil_ref   = vapply(seq_along(lab_levels), function(i) {
        mean(sil_ref[label_int == i, 3])
      }, numeric(1)),
      stringsAsFactors = FALSE
    )

    list(
      emb_mean = mean(sil_emb[, 3]),
      ref_mean = mean(sil_ref[, 3]),
      per_type = per_type
    )
  } else {
    message("  Package 'cluster' not available; skipping silhouette.")
    NULL
  }
}


# --------------------------------------------------------------------------
# Internal plot helpers
# --------------------------------------------------------------------------

#' @keywords internal
.plot_elbow <- function(pca_result, base_size = 12) {
  df <- data.frame(
    PC           = seq_along(pca_result$var_explained),
    VarExplained = pca_result$var_explained
  )
  knee <- pca_result$effective_dim

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC,
                                         y = .data$VarExplained)) +
    ggplot2::geom_line(color = "#2166AC", linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5, color = "#2166AC") +
    ggplot2::geom_vline(xintercept = knee, linetype = "dashed",
                         color = "red", linewidth = 0.5) +
    ggplot2::annotate("text", x = knee + 1,
                       y = max(df$VarExplained) * 0.9,
                       label = paste0("Knee = ", knee),
                       hjust = 0, color = "red", size = 4) +
    ggplot2::labs(
      x     = "PC of DNN Embedding",
      y     = "Proportion of Variance Explained",
      title = "Intrinsic Dimensionality (ElbowPlot)"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  p
}


#' @keywords internal
.plot_cumvar <- function(pca_result, base_size = 12) {
  df <- data.frame(
    PC     = seq_along(pca_result$cum_var),
    CumVar = pca_result$cum_var
  )
  knee <- pca_result$effective_dim

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC,
                                         y = .data$CumVar)) +
    ggplot2::geom_line(color = "#B2182B", linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5, color = "#B2182B") +
    ggplot2::geom_hline(yintercept = c(0.8, 0.9, 0.95),
                         linetype = "dotted", color = "grey50") +
    ggplot2::geom_vline(xintercept = knee, linetype = "dashed",
                         color = "red", linewidth = 0.5) +
    ggplot2::annotate("text", x = knee + 1, y = 0.5,
                       label = paste0("Knee = ", knee),
                       hjust = 0, color = "red", size = 4) +
    ggplot2::labs(
      x     = "PC of DNN Embedding",
      y     = "Cumulative Variance Explained",
      title = "Cumulative Variance"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  p
}


#' @keywords internal
.plot_knn_overlap <- function(knn_df, base_size = 12) {
  p <- ggplot2::ggplot(knn_df, ggplot2::aes(x = .data$k,
                                              y = .data$mean_jaccard)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = pmax(0, .data$mean_jaccard - .data$sd_jaccard),
        ymax = pmin(1, .data$mean_jaccard + .data$sd_jaccard)
      ),
      fill = "#2166AC", alpha = 0.2
    ) +
    ggplot2::geom_line(color = "#2166AC", linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5, color = "#2166AC") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$mean_jaccard)),
      vjust = -1.2, size = 3.5, color = "#2166AC"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x     = "k (Number of Neighbors)",
      y     = "Mean Jaccard Index",
      title = "KNN Neighborhood Overlap (DNN vs PCA)"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  p
}


#' @keywords internal
.plot_distcor <- function(dist_result, base_size = 12) {
  n_plot <- min(length(dist_result$d_emb), 2000L)
  idx <- sample(length(dist_result$d_emb), n_plot)
  df <- data.frame(
    d_ref = dist_result$d_ref[idx],
    d_emb = dist_result$d_emb[idx]
  )
  rho_label <- sprintf("rho == %.3f", dist_result$correlation)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$d_ref,
                                         y = .data$d_emb)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.3, color = "grey30") +
    ggplot2::geom_smooth(method = "lm", se = FALSE,
                          color = "#B2182B", linewidth = 0.8) +
    ggplot2::annotate("text", x = Inf, y = Inf,
                       label = rho_label, parse = TRUE,
                       hjust = 1.1, vjust = 1.5, size = 5,
                       color = "#B2182B") +
    ggplot2::labs(
      x     = "Pairwise Distance (PCA)",
      y     = "Pairwise Distance (DNN Embedding)",
      title = "Distance Correlation"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  p
}


#' @keywords internal
.plot_silhouette <- function(sil_result, base_size = 12) {
  if (is.null(sil_result)) return(NULL)
  df <- data.frame(
    cell_type  = rep(sil_result$per_type$cell_type, 2),
    silhouette = c(sil_result$per_type$sil_emb,
                   sil_result$per_type$sil_ref),
    space      = rep(c("DNN Embedding", "PCA"),
                     each = nrow(sil_result$per_type)),
    stringsAsFactors = FALSE
  )

  # Order by DNN silhouette
  emb_order <- sil_result$per_type$cell_type[
    order(sil_result$per_type$sil_emb)
  ]
  df$cell_type <- factor(df$cell_type, levels = emb_order)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x    = .data$cell_type,
    y    = .data$silhouette,
    fill = .data$space
  )) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::scale_fill_manual(values = c("DNN Embedding" = "#2166AC",
                                           "PCA"           = "#B2182B")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x     = NULL,
      y     = "Mean Silhouette Width",
      title = "Silhouette Comparison by Cell Type",
      fill  = "Space"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text      = ggplot2::element_text(color = "black"),
      plot.title     = ggplot2::element_text(hjust = 0.5),
      legend.position = "top"
    )
  p
}


# ============================================================================
# EvaluateEmbedding
# ============================================================================

#' Evaluate DNN Embedding Quality
#'
#' Assess the information content of the DNN shared embedding and compare
#' it with PCA space. Returns intrinsic dimensionality, KNN overlap,
#' distance correlation, RV coefficient, and silhouette comparison.
#'
#' @param embedding Numeric matrix (n_cells x d), typically
#'   \code{pred$shared_embedding} from
#'   \code{DNN_predict(return_embedding = TRUE)}.
#'   Rownames must be cell IDs.
#' @param seurat_obj A Seurat object (optional). If provided, PCA
#'   embeddings and cell type labels are extracted automatically.
#'   At least one of \code{seurat_obj} or \code{pca_ref} is needed
#'   for consistency metrics.
#' @param pca_ref Numeric matrix (n_cells x n_pcs), alternative to
#'   \code{seurat_obj}. Rownames must be cell IDs.
#' @param cell_type_col Character. Column name in Seurat metadata or
#'   a named character vector of cell type labels.
#'   Used for silhouette comparison. Default \code{"cell_type"}.
#' @param n_pcs Integer. Number of PCs to compute on the DNN embedding
#'   for intrinsic dimensionality assessment (default 50).
#' @param k_values Integer vector. Values of k for KNN overlap curve
#'   (default \code{c(10, 20, 50, 100, 200)}).
#' @param n_pairs Integer. Number of cell pairs for distance correlation
#'   (default 5000).
#' @param max_cells Integer. Downsample to this many cells for
#'   performance (default 20000).
#' @param seed Integer. Random seed (default 42).
#' @param verbose Logical. Print progress messages (default TRUE).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{pca_of_embedding}{List: sdev, var_explained, cum_var,
#'       effective_dim}
#'     \item{consistency}{List: knn_overlap (data.frame), dist_cor (list),
#'       rv_coef (numeric), silhouette (list). NULL if no PCA reference.}
#'     \item{summary}{Character string summarising key findings.}
#'     \item{params}{Parameters used.}
#'   }
#'
#' @examples
#' \dontrun{
#' pred <- DNN_predict(query, model_path, return_embedding = TRUE)
#' eval_res <- EvaluateEmbedding(pred$shared_embedding,
#'                                seurat_obj = seu,
#'                                cell_type_col = "cell_type")
#' print(eval_res$summary)
#' PlotEmbeddingEval(eval_res)
#' PlotEmbeddingEval(eval_res, which = "elbow")
#' }
#'
#' @export
EvaluateEmbedding <- function(embedding,
                               seurat_obj    = NULL,
                               pca_ref       = NULL,
                               cell_type_col = "cell_type",
                               n_pcs         = 50L,
                               k_values      = c(10L, 20L, 50L, 100L, 200L),
                               n_pairs       = 5000L,
                               max_cells     = 20000L,
                               seed          = 42L,
                               verbose       = TRUE) {

  set.seed(seed)

  # -- Validate embedding -----------------------------------------------------
  if (!is.matrix(embedding) || !is.numeric(embedding))
    stop("'embedding' must be a numeric matrix (e.g. pred$shared_embedding).")
  if (is.null(rownames(embedding)))
    stop("'embedding' must have rownames (cell IDs).")

  # -- Extract PCA reference & labels -----------------------------------------
  labels <- NULL

  # cell_type_col can be a column name (string) or a named character vector
  if (is.character(cell_type_col) && length(cell_type_col) > 1) {
    # Named vector of labels passed directly
    labels <- cell_type_col
  }

  if (!is.null(seurat_obj)) {
    if (!requireNamespace("SeuratObject", quietly = TRUE))
      stop("Package 'SeuratObject' is required when providing seurat_obj.")

    # PCA
    if (is.null(pca_ref)) {
      pca_red <- SeuratObject::Reductions(seurat_obj, "pca")
      if (is.null(pca_red))
        stop("No PCA reduction found in seurat_obj. ",
             "Run Seurat::RunPCA() first or provide 'pca_ref'.")
      pca_ref <- SeuratObject::Embeddings(pca_red)
    }

    # Cell type labels (from Seurat metadata, if not already set)
    if (is.null(labels) && is.character(cell_type_col) &&
        length(cell_type_col) == 1) {
      if (cell_type_col %in% colnames(seurat_obj@meta.data)) {
        labels <- as.character(seurat_obj@meta.data[[cell_type_col]])
        names(labels) <- rownames(seurat_obj@meta.data)
      }
    }
  }

  has_ref <- !is.null(pca_ref)

  # -- Intersect cell IDs -----------------------------------------------------
  if (has_ref) {
    if (is.null(rownames(pca_ref)))
      stop("'pca_ref' must have rownames (cell IDs).")
    common <- intersect(rownames(embedding), rownames(pca_ref))
    if (length(common) == 0)
      stop("No common cell IDs between embedding and pca_ref.")
    if (verbose && length(common) < nrow(embedding))
      message(sprintf("  Using %d / %d cells (intersection).",
                      length(common), nrow(embedding)))
    embedding <- embedding[common, , drop = FALSE]
    pca_ref   <- pca_ref[common, , drop = FALSE]
    if (!is.null(labels)) labels <- labels[common]
  }

  # -- Downsample for performance ----------------------------------------------
  n_cells <- nrow(embedding)
  if (n_cells > max_cells) {
    if (verbose)
      message(sprintf("  Downsampling: %d -> %d cells", n_cells, max_cells))
    idx <- sample(n_cells, max_cells)
    embedding <- embedding[idx, , drop = FALSE]
    if (has_ref) pca_ref <- pca_ref[idx, , drop = FALSE]
    if (!is.null(labels)) labels <- labels[idx]
    n_cells <- max_cells
  }

  if (verbose) message(sprintf("  Evaluating embedding: %d cells x %d dims",
                                n_cells, ncol(embedding)))

  # ===== A. Intrinsic dimensionality ==========================================
  if (verbose) message("  [1/5] Computing intrinsic dimensionality (PCA on embedding) ...")

  n_pcs_actual <- min(n_pcs, ncol(embedding), n_cells - 1L)
  # Total variance = sum of column variances (all 512 dims)
  emb_centered <- scale(embedding, center = TRUE, scale = FALSE)
  total_var <- sum(colSums(emb_centered^2)) / (n_cells - 1)
  pca_emb <- stats::prcomp(embedding, center = TRUE, scale. = FALSE,
                            rank. = n_pcs_actual)
  var_explained <- pca_emb$sdev^2 / total_var
  cum_var       <- cumsum(var_explained)
  effective_dim <- .find_knee(var_explained)

  pca_of_embedding <- list(
    sdev          = pca_emb$sdev,
    var_explained = var_explained,
    cum_var       = cum_var,
    effective_dim = effective_dim
  )

  if (verbose)
    message(sprintf("    Effective dimensions: ~%d (knee), ",
                    effective_dim),
            sprintf("90%% variance at PC %d, 95%% at PC %d",
                    which(cum_var >= 0.90)[1],
                    which(cum_var >= 0.95)[1]))

  # ===== B. Consistency metrics ===============================================
  consistency <- NULL
  if (has_ref) {
    if (!requireNamespace("FNN", quietly = TRUE))
      stop("Package 'FNN' is required for consistency metrics.\n",
           "Install with: install.packages('FNN')")

    # -- B1. KNN overlap -------------------------------------------------------
    if (verbose) message("  [2/5] Computing KNN overlap ...")
    # Filter k_values to be valid
    k_values <- k_values[k_values < n_cells]
    knn_results <- lapply(k_values, function(k) {
      if (verbose) message(sprintf("    k = %d ...", k))
      jac <- .knn_jaccard(embedding, pca_ref, k)
      data.frame(k = k, mean_jaccard = jac["mean"],
                 sd_jaccard = jac["sd"],
                 stringsAsFactors = FALSE)
    })
    knn_overlap <- do.call(rbind, knn_results)
    rownames(knn_overlap) <- NULL

    # -- B2. Distance correlation ----------------------------------------------
    if (verbose) message("  [3/5] Computing distance correlation ...")
    dist_cor <- .dist_correlation(embedding, pca_ref, n_pairs, seed)
    if (verbose)
      message(sprintf("    Spearman rho = %.3f (p = %.2e)",
                      dist_cor$correlation, dist_cor$p_value))

    # -- B3. RV coefficient ----------------------------------------------------
    if (verbose) message("  [4/5] Computing RV coefficient ...")
    rv_coef <- .rv_coefficient_fast(embedding, pca_ref)
    if (verbose)
      message(sprintf("    RV coefficient = %.3f", rv_coef))

    # -- B4. Silhouette --------------------------------------------------------
    sil_result <- NULL
    if (!is.null(labels)) {
      if (verbose) message("  [5/5] Computing silhouette comparison ...")
      sil_result <- .silhouette_comparison(embedding, pca_ref, labels)
      if (!is.null(sil_result) && verbose)
        message(sprintf("    Silhouette: DNN = %.3f, PCA = %.3f",
                        sil_result$emb_mean, sil_result$ref_mean))
    } else {
      if (verbose) message("  [5/5] Skipping silhouette (no cell type labels)")
    }

    consistency <- list(
      knn_overlap = knn_overlap,
      dist_cor    = dist_cor,
      rv_coef     = rv_coef,
      silhouette  = sil_result
    )
  } else {
    if (verbose) message("  [2-5/5] Skipping consistency metrics (no PCA reference)")
  }

  # ===== Summary string =======================================================
  summary_parts <- sprintf(
    "DNN Embedding: %d cells x %d dims, effective dim ~%d (knee)",
    n_cells, ncol(embedding), effective_dim
  )
  if (has_ref) {
    summary_parts <- c(summary_parts,
      sprintf("KNN overlap (k=%d): Jaccard = %.3f",
              knn_overlap$k[which.min(abs(knn_overlap$k - 30))],
              knn_overlap$mean_jaccard[which.min(abs(knn_overlap$k - 30))]),
      sprintf("Distance correlation: rho = %.3f", dist_cor$correlation),
      sprintf("RV coefficient: %.3f", rv_coef)
    )
    if (!is.null(sil_result)) {
      summary_parts <- c(summary_parts,
        sprintf("Silhouette: DNN = %.3f, PCA = %.3f (%s)",
                sil_result$emb_mean, sil_result$ref_mean,
                if (sil_result$emb_mean > sil_result$ref_mean)
                  "DNN better" else "PCA better"))
    }
  }
  summary_str <- paste(summary_parts, collapse = "\n")
  if (verbose) message("\n  === Summary ===\n  ",
                        gsub("\n", "\n  ", summary_str))

  result <- list(
    pca_of_embedding = pca_of_embedding,
    consistency      = consistency,
    summary          = summary_str,
    params           = list(
      n_pcs     = n_pcs_actual,
      k_values  = k_values,
      n_pairs   = n_pairs,
      max_cells = max_cells,
      seed      = seed,
      n_cells   = n_cells,
      emb_dims  = ncol(embedding)
    )
  )

  if (verbose) message("  Done.")
  result
}


# ============================================================================
# PlotEmbeddingEval
# ============================================================================

#' Plot Embedding Evaluation Results
#'
#' Visualise results from \code{\link{EvaluateEmbedding}}.
#'
#' @param eval_result Output from \code{EvaluateEmbedding()}.
#' @param which Character. Which plot(s) to draw:
#'   \code{"all"} (default, returns a list of all available plots),
#'   \code{"elbow"}, \code{"cumvar"}, \code{"knn"}, \code{"distcor"},
#'   or \code{"silhouette"}.
#' @param base_size Numeric. Base font size (default 12).
#'
#' @return A single ggplot object or a named list of ggplot objects
#'   (when \code{which = "all"}).
#'
#' @examples
#' \dontrun{
#' eval_res <- EvaluateEmbedding(pred$shared_embedding, seurat_obj = seu)
#'
#' # All plots as a list
#' plots <- PlotEmbeddingEval(eval_res)
#' plots$elbow
#' plots$knn
#'
#' # Single plot
#' PlotEmbeddingEval(eval_res, which = "silhouette")
#'
#' # Combine with patchwork
#' library(patchwork)
#' plots$elbow + plots$cumvar + plots$knn + plots$distcor
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_col geom_ribbon
#'   geom_smooth geom_text geom_hline geom_vline annotate coord_flip
#'   scale_y_continuous scale_fill_manual labs theme_classic theme
#'   element_text
#' @importFrom rlang .data
#' @export
PlotEmbeddingEval <- function(eval_result,
                               which     = "all",
                               base_size = 12) {

  which <- match.arg(which, c("all", "elbow", "cumvar",
                               "knn", "distcor", "silhouette"))

  plots <- list()

  # Always available: intrinsic dimensionality
  if (which %in% c("all", "elbow")) {
    plots$elbow <- .plot_elbow(eval_result$pca_of_embedding, base_size)
  }
  if (which %in% c("all", "cumvar")) {
    plots$cumvar <- .plot_cumvar(eval_result$pca_of_embedding, base_size)
  }

  # Consistency metrics (only if available)
  cons <- eval_result$consistency
  if (!is.null(cons)) {
    if (which %in% c("all", "knn")) {
      plots$knn <- .plot_knn_overlap(cons$knn_overlap, base_size)
    }
    if (which %in% c("all", "distcor")) {
      plots$distcor <- .plot_distcor(cons$dist_cor, base_size)
    }
    if (which %in% c("all", "silhouette") && !is.null(cons$silhouette)) {
      plots$silhouette <- .plot_silhouette(cons$silhouette, base_size)
    }
  }

  if (length(plots) == 0)
    stop("No plots available for which = '", which, "'.")
  if (which != "all") return(plots[[1]])
  plots
}
