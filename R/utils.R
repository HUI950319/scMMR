# ══════════════════════════════════════════════════════════════════════════════
# utils.R — Internal helper functions
# ══════════════════════════════════════════════════════════════════════════════

#' Resolve input to AnnData (internal)
#'
#' Accepts either an h5ad file path or a Seurat object,
#' returns a Python AnnData object.
#'
#' @param input Character (h5ad path) or Seurat object.
#' @param embedding_key Embedding key (used for Seurat objects).
#' @return Python AnnData object.
#' @keywords internal
.resolve_input <- function(input, embedding_key = "umap") {
  if (is.character(input)) {
    # ── h5ad file path ─────────────────────────────────────────────────────
    input <- normalizePath(input, mustWork = TRUE)
    message("  Loading from h5ad: ", input)
    adata <- .scMMR_env$mt_load_query(input)
  } else if (inherits(input, "Seurat")) {
    # ── Seurat object → AnnData ────────────────────────────────────────────
    message("  Converting Seurat object to AnnData ...")
    adata <- .seurat_to_adata(input, embedding_key = embedding_key)
  } else {
    stop("'input' must be a character (h5ad path) or a Seurat object.")
  }
  n_obs  <- reticulate::py_to_r(adata$n_obs)
  n_vars <- reticulate::py_to_r(adata$n_vars)
  message(sprintf("  Data: %d cells x %d genes", n_obs, n_vars))
  adata
}


#' Convert Seurat object to Python AnnData (internal)
#'
#' Extracts counts matrix, metadata, and embeddings from a Seurat object
#' and passes them to the Python helper.
#'
#' @param srt Seurat object.
#' @param embedding_key Which reduction to include (default "umap").
#' @return Python AnnData object.
#' @keywords internal
.seurat_to_adata <- function(srt, embedding_key = "umap") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required to convert Seurat objects.\n",
         "Install with: install.packages('SeuratObject')")
  }

  # Extract counts (genes x cells, dgCMatrix)
  # Auto-detect assay: use default assay (handles "RNA", "originalexp", etc.)
  assay_name <- SeuratObject::DefaultAssay(srt)
  counts <- tryCatch(
    SeuratObject::GetAssayData(srt, assay = assay_name, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(counts) || ncol(counts) == 0) {
    counts <- srt[[assay_name]]@counts
  }

  # Metadata, names
  obs_df    <- srt@meta.data
  var_names <- rownames(counts)
  obs_names <- colnames(counts)

  # Embeddings
  obsm_dict <- list()
  for (red_name in names(srt@reductions)) {
    emb <- srt@reductions[[red_name]]@cell.embeddings
    obsm_dict[[red_name]] <- emb
  }

  # Convert for reticulate
  counts_py <- reticulate::r_to_py(counts)
  obs_py    <- reticulate::r_to_py(obs_df)
  obsm_py   <- if (length(obsm_dict) > 0) reticulate::r_to_py(obsm_dict) else NULL

  .scMMR_env$mt_srt_to_adata(
    counts_matrix = counts_py,
    obs_df        = obs_py,
    var_names     = var_names,
    obs_names     = obs_names,
    obsm_dict     = obsm_py
  )
}


#' Select highly variable genes with batch fallback (internal)
#'
#' Wraps \code{mt_select_hvg} with a \code{tryCatch} that falls back
#' to non-batch mode when the data is too small for batch-aware selection.
#'
#' @param adata     Python AnnData object.
#' @param n_top     Number of HVGs (default 6000).
#' @param batch_key obs column for batch-aware selection (NULL = no batch).
#' @return Character vector of HVG names.
#' @keywords internal
.select_hvg <- function(adata, n_top = 6000L, batch_key = NULL) {
  hvg <- tryCatch({
    .scMMR_env$mt_select_hvg(adata,
                              n_top_genes = as.integer(n_top),
                              batch_key   = batch_key)
  }, error = function(e) {
    if (!is.null(batch_key) &&
        grepl("reciprocal condition number", e$message)) {
      warning("Batch-aware HVG failed (data too small?), ",
              "retrying without batch_key ...", call. = FALSE)
      .scMMR_env$mt_select_hvg(adata,
                                n_top_genes = as.integer(n_top),
                                batch_key   = NULL)
    } else {
      stop(e)
    }
  })
  reticulate::py_to_r(hvg)
}


#' Prepare bulk expression input for deconvolution (internal)
#'
#' Handles multiple input formats (matrix, data.frame, CSV path, h5ad path,
#' Seurat object) and returns a standardised list.
#'
#' @param bulk_expr Input bulk expression data.
#' @return List with: \code{matrix} (genes x samples, numeric),
#'   \code{sample_names} (character), \code{gene_names} (character).
#' @keywords internal
.prepare_bulk_input <- function(bulk_expr) {

  if (is.character(bulk_expr) && length(bulk_expr) == 1 &&
      grepl("\\.csv$", bulk_expr, ignore.case = TRUE)) {
    # ── CSV file: first column = gene names, rest = samples ──────────────
    bulk_expr <- normalizePath(bulk_expr, mustWork = TRUE)
    df <- utils::read.csv(bulk_expr, row.names = 1, check.names = FALSE)
    list(matrix       = as.matrix(df),
         sample_names = colnames(df),
         gene_names   = rownames(df))

  } else if (is.character(bulk_expr) && length(bulk_expr) == 1 &&
             grepl("\\.h5ad$", bulk_expr, ignore.case = TRUE)) {
    # ── h5ad file ─────────────────────────────────────────────────────────
    bulk_expr <- normalizePath(bulk_expr, mustWork = TRUE)
    adata <- .scMMR_env$mt_load_query(bulk_expr)
    X <- reticulate::py_to_r(adata$X)
    if (inherits(X, "dgCMatrix") || inherits(X, "dgRMatrix")) {
      X <- as.matrix(X)
    }
    # AnnData: cells/samples x genes → transpose to genes x samples
    list(matrix       = t(X),
         sample_names = as.character(reticulate::py_to_r(adata$obs_names)),
         gene_names   = as.character(reticulate::py_to_r(adata$var_names)))

  } else if (inherits(bulk_expr, "Seurat")) {
    # ── Seurat object: use RNA assay counts ──────────────────────────────
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Package 'SeuratObject' is required for Seurat input.\n",
           "Install with: install.packages('SeuratObject')")
    }
    counts <- SeuratObject::GetAssayData(bulk_expr, assay = "RNA",
                                          layer = "counts")
    if (is.null(counts) || ncol(counts) == 0) {
      counts <- bulk_expr[["RNA"]]@counts
    }
    list(matrix       = as.matrix(counts),
         sample_names = colnames(counts),
         gene_names   = rownames(counts))

  } else if (is.matrix(bulk_expr) || is.data.frame(bulk_expr)) {
    # ── Matrix / data.frame: assume genes x samples ──────────────────────
    mat <- as.matrix(bulk_expr)
    sn  <- colnames(mat)
    if (is.null(sn)) sn <- paste0("sample_", seq_len(ncol(mat)))
    gn  <- rownames(mat)
    if (is.null(gn)) stop("bulk_expr must have rownames (gene names).")
    list(matrix       = mat,
         sample_names = sn,
         gene_names   = gn)

  } else {
    stop("'bulk_expr' must be a matrix, data.frame, CSV path, ",
         "h5ad path, or Seurat object.")
  }
}
