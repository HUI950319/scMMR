#' Standard single-cell analysis pipeline
#'
#' Run a standard Seurat workflow: Normalize, FindVariableFeatures, ScaleData,
#' PCA, (optional Harmony batch correction), UMAP, FindNeighbors, FindClusters.
#'
#' @param seu Seurat object
#' @param nfeatures Number of highly variable features. Default: 2000
#' @param npcs Number of PCA dimensions to compute. Default: 50
#' @param dims Dimensions to use for downstream steps (UMAP, neighbors). Default: 1:30
#' @param use.harmony Whether to run Harmony batch correction. Default: FALSE
#' @param harmony.group.by Character vector of metadata columns for Harmony grouping.
#'   Required when \code{use.harmony = TRUE}. Default: NULL
#' @param harmony.dims Dimensions passed to \code{harmony::RunHarmony()}. Default: 1:50
#' @param resolution Clustering resolution for \code{FindClusters}. Default: 0.6
#' @param k.param Number of nearest neighbors for \code{FindNeighbors}. Default: 20
#' @param vars.to.regress Variables to regress out during \code{ScaleData}. Default: NULL
#' @param verbose Whether to print Seurat messages. Default: FALSE
#' @param seed Random seed. Default: 42
#'
#' @return Seurat object with PCA, (Harmony), UMAP reductions and cluster assignments.
#' @export
StandardPipeline <- function(seu,
                             nfeatures = 2000,
                             npcs = 50,
                             dims = 1:30,
                             use.harmony = FALSE,
                             harmony.group.by = NULL,
                             harmony.dims = 1:50,
                             resolution = 0.6,
                             k.param = 20,
                             vars.to.regress = NULL,
                             verbose = FALSE,
                             seed = 42) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }

  set.seed(seed)

  # 1. Normalize
  message(">>> NormalizeData ...")
  seu <- Seurat::NormalizeData(seu, verbose = verbose)

  # 2. FindVariableFeatures
  message(">>> FindVariableFeatures (n = ", nfeatures, ") ...")
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = nfeatures, verbose = verbose)

  # 3. ScaleData
  message(">>> ScaleData ...")
  seu <- Seurat::ScaleData(seu, vars.to.regress = vars.to.regress, verbose = verbose)

  # 4. PCA
  message(">>> RunPCA (npcs = ", npcs, ") ...")
  seu <- Seurat::RunPCA(seu, npcs = npcs, verbose = verbose)

  # 5. Harmony (optional)
  reduction.use <- "pca"
  if (use.harmony) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("Package 'harmony' is required. Install it with install.packages('harmony').")
    }
    if (is.null(harmony.group.by)) {
      stop("harmony.group.by must be specified when use.harmony = TRUE.")
    }
    message(">>> RunHarmony (group.by = ", paste(harmony.group.by, collapse = ", "), ") ...")
    seu <- harmony::RunHarmony(seu,
                               group.by.vars = harmony.group.by,
                               reduction.use = "pca",
                               dims.use = harmony.dims,
                               verbose = verbose)
    reduction.use <- "harmony"
  }

  # 6. UMAP
  message(">>> RunUMAP (dims = ", min(dims), ":", max(dims), ") ...")
  seu <- Seurat::RunUMAP(seu, reduction = reduction.use, dims = dims, verbose = verbose)

  # 7. FindNeighbors
  message(">>> FindNeighbors (k = ", k.param, ") ...")
  seu <- Seurat::FindNeighbors(seu, reduction = reduction.use, dims = dims,
                               k.param = k.param, verbose = verbose)

  # 8. FindClusters
  message(">>> FindClusters (resolution = ", resolution, ") ...")
  seu <- Seurat::FindClusters(seu, resolution = resolution, verbose = verbose)

  message(">>> Done. Clusters stored in 'seurat_clusters'.")
  return(seu)
}
