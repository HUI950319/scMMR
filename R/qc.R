# QC functions: doublet detection and ambient RNA removal

#' Mark doublets via DoubletFinder
#'
#' Detect doublets in a Seurat object using DoubletFinder.
#' Compatible with Seurat V5 and DoubletFinder >= 2.06.
#'
#' @param seu Seurat object
#' @param PCs Vectors indicating used principal components. Default: 1:10
#' @param split.by Name of a metadata column to split by (run DoubletFinder per group). Default: NULL
#' @param num.cores Threads for calculation. Default: 1
#' @return Seurat object with `DF.classifications` column added to metadata.
#'   Values: "Doublet.hc" (high confidence), "Doublet.lc" (low confidence), "Singlet".
#' @export
ComputeDoublets <- function(seu, PCs = 1:10, split.by = NULL, num.cores = 1) {

  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    stop("Package 'DoubletFinder' is required. Install it with remotes::install_github('chris-mcginnis-ucsf/DoubletFinder').")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    stop("Package 'pbapply' is required.")
  }

  PreprocessSeurat <- function(seu, PCs = 1:10) {
    message("  Normalize ...")
    seu <- Seurat::NormalizeData(seu, verbose = FALSE)
    message("  Find variable genes ...")
    seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    message("  Scale data ...")
    seu <- Seurat::ScaleData(seu, verbose = FALSE)
    message("  Run PCA ...")
    seu <- Seurat::RunPCA(seu, verbose = FALSE)
    message("  Find neighbors ...")
    k <- round(ncol(seu) * 0.02)
    k <- ifelse(k < 20, 20, k)
    seu <- Seurat::FindNeighbors(seu, dims = 1:10, reduction = "pca", k.param = k, verbose = FALSE)
    message("  Find clusters ...")
    seu <- Seurat::FindClusters(seu, resolution = 0.4, verbose = FALSE)
    seu
  }

  FindOptimalpK <- function(seu, PCs = 1:10, num.cores = 1) {
    sweep.res.list <- DoubletFinder::paramSweep(seu, PCs = PCs, sct = FALSE, num.cores = num.cores)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.stats)
    pK <- bcmvn[which.max(bcmvn$BCmetric), ]$pK
    as.numeric(as.character(pK))
  }

  DF <- function(seu, PCs = 1:10, auto.pK = TRUE, auto.cluster = TRUE, num.cores = 1) {
    message("1. Preprocessing ...")
    if (auto.cluster) {
      seu <- PreprocessSeurat(seu, PCs = PCs)
    }

    message("2. Find optimal pK ...")
    if (auto.pK) {
      optimal.pK <- FindOptimalpK(seu, PCs = PCs, num.cores = num.cores)
      message(paste0("   Optimal pK = ", optimal.pK))
    } else {
      optimal.pK <- 0.09
    }

    message("3. Mark doublets ...")
    homotypic.prop <- DoubletFinder::modelHomotypic(seu$seurat_clusters)
    nExp_poi <- round(0.075 * ncol(seu))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    pN <- 0.25
    pANN_name <- paste("pANN", pN, optimal.pK, nExp_poi, sep = "_")
    DF_name <- paste("DF.classifications", pN, optimal.pK, nExp_poi, sep = "_")
    DF_name.adj <- paste("DF.classifications", pN, optimal.pK, nExp_poi.adj, sep = "_")

    seu <- DoubletFinder::doubletFinder(seu, PCs = PCs, pN = pN, pK = optimal.pK,
                                        nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    seu <- DoubletFinder::doubletFinder(seu, PCs = PCs, pN = pN, pK = optimal.pK,
                                        nExp = nExp_poi.adj, reuse.pANN = pANN_name, sct = FALSE)

    seu$DF.classifications <- ifelse(seu@meta.data[[DF_name.adj]] == "Doublet", "Doublet.hc",
                                     ifelse(seu@meta.data[[DF_name]] == "Doublet", "Doublet.lc", "Singlet"))
    seu@meta.data[[pANN_name]] <- NULL
    seu@meta.data[[DF_name]] <- NULL
    seu@meta.data[[DF_name.adj]] <- NULL
    seu
  }

  if (is.null(split.by)) {
    seu.list <- list(seu)
  } else {
    seu.list <- Seurat::SplitObject(seu, split.by = split.by)
  }
  names(seu.list) <- NULL
  metadata.new <- pbapply::pblapply(seu.list, function(xx) {
    seu.tmp <- DF(xx, PCs = PCs, auto.pK = TRUE, auto.cluster = TRUE, num.cores = num.cores)
    seu.tmp@meta.data
  })
  metadata.new <- do.call(rbind, metadata.new)
  seu$DF.classifications <- metadata.new[rownames(seu@meta.data), ]$DF.classifications
  return(seu)
}


#' Remove ambient RNA contamination via decontX
#'
#' Estimate and remove ambient RNA contamination from a Seurat object
#' using celda::decontX. The corrected counts are stored in a new assay "decontX".
#'
#' @param seu Seurat object
#' @param split.by The grouping variable used to split the input object (batch correction). Default: NULL.
#' @param cluster.name Name of cluster field in Seurat object metadata. Must be pre-defined.
#' @return Seurat object with an additional "decontX" assay (corrected counts)
#'   and a `decontX_contamination` column in metadata.
#' @export
ComputeAmbientRNA <- function(seu, split.by = NULL, cluster.name = NULL) {

  if (!requireNamespace("sceasy", quietly = TRUE)) {
    stop("Package 'sceasy' is required. Install it with devtools::install_github('cellgeni/sceasy').")
  }
  if (!requireNamespace("celda", quietly = TRUE)) {
    stop("Package 'celda' is required.")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }

  if (is.null(cluster.name)) {
    stop("The clusters must be pre-defined. Please provide cluster.name.")
  }

  sce <- sceasy::convertFormat(seu, from = "seurat", to = "sce")
  if (is.null(split.by)) {
    sce <- celda::decontX(sce, z = sce[[cluster.name]])
  } else {
    sce <- celda::decontX(sce, z = sce[[cluster.name]],
                          batch = sce[[split.by]])
  }
  seu[["decontX"]] <- Seurat::CreateAssayObject(
    counts = round(celda::decontXcounts(sce), 0)
  )
  seu$decontX_contamination <- SingleCellExperiment::colData(sce)$decontX_contamination
  return(seu)
}
