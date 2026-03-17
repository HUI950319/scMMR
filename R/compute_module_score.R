#' Read a GMT (Gene Matrix Transposed) file
#'
#' Parse a GMT (Gene Matrix Transposed) file and return gene sets as a named
#' list. Each element contains the gene symbols belonging to that pathway.
#'
#' @param path Path to the GMT file.
#' @return A named list of character vectors (pathway name → gene symbols).
#' @examples
#' \dontrun{
#' gmt <- read_gmt("h.all.v2022.1.Hs.symbols.gmt")
#' names(gmt)[1:5]
#' }
#' @export
read_gmt <- function(path) {
  lines <- readLines(path)
  gene_sets <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      name <- parts[1]
      genes <- parts[-(1:2)]
      genes <- genes[nzchar(genes)]
      gene_sets[[name]] <- genes
    }
  }
  gene_sets
}


#' Parse gene sets from various input formats
#'
#' Converts gene sets from multiple input formats into a standardized named
#' list. Supports named list (pass-through), data.frame (term-to-gene mapping),
#' or a file path to a GMT file.
#'
#' @param gene.sets Gene sets in one of three formats:
#'   \itemize{
#'     \item A named list of character vectors (returned as-is).
#'     \item A data.frame with at least two columns: term (1st) and gene (2nd).
#'     \item A file path to a GMT file (parsed via \code{\link{read_gmt}}).
#'   }
#' @return A named list of character vectors (pathway name → gene symbols).
#' @examples
#' \dontrun{
#' # From GMT file
#' gs <- parse_gene_sets("hallmark.gmt")
#'
#' # From named list
#' gs <- parse_gene_sets(list(PathA = c("TP53", "MDM2"), PathB = c("EGFR", "ERBB2")))
#'
#' # From data.frame
#' df <- data.frame(term = c("PathA", "PathA"), gene = c("TP53", "MDM2"))
#' gs <- parse_gene_sets(df)
#' }
#' @seealso \code{\link{read_gmt}}, \code{\link{ComputeModuleScore}},
#'   \code{\link{RunPathwayAnalysis}}
#' @export
parse_gene_sets <- function(gene.sets) {
  # GMT file path
  if (is.character(gene.sets) && length(gene.sets) == 1) {
    if (!file.exists(gene.sets)) {
      stop("File not found: ", gene.sets)
    }
    return(read_gmt(gene.sets))
  }
  # data.frame: term-to-gene mapping (first col = term, second col = gene)
  if (is.data.frame(gene.sets)) {
    if (ncol(gene.sets) < 2) {
      stop("data.frame 'gene.sets' must have at least 2 columns (term, gene).")
    }
    gs <- split(as.character(gene.sets[[2]]), as.character(gene.sets[[1]]))
    return(lapply(gs, unique))
  }
  # named list
  if (is.list(gene.sets)) {
    if (is.null(names(gene.sets))) {
      stop("'gene.sets' list must be named.")
    }
    return(gene.sets)
  }
  stop("'gene.sets' must be a named list, a data.frame (term, gene), or a GMT file path.")
}


#' Calculate gene module scores
#'
#' Compute gene module activity scores using AUCell, Seurat (AddModuleScore), or UCell methods.
#' Supports multiple input formats for gene sets: named list, data.frame, or GMT file.
#'
#' @param x A gene expression matrix (genes x cells) or a Seurat object.
#' @param ... Additional arguments passed to scoring methods.
#' @return A score matrix (gene sets x cells) or a Seurat object with scores stored as an assay.
#'
#' @export
ComputeModuleScore <- function(x, ...) UseMethod("ComputeModuleScore")


#' @param gene.sets Gene sets in one of three formats:
#'   \itemize{
#'     \item A named list of character vectors (gene names).
#'     \item A data.frame with at least two columns: term (gene set name) and gene (gene symbol).
#'     \item A file path to a GMT file.
#'   }
#' @param method Scoring method: \code{"AUCell"} (default), \code{"Seurat"}, or \code{"UCell"}.
#' @param min.size Minimum number of genes (after filtering) for a gene set to be scored. Default: 20.
#' @param batch.size Number of cells per batch for AUCell method to reduce memory. Default: 500.
#' @param nbin Number of expression bins for Seurat method control gene selection. Default: 24.
#' @param ctrl Number of control genes per feature gene for Seurat method. Default: 100.
#' @param cores Number of parallel cores. Default: 1.
#' @param seed Random seed for reproducibility. Default: 11.
#'
#' @rdname ComputeModuleScore
#' @export
ComputeModuleScore.default <- function(x, gene.sets,
                                       method = c("AUCell", "Seurat", "UCell"),
                                       min.size = 20, batch.size = 500,
                                       nbin = 24, ctrl = 100,
                                       cores = 1, seed = 11, ...) {
  set.seed(seed)
  method <- match.arg(method)

  # Parse gene sets from list / data.frame / GMT
  gene.sets <- parse_gene_sets(gene.sets)

  # Filter: keep only genes present in the matrix, then apply min.size
  gene.sets <- lapply(gene.sets, function(g) intersect(g, rownames(x)))
  gene.sets <- gene.sets[vapply(gene.sets, length, integer(1)) >= min.size]

  if (length(gene.sets) == 0) {
    stop("No gene sets have >= ", min.size, " genes present in the expression matrix.")
  }

  message("Scoring ", length(gene.sets), " gene sets with method: ", method)

  score_mat <- switch(method,
    AUCell = .score_aucell(x, gene.sets, batch.size, cores),
    Seurat = .score_seurat(x, gene.sets, nbin, ctrl, cores),
    UCell  = .score_ucell(x, gene.sets, cores)
  )

  return(score_mat)
}


#' @param assay Name of the Seurat assay to use. Default: \code{DefaultAssay(x)}.
#' @param layer Data layer to extract: \code{"counts"} or \code{"data"}.
#'   Defaults to \code{"counts"} for AUCell/UCell and \code{"data"} for Seurat method.
#' @param store Storage destination for scores:
#'   \code{"assay"} (default) stores as a new assay,
#'   \code{"metadata"} stores as columns in \code{meta.data}.
#' @param assay.name Name of the new assay when \code{store = "assay"}. Default: the method name.
#' @param prefix Column name prefix when \code{store = "metadata"}. Default: the method name
#'   followed by underscore (e.g. \code{"AUCell_setA"}).
#'
#' @rdname ComputeModuleScore
#' @export
ComputeModuleScore.Seurat <- function(x, gene.sets,
                                      method = c("AUCell", "Seurat", "UCell"),
                                      min.size = 20, batch.size = 500,
                                      nbin = 24, ctrl = 100,
                                      cores = 1, seed = 11,
                                      assay = Seurat::DefaultAssay(x),
                                      layer = NULL,
                                      store = c("assay", "metadata"),
                                      assay.name = NULL,
                                      prefix = NULL, ...) {
  method <- match.arg(method)
  store <- match.arg(store)

  # Default layer: Seurat method uses normalized data, AUCell/UCell use raw counts
  if (is.null(layer)) {
    layer <- if (method == "Seurat") "data" else "counts"
  }

  dge <- Seurat::GetAssayData(x, assay = assay, layer = layer)

  score_mat <- ComputeModuleScore.default(
    x = dge, gene.sets = gene.sets, method = method,
    min.size = min.size, batch.size = batch.size,
    nbin = nbin, ctrl = ctrl, cores = cores, seed = seed, ...
  )

  if (store == "assay") {
    if (is.null(assay.name)) assay.name <- method
    x[[assay.name]] <- Seurat::CreateAssayObject(data = score_mat)
    message("Scores stored in assay: '", assay.name, "'")
  } else {
    # store == "metadata"
    if (is.null(prefix)) prefix <- method
    meta <- as.data.frame(t(score_mat))
    colnames(meta) <- paste0(prefix, "_", rownames(score_mat))
    x <- Seurat::AddMetaData(x, metadata = meta)
    message("Scores stored in meta.data: ", paste(colnames(meta), collapse = ", "))
  }

  return(x)
}


# ============================================================================
# Internal scoring functions
# ============================================================================

#' AUCell scoring with batch processing to reduce memory usage
#' @keywords internal
.score_aucell <- function(x, gene.sets, batch.size, cores) {
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    stop("Package 'AUCell' is required. Install with: BiocManager::install('AUCell')")
  }
  n.cells <- ncol(x)
  batches <- floor((seq_len(n.cells) - 1) / batch.size)
  batch.levels <- unique(batches)

  aucell_batch <- function(i) {
    dge.tmp <- x[, batches == i, drop = FALSE]
    cr <- AUCell::AUCell_buildRankings(dge.tmp, nCores = 1, plotStats = FALSE, verbose = FALSE)
    auc <- AUCell::AUCell_calcAUC(gene.sets, cr, nCores = 1, verbose = FALSE)
    AUCell::getAUC(auc)
  }

  auc_scores <- parallel::mclapply(batch.levels, aucell_batch, mc.cores = cores)
  do.call(cbind, auc_scores)
}


#' Seurat-style module scoring (AddModuleScore-like)
#'
#' Bins genes by mean expression, selects expression-matched control genes,
#' and computes score = mean(feature genes) - mean(control genes).
#' @keywords internal
.score_seurat <- function(x, gene.sets, nbin, ctrl, cores) {
  # Compute mean expression per gene and bin
  data_avg <- Matrix::rowMeans(x)
  data_avg <- sort(data_avg)
  data_cut <- ggplot2::cut_number(
    data_avg + stats::rnorm(length(data_avg)) / 1e30,
    n = nbin, labels = FALSE, right = FALSE
  )
  names(data_cut) <- names(data_avg)

  score_fun <- function(i) {
    features_use <- gene.sets[[i]]
    # Collect control gene candidates from the same expression bins
    ctrl_pool <- unlist(lapply(features_use, function(g) {
      names(data_cut[data_cut == data_cut[g]])
    }))
    n_sample <- min(ctrl * length(features_use), length(ctrl_pool))
    ctrl_use <- sample(ctrl_pool, size = n_sample, replace = FALSE)

    feat_score <- Matrix::colMeans(x[features_use, , drop = FALSE])
    ctrl_score <- Matrix::colMeans(x[ctrl_use, , drop = FALSE])
    feat_score - ctrl_score
  }

  scores <- parallel::mclapply(seq_along(gene.sets), score_fun, mc.cores = cores)
  score_mat <- do.call(rbind, scores)
  rownames(score_mat) <- names(gene.sets)
  score_mat
}


#' UCell scoring
#' @keywords internal
.score_ucell <- function(x, gene.sets, cores) {
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("Package 'UCell' is required. Install with: BiocManager::install('UCell')")
  }
  scores <- UCell::ScoreSignatures_UCell(x, features = gene.sets, ncores = cores)
  score_mat <- t(as.matrix(scores))
  # Remove "_UCell" suffix added by UCell
  rownames(score_mat) <- sub("_UCell$", "", rownames(score_mat))
  score_mat
}
