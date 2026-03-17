# ============================================================================
# Pathway Analysis Functions
# ============================================================================
#
# Refactored from: GRN/downstream/2-core_questions_of_perturbed_dataset_II.Rmd
#
# Exported functions:
#   - RunPathwayAnalysis()  : GSEA or SCPA across cell types
#   - PlotPathwayBubble()   : bubble plot for pathway results
#
# Two analysis methods:
#   - GSEA:  FindMarkers -> logFC ranking -> clusterProfiler::GSEA
#   - SCPA:  seurat_extract -> compare_pathways (distribution-based)
#
# Dependencies:
#   - Seurat, dplyr, ggplot2
#   - clusterProfiler (for GSEA method)
#   - SCPA (for SCPA method)
#   - parse_gene_sets() from compute_module_score.R or scMMR package
# ============================================================================


# ── Helper: clean pathway names ─────────────────────────────────────────────

#' Clean pathway names by removing common prefixes and converting to Title Case
#'
#' Removes database-specific prefixes such as \code{HALLMARK_}, \code{KEGG_},
#' \code{REACTOME_}, etc., replaces underscores with spaces, and converts the
#' result to Title Case for cleaner visualization.
#'
#' @param names_vec Character vector of pathway names.
#' @param prefix Optional custom prefix to remove (regex pattern).
#'   If \code{NULL} (default), common prefixes are removed automatically.
#' @return Cleaned character vector.
#'
#' @examples
#' \dontrun{
#' .clean_pathway_names(c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "KEGG_GLYCOLYSIS"))
#' # [1] "Tnfa Signaling Via Nfkb" "Glycolysis"
#' }
#'
#' @keywords internal
.clean_pathway_names <- function(names_vec, prefix = NULL) {
  if (!is.null(prefix)) {
    names_vec <- sub(prefix, "", names_vec)
  } else {
    names_vec <- sub("^HALLMARK_", "", names_vec)
    names_vec <- sub("^KEGG_", "", names_vec)
    names_vec <- sub("^REACTOME_", "", names_vec)
    names_vec <- sub("^GOBP_", "", names_vec)
    names_vec <- sub("^GOCC_", "", names_vec)
    names_vec <- sub("^GOMF_", "", names_vec)
    names_vec <- sub("^WP_", "", names_vec)
    names_vec <- sub("^BIOCARTA_", "", names_vec)
    names_vec <- sub("^PID_", "", names_vec)
  }
  names_vec <- gsub("_", " ", names_vec)
  names_vec <- stringr::str_to_title(tolower(names_vec))
  names_vec
}


# ── Helper: convert gene set formats ────────────────────────────────────────

#' Convert named list of gene sets to TERM2GENE data.frame
#'
#' Creates a two-column data.frame (term, gene) required by
#' \code{clusterProfiler::GSEA()}.
#'
#' @param gene.sets Named list of character vectors.
#' @return A 2-column data.frame with columns \code{term} and \code{gene}.
#' @keywords internal
.gs_to_term2gene <- function(gene.sets) {
  data.frame(
    term = rep(names(gene.sets), vapply(gene.sets, length, integer(1))),
    gene = unlist(gene.sets, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}


#' Convert named list of gene sets to SCPA pathway format
#'
#' SCPA's \code{compare_pathways()} expects pathways as a list of data.frames,
#' each with columns \code{Pathway} and \code{Genes}.
#'
#' @param gene.sets Named list of character vectors.
#' @return A list of data.frames suitable for \code{SCPA::compare_pathways()}.
#' @keywords internal
.gs_to_scpa_format <- function(gene.sets) {
  lapply(names(gene.sets), function(name) {
    data.frame(Pathway = name, Genes = gene.sets[[name]], stringsAsFactors = FALSE)
  })
}


# ── Internal: GSEA per cell type ────────────────────────────────────────────

#' Run GSEA for each cell type via FindMarkers + clusterProfiler::GSEA
#'
#' For each cell type, calls \code{FindMarkers()} to compute log2 fold changes,
#' ranks genes by logFC, and runs \code{clusterProfiler::GSEA()}.
#'
#' @param seu Seurat object with Idents set to composite group.
#' @param gene.sets Named list of gene sets.
#' @param celltypes Character vector of cell types to process.
#' @param ident.1,ident.2 Condition labels.
#' @param pvalue.cutoff Passed to \code{clusterProfiler::GSEA()}.
#' @param logfc.threshold Passed to \code{FindMarkers()}.
#' @param eps Passed to \code{clusterProfiler::GSEA()}.
#' @param minGSSize,maxGSSize Gene set size limits for GSEA.
#' @param verbose Logical; print progress.
#' @return A data.frame with standardized columns.
#' @keywords internal
.run_gsea_across_celltypes <- function(seu, gene.sets, celltypes,
                                       ident.1, ident.2,
                                       pvalue.cutoff, logfc.threshold, eps,
                                       minGSSize, maxGSSize, verbose) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.\n",
         "Install with: BiocManager::install('clusterProfiler')")
  }

  term2gene <- .gs_to_term2gene(gene.sets)

  run_one <- function(ct) {
    id1 <- paste0(ct, "_", ident.1)
    id2 <- paste0(ct, "_", ident.2)

    tryCatch({
      all_ids <- unique(as.character(Seurat::Idents(seu)))
      if (!id1 %in% all_ids) {
        if (verbose) message("  Skipping ", ct, ": '", id1, "' not found.")
        return(NULL)
      }
      if (!id2 %in% all_ids) {
        if (verbose) message("  Skipping ", ct, ": '", id2, "' not found.")
        return(NULL)
      }

      if (verbose) message("  [GSEA] ", ct, ": ", id1, " vs ", id2)

      fc <- Seurat::FindMarkers(seu, ident.1 = id1, ident.2 = id2,
                                logfc.threshold = logfc.threshold,
                                verbose = FALSE)

      genes.fc <- stats::setNames(fc$avg_log2FC, rownames(fc))
      genes.fc <- sort(genes.fc, decreasing = TRUE)

      gsea.res <- clusterProfiler::GSEA(
        genes.fc, TERM2GENE = term2gene,
        pvalueCutoff = pvalue.cutoff, eps = eps,
        minGSSize = minGSSize, maxGSSize = maxGSSize,
        verbose = FALSE
      )
      gsea.df <- as.data.frame(gsea.res)

      if (nrow(gsea.df) == 0) {
        if (verbose) message("    No significant pathways.")
        return(NULL)
      }

      data.frame(
        Celltype         = ct,
        Pathway          = gsea.df$Description,
        Score            = gsea.df$NES,
        PValue           = gsea.df$pvalue,
        AdjPValue        = gsea.df$p.adjust,
        QValue           = gsea.df$qvalue,
        Sign             = ifelse(gsea.df$NES > 0, "Activated", "Repressed"),
        SetSize          = gsea.df$setSize,
        EnrichmentScore  = gsea.df$enrichmentScore,
        CoreEnrichment   = gsea.df$core_enrichment,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  Error in ", ct, ": ", conditionMessage(e))
      return(NULL)
    })
  }

  if (requireNamespace("pbapply", quietly = TRUE)) {
    results <- pbapply::pblapply(celltypes, run_one)
  } else {
    results <- lapply(celltypes, run_one)
  }

  result_df <- do.call(rbind, results)
  if (is.null(result_df) || nrow(result_df) == 0) {
    warning("No significant pathways found for any cell type.")
    return(data.frame())
  }
  result_df
}


# ── Internal: SCPA per cell type ────────────────────────────────────────────

#' Run SCPA for each cell type via seurat_extract + compare_pathways
#'
#' For each cell type, extracts expression matrices for the two conditions
#' and runs \code{SCPA::compare_pathways()} for distribution-based comparison.
#'
#' @param seu Seurat object with composite group column \code{..pw_group..}.
#' @param gene.sets Named list of gene sets.
#' @param celltypes Character vector of cell types to process.
#' @param ident.1,ident.2 Condition labels.
#' @param cores Number of parallel cores for SCPA.
#' @param verbose Logical; print progress.
#' @return A data.frame with standardized columns.
#' @keywords internal
.run_scpa_across_celltypes <- function(seu, gene.sets, celltypes,
                                       ident.1, ident.2,
                                       cores, verbose) {
  if (!requireNamespace("SCPA", quietly = TRUE)) {
    stop("Package 'SCPA' is required.\n",
         "Install with: devtools::install_github('jackbibby1/SCPA')")
  }

  scpa_pathways <- .gs_to_scpa_format(gene.sets)

  run_one <- function(ct) {
    id1 <- paste0(ct, "_", ident.1)
    id2 <- paste0(ct, "_", ident.2)

    tryCatch({
      if (verbose) message("  [SCPA] ", ct, ": ", id1, " vs ", id2)

      mat.1 <- SCPA::seurat_extract(seu, meta1 = "..pw_group..",
                                    value_meta1 = id1)
      mat.2 <- SCPA::seurat_extract(seu, meta1 = "..pw_group..",
                                    value_meta1 = id2)

      scpa_out <- SCPA::compare_pathways(
        samples  = list(mat.1, mat.2),
        pathways = scpa_pathways,
        parallel = (cores > 1)
      )

      if (nrow(scpa_out) == 0) return(NULL)

      data.frame(
        Celltype         = ct,
        Pathway          = scpa_out$Pathway,
        Score            = scpa_out$FC,
        PValue           = scpa_out$Pval,
        AdjPValue        = scpa_out$adjPval,
        QValue           = scpa_out$qval,
        Sign             = ifelse(scpa_out$FC > 0, "Activated", "Repressed"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  Error in ", ct, ": ", conditionMessage(e))
      return(NULL)
    })
  }

  if (requireNamespace("pbapply", quietly = TRUE)) {
    results <- pbapply::pblapply(celltypes, run_one)
  } else {
    results <- lapply(celltypes, run_one)
  }

  result_df <- do.call(rbind, results)
  if (is.null(result_df) || nrow(result_df) == 0) {
    warning("No pathways found for any cell type.")
    return(data.frame())
  }
  result_df
}


# ============================================================================
# Main Function: RunPathwayAnalysis
# ============================================================================

#' Run Pathway Analysis across Cell Types
#'
#' Perform GSEA or SCPA pathway analysis for each cell type in a Seurat object.
#' Compares two conditions (\code{ident.1} vs \code{ident.2}) within each cell
#' type, and returns a standardized data.frame suitable for visualization with
#' \code{\link{PlotPathwayBubble}}.
#'
#' @param seu A Seurat object with cell type and condition labels in
#'   \code{meta.data}.
#' @param gene.sets Gene sets in one of three formats:
#'   \itemize{
#'     \item A named list of character vectors (gene names).
#'     \item A data.frame with at least two columns: term and gene.
#'     \item A file path to a GMT file.
#'   }
#'   Internally parsed by \code{\link{parse_gene_sets}}.
#' @param method Pathway analysis method: \code{"GSEA"} or \code{"SCPA"}.
#'   \describe{
#'     \item{\code{"GSEA"}}{Uses \code{FindMarkers()} to obtain log2 fold
#'       changes, ranks genes, and runs \code{clusterProfiler::GSEA()}.
#'       Requires \pkg{clusterProfiler}.}
#'     \item{\code{"SCPA"}}{Uses \code{SCPA::seurat_extract()} and
#'       \code{SCPA::compare_pathways()} for distribution-based pathway
#'       comparison. Requires \pkg{SCPA}.}
#'   }
#' @param celltype.col Column name in \code{meta.data} containing cell type
#'   labels. Default: \code{"celltype"}.
#' @param group.col Column name in \code{meta.data} containing
#'   condition/group labels. Default: \code{"group"}.
#' @param ident.1 Treatment/test condition label (a value in \code{group.col}).
#'   Fold change direction: \code{ident.1 / ident.2}.
#' @param ident.2 Control/reference condition label (a value in
#'   \code{group.col}).
#' @param celltypes Character vector specifying which cell types to analyze.
#'   Default: \code{NULL} (all cell types in \code{celltype.col}).
#' @param pvalue.cutoff P-value cutoff for GSEA enrichment results.
#'   Default: 0.05.
#' @param logfc.threshold Minimum log2 fold-change threshold for
#'   \code{FindMarkers()} (GSEA only). Default: 0 (include all genes).
#' @param eps Boundary correction for GSEA p-value calculation.
#'   Default: 0 (exact p-values).
#' @param minGSSize Minimum number of genes in a gene set for GSEA.
#'   Default: 10.
#' @param maxGSSize Maximum number of genes in a gene set for GSEA.
#'   Default: 500.
#' @param clean.names Logical; if \code{TRUE}, remove common database prefixes
#'   (e.g. \code{HALLMARK_}, \code{KEGG_}) and convert to Title Case.
#'   Default: \code{TRUE}.
#' @param prefix Custom regex prefix to remove when \code{clean.names = TRUE}.
#'   Default: \code{NULL} (auto-detect common prefixes).
#' @param cores Number of parallel cores for SCPA. Default: 1.
#' @param verbose Logical; print progress messages. Default: \code{TRUE}.
#'
#' @return A data.frame with standardized columns:
#'   \describe{
#'     \item{\code{Celltype}}{Cell type name.}
#'     \item{\code{Pathway}}{Pathway name (cleaned if
#'       \code{clean.names = TRUE}).}
#'     \item{\code{Score}}{Enrichment score: NES (GSEA) or FC (SCPA).}
#'     \item{\code{PValue}}{Raw p-value.}
#'     \item{\code{AdjPValue}}{Adjusted p-value (BH).}
#'     \item{\code{QValue}}{Q-value.}
#'     \item{\code{Sign}}{\code{"Activated"} (Score > 0) or
#'       \code{"Repressed"} (Score < 0).}
#'     \item{\code{Method}}{\code{"GSEA"} or \code{"SCPA"}.}
#'   }
#'   GSEA results additionally include: \code{SetSize},
#'   \code{EnrichmentScore}, \code{CoreEnrichment}.
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#' library(Seurat)
#'
#' # Gene sets from GMT file
#' gmt_file <- system.file("extdata", "gmt",
#'   "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")
#'
#' # --- GSEA ---
#' gsea_res <- RunPathwayAnalysis(
#'   seu, gene.sets = gmt_file, method = "GSEA",
#'   celltype.col = "celltype", group.col = "group",
#'   ident.1 = "PH", ident.2 = "PT"
#' )
#' PlotPathwayBubble(gsea_res)
#'
#' # --- SCPA (specific cell types) ---
#' scpa_res <- RunPathwayAnalysis(
#'   seu, gene.sets = gmt_file, method = "SCPA",
#'   celltype.col = "celltype", group.col = "group",
#'   ident.1 = "PH", ident.2 = "PT",
#'   celltypes = c("Parathyroid cells", "Fibroblasts")
#' )
#' PlotPathwayBubble(scpa_res, size.by = "qvalue")
#' }
#'
#' @seealso \code{\link{PlotPathwayBubble}}, \code{\link{parse_gene_sets}}
#' @export
RunPathwayAnalysis <- function(seu, gene.sets,
                               method = c("GSEA", "SCPA"),
                               celltype.col = "celltype",
                               group.col = "group",
                               ident.1 = NULL,
                               ident.2 = NULL,
                               celltypes = NULL,
                               pvalue.cutoff = 0.05,
                               logfc.threshold = 0,
                               eps = 0,
                               minGSSize = 10,
                               maxGSSize = 500,
                               clean.names = TRUE,
                               prefix = NULL,
                               cores = 1,
                               verbose = TRUE) {
  method <- match.arg(method)

  # ── Input validation ──
  if (!inherits(seu, "Seurat")) {
    stop("'seu' must be a Seurat object.")
  }
  if (!celltype.col %in% colnames(seu@meta.data)) {
    stop("Column '", celltype.col, "' not found in meta.data.")
  }
  if (!group.col %in% colnames(seu@meta.data)) {
    stop("Column '", group.col, "' not found in meta.data.")
  }
  if (is.null(ident.1) || is.null(ident.2)) {
    stop("Both 'ident.1' and 'ident.2' must be specified.\n",
         "  Available groups: ",
         paste(unique(seu@meta.data[[group.col]]), collapse = ", "))
  }

  # ── Determine cell types ──
  if (is.null(celltypes)) {
    ct_vals <- seu@meta.data[[celltype.col]]
    celltypes <- if (is.factor(ct_vals)) levels(ct_vals) else sort(unique(ct_vals))
  }

  if (verbose) {
    message("=== RunPathwayAnalysis ===")
    message("Method: ", method)
    message("Comparison: ", ident.1, " vs ", ident.2)
    message("Cell types (", length(celltypes), "): ",
            paste(celltypes, collapse = ", "))
  }

  # ── Parse gene sets ──
  if (!exists("parse_gene_sets", mode = "function")) {
    stop("Function 'parse_gene_sets' not found.\n",
         "  Please run: source('compute_module_score.R') or library(scMMR)")
  }
  gene.sets.list <- parse_gene_sets(gene.sets)

  if (verbose) message("Gene sets: ", length(gene.sets.list), " pathways loaded.")

  # Clean pathway names
  if (clean.names) {
    names(gene.sets.list) <- .clean_pathway_names(names(gene.sets.list),
                                                  prefix = prefix)
  }

  # ── Create composite group column ──
  orig_idents <- Seurat::Idents(seu)
  pw_vec <- paste(seu@meta.data[[celltype.col]],
                  seu@meta.data[[group.col]], sep = "_")
  seu <- Seurat::AddMetaData(seu, metadata = pw_vec, col.name = "..pw_group..")
  Seurat::Idents(seu) <- "..pw_group.."

  # ── Run analysis ──
  result <- switch(method,
    GSEA = .run_gsea_across_celltypes(
      seu, gene.sets.list, celltypes, ident.1, ident.2,
      pvalue.cutoff, logfc.threshold, eps, minGSSize, maxGSSize, verbose
    ),
    SCPA = .run_scpa_across_celltypes(
      seu, gene.sets.list, celltypes, ident.1, ident.2,
      cores, verbose
    )
  )

  if (nrow(result) > 0) {
    result$Method <- method
  }

  if (verbose) {
    n_sig <- sum(result$AdjPValue < 0.05, na.rm = TRUE)
    message("Done. Total results: ", nrow(result),
            " (", n_sig, " with adjP < 0.05)")
  }

  return(result)
}


# ============================================================================
# Plotting Function: PlotPathwayBubble
# ============================================================================

#' Pathway Bubble Plot
#'
#' Create a bubble plot visualizing pathway analysis results across cell types.
#' Automatically adapts visual encoding to GSEA or SCPA results based on the
#' \code{Method} column.
#'
#' @section Visual encoding:
#' \describe{
#'   \item{GSEA}{
#'     \itemize{
#'       \item \strong{fill}: NES value (continuous blue-white-red gradient).
#'       \item \strong{size}: determined by \code{size.by} (default:
#'         \code{-log10(q-value)}).
#'       \item \strong{facet}: Activated / Repressed (by default).
#'     }
#'   }
#'   \item{SCPA}{
#'     \itemize{
#'       \item \strong{fill}: Direction (Activated = red, Repressed = blue).
#'       \item \strong{size}: determined by \code{size.by} (default:
#'         \code{|FC|}).
#'       \item \strong{no facet} by default.
#'     }
#'   }
#' }
#'
#' @param data A data.frame returned by \code{\link{RunPathwayAnalysis}}.
#'   Must contain columns: \code{Celltype}, \code{Pathway}, \code{Score},
#'   \code{Sign}, \code{Method}, and at least one of \code{QValue} or
#'   \code{AdjPValue}.
#' @param top.n Number of top pathways to display. Default: 10.
#' @param size.by What to map to point size:
#'   \describe{
#'     \item{\code{"pvalue"}}{Point size = \code{-log10(pvalue.col)}.
#'       Default for GSEA.}
#'     \item{\code{"score"}}{Point size = \code{|Score|} (i.e. |NES| or |FC|).
#'       Default for SCPA.}
#'     \item{\code{"qvalue"}}{Point size = raw \code{QValue}.}
#'   }
#'   Default: \code{NULL} (auto-selected based on method).
#' @param pvalue.col Column name for significance filtering and
#'   \code{size.by = "pvalue"} mapping. One of \code{"QValue"} or
#'   \code{"AdjPValue"}. Default: auto-detected (\code{"QValue"} for GSEA,
#'   \code{"AdjPValue"} for SCPA).
#' @param pvalue.cutoff Significance threshold for selecting top pathways.
#'   Default: 0.05.
#' @param select.by Strategy for choosing top pathways:
#'   \describe{
#'     \item{\code{"frequency"}}{Pathways most frequently significant across
#'       cell types (default).}
#'     \item{\code{"score"}}{Pathways with highest average |Score|.}
#'   }
#' @param pathways Character vector of specific pathway names to display.
#'   Overrides \code{top.n} and \code{select.by} when provided.
#' @param order.by How to order pathways on the y-axis:
#'   \describe{
#'     \item{\code{"frequency"}}{By count of significant appearances
#'       (default).}
#'     \item{\code{"score"}}{By mean Score value (low to high).}
#'     \item{\code{"sign"}}{By mean Score (high to low, Activated first).}
#'   }
#' @param base.size Base font size for the plot. Default: 15.
#' @param colors.gsea Color gradient vector for GSEA NES fill.
#'   Default: \code{c("blue", "white", "red")}.
#' @param colors.scpa Named vector for SCPA Sign fill.
#'   Default: \code{c(Activated = "red", Repressed = "blue")}.
#' @param facet.by.sign Logical; whether to facet by Activated/Repressed.
#'   Default: \code{TRUE} for GSEA, \code{FALSE} for SCPA. Set explicitly
#'   to override.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # GSEA: default (fill = NES, size = -log10(q))
#' PlotPathwayBubble(gsea_res)
#'
#' # GSEA: size by |NES|
#' PlotPathwayBubble(gsea_res, size.by = "score")
#'
#' # SCPA: default (fill = Sign, size = |FC|)
#' PlotPathwayBubble(scpa_res)
#'
#' # SCPA: size by q-value
#' PlotPathwayBubble(scpa_res, size.by = "qvalue")
#'
#' # Show specific pathways
#' PlotPathwayBubble(gsea_res, pathways = c("Oxidative Phosphorylation",
#'                                           "Interferon Alpha Response"))
#'
#' # Top 15, selected by score, ordered by score
#' PlotPathwayBubble(gsea_res, top.n = 15, select.by = "score",
#'                   order.by = "score")
#' }
#'
#' @seealso \code{\link{RunPathwayAnalysis}}
#' @import ggplot2
#' @importFrom dplyr filter group_by summarise arrange slice_head pull mutate
#' @importFrom rlang .data
#' @export
PlotPathwayBubble <- function(data,
                              top.n = 10,
                              size.by = NULL,
                              pvalue.col = NULL,
                              pvalue.cutoff = 0.05,
                              select.by = c("frequency", "score"),
                              pathways = NULL,
                              order.by = c("frequency", "score", "sign"),
                              base.size = 15,
                              colors.gsea = c("blue", "white", "red"),
                              colors.scpa = c(Activated = "red", Repressed = "blue"),
                              facet.by.sign = NULL) {
  select.by <- match.arg(select.by)
  order.by  <- match.arg(order.by)

  if (nrow(data) == 0) stop("Input data is empty.")

  method <- unique(data$Method)[1]

  # ── Default settings based on method ──
  if (is.null(pvalue.col)) {
    pvalue.col <- if (method == "GSEA" && "QValue" %in% colnames(data)) {
      "QValue"
    } else {
      "AdjPValue"
    }
  }
  if (is.null(facet.by.sign)) {
    facet.by.sign <- (method == "GSEA")
  }
  if (is.null(size.by)) {
    size.by <- if (method == "GSEA") "pvalue" else "score"
  }
  size.by <- match.arg(size.by, c("pvalue", "score", "qvalue"))

  # ── Validate size.by ──
  if (size.by == "qvalue" && !"QValue" %in% colnames(data)) {
    warning("'QValue' column not found, falling back to size.by = 'pvalue'.")
    size.by <- "pvalue"
  }

  # ── Select top pathways ──
  if (is.null(pathways)) {
    sig_data <- data[data[[pvalue.col]] < pvalue.cutoff, , drop = FALSE]
    if (nrow(sig_data) == 0) {
      warning("No pathways with ", pvalue.col, " < ", pvalue.cutoff,
              ". Try increasing pvalue.cutoff.")
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::ggtitle("No significant pathways"))
    }

    if (select.by == "frequency") {
      path_rank <- sig_data %>%
        dplyr::group_by(.data$Pathway) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(.data$n)) %>%
        dplyr::slice_head(n = top.n)
    } else {
      path_rank <- sig_data %>%
        dplyr::group_by(.data$Pathway) %>%
        dplyr::summarise(avg_score = mean(abs(.data$Score)), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(.data$avg_score)) %>%
        dplyr::slice_head(n = top.n)
    }
    top_paths <- path_rank$Pathway
  } else {
    top_paths <- pathways
  }

  # ── Subset and order ──
  data.plot <- data %>%
    dplyr::filter(.data$Pathway %in% top_paths)

  if (nrow(data.plot) == 0) {
    warning("No data to plot after filtering.")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::ggtitle("No pathways to display"))
  }

  # Determine pathway order on y-axis
  if (order.by == "frequency") {
    path_order <- data.plot %>%
      dplyr::filter(.data[[pvalue.col]] < pvalue.cutoff) %>%
      dplyr::group_by(.data$Pathway) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$n)) %>%
      dplyr::pull(.data$Pathway)
  } else if (order.by == "score") {
    path_order <- data.plot %>%
      dplyr::group_by(.data$Pathway) %>%
      dplyr::summarise(avg = mean(.data$Score), .groups = "drop") %>%
      dplyr::arrange(.data$avg) %>%
      dplyr::pull(.data$Pathway)
  } else {
    # order.by == "sign": Activated (high score) first
    path_order <- data.plot %>%
      dplyr::group_by(.data$Pathway) %>%
      dplyr::summarise(avg = mean(.data$Score), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$avg)) %>%
      dplyr::pull(.data$Pathway)
  }
  data.plot$Pathway <- factor(data.plot$Pathway, levels = rev(path_order))

  # ── Compute size variable ──
  if (size.by == "pvalue") {
    data.plot$.size_var. <- -log10(data.plot[[pvalue.col]])
    size_label <- bquote(-log[10](italic(q)))
  } else if (size.by == "score") {
    data.plot$.size_var. <- abs(data.plot$Score)
    size_label <- if (method == "GSEA") "|NES|" else "|FC|"
  } else {
    # size.by == "qvalue"
    data.plot$.size_var. <- data.plot$QValue
    size_label <- "q value"
  }

  # ── Build plot based on method ──
  if (method == "GSEA") {
    p <- ggplot2::ggplot(
      data.plot,
      ggplot2::aes(x = .data$Celltype, y = .data$Pathway,
                   fill = .data$Score,
                   size = .data$.size_var.)
    ) +
      ggplot2::geom_point(shape = 22, color = "black", stroke = 1) +
      ggplot2::scale_fill_gradientn(colors = colors.gsea, name = "NES") +
      ggplot2::scale_size_continuous(name = size_label)

  } else {
    # SCPA
    p <- ggplot2::ggplot(
      data.plot,
      ggplot2::aes(x = .data$Celltype, y = .data$Pathway,
                   fill = .data$Sign,
                   size = .data$.size_var.)
    ) +
      ggplot2::geom_point(shape = 22, color = "black", stroke = 1) +
      ggplot2::scale_fill_manual(values = colors.scpa, name = "") +
      ggplot2::scale_size_continuous(name = size_label) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(override.aes = list(size = 5))
      )
  }

  # Facet
  if (facet.by.sign) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$Sign),
                                 scales = "free_y")
  }

  # ── Common theme ──
  p <- p +
    ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
    ggplot2::theme_bw(base_size = base.size) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(size = base.size, face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1,
                                                color = "black"),
      axis.text.y      = ggplot2::element_text(color = "black"),
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}
