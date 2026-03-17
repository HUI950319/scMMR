# ============================================================================
# Trajectory analysis: pseudotime-associated genes and GSEA
# ============================================================================
#
# Exported functions:
#   - RunTraceGene()  : Find pseudotime-associated genes (Spearman / GAM)
#   - RunTraceGSEA()  : Spearman correlation ranking -> GSEA per lineage
#
# Dependencies:
#   - Seurat, SeuratObject
#   - clusterProfiler (for RunTraceGSEA)
#   - mgcv (for GAM in RunTraceGene)
#   - .fit_dynamic_feature() from dynamic_plot.R
#   - parse_gene_sets(), .gs_to_term2gene(), .clean_pathway_names() from scMMR
# ============================================================================


#' Run Trajectory GSEA Along Pseudotime
#'
#' Perform Gene Set Enrichment Analysis (GSEA) along pseudotime trajectories.
#' For each lineage, computes Spearman correlation between gene expression and
#' pseudotime, uses the correlation-based ranking to run
#' \code{clusterProfiler::GSEA()}, and identifies pathways activated or
#' repressed along the trajectory.
#'
#' @details
#' The workflow for each lineage is:
#' \enumerate{
#'   \item Subset cells with finite pseudotime values (\code{is.finite()}).
#'   \item Use all genes by default (recommended for GSEA); or use HVG
#'     if \code{features = "HVG"}.
#'   \item Filter genes by minimum expression percentage (\code{min.pct}).
#'   \item Compute Spearman correlation (\code{cor.test(method = "spearman")})
#'     between each gene's expression and pseudotime.
#'   \item Build a gene ranking from the correlation statistics.
#'   \item Run \code{clusterProfiler::GSEA()} with the ranking.
#' }
#'
#' Two ranking metrics are available:
#' \describe{
#'   \item{\code{"rho"}}{Uses Spearman rho directly. Genes positively
#'     correlated with pseudotime rank high (expression increases along
#'     trajectory); negatively correlated rank low. This mode uses
#'     vectorized \code{stats::cor()} for speed.}
#'   \item{\code{"signed_logp"}}{Uses \code{-log10(p) * sign(rho)} for
#'     sharper separation of highly significant correlations. This mode
#'     calls \code{cor.test()} per gene and is slower but provides
#'     stronger discrimination.}
#' }
#'
#' @param seu A Seurat object with pseudotime stored in \code{meta.data}.
#' @param lineages Character vector of column names in \code{seu@@meta.data}
#'   containing pseudotime values. Cells not on a lineage should have
#'   \code{NA} or \code{Inf}. Multiple lineages are analyzed independently.
#' @param gene.sets Gene sets in one of three formats:
#'   \itemize{
#'     \item A named list of character vectors (gene names).
#'     \item A data.frame with columns term and gene.
#'     \item A file path to a GMT file.
#'   }
#'   Parsed internally by \code{\link{parse_gene_sets}}.
#' @param features Character vector of gene names to include in the
#'   correlation analysis. Default: \code{NULL} (use all genes in the assay,
#'   which is recommended for GSEA).
#' @param n_candidates Not used when \code{features = NULL} (all genes are
#'   used by default). Only applies when you want to restrict to top variable
#'   features by setting \code{features = "HVG"}. Default: 2000.
#' @param min.pct Minimum fraction of lineage cells expressing a gene
#'   (expression > 0) for the gene to be included. Default: 0.1.
#' @param ranking.metric Ranking method for GSEA:
#'   \code{"rho"} (Spearman rho, fast) or \code{"signed_logp"}
#'   (\code{-log10(p) * sign(rho)}, more discriminative). Default: \code{"rho"}.
#' @param assay Seurat assay to use. Default: \code{DefaultAssay(seu)}.
#' @param layer Data layer to extract. Default: \code{"data"}
#'   (log-normalized expression recommended for correlation analysis).
#' @param pvalue.cutoff P-value cutoff for \code{clusterProfiler::GSEA()}.
#'   Default: 1 (return all results; filter post-hoc).
#' @param minGSSize Minimum number of genes in a gene set. Default: 10.
#' @param maxGSSize Maximum number of genes in a gene set. Default: 500.
#' @param eps Boundary correction for GSEA p-value calculation.
#'   Default: 0 (exact p-values).
#' @param clean.names Logical; if \code{TRUE}, remove common database prefixes
#'   (e.g. \code{HALLMARK_}, \code{KEGG_}) and convert to Title Case.
#'   Default: \code{TRUE}.
#' @param prefix Custom regex prefix to remove when \code{clean.names = TRUE}.
#'   Default: \code{NULL} (auto-detect common prefixes).
#' @param verbose Logical; print progress messages. Default: \code{TRUE}.
#'
#' @return A data.frame compatible with \code{\link{PlotPathwayBubble}},
#'   with columns:
#'   \describe{
#'     \item{\code{Celltype}}{Lineage name (for \code{PlotPathwayBubble}
#'       x-axis compatibility).}
#'     \item{\code{Lineage}}{Lineage name.}
#'     \item{\code{Pathway}}{Pathway name (cleaned if
#'       \code{clean.names = TRUE}).}
#'     \item{\code{Score}}{NES (Normalized Enrichment Score). Positive =
#'       pathway activated along trajectory; negative = repressed.}
#'     \item{\code{PValue}}{Raw p-value.}
#'     \item{\code{AdjPValue}}{BH-adjusted p-value.}
#'     \item{\code{QValue}}{Q-value.}
#'     \item{\code{Sign}}{\code{"Activated"} (NES > 0) or
#'       \code{"Repressed"} (NES < 0).}
#'     \item{\code{SetSize}}{Number of genes in set.}
#'     \item{\code{EnrichmentScore}}{Raw enrichment score.}
#'     \item{\code{CoreEnrichment}}{Core enrichment genes
#'       (slash-separated).}
#'     \item{\code{RankingMetric}}{Ranking method used
#'       (\code{"rho"} or \code{"signed_logp"}).}
#'     \item{\code{NCells}}{Number of cells on this lineage.}
#'     \item{\code{Method}}{\code{"GSEA"}.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#'
#' # Hallmark gene sets
#' gmt <- system.file("extdata", "gmt",
#'   "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")
#'
#' # Single lineage
#' res <- RunTraceGSEA(seu, lineages = "Lineage1", gene.sets = gmt)
#' PlotPathwayBubble(res)
#'
#' # Multiple lineages with signed_logp ranking
#' res <- RunTraceGSEA(
#'   seu,
#'   lineages = c("Lineage1", "Lineage2", "Lineage3"),
#'   gene.sets = gmt,
#'   ranking.metric = "signed_logp"
#' )
#' PlotPathwayBubble(res, top.n = 15, facet.by.sign = TRUE)
#'
#' # KEGG pathways with custom gene list
#' kegg_gmt <- system.file("extdata", "gmt",
#'   "c2.cp.kegg.v2022.1.Hs.symbols.gmt", package = "scMMR")
#' res <- RunTraceGSEA(
#'   seu,
#'   lineages = "Lineage1",
#'   gene.sets = kegg_gmt,
#'   features = Seurat::VariableFeatures(seu),
#'   min.pct = 0.05,
#'   ranking.metric = "signed_logp"
#' )
#' }
#'
#' @seealso \code{\link{PlotPathwayBubble}}, \code{\link{RunPathwayAnalysis}},
#'   \code{\link{PlotDynamicFeatures}}, \code{\link{parse_gene_sets}}
#' @export
RunTraceGSEA <- function(seu,
                           lineages,
                           gene.sets,
                           features       = NULL,
                           n_candidates   = 2000,
                           min.pct        = 0.1,
                           ranking.metric = c("rho", "signed_logp"),
                           assay          = NULL,
                           layer          = "data",
                           pvalue.cutoff  = 1,
                           minGSSize      = 10,
                           maxGSSize      = 500,
                           eps            = 0,
                           clean.names    = TRUE,
                           prefix         = NULL,
                           verbose        = TRUE) {

  ranking.metric <- match.arg(ranking.metric)

  # ── Input validation ──
  if (!inherits(seu, "Seurat")) {
    stop("'seu' must be a Seurat object.")
  }

  missing_lin <- lineages[!lineages %in% colnames(seu@meta.data)]
  if (length(missing_lin) > 0) {
    stop("Lineage column(s) not found in meta.data: ",
         paste(missing_lin, collapse = ", "))
  }

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.\n",
         "Install with: BiocManager::install('clusterProfiler')")
  }

  # ── Resolve assay ──
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seu)
  }
  all_genes <- rownames(seu[[assay]])

  # ── Determine feature universe ──
  if (is.null(features)) {
    # Default: use ALL genes (standard for GSEA)
    features <- all_genes
    if (verbose) message("Using all ", length(features), " genes (recommended for GSEA)")
  } else if (length(features) == 1 && features == "HVG") {
    # Shortcut: use top n_candidates variable features
    vf <- SeuratObject::VariableFeatures(seu, assay = assay)
    if (length(vf) < n_candidates) {
      if (verbose) message("Finding top ", n_candidates, " variable features...")
      seu_tmp <- Seurat::FindVariableFeatures(seu, nfeatures = n_candidates,
                                               assay = assay, verbose = FALSE)
      features <- SeuratObject::VariableFeatures(seu_tmp, assay = assay)
    } else {
      features <- vf[seq_len(n_candidates)]
    }
  }
  features <- intersect(features, all_genes)
  if (length(features) == 0) {
    stop("No valid features found in assay '", assay, "'.")
  }

  # ── Parse gene sets ──
  gene.sets.list <- parse_gene_sets(gene.sets)

  if (clean.names) {
    names(gene.sets.list) <- .clean_pathway_names(names(gene.sets.list),
                                                   prefix = prefix)
  }

  term2gene <- .gs_to_term2gene(gene.sets.list)

  if (verbose) {
    message("=== RunTraceGSEA ===")
    message("Lineages: ", paste(lineages, collapse = ", "))
    message("Ranking metric: ", ranking.metric)
    message("Feature universe: ", length(features), " genes")
    message("Gene sets: ", length(gene.sets.list), " pathways loaded")
  }

  # ── Internal function: GSEA for one lineage ──
  .run_one_lineage <- function(lin) {
    tryCatch({
      # Subset cells with finite pseudotime
      pt_vals <- seu@meta.data[[lin]]
      valid_idx <- which(is.finite(pt_vals))

      if (length(valid_idx) < 30) {
        if (verbose) message("  Skipping ", lin, ": only ",
                             length(valid_idx), " cells.")
        return(NULL)
      }

      valid_cells <- colnames(seu)[valid_idx]
      pseudotime  <- pt_vals[valid_idx]

      if (verbose) message("  [", lin, "] ", length(valid_cells), " cells")

      # Extract expression matrix
      expr_mat <- as.matrix(
        SeuratObject::GetAssayData(seu, assay = assay, layer = layer)[
          features, valid_cells, drop = FALSE
        ]
      )

      # Filter genes by min.pct
      pct_expr <- rowMeans(expr_mat > 0)
      genes_keep <- names(pct_expr[pct_expr >= min.pct])

      if (length(genes_keep) == 0) {
        if (verbose) message("    No genes pass min.pct filter.")
        return(NULL)
      }

      expr_mat <- expr_mat[genes_keep, , drop = FALSE]
      if (verbose) message("    ", length(genes_keep),
                           " genes after min.pct filter")

      # Compute gene ranking
      if (ranking.metric == "rho") {
        # Fast vectorized Spearman correlation
        rho_vec <- stats::cor(
          t(expr_mat), pseudotime, method = "spearman"
        )[, 1]
        rho_vec <- rho_vec[is.finite(rho_vec)]

        if (length(rho_vec) == 0) {
          if (verbose) message("    No genes with finite correlation.")
          return(NULL)
        }

        gene_rank <- sort(rho_vec, decreasing = TRUE)

      } else {
        # signed_logp: need per-gene cor.test for p-values
        cor_fun <- function(gene_expr) {
          ct <- suppressWarnings(
            stats::cor.test(gene_expr, pseudotime, method = "spearman",
                            exact = FALSE)
          )
          c(rho = unname(ct$estimate), pvalue = ct$p.value)
        }

        if (requireNamespace("pbapply", quietly = TRUE) && verbose) {
          cor_results <- pbapply::pbapply(expr_mat, 1, cor_fun)
        } else {
          cor_results <- apply(expr_mat, 1, cor_fun)
        }

        rho_vals <- cor_results["rho", ]
        p_vals   <- cor_results["pvalue", ]

        # Keep finite values
        keep <- is.finite(rho_vals) & is.finite(p_vals)
        rho_vals <- rho_vals[keep]
        p_vals   <- p_vals[keep]

        if (length(rho_vals) == 0) {
          if (verbose) message("    No genes with finite correlation.")
          return(NULL)
        }

        # -log10(p) * sign(rho), cap p at machine minimum
        safe_p <- pmax(p_vals, .Machine$double.xmin)
        signed_stat <- -log10(safe_p) * sign(rho_vals)
        gene_rank <- sort(signed_stat, decreasing = TRUE)
      }

      if (verbose) message("    Running GSEA (", ranking.metric,
                           " ranking)...")

      # Run clusterProfiler::GSEA
      gsea_res <- clusterProfiler::GSEA(
        gene_rank,
        TERM2GENE     = term2gene,
        pvalueCutoff  = pvalue.cutoff,
        eps           = eps,
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        verbose       = FALSE
      )
      gsea_df <- as.data.frame(gsea_res)

      if (nrow(gsea_df) == 0) {
        if (verbose) message("    No pathways enriched.")
        return(NULL)
      }

      # Standardized output
      data.frame(
        Celltype         = lin,
        Lineage          = lin,
        Pathway          = gsea_df$Description,
        Score            = gsea_df$NES,
        PValue           = gsea_df$pvalue,
        AdjPValue        = gsea_df$p.adjust,
        QValue           = gsea_df$qvalue,
        Sign             = ifelse(gsea_df$NES > 0, "Activated", "Repressed"),
        SetSize          = gsea_df$setSize,
        EnrichmentScore  = gsea_df$enrichmentScore,
        CoreEnrichment   = gsea_df$core_enrichment,
        RankingMetric    = ranking.metric,
        NCells           = length(valid_cells),
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      if (verbose) message("  Error in ", lin, ": ", conditionMessage(e))
      return(NULL)
    })
  }

  # ── Execute across lineages ──
  results <- lapply(lineages, .run_one_lineage)
  result_df <- do.call(rbind, results)

  if (is.null(result_df) || nrow(result_df) == 0) {
    warning("No pathways found for any lineage.")
    return(data.frame())
  }

  result_df$Method <- "GSEA"
  rownames(result_df) <- NULL

  if (verbose) {
    n_sig <- sum(result_df$AdjPValue < 0.05, na.rm = TRUE)
    message("Done. Total results: ", nrow(result_df),
            " (", n_sig, " with adjP < 0.05)")
  }

  return(result_df)
}


# ============================================================================
# RunTraceGene: Find pseudotime-associated genes
# ============================================================================

#' Find Pseudotime-Associated Genes Along Trajectories
#'
#' Identify genes whose expression significantly changes along pseudotime
#' trajectories. Supports two testing methods: fast Spearman correlation
#' and more rigorous GAM (Generalized Additive Model) fitting.
#'
#' @details
#' For each lineage, the function:
#' \enumerate{
#'   \item Subsets cells with finite pseudotime (\code{is.finite()}).
#'   \item Selects candidate genes (HVG by default, or all / user-specified).
#'   \item Filters genes by minimum expression percentage (\code{min.pct}).
#'   \item Tests each gene for association with pseudotime.
#'   \item Computes effect sizes, peak/valley times, and adjusted p-values.
#' }
#'
#' Two test methods are available:
#' \describe{
#'   \item{\code{"spearman"}}{Computes Spearman correlation between gene
#'     expression and pseudotime. Fast (vectorized rho + per-gene p-value).
#'     Output includes \code{rho} and \code{pvalue}.}
#'   \item{\code{"gam"}}{Fits a GAM smooth spline via \code{mgcv::gam()}
#'     and tests whether the smooth term is significant (F-test).
#'     Output includes \code{r_sq} (R-squared), \code{dev_expl}
#'     (deviance explained), and \code{pvalue}. Additionally computes
#'     Spearman rho for direction. Requires \pkg{mgcv}.}
#' }
#'
#' @param seu A Seurat object with pseudotime stored in \code{meta.data}.
#' @param lineages Character vector of column names in \code{seu@@meta.data}
#'   containing pseudotime values. Cells not on a lineage should have
#'   \code{NA}. Multiple lineages are analyzed independently.
#' @param features Character vector of gene names to test. Default:
#'   \code{NULL} (use top \code{n_candidates} variable features).
#'   Set to \code{"all"} to test all genes.
#' @param n_candidates Number of variable features to select when
#'   \code{features = NULL}. Default: 2000.
#' @param min.pct Minimum fraction of lineage cells expressing a gene
#'   (expression > 0). Default: 0.05.
#' @param test.method Testing method: \code{"spearman"} (fast) or
#'   \code{"gam"} (GAM smooth, more rigorous). Default: \code{"spearman"}.
#' @param assay Seurat assay to use. Default: \code{DefaultAssay(seu)}.
#' @param layer Data layer to extract. Default: \code{"data"}.
#' @param smooth_k Number of basis dimensions for GAM smooth term
#'   (\code{mgcv::s(k = ...)}). Only used when \code{test.method = "gam"}.
#'   Smaller values = smoother curves. Default: 10.
#' @param family Distribution family for GAM. Default: \code{"gaussian"}.
#'   Use \code{"nb"} for raw count data (requires \code{layer = "counts"}).
#' @param padjust.method Method for p-value adjustment. Default: \code{"fdr"}.
#' @param verbose Logical; print progress messages. Default: \code{TRUE}.
#'
#' @return A data.frame with one row per gene-lineage combination:
#'   \describe{
#'     \item{\code{gene}}{Gene name.}
#'     \item{\code{lineage}}{Lineage name.}
#'     \item{\code{rho}}{Spearman correlation coefficient (positive =
#'       expression increases along pseudotime).}
#'     \item{\code{pvalue}}{Raw p-value (Spearman test or GAM F-test).}
#'     \item{\code{padjust}}{Adjusted p-value.}
#'     \item{\code{r_sq}}{R-squared (GAM only; NA for Spearman).}
#'     \item{\code{dev_expl}}{Deviance explained (GAM only; NA for Spearman).}
#'     \item{\code{peaktime}}{Pseudotime of peak expression (median of
#'       top 1\% fitted/raw values).}
#'     \item{\code{valleytime}}{Pseudotime of valley expression (median of
#'       bottom 1\% fitted/raw values).}
#'     \item{\code{exp_ncells}}{Number of cells expressing the gene.}
#'     \item{\code{direction}}{\code{"up"} (rho > 0) or \code{"down"}
#'       (rho < 0).}
#'     \item{\code{test_method}}{Testing method used.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#'
#' # Fast Spearman correlation (default)
#' res <- RunTraceGene(seu, lineages = c("Lineage1", "Lineage2"))
#' sig <- res[res$padjust < 0.05, ]
#' head(sig[order(sig$padjust), ], 20)
#'
#' # GAM fitting (more rigorous)
#' res_gam <- RunTraceGene(seu, lineages = "Lineage1",
#'                          test.method = "gam")
#' sig_gam <- res_gam[res_gam$padjust < 0.05, ]
#'
#' # Visualize top up-regulated genes
#' up_genes <- head(sig$gene[sig$direction == "up"], 10)
#' PlotDynamicFeatures(seu, pseudotime = "Lineage1", features = up_genes)
#'
#' # Test all genes
#' res_all <- RunTraceGene(seu, lineages = "Lineage1", features = "all")
#' }
#'
#' @seealso \code{\link{RunTraceGSEA}}, \code{\link{PlotDynamicFeatures}}
#' @export
RunTraceGene <- function(seu,
                         lineages,
                         features       = NULL,
                         n_candidates   = 2000,
                         min.pct        = 0.05,
                         test.method    = c("spearman", "gam"),
                         assay          = NULL,
                         layer          = "data",
                         smooth_k       = 10,
                         family         = "gaussian",
                         padjust.method = "fdr",
                         verbose        = TRUE) {

  test.method <- match.arg(test.method)

  # ── Input validation ──
  if (!inherits(seu, "Seurat")) {
    stop("'seu' must be a Seurat object.")
  }

  missing_lin <- lineages[!lineages %in% colnames(seu@meta.data)]
  if (length(missing_lin) > 0) {
    stop("Lineage column(s) not found in meta.data: ",
         paste(missing_lin, collapse = ", "))
  }

  if (test.method == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for test.method = 'gam'.\n",
           "Install with: install.packages('mgcv')")
    }
  }

  # ── Resolve assay ──
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seu)
  }
  all_genes <- rownames(seu[[assay]])

  # ── Determine feature universe ──
  if (is.null(features)) {
    # Default: use HVG (gene discovery benefits from focused set)
    vf <- SeuratObject::VariableFeatures(seu, assay = assay)
    if (length(vf) < n_candidates) {
      if (verbose) message("Finding top ", n_candidates, " variable features...")
      seu_tmp <- Seurat::FindVariableFeatures(seu, nfeatures = n_candidates,
                                               assay = assay, verbose = FALSE)
      features <- SeuratObject::VariableFeatures(seu_tmp, assay = assay)
    } else {
      features <- vf[seq_len(n_candidates)]
    }
  } else if (length(features) == 1 && features == "all") {
    features <- all_genes
  }
  features <- intersect(features, all_genes)
  if (length(features) == 0) {
    stop("No valid features found in assay '", assay, "'.")
  }

  if (verbose) {
    message("=== RunTraceGene ===")
    message("Lineages: ", paste(lineages, collapse = ", "))
    message("Test method: ", test.method)
    message("Candidate features: ", length(features), " genes")
  }

  # ── Internal: test one lineage ──
  .test_one_lineage <- function(lin) {
    tryCatch({
      # Subset cells with finite pseudotime
      pt_vals <- seu@meta.data[[lin]]
      valid_idx <- which(is.finite(pt_vals))

      if (length(valid_idx) < 30) {
        if (verbose) message("  Skipping ", lin, ": only ",
                             length(valid_idx), " cells.")
        return(NULL)
      }

      valid_cells <- colnames(seu)[valid_idx]
      pseudotime  <- pt_vals[valid_idx]

      if (verbose) message("  [", lin, "] ", length(valid_cells), " cells")

      # Extract expression matrix
      expr_mat <- as.matrix(
        SeuratObject::GetAssayData(seu, assay = assay, layer = layer)[
          features, valid_cells, drop = FALSE
        ]
      )

      # Filter genes by min.pct
      pct_expr <- rowMeans(expr_mat > 0)
      genes_keep <- names(pct_expr[pct_expr >= min.pct])

      if (length(genes_keep) == 0) {
        if (verbose) message("    No genes pass min.pct filter.")
        return(NULL)
      }

      expr_mat <- expr_mat[genes_keep, , drop = FALSE]
      if (verbose) message("    ", length(genes_keep),
                           " genes after min.pct filter")

      # ── Spearman method ──
      if (test.method == "spearman") {
        if (verbose) message("    Computing Spearman correlations...")

        # Vectorized rho
        rho_vec <- stats::cor(
          t(expr_mat), pseudotime, method = "spearman"
        )[, 1]

        # Per-gene p-values
        cor_fun <- function(gene_expr) {
          ct <- suppressWarnings(
            stats::cor.test(gene_expr, pseudotime, method = "spearman",
                            exact = FALSE)
          )
          ct$p.value
        }

        if (requireNamespace("pbapply", quietly = TRUE) && verbose) {
          p_vec <- pbapply::pbapply(expr_mat, 1, cor_fun)
        } else {
          p_vec <- apply(expr_mat, 1, cor_fun)
        }

        # Peak/valley from raw expression
        peak_valley <- t(apply(expr_mat, 1, function(y) {
          peaktime <- stats::median(
            pseudotime[y >= stats::quantile(y, 0.99, na.rm = TRUE)]
          )
          valleytime <- stats::median(
            pseudotime[y <= stats::quantile(y, 0.01, na.rm = TRUE)]
          )
          c(peaktime = peaktime, valleytime = valleytime)
        }))

        # Build result
        res_df <- data.frame(
          gene       = genes_keep,
          lineage    = lin,
          rho        = rho_vec[genes_keep],
          pvalue     = p_vec[genes_keep],
          r_sq       = NA_real_,
          dev_expl   = NA_real_,
          peaktime   = peak_valley[, "peaktime"],
          valleytime = peak_valley[, "valleytime"],
          exp_ncells = as.integer(rowSums(expr_mat > 0)),
          stringsAsFactors = FALSE
        )

      } else {
        # ── GAM method ──
        if (verbose) message("    Fitting GAM for each gene...")

        # Vectorized rho for direction
        rho_vec <- stats::cor(
          t(expr_mat), pseudotime, method = "spearman"
        )[, 1]

        # Per-gene GAM fitting
        gam_fun <- function(gene_name) {
          y <- expr_mat[gene_name, ]
          fit <- tryCatch(
            .fit_dynamic_feature(
              y = y, x = pseudotime,
              fit_method = "gam",
              smooth_k = smooth_k,
              loess_span = 0.75,
              bspline_knot = 3,
              family = family
            ),
            error = function(e) NULL
          )

          if (is.null(fit)) {
            return(c(pvalue = NA_real_, r_sq = NA_real_, dev_expl = NA_real_,
                     peaktime = NA_real_, valleytime = NA_real_))
          }

          # dev.expl from GAM summary (refit to get it)
          dev_expl <- tryCatch({
            k_use <- min(smooth_k, length(y) - 1)
            mod <- mgcv::gam(
              y_val ~ s(x_val, bs = "cs", k = k_use),
              family = family,
              data = data.frame(y_val = y[order(pseudotime)],
                                x_val = sort(pseudotime))
            )
            summary(mod)$dev.expl
          }, error = function(e) NA_real_)

          # Peak/valley from fitted values
          fitted_vals <- fit$curve$fitted
          pt_ordered  <- fit$curve$pseudotime
          peaktime <- stats::median(
            pt_ordered[fitted_vals >= stats::quantile(fitted_vals, 0.99,
                                                      na.rm = TRUE)]
          )
          valleytime <- stats::median(
            pt_ordered[fitted_vals <= stats::quantile(fitted_vals, 0.01,
                                                      na.rm = TRUE)]
          )

          c(pvalue = fit$pvalue, r_sq = fit$r_sq, dev_expl = dev_expl,
            peaktime = peaktime, valleytime = valleytime)
        }

        if (requireNamespace("pbapply", quietly = TRUE) && verbose) {
          gam_results <- pbapply::pbsapply(genes_keep, gam_fun)
        } else {
          gam_results <- sapply(genes_keep, gam_fun)
        }

        # Build result
        res_df <- data.frame(
          gene       = genes_keep,
          lineage    = lin,
          rho        = rho_vec[genes_keep],
          pvalue     = gam_results["pvalue", ],
          r_sq       = gam_results["r_sq", ],
          dev_expl   = gam_results["dev_expl", ],
          peaktime   = gam_results["peaktime", ],
          valleytime = gam_results["valleytime", ],
          exp_ncells = as.integer(rowSums(expr_mat > 0)),
          stringsAsFactors = FALSE
        )
      }

      # Filter out genes with NA p-values
      res_df <- res_df[is.finite(res_df$pvalue), , drop = FALSE]

      if (nrow(res_df) == 0) return(NULL)

      # Direction and test method
      res_df$direction   <- ifelse(res_df$rho > 0, "up", "down")
      res_df$test_method <- test.method

      res_df

    }, error = function(e) {
      if (verbose) message("  Error in ", lin, ": ", conditionMessage(e))
      return(NULL)
    })
  }

  # ── Execute across lineages ──
  results <- lapply(lineages, .test_one_lineage)
  result_df <- do.call(rbind, results)

  if (is.null(result_df) || nrow(result_df) == 0) {
    warning("No significant genes found for any lineage.")
    return(data.frame())
  }

  # Adjust p-values per lineage
  result_df$padjust <- ave(
    result_df$pvalue,
    result_df$lineage,
    FUN = function(p) stats::p.adjust(p, method = padjust.method)
  )

  # Sort by p-value
  result_df <- result_df[order(result_df$lineage, result_df$pvalue), ]
  rownames(result_df) <- NULL

  # Reorder columns
  col_order <- c("gene", "lineage", "rho", "pvalue", "padjust",
                  "r_sq", "dev_expl", "peaktime", "valleytime",
                  "exp_ncells", "direction", "test_method")
  result_df <- result_df[, col_order]

  if (verbose) {
    n_sig <- sum(result_df$padjust < 0.05, na.rm = TRUE)
    message("Done. Total genes tested: ", nrow(result_df),
            " (", n_sig, " with padjust < 0.05)")
  }

  return(result_df)
}
