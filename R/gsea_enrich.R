# ============================================================================
# RunGseaEnrich: compareCluster-based enricher + GSEA analysis
# ============================================================================


#' Parse geneset input for RunGseaEnrich
#'
#' Converts various geneset input formats into a TERM2GENE data.frame
#' with columns \code{gs_name} and \code{gene_symbol}.
#'
#' @param geneset Gene set specification. One of:
#'   \describe{
#'     \item{msigdbr args}{A list with msigdbr parameters (e.g.
#'       \code{list(species = "Homo sapiens", collection = "H")}).
#'       Recognized keys: \code{species}, \code{collection},
#'       \code{subcollection}, \code{category}, \code{subcategory}.}
#'     \item{data.frame}{A data.frame with columns \code{gs_name} and
#'       \code{gene_symbol} (TERM2GENE format).}
#'     \item{named list}{A named list of character vectors, e.g.
#'       \code{list(PathA = c("TP53","BRCA1"))}.}
#'     \item{GMT file path}{A character string path to a GMT file.}
#'   }
#' @return A data.frame with columns \code{gs_name} and \code{gene_symbol}.
#' @keywords internal
.parse_geneset_to_term2gene <- function(geneset) {
  # 1) data.frame with gs_name + gene_symbol

if (is.data.frame(geneset)) {
    if (!all(c("gs_name", "gene_symbol") %in% colnames(geneset))) {
      # Try first two columns as fallback (term, gene)
      if (ncol(geneset) >= 2) {
        geneset <- data.frame(
          gs_name = as.character(geneset[[1]]),
          gene_symbol = as.character(geneset[[2]]),
          stringsAsFactors = FALSE
        )
      } else {
        stop("geneset data.frame must contain 'gs_name' and 'gene_symbol' columns,",
             " or at least 2 columns (term, gene).")
      }
    }
    return(dplyr::distinct(geneset, .data$gs_name, .data$gene_symbol))
  }

  # 2) Character string -> GMT file path
  if (is.character(geneset) && length(geneset) == 1) {
    gs_list <- read_gmt(geneset)
    return(data.frame(
      gs_name = rep(names(gs_list), vapply(gs_list, length, integer(1))),
      gene_symbol = unlist(gs_list, use.names = FALSE),
      stringsAsFactors = FALSE
    ))
  }

  # 3) Named list
  if (is.list(geneset) && !is.null(names(geneset))) {
    msigdbr_keys <- c("species", "collection", "subcollection",
                      "category", "subcategory")
    if (any(names(geneset) %in% msigdbr_keys)) {
      # msigdbr args
      if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("Package 'msigdbr' is required when geneset is msigdbr args.\n",
             "Install with: install.packages('msigdbr')")
      }
      gs_df <- rlang::exec(msigdbr::msigdbr, !!!geneset)
      return(dplyr::distinct(gs_df, .data$gs_name, .data$gene_symbol))
    } else {
      # Named list of gene vectors
      return(data.frame(
        gs_name = rep(names(geneset),
                      vapply(geneset, length, integer(1))),
        gene_symbol = unlist(geneset, use.names = FALSE),
        stringsAsFactors = FALSE
      ))
    }
  }

  stop("geneset must be one of:\n",
       "  1) list(species=..., collection=...) for msigdbr\n",
       "  2) data.frame with gs_name + gene_symbol columns\n",
       "  3) named list, e.g. list(pathway = c('gene1','gene2'))\n",
       "  4) character path to a GMT file")
}


#' Run compareCluster-based Enricher and GSEA Analysis
#'
#' Perform enricher (ORA) and GSEA enrichment analysis across cell types and
#' groups using \code{clusterProfiler::compareCluster()}. Requires
#' \code{FindAllMarkers()} results grouped by cell type and condition.
#'
#' @param scobj A Seurat object.
#' @param celltype Column name in \code{meta.data} for cell type labels.
#'   Default: \code{"celltype"}.
#' @param group Column name in \code{meta.data} for condition/group labels.
#'   Default: \code{"group"}.
#' @param allMarkers.args List of arguments passed to
#'   \code{Seurat::FindAllMarkers()}.
#' @param geneset Gene set specification. Supports multiple formats:
#'   \describe{
#'     \item{msigdbr args}{A list with msigdbr parameters, e.g.
#'       \code{list(species = "Homo sapiens", collection = "H")}.}
#'     \item{data.frame}{With columns \code{gs_name} and \code{gene_symbol}.}
#'     \item{named list}{E.g. \code{list(PathA = c("TP53", "BRCA1"))}.}
#'     \item{GMT file}{Character path to a \code{.gmt} file.}
#'   }
#' @param padj_cutoff P-value adjusted cutoff for enricher ORA analysis.
#'   Default: 0.001.
#' @param plot Logical; whether to generate dotplots. Default: \code{TRUE}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{markers}}{Data.frame from \code{FindAllMarkers()}, with
#'       added \code{celltype} and \code{group} columns.}
#'     \item{\code{enricher}}{\code{compareClusterResult} from enricher.}
#'     \item{\code{gsea}}{\code{compareClusterResult} from GSEA.}
#'     \item{\code{plots}}{(if \code{plot = TRUE}) List of dotplot lists.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(scMMR)
#'
#' # Using msigdbr Hallmark gene sets (default)
#' res <- RunGseaEnrich(seu)
#'
#' # Using custom column names
#' res <- RunGseaEnrich(seu, celltype = "cell_type", group = "condition")
#'
#' # Using a data.frame as geneset
#' my_gs <- data.frame(gs_name = c("PathA","PathA","PathB"),
#'                     gene_symbol = c("TP53","BRCA1","EGFR"))
#' res <- RunGseaEnrich(seu, geneset = my_gs)
#'
#' # Using a named list
#' res <- RunGseaEnrich(seu,
#'   geneset = list(PathA = c("TP53","BRCA1"), PathB = c("EGFR","KRAS")))
#'
#' # Using a GMT file
#' gmt <- system.file("extdata", "gmt",
#'   "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")
#' res <- RunGseaEnrich(seu, geneset = gmt)
#' }
#'
#' @export
RunGseaEnrich <- function(
    scobj,
    celltype = "celltype",
    group = "group",
    allMarkers.args = list(only.pos = FALSE, logfc.threshold = 0.1,
                           min.pct = 0.25, thresh.use = 0.25),
    geneset = list(species = "Homo sapiens", collection = "H",
                   subcollection = NULL),
    padj_cutoff = 0.001,
    plot = TRUE) {

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.\n",
         "Install with: BiocManager::install('clusterProfiler')")
  }

  # ── Validate inputs ──
  if (!inherits(scobj, "Seurat")) {
    stop("'scobj' must be a Seurat object.")
  }
  if (!celltype %in% colnames(scobj@meta.data)) {
    stop(sprintf("Column '%s' not found in meta.data.", celltype))
  }
  if (!group %in% colnames(scobj@meta.data)) {
    stop(sprintf("Column '%s' not found in meta.data.", group))
  }

  # ── Create composite identity ──
  scobj$..enrich_ct_grp.. <- paste(scobj@meta.data[[celltype]],
                                   scobj@meta.data[[group]], sep = "_")
  Seurat::Idents(scobj) <- "..enrich_ct_grp.."

  # ── 1. FindAllMarkers ──
  message("Step 1/3: Running FindAllMarkers...")
  allMarkers <- rlang::exec(Seurat::FindAllMarkers,
                            object = scobj,
                            !!!allMarkers.args) %>%
    tidyr::separate(.data$cluster, into = c("celltype", "group"),
                    sep = "_", remove = FALSE)

  # ── 2. Parse gene sets ──
  message("Step 2/3: Parsing gene sets...")
  geneSet <- .parse_geneset_to_term2gene(geneset)
  message(sprintf("  -> %d term-gene pairs, %d unique terms.",
                  nrow(geneSet), length(unique(geneSet$gs_name))))

  # ── 3. Enrichment analysis ──
  message("Step 3/3: Running enrichment analysis...")

  # Enricher (ORA)
  message("  Running enricher (ORA)...")
  enricher_data <- allMarkers %>% dplyr::filter(.data$p_val_adj < padj_cutoff)
  if (nrow(enricher_data) == 0) {
    warning("No markers with p_val_adj < ", padj_cutoff,
            ". Enricher will be skipped. Consider increasing padj_cutoff.")
    enricher_result <- NULL
  } else {
    enricher_result <- clusterProfiler::compareCluster(
      gene ~ celltype + group,
      data = enricher_data,
      fun = "enricher",
      TERM2GENE = geneSet
    )
  }

  # GSEA
  message("  Running GSEA...")
  gsea_result <- clusterProfiler::compareCluster(
    gene | avg_log2FC ~ celltype + group,
    data = allMarkers,
    fun = "GSEA",
    TERM2GENE = geneSet
  )

  # ── 4. Plots ──
  plots <- NULL
  if (plot) {
    message("Generating plots...")
    plots <- .build_gsea_enrich_plots(enricher_result, gsea_result)
  }

  message("Done.")

  result <- list(
    markers = allMarkers,
    enricher = enricher_result,
    gsea = gsea_result
  )
  if (!is.null(plots)) result$plots <- plots

  return(result)
}


# ============================================================================
# FindRegulonDrivers: identify upstream TF regulators of pathways
# ============================================================================

#' Find Upstream Regulon Drivers of Pathways
#'
#' Identify transcription factor regulons that may drive pathway activity by
#' integrating two complementary metrics:
#' \enumerate{
#'   \item \strong{Correlation}: Pearson/Spearman correlation between regulon
#'     activity scores and pathway scores at single-cell level.
#'   \item \strong{Gene overlap}: intersection of regulon target genes with
#'     pathway gene sets, plus fraction of hits (\code{n_overlap / regulon_size}).
#' }
#'
#' Internally uses \code{\link{ComputeModuleScore}} to score both regulon and
#' pathway gene sets in the Seurat object.
#'
#' @param seu A Seurat object.
#' @param pathway.genes Pathway gene sets. Accepts:
#'   \itemize{
#'     \item Named list of character vectors.
#'     \item GMT file path.
#'     \item data.frame with columns (term, gene).
#'   }
#' @param regulons Regulon target gene sets. Same formats as
#'   \code{pathway.genes} (typically from SCENIC \code{regulons.rds}).
#' @param pathways Character vector of specific pathways to analyze.
#'   Default: \code{NULL} (all pathways).
#' @param score.method Scoring method passed to \code{ComputeModuleScore}:
#'   \code{"AUCell"}, \code{"Seurat"}, or \code{"UCell"}.
#'   Default: \code{"AUCell"}.
#' @param min.size Minimum gene set size for scoring. Default: 5.
#' @param cor.method Correlation method: \code{"pearson"} or \code{"spearman"}.
#'   Default: \code{"pearson"}.
#' @param min.overlap Minimum number of overlapping genes to report.
#'   Default: 1.
#' @param label.n Number of top/bottom regulons to label in PCC rank bar plot.
#'   Default: 5.
#' @param label.overlap Overlap count threshold for labeling in scatter plot.
#'   Default: \code{NULL} (auto: top 10\% quantile).
#' @param plot Logical; generate plots. Default: \code{TRUE}.
#' @param cores Number of parallel cores for scoring. Default: 1.
#' @param verbose Logical; print progress. Default: \code{TRUE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{result}}{Data.frame with columns: \code{regulon}, \code{PCC},
#'       \code{rank}, \code{n_overlap}, \code{regulon_size},
#'       \code{frac_of_hits}, \code{pathway}. One row per regulon × pathway.}
#'     \item{\code{cor.mat}}{Correlation matrix (pathways × regulons).}
#'     \item{\code{plots}}{(if \code{plot = TRUE}) List of ggplot objects:
#'       \code{rank_bar} (PCC rank) and \code{scatter} (PCC vs overlap),
#'       one set per pathway.}
#'   }
#'
#' @examples
#' \dontrun{
#' regulons <- readRDS("regulons.rds")
#' hallmark <- system.file("extdata", "gmt",
#'   "h.all.v2022.1.Hs.symbols.gmt", package = "scMMR")
#'
#' res <- FindRegulonDrivers(seu, pathway.genes = hallmark, regulons = regulons)
#' head(res$result)
#' res$plots[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]$scatter
#'
#' # Focus on specific pathways
#' res2 <- FindRegulonDrivers(seu,
#'   pathway.genes = hallmark, regulons = regulons,
#'   pathways = c("HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"))
#' }
#'
#' @seealso \code{\link{ComputeModuleScore}}, \code{\link{CrossEnrichOverlap}}
#' @importFrom dplyr arrange desc mutate
#' @importFrom rlang .data
#' @import ggplot2
#' @export
FindRegulonDrivers <- function(seu,
                               pathway.genes,
                               regulons,
                               pathways = NULL,
                               score.method = c("AUCell", "Seurat", "UCell"),
                               min.size = 5,
                               cor.method = c("pearson", "spearman"),
                               min.overlap = 1,
                               label.n = 5,
                               label.overlap = NULL,
                               plot = TRUE,
                               cores = 1,
                               verbose = TRUE) {

  score.method <- match.arg(score.method)
  cor.method <- match.arg(cor.method)

  if (!inherits(seu, "Seurat")) stop("'seu' must be a Seurat object.")

  # ── Parse gene sets ──
  pw_list <- parse_gene_sets(pathway.genes)
  reg_list <- parse_gene_sets(regulons)

  # Clean regulon names: remove "(+)" suffix from SCENIC
  names(reg_list) <- sub("\\(\\+\\)$", "", names(reg_list))

  if (!is.null(pathways)) {
    found <- pathways %in% names(pw_list)
    if (!any(found)) stop("None of the specified pathways found in pathway.genes.")
    if (any(!found) && verbose) {
      message("Warning: pathways not found: ",
              paste(pathways[!found], collapse = ", "))
    }
    pw_list <- pw_list[pathways[found]]
  }

  if (verbose) {
    message("=== FindRegulonDrivers ===")
    message("Pathways: ", length(pw_list))
    message("Regulons: ", length(reg_list))
    message("Score method: ", score.method)
  }

  # ── Score pathways ──
  if (verbose) message("Step 1/3: Scoring pathway gene sets...")
  dge <- Seurat::GetAssayData(seu, layer = if (score.method == "Seurat") "data" else "counts")
  pw_mat <- ComputeModuleScore(dge, gene.sets = pw_list,
                               method = score.method, min.size = min.size,
                               cores = cores)
  # pw_mat: gene_sets × cells -> transpose to cells × gene_sets
  pw_mat <- t(pw_mat)

  # ── Score regulons ──
  if (verbose) message("Step 2/3: Scoring regulon gene sets...")
  reg_mat <- ComputeModuleScore(dge, gene.sets = reg_list,
                                method = score.method, min.size = min.size,
                                cores = cores)
  reg_mat <- t(reg_mat)

  # Update lists to only scored sets
  pw_list <- pw_list[colnames(pw_mat)]
  reg_list <- reg_list[colnames(reg_mat)]

  if (verbose) {
    message(sprintf("  Scored: %d pathways, %d regulons",
                    ncol(pw_mat), ncol(reg_mat)))
  }

  # ── Compute correlation matrix ──
  if (verbose) message("Step 3/3: Computing correlation & overlap...")
  # Ensure same cells
  common_cells <- intersect(rownames(pw_mat), rownames(reg_mat))
  cor_mat <- stats::cor(pw_mat[common_cells, , drop = FALSE],
                        reg_mat[common_cells, , drop = FALSE],
                        method = cor.method)
  # cor_mat: pathways × regulons

  # ── Compute overlap for each pathway × regulon ──
  records <- list()
  for (pw_name in rownames(cor_mat)) {
    pw_genes_vec <- pw_list[[pw_name]]
    # PCC values for this pathway
    pcc_vec <- cor_mat[pw_name, ]
    pcc_ranked <- sort(pcc_vec, decreasing = TRUE)
    rank_map <- stats::setNames(seq_along(pcc_ranked), names(pcc_ranked))

    for (reg_name in colnames(cor_mat)) {
      reg_genes_vec <- reg_list[[reg_name]]
      common <- intersect(reg_genes_vec, pw_genes_vec)
      n_ov <- length(common)

      if (n_ov >= min.overlap || TRUE) {
        # Always record (for correlation), filter later for overlap
        records[[length(records) + 1L]] <- data.frame(
          regulon = reg_name,
          pathway = pw_name,
          PCC = round(pcc_vec[reg_name], 4),
          rank = as.integer(rank_map[reg_name]),
          n_overlap = n_ov,
          regulon_size = length(reg_genes_vec),
          pathway_size = length(pw_genes_vec),
          frac_of_hits = round(n_ov / length(reg_genes_vec), 4),
          overlap_genes = if (n_ov > 0) paste(common, collapse = ",") else "",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  result_df <- do.call(rbind, records) %>%
    dplyr::arrange(.data$pathway, dplyr::desc(.data$PCC))

  if (verbose) {
    n_with_ov <- sum(result_df$n_overlap >= min.overlap)
    message(sprintf("Done. %d regulon-pathway pairs total, %d with >= %d overlap.",
                    nrow(result_df), n_with_ov, min.overlap))
  }

  out <- list(
    result = result_df,
    cor.mat = cor_mat
  )

  # ── Plots (per pathway) ──
  if (plot) {
    out$plots <- list()
    for (pw_name in unique(result_df$pathway)) {
      df_pw <- result_df %>% dplyr::filter(.data$pathway == pw_name)
      out$plots[[pw_name]] <- .plot_regulon_drivers(
        df_pw, pw_name, label.n, label.overlap
      )
    }
  }

  return(out)
}


#' Generate plots for one pathway in FindRegulonDrivers
#'
#' @param df Data.frame for one pathway (from result_df).
#' @param pw_name Pathway name.
#' @param label.n Top N for rank barplot.
#' @param label.overlap Overlap threshold for scatter labeling.
#' @return List with rank_bar and scatter plots.
#' @keywords internal
.plot_regulon_drivers <- function(df, pw_name, label.n = 5,
                                  label.overlap = NULL) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    warning("Package 'ggrepel' not found. Scatter labels will use geom_text.")
    has_repel <- FALSE
  } else {
    has_repel <- TRUE
  }

  plots <- list()

  # ── 1. PCC rank barplot (top N + bottom N) ──
  df_sorted <- df %>% dplyr::arrange(dplyr::desc(.data$PCC))
  df_sorted$rank <- seq_len(nrow(df_sorted))

  top_bottom <- rbind(
    utils::head(df_sorted, label.n),
    utils::tail(df_sorted, label.n)
  ) %>% dplyr::distinct()

  top_bottom <- top_bottom %>%
    dplyr::arrange(.data$rank) %>%
    dplyr::mutate(regulon = factor(.data$regulon, levels = rev(.data$regulon)))

  plots$rank_bar <- ggplot2::ggplot(
    top_bottom,
    ggplot2::aes(x = .data$regulon, y = .data$PCC)
  ) +
    ggplot2::geom_col(color = "black", fill = "lightblue") +
    ggplot2::coord_flip() +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("PCC") +
    ggplot2::ggtitle(pw_name) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  # ── 2. PCC vs n_overlap scatter (size = frac_of_hits) ──
  df_ov <- df %>% dplyr::filter(.data$n_overlap > 0)

  if (nrow(df_ov) > 0) {
    # Auto threshold for labeling
    if (is.null(label.overlap)) {
      label.overlap <- stats::quantile(df_ov$n_overlap, 0.90, na.rm = TRUE)
      label.overlap <- max(label.overlap, 1)
    }

    df_label <- df_ov %>% dplyr::filter(.data$n_overlap >= label.overlap)

    p_scatter <- ggplot2::ggplot(
      df_ov,
      ggplot2::aes(x = .data$PCC, y = .data$n_overlap,
                   size = .data$frac_of_hits)
    ) +
      ggplot2::geom_point(shape = 21, color = "black",
                          fill = "lightblue", stroke = 0.5) +
      ggplot2::scale_size_continuous(name = "Frac. of hits",
                                     range = c(1, 6)) +
      ggplot2::labs(
        x = "PCC",
        y = "Number of overlapped targets",
        title = pw_name
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

    if (nrow(df_label) > 0) {
      if (has_repel) {
        p_scatter <- p_scatter +
          ggrepel::geom_text_repel(
            data = df_label,
            ggplot2::aes(label = .data$regulon),
            size = 3.5, show.legend = FALSE
          )
      } else {
        p_scatter <- p_scatter +
          ggplot2::geom_text(
            data = df_label,
            ggplot2::aes(label = .data$regulon),
            size = 3.5, vjust = -1, show.legend = FALSE
          )
      }
    }

    plots$scatter <- p_scatter
  }

  plots
}


#' Build dotplots for RunGseaEnrich results
#'
#' @param enricher_result compareCluster enricher result (or NULL).
#' @param gsea_result compareCluster GSEA result.
#' @return A list of plot lists.
#' @keywords internal
.build_gsea_enrich_plots <- function(enricher_result, gsea_result) {
  dot_theme <- ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

  plots <- list()

  # Enricher plots
  if (!is.null(enricher_result)) {
    plots$enricher <- list(
      by_group = clusterProfiler::dotplot(
        enricher_result, label_format = 60, font.size = 8, x = "group"
      ) + ggplot2::facet_grid(~celltype) + dot_theme,
      by_celltype = clusterProfiler::dotplot(
        enricher_result, label_format = 60, font.size = 8, x = "celltype"
      ) + ggplot2::facet_grid(~group) + dot_theme
    )
  }

  # GSEA basic
  plots$gsea_basic <- list(
    by_group = clusterProfiler::dotplot(
      gsea_result, label_format = 60, font.size = 8, x = "group"
    ) + ggplot2::facet_grid(~celltype) + dot_theme,
    by_celltype = clusterProfiler::dotplot(
      gsea_result, label_format = 60, font.size = 8, x = "celltype"
    ) + ggplot2::facet_grid(~group) + dot_theme
  )

  # GSEA by sign
  plots$gsea_by_sign <- list(
    by_celltype_sign = clusterProfiler::dotplot(
      gsea_result, label_format = 60, font.size = 8, x = "celltype"
    ) + ggplot2::facet_grid(.sign ~ group) + dot_theme,
    by_group_sign = clusterProfiler::dotplot(
      gsea_result, label_format = 60, font.size = 8, x = "group"
    ) + ggplot2::facet_grid(.sign ~ celltype) + dot_theme
  )

  # GSEA combined
  plots$gsea_combined <- list(
    by_celltype_sign_group = clusterProfiler::dotplot(
      gsea_result, label_format = 60, font.size = 8, x = "celltype"
    ) + ggplot2::facet_grid(~ .sign + group) + dot_theme,
    by_group_sign_celltype = clusterProfiler::dotplot(
      gsea_result, label_format = 60, font.size = 8, x = "group"
    ) + ggplot2::facet_grid(~ .sign + celltype) + dot_theme
  )

  plots
}


# ============================================================================
# CrossEnrichOverlap: TF-Pathway core_enrichment overlap
# ============================================================================

#' Extract core_enrichment genes from RunGseaEnrich result
#'
#' @param res A \code{RunGseaEnrich} result list or a
#'   \code{compareClusterResult} object.
#' @param type Which result to extract: \code{"gsea"} or \code{"enricher"}.
#' @param p.cutoff Filter terms by p.adjust cutoff. Default: 0.05.
#' @return A named list: term name -> character vector of core_enrichment genes.
#' @keywords internal
.extract_core_genes <- function(res, type = "gsea", p.cutoff = 0.05) {
  # Accept either RunGseaEnrich list or compareClusterResult directly
  if (is.list(res) && !isS4(res) && type %in% names(res)) {
    obj <- res[[type]]
  } else if (isS4(res)) {
    obj <- res
  } else {
    stop("'res' must be a RunGseaEnrich result list or a compareClusterResult.")
  }

  if (is.null(obj)) stop("No '", type, "' result found.")

  df <- as.data.frame(obj)

  # Determine gene column name (core_enrichment for GSEA, geneID for enricher)
  gene_col <- if ("core_enrichment" %in% colnames(df)) {
    "core_enrichment"
  } else if ("geneID" %in% colnames(df)) {
    "geneID"
  } else {
    stop("Cannot find 'core_enrichment' or 'geneID' column.")
  }

  # Filter by p.adjust
  if ("p.adjust" %in% colnames(df)) {
    df <- df[df$p.adjust < p.cutoff, , drop = FALSE]
  }
  if (nrow(df) == 0) return(list())

  # Build unique key: "Description (celltype|group)" or just "Description"
  has_ct <- "celltype" %in% colnames(df)
  has_grp <- "group" %in% colnames(df)
  if (has_ct && has_grp) {
    df$..key.. <- paste0(df$Description, " [", df$celltype, "|", df$group, "]")
  } else if (has_ct) {
    df$..key.. <- paste0(df$Description, " [", df$celltype, "]")
  } else {
    df$..key.. <- df$Description
  }

  stats::setNames(
    lapply(strsplit(df[[gene_col]], "/"), function(x) x[nzchar(x)]),
    df$..key..
  )
}


#' Cross-comparison of Core Enrichment Genes between Two GSEA Results
#'
#' Compare core enrichment genes between two \code{RunGseaEnrich} results
#' (e.g., TF regulons vs pathways) to identify shared regulatory genes.
#' Returns an overlap summary table and optionally a heatmap or tile plot.
#'
#' @param res_tf TF/GRN \code{RunGseaEnrich} result (or \code{compareClusterResult}).
#' @param res_pathway Pathway \code{RunGseaEnrich} result (or \code{compareClusterResult}).
#' @param type Which sub-result to use: \code{"gsea"} or \code{"enricher"}.
#'   Default: \code{"gsea"}.
#' @param p.cutoff P.adjust cutoff for filtering enriched terms. Default: 0.05.
#' @param min.overlap Minimum number of overlapping genes to report. Default: 1.
#' @param show.genes Logical; include overlapping gene names in result.
#'   Default: \code{TRUE}.
#' @param plot Logical; generate a tile heatmap of overlap counts.
#'   Default: \code{TRUE}.
#' @param top.tf Maximum number of TF terms to display (by overlap count).
#'   Default: 20.
#' @param top.pathway Maximum number of pathway terms to display.
#'   Default: 20.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{overlap}}{Data.frame with columns: \code{TF}, \code{Pathway},
#'       \code{n_overlap}, \code{overlap_genes}, \code{n_tf_genes},
#'       \code{n_pathway_genes}, \code{jaccard}.}
#'     \item{\code{tf_genes}}{Named list of core_enrichment genes per TF term.}
#'     \item{\code{pathway_genes}}{Named list of core_enrichment genes per
#'       pathway term.}
#'     \item{\code{plot}}{(if \code{plot = TRUE}) A ggplot tile plot of overlap
#'       counts.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compare TF regulon enrichment with pathway enrichment
#' ov <- CrossEnrichOverlap(res_PH_TF, res_PH)
#' head(ov$overlap)
#' ov$plot
#'
#' # Only enricher results, stricter cutoff
#' ov2 <- CrossEnrichOverlap(res_TF, res_pathway,
#'   type = "enricher", p.cutoff = 0.01, min.overlap = 3)
#' }
#'
#' @importFrom dplyr filter arrange desc
#' @import ggplot2
#' @export
CrossEnrichOverlap <- function(res_tf,
                               res_pathway,
                               type = c("gsea", "enricher"),
                               p.cutoff = 0.05,
                               min.overlap = 1,
                               show.genes = TRUE,
                               plot = TRUE,
                               top.tf = 20,
                               top.pathway = 20) {
  type <- match.arg(type)

  # Extract core_enrichment gene lists
  tf_genes <- .extract_core_genes(res_tf, type = type, p.cutoff = p.cutoff)
  pw_genes <- .extract_core_genes(res_pathway, type = type, p.cutoff = p.cutoff)

  if (length(tf_genes) == 0) stop("No enriched TF terms found (p.adjust < ", p.cutoff, ").")
  if (length(pw_genes) == 0) stop("No enriched pathway terms found (p.adjust < ", p.cutoff, ").")

  message(sprintf("TF terms: %d, Pathway terms: %d", length(tf_genes), length(pw_genes)))

  # Compute pairwise overlap
  records <- list()
  for (tf_name in names(tf_genes)) {
    tg <- tf_genes[[tf_name]]
    for (pw_name in names(pw_genes)) {
      pg <- pw_genes[[pw_name]]
      common <- intersect(tg, pg)
      if (length(common) >= min.overlap) {
        union_size <- length(union(tg, pg))
        records[[length(records) + 1L]] <- data.frame(
          TF = tf_name,
          Pathway = pw_name,
          n_overlap = length(common),
          overlap_genes = if (show.genes) paste(common, collapse = ",") else NA_character_,
          n_tf_genes = length(tg),
          n_pathway_genes = length(pg),
          jaccard = round(length(common) / union_size, 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(records) == 0) {
    message("No overlapping genes found between TF and pathway core_enrichment.")
    return(list(
      overlap = data.frame(),
      tf_genes = tf_genes,
      pathway_genes = pw_genes
    ))
  }

  overlap_df <- do.call(rbind, records) %>%
    dplyr::arrange(dplyr::desc(.data$n_overlap))

  message(sprintf("Found %d TF-Pathway pairs with >= %d overlapping genes.",
                  nrow(overlap_df), min.overlap))

  result <- list(
    overlap = overlap_df,
    tf_genes = tf_genes,
    pathway_genes = pw_genes
  )

  # Plot
  if (plot && nrow(overlap_df) > 0) {
    # Select top TFs and pathways by total overlap count
    top_tfs <- overlap_df %>%
      dplyr::group_by(.data$TF) %>%
      dplyr::summarise(total = sum(.data$n_overlap), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$total)) %>%
      dplyr::slice_head(n = top.tf) %>%
      dplyr::pull(.data$TF)

    top_pws <- overlap_df %>%
      dplyr::group_by(.data$Pathway) %>%
      dplyr::summarise(total = sum(.data$n_overlap), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$total)) %>%
      dplyr::slice_head(n = top.pathway) %>%
      dplyr::pull(.data$Pathway)

    plot_df <- overlap_df %>%
      dplyr::filter(.data$TF %in% top_tfs, .data$Pathway %in% top_pws)

    if (nrow(plot_df) > 0) {
      # Order by total overlap
      tf_order <- plot_df %>%
        dplyr::group_by(.data$TF) %>%
        dplyr::summarise(total = sum(.data$n_overlap), .groups = "drop") %>%
        dplyr::arrange(.data$total) %>%
        dplyr::pull(.data$TF)
      pw_order <- plot_df %>%
        dplyr::group_by(.data$Pathway) %>%
        dplyr::summarise(total = sum(.data$n_overlap), .groups = "drop") %>%
        dplyr::arrange(.data$total) %>%
        dplyr::pull(.data$Pathway)

      plot_df$TF <- factor(plot_df$TF, levels = tf_order)
      plot_df$Pathway <- factor(plot_df$Pathway, levels = pw_order)

      result$plot <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$TF, y = .data$Pathway,
                     fill = .data$n_overlap, size = .data$jaccard)
      ) +
        ggplot2::geom_point(shape = 21, color = "black", stroke = 0.3) +
        ggplot2::scale_fill_gradient(low = "lightyellow", high = "red",
                                     name = "Overlap\ncount") +
        ggplot2::scale_size_continuous(name = "Jaccard", range = c(1, 8)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::xlab("TF / Regulon") +
        ggplot2::ylab("Pathway") +
        ggplot2::ggtitle("TF-Pathway Core Enrichment Overlap")
    }
  }

  return(result)
}
