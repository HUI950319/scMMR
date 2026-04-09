# ============================================================================
# Monocle3 wrapper — pseudotime trajectory inference
# ============================================================================


#' Run Monocle3 trajectory inference on a Seurat object
#'
#' Thin wrapper around the public \pkg{monocle3} API that converts a Seurat
#' object to a \code{cell_data_set}, runs the standard
#' \code{cluster_cells -> learn_graph -> order_cells} pipeline, and writes
#' clusters, partitions and pseudotime back into the Seurat object.
#'
#' This wrapper deliberately uses only public Monocle3 functions (no
#' \code{monocle3:::} internals) and is non-interactive: roots must be
#' supplied as parameters.
#'
#' @param srt A Seurat object.
#' @param assay Assay to pull counts from. Default \code{NULL} (active assay).
#' @param layer Layer/slot to pull counts from. Default \code{"counts"}.
#' @param reduction Name of the reduction in \code{srt} to use as the
#'   2D embedding for Monocle3. Default \code{"umap"}.
#' @param root_cells Character vector of cell barcodes to use as the
#'   trajectory root. Either \code{root_cells}, \code{root_pr_nodes}, or
#'   \code{root_group} must be supplied.
#' @param root_pr_nodes Character vector of principal-graph node names
#'   (e.g. \code{"Y_1"}) to use as the root.
#' @param root_group A length-1 character: name of a metadata column. The
#'   most-stem-like value in that column will be used to pick root cells
#'   (cells whose values match \code{root_group_value}).
#' @param root_group_value The level of \code{root_group} to use as root
#'   (e.g. \code{"Stem"}).
#' @param use_partition Passed to \code{monocle3::learn_graph()}. If
#'   \code{TRUE} (default), each partition gets its own disjoint graph.
#' @param close_loop Passed to \code{monocle3::learn_graph()}. Default
#'   \code{TRUE}.
#' @param k Number of nearest neighbours for \code{cluster_cells}.
#'   Default 50.
#' @param cluster_method Clustering method for \code{cluster_cells}
#'   (\code{"louvain"} or \code{"leiden"}). Default \code{"louvain"}.
#' @param resolution Resolution parameter for \code{cluster_cells}.
#'   Default \code{NULL} (Monocle3 auto-selects).
#' @param partition_qval q-value threshold for partitioning. Default 0.05.
#' @param seed Random seed. Default 11.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return The input Seurat object with:
#' \itemize{
#'   \item Metadata columns added: \code{Monocle3_clusters},
#'         \code{Monocle3_partitions}, \code{Monocle3_Pseudotime}.
#'   \item \code{srt@tools$Monocle3} list containing:
#'         \code{cds} (the fitted \code{cell_data_set}),
#'         \code{edge_df} (principal-graph edges in 2D),
#'         \code{node_df} (branch / leaf / internal nodes in 2D),
#'         \code{trajectory} (a list of ggplot2 layers — \code{geom_segment}
#'         for the principal graph; just add to any ggplot),
#'         \code{milestones} (a list of ggplot2 layers — leaf nodes in
#'         black, branch nodes in red),
#'         \code{root_cells} or \code{root_pr_nodes} (whichever was used).
#' }
#'
#' @details
#' \strong{Choosing roots.} You must give the function a way to pick the
#' starting point of pseudotime. The three options, in order of preference:
#' \enumerate{
#'   \item \code{root_cells}: explicit cell barcodes you trust as the start.
#'   \item \code{root_pr_nodes}: principal-graph node names. Run the function
#'         once without a root to see the available nodes in the returned
#'         \code{node_df}, then re-run with the chosen node(s).
#'   \item \code{root_group} + \code{root_group_value}: pick all cells whose
#'         metadata column \code{root_group} equals \code{root_group_value}.
#' }
#' If none is supplied, the function will run everything up to
#' \code{learn_graph} and return without ordering, so you can inspect
#' \code{node_df} and decide.
#'
#' \strong{Plotting the trajectory.} The returned \code{edge_df} and
#' \code{node_df} are plain data frames in the coordinates of \code{reduction},
#' so you can overlay them on any \code{ggplot2} scatter:
#' \preformatted{
#'   library(ggplot2)
#'   d <- cbind(srt@meta.data, srt[["umap"]]@cell.embeddings)
#'   ggplot(d, aes(umap_1, umap_2)) +
#'     geom_point(aes(colour = Monocle3_Pseudotime), size = .3) +
#'     geom_segment(data = srt@tools$Monocle3$edge_df,
#'                  aes(x = x, y = y, xend = xend, yend = yend),
#'                  inherit.aes = FALSE) +
#'     geom_point(data = srt@tools$Monocle3$node_df,
#'                aes(x, y), inherit.aes = FALSE,
#'                shape = 21, fill = "red", size = 2)
#' }
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scMMR)
#'
#' seu <- qs::qread(system.file("extdata", "toy_test.qs", package = "scMMR"))
#'
#' # First pass: no root, just to see available nodes
#' seu <- RunMonocle3(seu, reduction = "umap")
#' head(seu@tools$Monocle3$node_df)
#'
#' # Second pass: pick a root by metadata
#' seu <- RunMonocle3(
#'   seu,
#'   reduction       = "umap",
#'   root_group      = "celltype",
#'   root_group_value = "B cells"
#' )
#' head(seu$Monocle3_Pseudotime)
#' }
#'
#' @export
RunMonocle3 <- function(srt,
                        assay            = NULL,
                        layer            = "counts",
                        reduction        = "umap",
                        root_cells       = NULL,
                        root_pr_nodes    = NULL,
                        root_group       = NULL,
                        root_group_value = NULL,
                        use_partition    = TRUE,
                        close_loop       = TRUE,
                        k                = 50,
                        cluster_method   = c("louvain", "leiden"),
                        resolution       = NULL,
                        partition_qval   = 0.05,
                        seed             = 11,
                        verbose          = TRUE) {

  # ── 0. checks ────────────────────────────────────────────────────────────
  if (!inherits(srt, "Seurat")) {
    stop("`srt` must be a Seurat object.", call. = FALSE)
  }
  for (pkg in c("monocle3", "SingleCellExperiment", "SeuratObject", "igraph")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package `", pkg, "` is required. Install with:\n",
        if (pkg == "monocle3")
          "  remotes::install_github('cole-trapnell-lab/monocle3')"
        else paste0("  install.packages('", pkg, "')"),
        call. = FALSE
      )
    }
  }
  cluster_method <- match.arg(cluster_method)

  set.seed(seed)
  msg <- function(...) if (isTRUE(verbose)) message("[RunMonocle3] ", ...)

  # ── 1. extract counts ────────────────────────────────────────────────────
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!reduction %in% names(srt@reductions)) {
    stop("Reduction `", reduction, "` not found in the Seurat object.",
         call. = FALSE)
  }
  msg("Building cell_data_set from assay '", assay, "', layer '", layer, "' ...")

  expr <- SeuratObject::GetAssayData(srt, assay = assay, layer = layer)
  if (is.null(expr) || length(expr) == 0L) {
    stop("Empty expression matrix from assay='", assay,
         "', layer='", layer, "'.", call. = FALSE)
  }
  expr <- methods::as(expr, "CsparseMatrix")

  cell_md <- srt@meta.data
  gene_md <- data.frame(
    gene_short_name = rownames(expr),
    row.names       = rownames(expr)
  )

  cds <- monocle3::new_cell_data_set(
    expression_data = expr,
    cell_metadata   = cell_md,
    gene_metadata   = gene_md
  )

  # ── 2. inject pre-computed reduction ─────────────────────────────────────
  emb <- SeuratObject::Embeddings(srt[[reduction]])
  if (ncol(emb) < 2L) {
    stop("Reduction `", reduction, "` must have at least 2 dimensions.",
         call. = FALSE)
  }
  emb <- emb[colnames(cds), 1:2, drop = FALSE]
  SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- emb

  # ── 3. cluster_cells -> learn_graph ─────────────────────────────────────
  msg("Clustering cells (k=", k, ", method='", cluster_method, "') ...")
  cds <- monocle3::cluster_cells(
    cds              = cds,
    reduction_method = "UMAP",
    k                = k,
    cluster_method   = cluster_method,
    resolution       = resolution,
    partition_qval   = partition_qval,
    verbose          = verbose
  )

  msg("Learning principal graph (use_partition=", use_partition,
      ", close_loop=", close_loop, ") ...")
  cds <- monocle3::learn_graph(
    cds            = cds,
    use_partition  = use_partition,
    close_loop     = close_loop,
    verbose        = verbose
  )

  # ── 4. extract graph nodes/edges in 2D ──────────────────────────────────
  node_coords <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  pg          <- cds@principal_graph[["UMAP"]]
  edge_df     <- igraph::as_data_frame(pg, what = "edges")
  edge_df[, c("x",    "y")]    <- node_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- node_coords[edge_df[["to"]],   1:2]

  deg <- igraph::degree(pg)
  node_role <- ifelse(deg == 1L, "leaf",
               ifelse(deg >= 3L, "branch", "internal"))
  node_df <- data.frame(
    node = igraph::V(pg)$name,
    x    = node_coords[igraph::V(pg)$name, 1],
    y    = node_coords[igraph::V(pg)$name, 2],
    role = node_role,
    stringsAsFactors = FALSE
  )

  # ── 5. resolve root ──────────────────────────────────────────────────────
  if (is.null(root_cells) && is.null(root_pr_nodes) && !is.null(root_group)) {
    if (!root_group %in% colnames(srt@meta.data)) {
      stop("`root_group` column '", root_group, "' not in meta.data.",
           call. = FALSE)
    }
    if (is.null(root_group_value)) {
      stop("`root_group_value` must be supplied with `root_group`.",
           call. = FALSE)
    }
    root_cells <- rownames(srt@meta.data)[
      srt@meta.data[[root_group]] == root_group_value
    ]
    root_cells <- intersect(root_cells, colnames(cds))
    if (length(root_cells) == 0L) {
      stop("No cells matched ", root_group, " == ", root_group_value,
           call. = FALSE)
    }
    msg("Using ", length(root_cells), " cells from '",
        root_group, " == ", root_group_value, "' as roots.")
  }

  has_root <- !is.null(root_cells) || !is.null(root_pr_nodes)

  # ── 6. order_cells (only if root supplied) ──────────────────────────────
  if (has_root) {
    msg("Ordering cells along pseudotime ...")
    cds <- monocle3::order_cells(
      cds           = cds,
      root_pr_nodes = root_pr_nodes,
      root_cells    = root_cells
    )
    pseudotime <- monocle3::pseudotime(cds)
    pseudotime[is.infinite(pseudotime)] <- NA_real_
    srt[["Monocle3_Pseudotime"]] <- pseudotime[colnames(srt)]
  } else {
    msg("No root supplied — skipping order_cells. ",
        "Inspect `node_df` and re-run with `root_pr_nodes` or `root_cells`.")
    srt[["Monocle3_Pseudotime"]] <- NA_real_
  }

  # ── 7. write back to Seurat ─────────────────────────────────────────────
  cl_vec <- cds@clusters[["UMAP"]]$clusters
  pt_vec <- cds@clusters[["UMAP"]]$partitions
  srt[["Monocle3_clusters"]]   <- as.factor(cl_vec[colnames(srt)])
  srt[["Monocle3_partitions"]] <- as.factor(pt_vec[colnames(srt)])

  # ── 7b. ggplot2-ready trajectory & milestones layers ────────────────────
  # `trajectory` and `milestones` are lists of ggplot2 layers, so users can
  # add them directly to any scatter plot built on the same reduction:
  #     ggplot(d, aes(umap_1, umap_2)) + geom_point() + trajectory + milestones
  trajectory <- list(
    ggplot2::geom_segment(
      data        = edge_df,
      mapping     = ggplot2::aes(x = .data$x, y = .data$y,
                                 xend = .data$xend, yend = .data$yend),
      inherit.aes = FALSE,
      linewidth   = 0.5,
      colour      = "black"
    )
  )
  leaf_df   <- node_df[node_df$role == "leaf",   , drop = FALSE]
  branch_df <- node_df[node_df$role == "branch", , drop = FALSE]
  milestones <- list(
    ggplot2::geom_point(
      data        = leaf_df,
      mapping     = ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      shape       = 21, colour = "white", fill = "black",
      size        = 3, stroke = 1
    ),
    ggplot2::geom_point(
      data        = branch_df,
      mapping     = ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      shape       = 21, colour = "white", fill = "red",
      size        = 3, stroke = 1
    )
  )

  srt@tools$Monocle3 <- list(
    cds           = cds,
    edge_df       = edge_df,
    node_df       = node_df,
    trajectory    = trajectory,
    milestones    = milestones,
    root_cells    = root_cells,
    root_pr_nodes = root_pr_nodes,
    reduction     = reduction
  )

  msg("Done. Cells: ", ncol(cds),
      " | Clusters: ", nlevels(as.factor(cl_vec)),
      " | Partitions: ", nlevels(as.factor(pt_vec)),
      if (has_root) paste0(" | Pseudotime range: [",
                           sprintf("%.2f", min(pseudotime, na.rm = TRUE)), ", ",
                           sprintf("%.2f", max(pseudotime, na.rm = TRUE)), "]")
      else " | Pseudotime: <not computed>")

  invisible(srt)
}
