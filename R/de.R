# --------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------

#' Combine p-values using metap
#' @keywords internal
.combine_pvals <- function(p, method = "maximump") {
  if (!requireNamespace("metap", quietly = TRUE)) {
    stop("Package 'metap' is required for conserved markers.\n",
         "Install with: install.packages('metap')")
  }
  res <- do.call(method, args = list(p = p))
  res$p
}

#' Post-process marker results: add gene/group/p_val_adj/test_group columns
#' @keywords internal
.post_process_markers <- function(markers,
                                  group1_str,
                                  group2_str = "others",
                                  p.adjust.method = "bonferroni",
                                  group_levels = NULL) {
  if (is.null(markers) || nrow(markers) == 0) return(NULL)

  markers[["gene"]] <- rownames(markers)
  markers[["group1"]] <- group1_str
  markers[["group2"]] <- group2_str
  rownames(markers) <- NULL
  if (!is.null(group_levels)) {
    markers[["group1"]] <- factor(markers[["group1"]], levels = group_levels)
  } else {
    markers[["group1"]] <- factor(markers[["group1"]],
                                  levels = unique(markers[["group1"]]))
  }
  if ("p_val" %in% colnames(markers)) {
    markers[["p_val_adj"]] <- stats::p.adjust(markers[["p_val"]],
                                              method = p.adjust.method)
  }
  markers
}

#' Add test_group_number and test_group columns
#' @keywords internal
#' Split composite group1/group2 columns back into celltype + condition
#' @keywords internal
.split_composite_cols <- function(markers, group_name, split_name) {
  if (is.null(markers) || nrow(markers) == 0) return(markers)
  sep <- "_"
  for (col in c("group1", "group2")) {
    if (!col %in% colnames(markers)) next
    vals <- as.character(markers[[col]])
    # Split on the last "_" to handle group names containing "_"
    # Pattern: everything up to last "_" = group, after last "_" = split
    ct_col <- paste0(col, "_", group_name)
    sp_col <- paste0(col, "_", split_name)
    last_sep <- regexpr(paste0(sep, "[^", sep, "]*$"), vals)
    markers[[ct_col]] <- ifelse(last_sep > 0,
                                substr(vals, 1, last_sep - 1), vals)
    markers[[sp_col]] <- ifelse(last_sep > 0,
                                substr(vals, last_sep + 1, nchar(vals)), NA)
  }
  markers
}

.add_test_group_info <- function(markers) {
  if (is.null(markers) || nrow(markers) == 0) return(markers)
  markers[["test_group_number"]] <- as.integer(
    table(markers[["gene"]])[markers[["gene"]]]
  )
  markers_matrix <- as.data.frame.matrix(
    table(markers[, c("gene", "group1")])
  )
  markers[["test_group"]] <- apply(markers_matrix, 1, function(x) {
    paste0(colnames(markers_matrix)[x > 0], collapse = ";")
  })[markers[["gene"]]]
  markers
}

#' Post-filter DE results by adjusted p-value and/or log2 fold change
#' @keywords internal
.filter_de_results <- function(markers, p_val_cutoff = NULL,
                                logfc_cutoff = NULL, verbose = TRUE) {
  if (is.null(markers) || nrow(markers) == 0) return(markers)
  n_before <- nrow(markers)

  if (!is.null(p_val_cutoff)) {
    pcol <- if ("p_val_adj" %in% colnames(markers)) "p_val_adj" else "p_val"
    markers <- markers[!is.na(markers[[pcol]]) &
                         markers[[pcol]] < p_val_cutoff, , drop = FALSE]
  }
  if (!is.null(logfc_cutoff)) {
    if ("avg_log2FC" %in% colnames(markers)) {
      markers <- markers[abs(markers[["avg_log2FC"]]) > logfc_cutoff, ,
                          drop = FALSE]
    }
  }
  if (verbose && nrow(markers) < n_before) {
    message("Post-filter: ", n_before, " -> ", nrow(markers), " markers")
  }
  markers
}

#' Find conserved markers across grouping variable levels
#'
#' Splits cells by grouping.var, runs FindMarkers within each level,
#' intersects significant genes, and combines p-values.
#'
#' @keywords internal
.find_conserved_markers <- function(object,
                                    assay,
                                    layer,
                                    cells.1,
                                    cells.2,
                                    grouping.var,
                                    features = NULL,
                                    test.use = "wilcox",
                                    logfc.threshold = 0.25,
                                    base = 2,
                                    pseudocount.use = 1,
                                    mean.fxn = NULL,
                                    min.pct = 0.1,
                                    min.diff.pct = -Inf,
                                    max.cells.per.ident = Inf,
                                    latent.vars = NULL,
                                    only.pos = FALSE,
                                    min.cells.group = 3,
                                    min.cells.feature = 3,
                                    meta.method = "maximump",
                                    norm.method = "LogNormalize",
                                    verbose = TRUE,
                                    ...) {
  object.var <- SeuratObject::FetchData(object = object, vars = grouping.var)
  levels.split <- names(sort(table(object.var[, 1])))
  num.groups <- length(levels.split)

  # Split cells by grouping.var levels
  cells_by_level <- lapply(seq_len(num.groups), function(i) {
    rownames(object.var[object.var[, 1] == levels.split[i], , drop = FALSE])
  })

  marker.test <- list()
  for (i in seq_len(num.groups)) {
    level.use <- levels.split[i]
    cells.1.use <- intersect(cells_by_level[[i]], cells.1)
    if (length(cells.1.use) < min.cells.group) {
      if (verbose) message("  ", level.use, ": fewer than ", min.cells.group,
                           " cells in group1, skipping")
      next
    }
    if (is.null(cells.2)) {
      cells.2.use <- setdiff(cells_by_level[[i]], cells.1.use)
    } else {
      cells.2.use <- intersect(cells_by_level[[i]], cells.2)
    }
    if (length(cells.2.use) < min.cells.group) {
      if (verbose) message("  ", level.use, ": fewer than ", min.cells.group,
                           " cells in group2, skipping")
      next
    }
    if (verbose) message("  Testing group ", level.use)
    marker.test[[level.use]] <- Seurat::FindMarkers(
      object = Seurat::GetAssay(object, assay),
      layer = layer,
      cells.1 = cells.1.use,
      cells.2 = cells.2.use,
      features = features,
      test.use = test.use,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct,
      min.diff.pct = min.diff.pct,
      max.cells.per.ident = max.cells.per.ident,
      min.cells.group = min.cells.group,
      min.cells.feature = min.cells.feature,
      norm.method = norm.method,
      base = base,
      pseudocount.use = pseudocount.use,
      mean.fxn = mean.fxn,
      latent.vars = latent.vars,
      only.pos = only.pos,
      verbose = FALSE,
      ...
    )
  }

  marker.test <- marker.test[!vapply(marker.test, is.null, logical(1))]
  if (length(marker.test) == 0) {
    if (verbose) warning("No group was tested")
    return(NULL)
  }

  # Intersect genes across all levels
  genes.conserved <- Reduce(
    f = intersect,
    x = lapply(marker.test, rownames)
  )
  if (length(genes.conserved) == 0) return(NULL)

  # Combine results with level-prefixed column names
  markers.conserved <- lapply(seq_along(marker.test), function(i) {
    df <- marker.test[[i]][genes.conserved, , drop = FALSE]
    colnames(df) <- paste(names(marker.test)[i], colnames(df), sep = "_")
    df
  })
  markers.combined <- Reduce(cbind, markers.conserved)

  # Overall fold change
  fc <- Seurat::FoldChange(
    Seurat::GetAssay(object, assay),
    layer = layer,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = genes.conserved,
    norm.method = norm.method,
    base = base,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn
  )
  markers.combined <- cbind(markers.combined, fc[genes.conserved, , drop = FALSE])

  # Filter by fold change direction
  logFC.codes <- colnames(markers.combined)[grepl("avg_log.*FC$",
                                                  colnames(markers.combined))]
  if (isTRUE(only.pos)) {
    keep <- apply(markers.combined[, logFC.codes, drop = FALSE] > 0, 1, all)
  } else {
    keep <- apply(markers.combined[, logFC.codes, drop = FALSE] < 0, 1, all) |
      apply(markers.combined[, logFC.codes, drop = FALSE] > 0, 1, all)
  }
  markers.combined <- markers.combined[keep, , drop = FALSE]
  if (nrow(markers.combined) == 0) return(NULL)

  # Combine p-values
  pval.codes <- colnames(markers.combined)[grepl("_p_val$",
                                                 colnames(markers.combined))]
  if (length(pval.codes) > 1) {
    markers.combined[["max_pval"]] <- apply(
      markers.combined[, pval.codes, drop = FALSE], 1, max
    )
    combined.pval <- apply(
      markers.combined[, pval.codes, drop = FALSE], 1,
      function(x) .combine_pvals(x, method = meta.method)
    )
    markers.combined[[paste0(meta.method, "_p_val")]] <- combined.pval
    markers.combined[["p_val"]] <- combined.pval
    markers.combined <- markers.combined[order(combined.pval), , drop = FALSE]
  } else {
    markers.combined[["max_pval"]] <- markers.combined[["p_val"]] <-
      markers.combined[, pval.codes]
    if (verbose) warning("Only a single group was tested")
  }

  markers.combined
}

# --------------------------------------------------------------------------
# Main exported function
# --------------------------------------------------------------------------

#' Run Differential Expression Test
#'
#' Perform differential expression (DE) testing across cell groups in a
#' Seurat object. Supports four marker types: one-vs-rest (\code{"all"}),
#' pairwise (\code{"paired"}), conserved across conditions
#' (\code{"conserved"}), and condition-specific (\code{"disturbed"}).
#' Optionally downsamples cells before testing for balanced comparisons.
#'
#' @param seu A Seurat object.
#' @param group.by Column name in \code{meta.data} for grouping. If
#'   \code{NULL}, uses \code{Idents(seu)}.
#' @param split.by Optional column name for a second grouping variable
#'   (e.g. condition). When provided, cells are grouped by
#'   \code{group.by x split.by} composite identity for DE testing.
#'   Result columns \code{celltype} and \code{condition} are added.
#'   Default: NULL.
#' @param group1 Character vector of group labels for the first group.
#' @param group2 Character vector of group labels for the second group.
#' @param cells1 Character vector of cell barcodes for the first group.
#'   Overrides \code{group1}.
#' @param cells2 Character vector of cell barcodes for the second group.
#' @param features Character vector of genes to test. Default: all genes.
#' @param markers_type Type of marker analysis. One of \code{"all"}
#'   (one-vs-rest), \code{"paired"} (all pairwise), \code{"conserved"}
#'   (across conditions), or \code{"disturbed"} (condition-specific).
#' @param grouping.var Column for condition variable (required for
#'   conserved/disturbed).
#' @param meta.method Method for combining p-values in conserved markers.
#'   One of \code{"maximump"}, \code{"minimump"}, \code{"wilkinsonp"},
#'   \code{"meanp"}, \code{"sump"}, \code{"votep"}.
#' @param test.use DE test method. Passed to \code{Seurat::FindMarkers}.
#' @param only.pos Logical; return only positive markers. Default: TRUE.
#' @param fc.threshold Fold change threshold (linear scale, >= 1).
#'   Pre-filter before testing: genes below this threshold are skipped.
#'   Default: 1.5.
#' @param p_val_cutoff Numeric or NULL. Post-filter: keep only results with
#'   \code{p_val_adj < p_val_cutoff}. Default: NULL (no filtering).
#' @param logfc_cutoff Numeric or NULL. Post-filter: keep only results with
#'   \code{|avg_log2FC| > logfc_cutoff}. Default: NULL (no filtering).
#' @param base Log base for fold change. Default: 2.
#' @param pseudocount.use Pseudocount for log fold change. Default: 1.
#' @param mean.fxn Custom mean function for fold change.
#' @param min.pct Minimum fraction of cells expressing the gene. Default: 0.1.
#' @param min.diff.pct Minimum difference in fraction between groups.
#' @param max.cells.per.ident Maximum cells per identity for testing.
#' @param latent.vars Latent variables for certain tests (MAST, LR, etc.).
#' @param min.cells.feature Minimum cells expressing feature. Default: 3.
#' @param min.cells.group Minimum cells per group. Default: 3.
#' @param norm.method Normalization method. Default: "LogNormalize".
#' @param p.adjust.method P-value adjustment method. Default: "bonferroni".
#' @param downsample Integer or NULL. If specified, downsample each identity
#'   class to at most this many cells before DE testing. Uses
#'   \code{subset(seu, downsample = downsample)}. Default: NULL.
#' @param layer Data layer to use. Default: "data".
#' @param assay Assay to use. Default: \code{DefaultAssay(seu)}.
#' @param cores Integer. Number of parallel cores. Default:
#'   \code{parallel::detectCores() \%/\% 2} (half of available cores).
#'   Set to 1 to disable parallelism.
#' @param seed Random seed. Default: 11.
#' @param verbose Logical; print progress messages. Default: TRUE.
#'
#' @return A data.frame of DE results with columns: p_val, avg_log2FC,
#'   pct.1, pct.2, p_val_adj, gene, group1, group2, test_group_number,
#'   test_group.
#'
#' @examples
#' \dontrun{
#' # One-vs-rest DE
#' markers <- RunDE(seu, group.by = "celltype")
#'
#' # Pairwise with downsampling
#' markers <- RunDE(seu, group.by = "celltype",
#'                  markers_type = "paired", downsample = 200)
#'
#' # Conserved markers across conditions
#' markers <- RunDE(seu, group.by = "celltype",
#'                  grouping.var = "group",
#'                  markers_type = "conserved")
#'
#' # Joint celltype x condition DE (for downstream GSEA)
#' markers <- RunDE(seu, group.by = "celltype", split.by = "condition")
#'
#' # With post-filtering: padj < 0.05 and |log2FC| > 1
#' markers <- RunDE(seu, group.by = "celltype",
#'                  p_val_cutoff = 0.05, logfc_cutoff = 1)
#' }
#'
#' @seealso \code{\link[Seurat]{FindMarkers}},
#'   \code{\link{RunPathwayAnalysis}}
#' @export
RunDE <- function(seu,
                  group.by = NULL,
                  split.by = NULL,
                  group1 = NULL,
                  group2 = NULL,
                  cells1 = NULL,
                  cells2 = NULL,
                  features = NULL,
                  markers_type = c("all", "paired", "conserved", "disturbed"),
                  grouping.var = NULL,
                  meta.method = c("maximump", "minimump", "wilkinsonp",
                                  "meanp", "sump", "votep"),
                  test.use = "wilcox",
                  only.pos = TRUE,
                  fc.threshold = 1.5,
                  p_val_cutoff = NULL,
                  logfc_cutoff = NULL,
                  base = 2,
                  pseudocount.use = 1,
                  mean.fxn = NULL,
                  min.pct = 0.1,
                  min.diff.pct = -Inf,
                  max.cells.per.ident = Inf,
                  latent.vars = NULL,
                  min.cells.feature = 3,
                  min.cells.group = 3,
                  norm.method = "LogNormalize",
                  p.adjust.method = "bonferroni",
                  downsample = NULL,
                  layer = "data",
                  assay = NULL,
                  cores = max(1L, parallel::detectCores() %/% 2),
                  seed = 11,
                  verbose = TRUE) {
  # --- validation ---
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Install with: install.packages('Seurat')")
  }
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required.")
  }
  markers_type <- match.arg(markers_type)
  meta.method <- match.arg(meta.method)
  if (markers_type %in% c("conserved", "disturbed") && is.null(grouping.var)) {
    stop("'grouping.var' must be provided for conserved/disturbed markers")
  }
  if (fc.threshold < 1) {
    stop("'fc.threshold' must be >= 1")
  }
  set.seed(seed)
  assay <- assay %||% Seurat::DefaultAssay(seu)
  logfc.threshold <- log(fc.threshold, base = base)

  # --- split.by: create composite identity ---
  .split_by_used <- FALSE
  if (!is.null(split.by)) {
    if (is.null(group.by)) {
      stop("'group.by' must be provided when 'split.by' is used.", call. = FALSE)
    }
    if (!split.by %in% colnames(seu@meta.data)) {
      stop("'split.by' column '", split.by, "' not found in meta.data.",
           call. = FALSE)
    }
    .split_by_used <- TRUE
    .orig_group.by <- group.by
    .orig_split.by <- split.by
    composite_col <- paste0("..de_", group.by, "_", split.by, "..")
    seu <- Seurat::AddMetaData(
      seu,
      paste(seu@meta.data[[group.by]], seu@meta.data[[split.by]], sep = "_"),
      col.name = composite_col
    )
    group.by <- composite_col
    if (verbose) {
      n_combo <- length(unique(seu@meta.data[[composite_col]]))
      message("split.by: created ", n_combo, " composite groups (",
              .orig_group.by, " x ", .orig_split.by, ")")
    }
  }

  # --- downsample ---
  if (!is.null(downsample)) {
    n_before <- ncol(seu)
    if (!is.null(group.by)) {
      Seurat::Idents(seu) <- group.by
    }
    seu <- subset(seu, downsample = as.integer(downsample))
    if (verbose) {
      message("Downsampled: ", n_before, " -> ", ncol(seu), " cells")
    }
  }

  # Common FindMarkers args
  fm_args <- list(
    layer = layer,
    features = features,
    test.use = test.use,
    logfc.threshold = logfc.threshold,
    base = base,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    max.cells.per.ident = max.cells.per.ident,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    latent.vars = latent.vars,
    only.pos = only.pos,
    norm.method = norm.method,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn,
    verbose = FALSE
  )

  # ======================================================================
  # Branch 1: Custom cells (cells1/group1 provided)
  # ======================================================================
  if (!is.null(cells1) || !is.null(group1)) {
    if (is.null(cells1)) {
      if (is.null(group.by)) {
        stop("'group.by' must be provided when 'group1' is specified")
      }
      cells1 <- colnames(seu)[seu[[group.by, drop = TRUE]] %in% group1]
    }
    if (is.null(cells2) && !is.null(group2)) {
      cells2 <- colnames(seu)[seu[[group.by, drop = TRUE]] %in% group2]
    }
    if (!all(cells1 %in% colnames(seu))) {
      stop("Some cells in 'cells1' are not in the Seurat object")
    }
    if (is.null(cells2)) {
      cells2 <- setdiff(colnames(seu), cells1)
      group2 <- "others"
    }
    if (!all(cells2 %in% colnames(seu))) {
      stop("Some cells in 'cells2' are not in the Seurat object")
    }
    if (length(cells1) < 3 || length(cells2) < 3) {
      stop("Cell groups must have at least 3 cells each")
    }

    group1_str <- if (is.null(group1)) {
      "group1"
    } else if (length(group1) > 1) {
      paste(group1, collapse = ";")
    } else {
      as.character(group1)
    }
    group2_str <- if (is.null(group2)) {
      "group2"
    } else if (length(group2) > 1) {
      paste(group2, collapse = ";")
    } else {
      as.character(group2)
    }

    if (verbose) message("Finding ", markers_type, " markers (", test.use,
                         ") for custom cell groups...")

    # -- all --
    if (markers_type == "all") {
      fm_args[["object"]] <- Seurat::GetAssay(seu, assay)
      fm_args[["cells.1"]] <- cells1
      fm_args[["cells.2"]] <- cells2
      markers <- do.call(Seurat::FindMarkers, fm_args)
      markers <- .post_process_markers(markers, group1_str, group2_str,
                                       p.adjust.method)
      if (is.null(markers)) {
        warning("No markers found")
        return(data.frame())
      }
      markers <- .add_test_group_info(markers)
      markers <- .filter_de_results(markers, p_val_cutoff, logfc_cutoff,
                                     verbose)
      return(markers)
    }

    # -- conserved --
    if (markers_type == "conserved") {
      markers <- .find_conserved_markers(
        object = seu, assay = assay, layer = layer,
        cells.1 = cells1, cells.2 = cells2,
        grouping.var = grouping.var, features = features,
        test.use = test.use, logfc.threshold = logfc.threshold,
        base = base, pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn, min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        latent.vars = latent.vars, only.pos = only.pos,
        min.cells.group = min.cells.group,
        min.cells.feature = min.cells.feature,
        meta.method = meta.method, norm.method = norm.method,
        verbose = verbose
      )
      markers <- .post_process_markers(markers, group1_str, group2_str,
                                       p.adjust.method)
      if (is.null(markers)) {
        warning("No markers found")
        return(data.frame())
      }
      markers <- .add_test_group_info(markers)
      markers <- .filter_de_results(markers, p_val_cutoff, logfc_cutoff,
                                     verbose)
      return(markers)
    }

    # -- disturbed --
    if (markers_type == "disturbed") {
      seu_tmp <- seu
      gv <- seu_tmp[[grouping.var, drop = TRUE]]
      gv[setdiff(colnames(seu_tmp), cells1)] <- NA
      seu_tmp <- Seurat::AddMetaData(seu_tmp, gv, col.name = grouping.var)

      if (length(stats::na.omit(unique(gv))) < 2) {
        warning("Fewer than 2 levels in grouping.var within cells1")
        return(data.frame())
      }
      markers <- RunDE(
        seu = seu_tmp, assay = assay, layer = layer,
        group.by = grouping.var, markers_type = "all",
        features = features, test.use = test.use,
        fc.threshold = fc.threshold, base = base,
        min.pct = min.pct, min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars, only.pos = only.pos,
        norm.method = norm.method, p.adjust.method = p.adjust.method,
        pseudocount.use = pseudocount.use, mean.fxn = mean.fxn,
        cores = 1L, seed = seed, verbose = FALSE
      )
      if (!is.null(markers) && nrow(markers) > 0) {
        colnames(markers) <- gsub("group", "var", colnames(markers))
        markers[["group1"]] <- group1_str
      } else {
        warning("No markers found")
        return(data.frame())
      }
      markers <- .filter_de_results(markers, p_val_cutoff, logfc_cutoff,
                                     verbose)
      return(markers)
    }

    # paired not applicable for custom cells
    if (markers_type == "paired") {
      warning("'paired' markers_type is not applicable for custom cell groups. ",
              "Use 'all' instead.")
      return(data.frame())
    }
  }

  # ======================================================================
  # Branch 2: Group-based analysis
  # ======================================================================
  if (is.null(group.by)) {
    cell_group <- Seurat::Idents(seu)
    group.by <- "active.ident"
  } else {
    cell_group <- seu[[group.by, drop = TRUE]]
  }
  if (!is.factor(cell_group)) {
    cell_group <- factor(cell_group, levels = unique(cell_group))
  }
  names(cell_group) <- colnames(seu)

  # Subsample per max.cells.per.ident
  if (is.finite(max.cells.per.ident)) {
    cell_group_list <- lapply(levels(cell_group), function(x) {
      cells <- cell_group[cell_group == x]
      sample(cells, size = min(max.cells.per.ident, length(cells)),
             replace = FALSE)
    })
    cell_group <- stats::setNames(
      unlist(lapply(cell_group_list, function(x) x), use.names = FALSE),
      unlist(lapply(cell_group_list, names))
    )
    cell_group <- factor(cell_group, levels = levels(cell_group))
  }

  fm_args[["object"]] <- Seurat::GetAssay(seu, assay)
  grp_levels <- levels(cell_group)
  n_groups <- length(grp_levels)

  if (verbose) message("Finding ", markers_type, " markers (", test.use,
                       ") among ", n_groups, " groups...")

  # choose apply function: parallel::mclapply when cores > 1
  cores <- max(1L, as.integer(cores))
  if (cores > 1) {
    apply_fun <- function(X, FUN, ...) {
      parallel::mclapply(X, FUN, ..., mc.cores = cores)
    }
    if (verbose) message("Using ", cores, " cores")
  } else if (requireNamespace("pbapply", quietly = TRUE) && verbose) {
    apply_fun <- pbapply::pblapply
  } else {
    apply_fun <- lapply
  }

  # ---- all ----
  if (markers_type == "all") {
    result_list <- apply_fun(grp_levels, function(group) {
      cells.1 <- names(cell_group)[cell_group == group]
      cells.2 <- names(cell_group)[cell_group != group]
      if (length(cells.1) < 3 || length(cells.2) < 3) return(NULL)
      fm_args[["cells.1"]] <- cells.1
      fm_args[["cells.2"]] <- cells.2
      markers <- do.call(Seurat::FindMarkers, fm_args)
      .post_process_markers(markers, as.character(group), "others",
                            p.adjust.method, grp_levels)
    })
    AllMarkers <- do.call(rbind.data.frame, result_list)
    if (is.null(AllMarkers) || nrow(AllMarkers) == 0) {
      warning("No markers found")
      return(data.frame())
    }
    rownames(AllMarkers) <- NULL
    AllMarkers[["group1"]] <- factor(AllMarkers[["group1"]],
                                     levels = grp_levels)
    AllMarkers <- .add_test_group_info(AllMarkers)
    if (.split_by_used) AllMarkers <- .split_composite_cols(
      AllMarkers, .orig_group.by, .orig_split.by)
    AllMarkers <- .filter_de_results(AllMarkers, p_val_cutoff, logfc_cutoff,
                                      verbose)
    return(AllMarkers)
  }

  # ---- paired ----
  if (markers_type == "paired") {
    pair <- expand.grid(x = grp_levels, y = grp_levels,
                        stringsAsFactors = FALSE)
    pair <- pair[pair[, 1] != pair[, 2], , drop = FALSE]

    result_list <- apply_fun(seq_len(nrow(pair)), function(i) {
      cells.1 <- names(cell_group)[cell_group == pair[i, 1]]
      cells.2 <- names(cell_group)[cell_group == pair[i, 2]]
      if (length(cells.1) < 3 || length(cells.2) < 3) return(NULL)
      fm_args[["cells.1"]] <- cells.1
      fm_args[["cells.2"]] <- cells.2
      markers <- do.call(Seurat::FindMarkers, fm_args)
      .post_process_markers(markers, as.character(pair[i, 1]),
                            as.character(pair[i, 2]),
                            p.adjust.method, grp_levels)
    })
    PairedMarkers <- do.call(rbind.data.frame, result_list)
    if (is.null(PairedMarkers) || nrow(PairedMarkers) == 0) {
      warning("No markers found")
      return(data.frame())
    }
    rownames(PairedMarkers) <- NULL
    PairedMarkers[["group1"]] <- factor(PairedMarkers[["group1"]],
                                        levels = grp_levels)
    PairedMarkers <- .add_test_group_info(PairedMarkers)
    if (.split_by_used) PairedMarkers <- .split_composite_cols(
      PairedMarkers, .orig_group.by, .orig_split.by)
    PairedMarkers <- .filter_de_results(PairedMarkers, p_val_cutoff,
                                         logfc_cutoff, verbose)
    return(PairedMarkers)
  }

  # ---- conserved ----
  if (markers_type == "conserved") {
    result_list <- apply_fun(grp_levels, function(group) {
      cells.1 <- names(cell_group)[cell_group == group]
      cells.2 <- names(cell_group)[cell_group != group]
      if (length(cells.1) < 3 || length(cells.2) < 3) return(NULL)
      markers <- .find_conserved_markers(
        object = seu, assay = assay, layer = layer,
        cells.1 = cells.1, cells.2 = cells.2,
        grouping.var = grouping.var, features = features,
        test.use = test.use, logfc.threshold = logfc.threshold,
        base = base, pseudocount.use = pseudocount.use,
        mean.fxn = mean.fxn, min.pct = min.pct,
        min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        latent.vars = latent.vars, only.pos = only.pos,
        min.cells.group = min.cells.group,
        min.cells.feature = min.cells.feature,
        meta.method = meta.method, norm.method = norm.method,
        verbose = FALSE
      )
      if (is.null(markers) || nrow(markers) == 0) return(NULL)
      markers <- .post_process_markers(markers, as.character(group),
                                       "others", p.adjust.method,
                                       grp_levels)
      markers
    })
    # Standardise columns to shared subset before rbind
    shared_cols <- c("avg_log2FC", "pct.1", "pct.2", "max_pval",
                     "p_val", "p_val_adj", "gene", "group1", "group2")
    ConservedMarkers <- do.call(
      rbind.data.frame,
      lapply(result_list, function(x) {
        if (is.null(x)) return(NULL)
        cols <- intersect(shared_cols, colnames(x))
        x[, cols, drop = FALSE]
      })
    )
    if (is.null(ConservedMarkers) || nrow(ConservedMarkers) == 0) {
      warning("No markers found")
      return(data.frame())
    }
    rownames(ConservedMarkers) <- NULL
    ConservedMarkers[["group1"]] <- factor(ConservedMarkers[["group1"]],
                                           levels = grp_levels)
    ConservedMarkers <- .add_test_group_info(ConservedMarkers)
    if (.split_by_used) ConservedMarkers <- .split_composite_cols(
      ConservedMarkers, .orig_group.by, .orig_split.by)
    ConservedMarkers <- .filter_de_results(ConservedMarkers, p_val_cutoff,
                                            logfc_cutoff, verbose)
    return(ConservedMarkers)
  }

  # ---- disturbed ----
  if (markers_type == "disturbed") {
    result_list <- apply_fun(grp_levels, function(group) {
      cells.1 <- names(cell_group)[cell_group == group]
      seu_tmp <- seu
      gv <- seu_tmp[[grouping.var, drop = TRUE]]
      gv[setdiff(colnames(seu_tmp), cells.1)] <- NA
      seu_tmp <- Seurat::AddMetaData(seu_tmp, gv, col.name = grouping.var)

      if (length(stats::na.omit(unique(gv))) < 2) return(NULL)
      markers <- RunDE(
        seu = seu_tmp, assay = assay, layer = layer,
        group.by = grouping.var, markers_type = "all",
        features = features, test.use = test.use,
        fc.threshold = fc.threshold, base = base,
        min.pct = min.pct, min.diff.pct = min.diff.pct,
        max.cells.per.ident = max.cells.per.ident,
        min.cells.feature = min.cells.feature,
        min.cells.group = min.cells.group,
        latent.vars = latent.vars, only.pos = only.pos,
        norm.method = norm.method, p.adjust.method = p.adjust.method,
        pseudocount.use = pseudocount.use, mean.fxn = mean.fxn,
        cores = 1L, seed = seed, verbose = FALSE
      )
      if (!is.null(markers) && nrow(markers) > 0) {
        colnames(markers) <- gsub("group", "var", colnames(markers))
        markers[["group1"]] <- as.character(group)
        return(markers)
      }
      NULL
    })
    DisturbedMarkers <- do.call(rbind.data.frame, result_list)
    if (is.null(DisturbedMarkers) || nrow(DisturbedMarkers) == 0) {
      warning("No markers found")
      return(data.frame())
    }
    rownames(DisturbedMarkers) <- NULL
    DisturbedMarkers[["group1"]] <- factor(DisturbedMarkers[["group1"]],
                                           levels = grp_levels)
    # test_group_number for disturbed uses unique gene-group1 pairs
    DisturbedMarkers[["test_group_number"]] <- as.integer(
      table(unique(DisturbedMarkers[, c("gene", "group1")])[["gene"]])[
        DisturbedMarkers[["gene"]]
      ]
    )
    dm_matrix <- as.data.frame.matrix(
      table(DisturbedMarkers[, c("gene", "group1")])
    )
    DisturbedMarkers[["test_group"]] <- apply(dm_matrix, 1, function(x) {
      paste0(colnames(dm_matrix)[x > 0], collapse = ";")
    })[DisturbedMarkers[["gene"]]]
    if (.split_by_used) DisturbedMarkers <- .split_composite_cols(
      DisturbedMarkers, .orig_group.by, .orig_split.by)
    DisturbedMarkers <- .filter_de_results(DisturbedMarkers, p_val_cutoff,
                                            logfc_cutoff, verbose)
    return(DisturbedMarkers)
  }
}

# ==========================================================================
# PlotDE: Visualization for DE results
# ==========================================================================

#' Normalize DE result column names
#' @keywords internal
.get_de_data <- function(res) {
  df <- as.data.frame(res)
  # Auto-detect BulkS4/limma column names
  if (!"gene" %in% colnames(df)) {
    if ("gene_id" %in% colnames(df)) {
      df[["gene"]] <- df[["gene_id"]]
    } else if (!is.null(rownames(df)) && nrow(df) > 0) {
      df[["gene"]] <- rownames(df)
    } else {
      stop("Cannot find gene column. Provide 'gene' or 'gene_id' column, or set row names.")
    }
  }
  if (!"avg_log2FC" %in% colnames(df)) {
    if ("logFC" %in% colnames(df)) {
      df[["avg_log2FC"]] <- df[["logFC"]]
    } else {
      stop("Cannot find fold change column ('avg_log2FC' or 'logFC')")
    }
  }
  if (!"p_val_adj" %in% colnames(df)) {
    if ("adj.P.Val" %in% colnames(df)) {
      df[["p_val_adj"]] <- df[["adj.P.Val"]]
    } else if ("padj" %in% colnames(df)) {
      df[["p_val_adj"]] <- df[["padj"]]
    } else {
      stop("Cannot find adjusted p-value column ('p_val_adj', 'adj.P.Val', or 'padj')")
    }
  }
  if (!"group1" %in% colnames(df)) {
    if ("cluster" %in% colnames(df)) {
      df[["group1"]] <- df[["cluster"]]
    } else {
      df[["group1"]] <- "All"
    }
  }
  if (!is.factor(df[["group1"]])) {
    df[["group1"]] <- factor(df[["group1"]], levels = unique(df[["group1"]]))
  }
  # diff_pct
  if (!"diff_pct" %in% colnames(df)) {
    if ("pct.1" %in% colnames(df) && "pct.2" %in% colnames(df)) {
      df[["diff_pct"]] <- df[["pct.1"]] - df[["pct.2"]]
    } else {
      df[["diff_pct"]] <- 0
    }
  }
  df[["-log10padj"]] <- -log10(df[["p_val_adj"]])
  rownames(df) <- NULL

  df
}

#' Clip log2FC to symmetric percentile range
#' @keywords internal
.clip_log2fc_symmetric <- function(df, fc_col = "avg_log2FC") {
  fc <- df[[fc_col]][is.finite(df[[fc_col]])]
  if (length(fc) == 0) return(list(df = df, fc_lim = c(-1, 1)))
  x_upper <- stats::quantile(fc, c(0.99, 1))
  x_lower <- stats::quantile(fc, c(0.01, 0))
  x_upper <- ifelse(x_upper[1] > 0, x_upper[1], x_upper[2])
  x_lower <- ifelse(x_lower[1] < 0, x_lower[1], x_lower[2])
  if (x_upper > 0 && x_lower < 0) {
    value_range <- min(abs(c(x_upper, x_lower)), na.rm = TRUE)
    x_upper <- value_range
    x_lower <- -value_range
  }
  df[df[[fc_col]] > x_upper, fc_col] <- x_upper
  df[df[[fc_col]] < x_lower, fc_col] <- x_lower
  list(df = df, fc_lim = c(x_lower, x_upper))
}

#' Get top markers for labeling
#' @keywords internal
.get_top_markers_for_label <- function(df, cluster_levels, nlabel,
                                       features_label = NULL) {
  top_list <- lapply(cluster_levels, function(x) {
    tmp <- df[df[["group1"]] == x, , drop = FALSE]
    if (nrow(tmp) == 0) return(NULL)
    if (is.null(features_label)) {
      top_up <- utils::head(
        tmp[order(tmp[["avg_log2FC"]], decreasing = TRUE), , drop = FALSE],
        nlabel
      )
      top_dn <- utils::head(
        tmp[order(tmp[["avg_log2FC"]], decreasing = FALSE), , drop = FALSE],
        nlabel
      )
      rbind(top_up, top_dn)
    } else {
      tmp[tmp[["gene"]] %in% features_label, , drop = FALSE]
    }
  })
  do.call(rbind, top_list[!vapply(top_list, is.null, logical(1))])
}

# --------------------------------------------------------------------------
# PlotDE main function
# --------------------------------------------------------------------------

#' Plot Differential Expression Results
#'
#' Unified visualization for differential expression analysis. Supports
#' volcano, MA, pct (diff_pct volcano), manhattan, ring, and PCA +
#' hierarchical clustering plots.
#'
#' @param x For \code{type = "volcano"}, \code{"ma"}, \code{"pct"},
#'   \code{"manhattan"}, \code{"ring"}: a data.frame of DE results
#'   (e.g., from \code{RunDE}).
#'   For \code{type = "pca_hc"}: a Seurat object.
#' @param type Plot type. One of \code{"volcano"}, \code{"ma"}, \code{"pct"},
#'   \code{"pct_mirror"}, \code{"manhattan"}, \code{"ring"}, or
#'   \code{"pca_hc"}.
#' @param fc_threshold Log2 fold change threshold for significance lines
#'   (volcano/ma). Default: 1.
#' @param p_threshold Adjusted p-value threshold. Default: 0.05.
#' @param highlight_fc_threshold FC threshold for gene labels
#'   (volcano/ma). Default: 3.
#' @param point_size Point size. Default: 1.2.
#' @param highlight_size Size of highlighted points. Default: 3.
#' @param colors Color vector for Up/NS/Down/Highlight groups. Default:
#'   \code{c("darkred", "grey50", "royalblue4", "green2")}.
#' @param symmetry Logical; symmetric x-axis for volcano. Default: FALSE.
#' @param label_genes Character vector of gene names to label manually.
#'   Default: NULL (auto-select by FC).
#' @param nlabel Number of top genes to label per group (manhattan/ring).
#'   Default: 5.
#' @param jitter_width Jitter width for manhattan/ring. Default: 0.5.
#' @param jitter_height Jitter height for manhattan. Default: 0.4.
#' @param label_size Label text size. Default: 4.
#' @param group.by Column name for grouping (pca_hc). Default: NULL.
#' @param reduction Reduction to use for PCA (pca_hc). Default: "pca".
#' @param dims Dimensions to plot (pca_hc). Default: 1:2.
#' @param dist_method Distance method for dendrogram (pca_hc).
#'   Default: "euclidean".
#' @param hclust_method Clustering method (pca_hc). Default: "ward.D2".
#' @param k Number of clusters for dendrogram coloring. Default: NULL
#'   (auto = number of groups).
#' @param add_ellipses Logical; add confidence ellipses (pca_hc).
#'   Default: TRUE.
#' @param combine Logical; combine multi-group volcano plots.
#'   Default: TRUE.
#' @param ncol Number of columns for combined facets. Default: NULL.
#' @param y_max Numeric; cap -log10(p) y-axis at this value (volcano/pct).
#'   Default: NULL (auto: 99.5th percentile).
#'
#' @return A ggplot object (or combined plot for pca_hc).
#'
#' @examples
#' \dontrun{
#' markers <- RunDE(seu, group.by = "celltype")
#' PlotDE(markers, type = "volcano")
#' PlotDE(markers, type = "pct")
#' PlotDE(markers, type = "manhattan")
#' PlotDE(seu, type = "pca_hc", group.by = "celltype")
#' }
#'
#' @seealso \code{\link{RunDE}}
#' @export
PlotDE <- function(x,
                   type = c("volcano", "ma", "pct", "pct_mirror", "manhattan", "ring", "pca_hc"),
                   fc_threshold = 1,
                   p_threshold = 0.05,
                   highlight_fc_threshold = 3,
                   point_size = 1.2,
                   highlight_size = 3,
                   colors = c("darkred", "grey50", "royalblue4", "green2"),
                   symmetry = FALSE,
                   label_genes = NULL,
                   nlabel = 5,
                   jitter_width = 0.5,
                   jitter_height = 0.4,
                   label_size = 4,
                   group.by = NULL,
                   reduction = "pca",
                   dims = 1:2,
                   dist_method = "euclidean",
                   hclust_method = "ward.D2",
                   k = NULL,
                   add_ellipses = TRUE,
                   combine = TRUE,
                   ncol = NULL,
                   y_max = NULL) {
  type <- match.arg(type)

  switch(type,
    volcano = .plot_de_volcano(x, fc_threshold, p_threshold,
                               highlight_fc_threshold, point_size,
                               highlight_size, colors, symmetry,
                               label_genes, combine, ncol, y_max),
    ma      = .plot_de_ma(x, fc_threshold, p_threshold,
                          highlight_fc_threshold, point_size,
                          highlight_size, colors, label_genes),
    pct     = .plot_de_pct(x, fc_threshold, p_threshold, point_size,
                            nlabel, label_genes, label_size,
                            y_max, ncol, mirror = FALSE),
    pct_mirror = .plot_de_pct(x, fc_threshold, p_threshold, point_size,
                               nlabel, label_genes, label_size,
                               y_max, ncol, mirror = TRUE),
    manhattan = .plot_de_manhattan(x, nlabel, label_genes, point_size,
                                  jitter_width, jitter_height,
                                  label_size, colors),
    ring    = .plot_de_ring(x, nlabel, label_genes, point_size,
                            jitter_width, label_size, colors),
    pca_hc  = .plot_de_pca_hc(x, group.by, reduction, dims,
                              dist_method, hclust_method, k,
                              add_ellipses, point_size)
  )
}

# --------------------------------------------------------------------------
# Volcano (BulkS4-style)
# --------------------------------------------------------------------------
#' @keywords internal
.plot_de_volcano <- function(df, fc_threshold, p_threshold,
                             highlight_fc_threshold, point_size,
                             highlight_size, colors, symmetry,
                             label_genes, combine, ncol, y_max) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install with: install.packages('ggrepel')")
  }
  df <- .get_de_data(df)

  # Cap p-values: replace 0 with minimum non-zero value to avoid -log10(0)=Inf
  nonzero_p <- df[["p_val_adj"]][df[["p_val_adj"]] > 0]
  if (length(nonzero_p) > 0) {
    min_p <- min(nonzero_p, na.rm = TRUE)
    df[["p_val_adj"]][df[["p_val_adj"]] == 0] <- min_p
  }

  df[["DE_group"]] <- dplyr::case_when(
    df[["p_val_adj"]] < p_threshold & df[["avg_log2FC"]] > fc_threshold ~ "Up",
    df[["p_val_adj"]] < p_threshold & df[["avg_log2FC"]] < -fc_threshold ~ "Down",
    TRUE ~ "NS"
  )
  df[["DE_group"]] <- factor(df[["DE_group"]], levels = c("Up", "NS", "Down"))

  # Compute -log10(p) and cap at y_max (quantile-based)
  df[["neg_log10_p"]] <- -log10(df[["p_val_adj"]])
  if (is.null(y_max)) {
    y_max <- quantile(df[["neg_log10_p"]][is.finite(df[["neg_log10_p"]])],
                      0.995, na.rm = TRUE)
    y_max <- max(y_max, -log10(p_threshold) + 1)
  }
  df[["neg_log10_p"]] <- pmin(df[["neg_log10_p"]], y_max)

  # Highlight genes
  if (is.null(label_genes)) {
    highlight <- df[abs(df[["avg_log2FC"]]) >= highlight_fc_threshold &
                      df[["p_val_adj"]] < p_threshold, , drop = FALSE]
  } else {
    highlight <- df[df[["gene"]] %in% label_genes, , drop = FALSE]
  }

  up_n <- sum(df[["DE_group"]] == "Up", na.rm = TRUE)
  dn_n <- sum(df[["DE_group"]] == "Down", na.rm = TRUE)
  title <- paste0(dn_n, " down, ", up_n, " up")

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["avg_log2FC"]],
    y = .data[["neg_log10_p"]],
    color = .data[["DE_group"]]
  )) +
    ggplot2::geom_point(alpha = 0.8, size = point_size) +
    ggplot2::scale_color_manual(values = stats::setNames(colors[1:3],
                                                         c("Up", "NS", "Down"))) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold),
                        lty = 4, lwd = 0.6, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                        lty = 4, lwd = 0.6, alpha = 0.8) +
    ggplot2::labs(x = "log2(Fold Change)", y = "-log10(adj.P)",
                  title = title, color = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "top"
    )

  if (nrow(highlight) > 0) {
    p <- p +
      ggplot2::geom_point(data = highlight, alpha = 0.8,
                          size = highlight_size, color = colors[4]) +
      ggrepel::geom_text_repel(
        data = highlight,
        ggplot2::aes(label = .data[["gene"]]),
        color = "black", alpha = 0.8, max.overlaps = 20
      )
  }

  if (symmetry) {
    ce <- ceiling(max(abs(df[["avg_log2FC"]]), na.rm = TRUE))
    p <- p + ggplot2::scale_x_continuous(limits = c(-ce, ce),
                                         expand = c(0, 0))
  }

  # Facet by group if multiple groups
  grps <- levels(df[["group1"]])
  if (length(grps) > 1 && !all(grps == "All")) {
    p <- p + ggplot2::facet_wrap(~ group1, ncol = ncol)
  }

  p
}

# --------------------------------------------------------------------------
# MA (BulkS4-style)
# --------------------------------------------------------------------------
#' @keywords internal
.plot_de_ma <- function(df, fc_threshold, p_threshold,
                        highlight_fc_threshold, point_size,
                        highlight_size, colors, label_genes) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install with: install.packages('ggrepel')")
  }
  df <- .get_de_data(df)

  # Compute x-axis: average expression
  if ("AveExpr" %in% colnames(df)) {
    df[["ave_expr"]] <- df[["AveExpr"]]
  } else if ("baseMean" %in% colnames(df)) {
    df[["ave_expr"]] <- log2(df[["baseMean"]] + 1)
  } else if ("logCPM" %in% colnames(df)) {
    df[["ave_expr"]] <- df[["logCPM"]]
  } else if ("pct.1" %in% colnames(df) && "pct.2" %in% colnames(df)) {
    df[["ave_expr"]] <- (df[["pct.1"]] + df[["pct.2"]]) / 2
  } else {
    stop("Cannot find average expression column. ",
         "Provide AveExpr, baseMean, logCPM, or pct.1/pct.2")
  }

  df[["DE_group"]] <- dplyr::case_when(
    df[["p_val_adj"]] < p_threshold & df[["avg_log2FC"]] > fc_threshold ~ "Up",
    df[["p_val_adj"]] < p_threshold & df[["avg_log2FC"]] < -fc_threshold ~ "Down",
    TRUE ~ "NS"
  )
  df[["DE_group"]] <- factor(df[["DE_group"]], levels = c("Up", "NS", "Down"))

  if (is.null(label_genes)) {
    highlight <- df[abs(df[["avg_log2FC"]]) >= highlight_fc_threshold &
                      df[["p_val_adj"]] < p_threshold, , drop = FALSE]
  } else {
    highlight <- df[df[["gene"]] %in% label_genes, , drop = FALSE]
  }

  up_n <- sum(df[["DE_group"]] == "Up", na.rm = TRUE)
  dn_n <- sum(df[["DE_group"]] == "Down", na.rm = TRUE)
  title <- paste0(dn_n, " down, ", up_n, " up")

  x_lab <- if ("pct.1" %in% colnames(df) && !"AveExpr" %in% colnames(df) &&
               !"baseMean" %in% colnames(df) && !"logCPM" %in% colnames(df)) {
    "Average Percent Expressed"
  } else {
    "log2(Average Expression)"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["ave_expr"]],
    y = .data[["avg_log2FC"]],
    color = .data[["DE_group"]]
  )) +
    ggplot2::geom_point(alpha = 0.8, size = point_size) +
    ggplot2::scale_color_manual(values = stats::setNames(colors[1:3],
                                                         c("Up", "NS", "Down"))) +
    ggplot2::geom_hline(yintercept = c(fc_threshold, -fc_threshold),
                        lty = 2, lwd = 1) +
    ggplot2::geom_hline(yintercept = 0, lwd = 1.2) +
    ggplot2::labs(x = x_lab, y = "log2(Fold Change)",
                  title = title, color = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "top"
    )

  if (nrow(highlight) > 0) {
    p <- p +
      ggplot2::geom_point(data = highlight, alpha = 0.8,
                          size = highlight_size, color = colors[4]) +
      ggrepel::geom_text_repel(
        data = highlight,
        ggplot2::aes(label = .data[["gene"]]),
        color = "black", alpha = 0.8, max.overlaps = 20
      )
  }

  # Facet by group if multiple groups
  grps <- levels(df[["group1"]])
  if (length(grps) > 1 && !all(grps == "All")) {
    p <- p + ggplot2::facet_wrap(~ group1)
  }

  p
}

# --------------------------------------------------------------------------
# Pct volcano (diff_pct x-axis, log2FC color gradient)
# --------------------------------------------------------------------------
#' @keywords internal
.plot_de_pct <- function(df, fc_threshold, p_threshold, point_size,
                         nlabel, label_genes, label_size, y_max, ncol,
                         mirror) {
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Package 'ggrepel' is required. Install with: install.packages('ggrepel')")
  if (!requireNamespace("scales", quietly = TRUE))
    stop("Package 'scales' is required. Install with: install.packages('scales')")

  df <- .get_de_data(df)

  # Ensure diff_pct exists
  if (!"diff_pct" %in% colnames(df)) {
    if ("pct.1" %in% colnames(df) && "pct.2" %in% colnames(df)) {
      df[["diff_pct"]] <- df[["pct.1"]] - df[["pct.2"]]
    } else {
      stop("Cannot compute diff_pct: columns 'pct.1' and 'pct.2' not found.",
           call. = FALSE)
    }
  }

  # Cap p-values
  nonzero_p <- df[["p_val_adj"]][df[["p_val_adj"]] > 0]
  if (length(nonzero_p) > 0) {
    min_p <- min(nonzero_p, na.rm = TRUE)
    df[["p_val_adj"]][df[["p_val_adj"]] == 0] <- min_p
  }
  df[["neg_log10_p"]] <- -log10(df[["p_val_adj"]])
  if (is.null(y_max)) {
    y_max <- quantile(df[["neg_log10_p"]][is.finite(df[["neg_log10_p"]])],
                      0.995, na.rm = TRUE)
    y_max <- max(y_max, -log10(p_threshold) + 1)
  }
  df[["neg_log10_p"]] <- pmin(df[["neg_log10_p"]], y_max)

  # Mirror: flip y for negative diff_pct
  if (mirror) {
    df[["neg_log10_p"]] <- ifelse(df[["diff_pct"]] < 0,
                                  -df[["neg_log10_p"]],
                                  df[["neg_log10_p"]])
  }

  # Clamp log2FC for color scale
  fc_cap <- quantile(abs(df[["avg_log2FC"]]), 0.99, na.rm = TRUE)
  fc_cap <- max(fc_cap, fc_threshold)
  df[["fc_color"]] <- pmax(pmin(df[["avg_log2FC"]], fc_cap), -fc_cap)

  # Select label genes
  cluster_levels <- unique(df[["group1"]])
  if (is.null(label_genes)) {
    lab_df <- .get_top_markers_for_label(df, cluster_levels, nlabel = nlabel,
                                         features_label = NULL)
  } else {
    lab_df <- df[df[["gene"]] %in% label_genes, , drop = FALSE]
  }

  # Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["diff_pct"]],
    y = .data[["neg_log10_p"]],
    color = .data[["fc_color"]]
  )) +
    ggplot2::geom_point(alpha = 0.8, size = point_size) +
    ggplot2::scale_color_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-fc_cap, fc_cap),
      oob = scales::squish,
      name = "log2FC"
    ) +
    ggplot2::geom_hline(
      yintercept = if (mirror) c(-log10(p_threshold), log10(p_threshold))
                   else -log10(p_threshold),
      lty = 4, lwd = 0.6, alpha = 0.8
    ) +
    ggplot2::labs(
      x = "diff_pct",
      y = if (mirror) "sign(diff_pct) * -log10(p-adjust)"
          else "-log10(p-adjust)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.border = ggplot2::element_rect(colour = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  # Labels
  if (!is.null(lab_df) && nrow(lab_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = .data[["gene"]]),
      size = label_size, color = "black",
      max.overlaps = 20, segment.color = "grey50"
    )
  }

  # Facet by group
  grps <- unique(df[["group1"]])
  if (length(grps) > 1 && !all(grps == "All")) {
    p <- p + ggplot2::facet_wrap(~ group1, scales = "free", ncol = ncol)
  }

  p
}

# --------------------------------------------------------------------------
# Manhattan plot
# --------------------------------------------------------------------------
#' @keywords internal
.plot_de_manhattan <- function(df, nlabel, features_label, point_size,
                               jitter_width, jitter_height,
                               label_size, colors) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }
  df <- .get_de_data(df)
  clip_res <- .clip_log2fc_symmetric(df)
  df <- clip_res$df
  fc_lim <- clip_res$fc_lim

  cluster_levels <- levels(df[["group1"]])
  n_m <- nrow(df)
  df[["x_num"]] <- as.numeric(df[["group1"]])
  df[["x_plot"]] <- df[["x_num"]] + (stats::runif(n_m) - 0.5) * jitter_width
  df[["y_plot"]] <- df[["avg_log2FC"]] + (stats::runif(n_m) - 0.5) * jitter_height

  top_marker <- .get_top_markers_for_label(df, cluster_levels, nlabel,
                                            features_label)
  if (!is.null(top_marker) && nrow(top_marker) > 0) {
    top_marker[["x_plot"]] <- df[match(
      paste(top_marker[["gene"]], top_marker[["group1"]]),
      paste(df[["gene"]], df[["group1"]])
    ), "x_plot"]
    top_marker[["y_plot"]] <- df[match(
      paste(top_marker[["gene"]], top_marker[["group1"]]),
      paste(df[["gene"]], df[["group1"]])
    ), "y_plot"]
  }

  # Background bars
  back_list <- lapply(cluster_levels, function(x) {
    tmp <- df[df[["group1"]] == x, , drop = FALSE]
    if (nrow(tmp) == 0) return(NULL)
    data.frame(cluster = x,
               min = min(tmp[["avg_log2FC"]], na.rm = TRUE) - 0.2,
               max = max(tmp[["avg_log2FC"]], na.rm = TRUE) + 0.2)
  })
  back_data <- do.call(rbind, back_list[!vapply(back_list, is.null, logical(1))])
  back_data[["x_num"]] <- match(back_data[["cluster"]], cluster_levels)

  # Group palette
  n_grp <- length(cluster_levels)
  grp_colors <- if (requireNamespace("scales", quietly = TRUE)) {
    scales::hue_pal()(n_grp)
  } else {
    grDevices::rainbow(n_grp)
  }
  names(grp_colors) <- cluster_levels

  tile_data <- data.frame(
    group1 = cluster_levels,
    x = seq_along(cluster_levels),
    y = 0
  )

  # RdBu-like gradient
  gradient_colors <- c("royalblue4", "grey90", "darkred")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x_plot"]],
                                         y = .data[["y_plot"]])) +
    ggplot2::geom_col(
      data = back_data,
      ggplot2::aes(x = .data[["x_num"]], y = .data[["min"]]),
      fill = "white", inherit.aes = FALSE
    ) +
    ggplot2::geom_col(
      data = back_data,
      ggplot2::aes(x = .data[["x_num"]], y = .data[["max"]]),
      fill = "white", inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["avg_log2FC"]]),
      size = point_size, alpha = 0.8
    ) +
    ggplot2::scale_color_gradientn(
      name = "log2FC",
      colors = gradient_colors,
      values = scales::rescale(unique(c(fc_lim[1], 0, fc_lim[2]))),
      limits = fc_lim,
      guide = ggplot2::guide_colorbar(
        frame.colour = "black", ticks.colour = "black",
        title.hjust = 0, order = 1
      )
    ) +
    ggplot2::geom_tile(
      data = tile_data,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   fill = .data[["group1"]]),
      color = "black", height = 0.5, show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = grp_colors) +
    ggplot2::geom_text(
      data = tile_data,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   label = .data[["group1"]]),
      inherit.aes = FALSE, size = 3, color = "black"
    ) +
    ggplot2::scale_x_continuous(breaks = seq_along(cluster_levels),
                                labels = cluster_levels) +
    ggplot2::scale_y_continuous(n.breaks = 6) +
    ggplot2::labs(x = NULL, y = "Average log2FoldChange") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (!is.null(top_marker) && nrow(top_marker) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_marker,
      ggplot2::aes(x = .data[["x_plot"]], y = .data[["y_plot"]],
                   label = .data[["gene"]]),
      inherit.aes = FALSE, min.segment.length = 0,
      max.overlaps = 100, segment.colour = "grey40",
      color = "black", size = label_size, force = 20
    )
  }

  p
}

# --------------------------------------------------------------------------
# Ring plot
# --------------------------------------------------------------------------
#' @keywords internal
.plot_de_ring <- function(df, nlabel, features_label, point_size,
                          jitter_width, label_size, colors,
                          seed = 11) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }
  if (!requireNamespace("geomtextpath", quietly = TRUE)) {
    stop("Package 'geomtextpath' is required for ring plot. ",
         "Install with: install.packages('geomtextpath')")
  }
  df <- .get_de_data(df)
  clip_res <- .clip_log2fc_symmetric(df)
  df <- clip_res$df
  fc_lim <- clip_res$fc_lim

  cluster_levels <- levels(df[["group1"]])
  n_grp <- length(cluster_levels)
  set.seed(seed)

  ring_r0 <- 2
  max_abs_fc <- max(abs(df[["avg_log2FC"]]), na.rm = TRUE)
  ring_k <- if (max_abs_fc > 0) 1.2 / max_abs_fc else 0.2
  tile_height <- 0.3
  tile_gap <- 0.1

  n_m <- nrow(df)
  df[["x_num"]] <- as.numeric(df[["group1"]])
  df[["x_angle"]] <- df[["x_num"]] + (stats::runif(n_m) - 0.5) * jitter_width
  df[["y_radius"]] <- ring_r0 + ring_k * df[["avg_log2FC"]]

  band_lo <- ring_r0 - tile_height / 2
  band_hi <- ring_r0 + tile_height / 2
  in_band <- df[["y_radius"]] >= band_lo & df[["y_radius"]] <= band_hi
  df[in_band & df[["y_radius"]] < ring_r0, "y_radius"] <- band_lo - tile_gap
  df[in_band & df[["y_radius"]] >= ring_r0, "y_radius"] <- band_hi + tile_gap

  top_marker <- .get_top_markers_for_label(df, cluster_levels, nlabel,
                                            features_label)
  if (!is.null(top_marker) && nrow(top_marker) > 0) {
    idx <- match(
      paste(top_marker[["gene"]], top_marker[["group1"]]),
      paste(df[["gene"]], df[["group1"]])
    )
    top_marker[["x_angle"]] <- df[idx, "x_angle"]
    top_marker[["y_radius"]] <- df[idx, "y_radius"]
  }

  grp_colors <- if (requireNamespace("scales", quietly = TRUE)) {
    scales::hue_pal()(n_grp)
  } else {
    grDevices::rainbow(n_grp)
  }
  names(grp_colors) <- cluster_levels

  tile_data <- data.frame(
    x = seq_len(n_grp), y = ring_r0,
    group1 = cluster_levels, label = cluster_levels
  )

  npt <- 40
  path_margin <- 0.04
  path_df <- do.call(rbind, lapply(seq_len(n_grp), function(i) {
    x_start <- i - 0.5 + path_margin
    x_end <- i + 0.5 - path_margin
    data.frame(
      x = x_start + (0:(npt - 1)) / (npt - 1) * (x_end - x_start),
      y = ring_r0, label = cluster_levels[i], group = i
    )
  }))

  gradient_colors <- c("royalblue4", "grey90", "darkred")

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["x_angle"]], y = .data[["y_radius"]]
  )) +
    ggplot2::geom_vline(
      xintercept = seq(0.5, n_grp - 0.5, by = 1),
      color = "grey85", linewidth = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = ring_r0, linetype = 2,
      color = "grey40", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[["avg_log2FC"]]),
      size = point_size, alpha = 0.8
    ) +
    ggplot2::scale_color_gradientn(
      name = "log2FC", colors = gradient_colors,
      values = scales::rescale(unique(c(fc_lim[1], 0, fc_lim[2]))),
      limits = fc_lim,
      guide = ggplot2::guide_colorbar(
        frame.colour = "black", ticks.colour = "black",
        title.hjust = 0, order = 1
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0.5, n_grp + 0.5),
      breaks = seq_len(n_grp), labels = NULL
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA), n.breaks = 5) +
    ggplot2::coord_polar(theta = "x", start = -pi / 2, direction = 1) +
    ggplot2::geom_tile(
      data = tile_data,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   fill = .data[["group1"]]),
      color = "black", height = tile_height,
      show.legend = FALSE, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = grp_colors) +
    geomtextpath::geom_textpath(
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   label = .data[["label"]], group = .data[["group"]]),
      data = path_df, inherit.aes = FALSE,
      size = 3, color = "black", upright = TRUE
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1,
      legend.position = "right",
      legend.background = ggplot2::element_blank()
    )

  if (!is.null(top_marker) && nrow(top_marker) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_marker,
      ggplot2::aes(x = .data[["x_angle"]], y = .data[["y_radius"]],
                   label = .data[["gene"]]),
      inherit.aes = FALSE, min.segment.length = 0,
      max.overlaps = 100, segment.colour = "grey40",
      color = "black", size = label_size, force = 5
    )
  }

  p
}

# --------------------------------------------------------------------------
# PCA + Hierarchical Clustering (BulkS4-style, adapted for Seurat)
# --------------------------------------------------------------------------
#' @keywords internal
.plot_de_pca_hc <- function(seu, group.by, reduction, dims,
                            dist_method, hclust_method, k,
                            add_ellipses, point_size) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  if (!requireNamespace("factoextra", quietly = TRUE)) {
    stop("Package 'factoextra' is required. Install with: install.packages('factoextra')")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  if (is.null(group.by)) {
    stop("'group.by' must be specified for pca_hc")
  }

  # PCA plot from Seurat embeddings
  embeddings <- Seurat::Embeddings(seu, reduction = reduction)[, dims]
  meta <- seu@meta.data
  pca_df <- data.frame(
    PC1 = embeddings[, 1],
    PC2 = embeddings[, 2],
    group = meta[[group.by]]
  )
  pc_labels <- colnames(embeddings)

  p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(
    x = .data[["PC1"]], y = .data[["PC2"]],
    color = .data[["group"]]
  )) +
    ggplot2::geom_point(size = point_size, alpha = 0.6) +
    ggplot2::labs(x = pc_labels[1], y = pc_labels[2], color = group.by) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom")

  if (add_ellipses) {
    p_pca <- p_pca + ggplot2::stat_ellipse(level = 0.95, linewidth = 0.8)
  }

  # Dendrogram from pseudo-bulk
  assay <- Seurat::DefaultAssay(seu)
  expr <- Seurat::GetAssayData(seu, layer = "data", assay = assay)
  groups <- meta[[group.by]]
  group_levels <- unique(groups)

  # Average expression per group
  avg_expr <- do.call(cbind, lapply(group_levels, function(g) {
    cells <- which(groups == g)
    Matrix::rowMeans(expr[, cells, drop = FALSE])
  }))
  colnames(avg_expr) <- group_levels

  dist_mat <- stats::dist(t(avg_expr), method = dist_method)
  hc <- stats::hclust(dist_mat, method = hclust_method)

  if (is.null(k)) k <- length(group_levels)

  p_dend <- factoextra::fviz_dend(
    hc, k = k, cex = 0.8,
    horiz = FALSE, type = "rectangle",
    color_labels_by_k = TRUE,
    rect = TRUE, rect_fill = TRUE
  )

  patchwork::wrap_plots(p_dend, p_pca, ncol = 2)
}


# ============================================================================
# RunGsea: GSEA analysis on RunDE results
# ============================================================================

#' Run Gene Set Enrichment Analysis on DE Results
#'
#' Perform GSEA for each group in a \code{\link{RunDE}} result data.frame.
#' Each group is analyzed independently using \code{clusterProfiler::GSEA()},
#' and the individual \code{gseaResult} objects are preserved for downstream
#' visualization (e.g., classic GSEA running score plots via
#' \code{enrichplot::gseaplot2()}).
#'
#' @param de_result A data.frame from \code{\link{RunDE}}, containing at
#'   least \code{gene} and \code{score.by} columns.
#' @param geneset Gene set specification. Supports multiple formats:
#'   \describe{
#'     \item{msigdbr args}{A list with msigdbr parameters, e.g.
#'       \code{list(species = "Homo sapiens", collection = "H")}.}
#'     \item{data.frame}{With columns \code{gs_name} and \code{gene_symbol}.}
#'     \item{named list}{E.g. \code{list(PathA = c("TP53", "BRCA1"))}.}
#'     \item{GMT file}{Character path to a \code{.gmt} file.}
#'   }
#' @param group.by Column in \code{de_result} to split groups. Default:
#'   \code{"group1"}.
#' @param score.by Column used as gene ranking score. Default:
#'   \code{"avg_log2FC"}.
#' @param de_threshold Optional character string expression to filter
#'   \code{de_result} before GSEA (e.g. \code{"p_val_adj < 0.05"}).
#'   Default: \code{NULL} (use all genes, standard GSEA practice).
#' @param scoreType GSEA score type: \code{"std"}, \code{"pos"}, or
#'   \code{"neg"}. Default: \code{"std"}.
#' @param minGSSize Minimum gene set size. Default: 10.
#' @param maxGSSize Maximum gene set size. Default: 500.
#' @param pvalueCutoff P-value cutoff for GSEA results. Default: 1
#'   (no filtering).
#' @param p.adjust.method P-value adjustment method. Default: \code{"BH"}.
#' @param clean.names Logical; clean pathway names by removing common
#'   prefixes (e.g. HALLMARK_). Default: TRUE.
#' @param prefix Optional prefix to remove from pathway names.
#' @param cores Integer; number of parallel cores. Default: 1.
#' @param seed Random seed. Default: 11.
#' @param verbose Logical; print progress. Default: TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{enrichment}}{Data.frame of all GSEA results across groups.}
#'     \item{\code{results}}{Named list of individual \code{gseaResult} objects
#'       (one per group), usable with \code{enrichplot::gseaplot2()}.}
#'     \item{\code{input}}{Data.frame of the (filtered) input gene lists.}
#'   }
#'
#' @examples
#' \dontrun{
#' markers <- RunDE(seu, group.by = "celltype")
#' res <- RunGsea(markers,
#'                geneset = list(species = "Homo sapiens", collection = "H"))
#'
#' # Summary table
#' head(res$enrichment)
#'
#' # Classic GSEA running score plot
#' enrichplot::gseaplot2(res$results[["T_cell"]], geneSetID = 1)
#'
#' # Custom GMT file
#' res2 <- RunGsea(markers, geneset = "path/to/custom.gmt")
#' }
#'
#' @seealso \code{\link{RunDE}}, \code{\link{RunGseaEnrich}}
#' @export
RunGsea <- function(de_result,
                    geneset = list(species = "Homo sapiens", collection = "H"),
                    group.by = "group1",
                    score.by = "avg_log2FC",
                    de_threshold = NULL,
                    scoreType = c("std", "pos", "neg"),
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    p.adjust.method = "BH",
                    clean.names = TRUE,
                    prefix = NULL,
                    cores = 1,
                    seed = 11,
                    verbose = TRUE) {

  scoreType <- match.arg(scoreType)
  set.seed(seed)

  # --- dependency check ---
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Package 'clusterProfiler' is required.\n",
         "  BiocManager::install('clusterProfiler')")

  # --- validate input ---
  if (!is.data.frame(de_result))
    stop("'de_result' must be a data.frame (e.g., from RunDE).", call. = FALSE)
  if (!"gene" %in% colnames(de_result))
    stop("'de_result' must contain a 'gene' column.", call. = FALSE)
  if (!score.by %in% colnames(de_result))
    stop("Column '", score.by, "' not found in de_result.", call. = FALSE)
  if (!group.by %in% colnames(de_result))
    stop("Column '", group.by, "' not found in de_result.", call. = FALSE)

  # --- optional DE filtering ---
  if (!is.null(de_threshold)) {
    n_before <- nrow(de_result)
    de_result <- de_result[
      with(de_result, eval(rlang::parse_expr(de_threshold))), , drop = FALSE
    ]
    if (verbose) message("DE filter (", de_threshold, "): ",
                         n_before, " -> ", nrow(de_result), " genes")
  }

  if (nrow(de_result) == 0) {
    warning("No genes remaining after filtering.")
    return(list(enrichment = data.frame(), results = list(),
                input = data.frame()))
  }

  # --- parse gene sets ---
  if (verbose) message("Parsing gene sets...")
  TERM2GENE <- .parse_geneset_to_term2gene(geneset)
  if (verbose) message("  ", length(unique(TERM2GENE[["gs_name"]])),
                        " gene sets, ",
                        nrow(TERM2GENE), " term-gene pairs")

  # --- handle Inf scores ---
  score_vec <- de_result[[score.by]]
  if (any(is.infinite(score_vec))) {
    finite_scores <- score_vec[is.finite(score_vec)]
    score_vec[is.infinite(score_vec) & score_vec > 0] <- max(finite_scores)
    score_vec[is.infinite(score_vec) & score_vec < 0] <- min(finite_scores)
    de_result[[score.by]] <- score_vec
  }

  # --- save input ---
  input_df <- de_result[, intersect(
    c("gene", score.by, group.by, "p_val_adj"),
    colnames(de_result)
  ), drop = FALSE]

  # --- group levels ---
  groups <- if (is.factor(de_result[[group.by]])) {
    levels(de_result[[group.by]])
  } else {
    unique(as.character(de_result[[group.by]]))
  }
  groups <- groups[groups %in% de_result[[group.by]]]

  if (verbose) message("Running GSEA for ", length(groups), " groups...")

  # --- GSEA per group ---
  .run_one_group <- function(grp) {
    df_grp <- de_result[as.character(de_result[[group.by]]) == grp, ,
                        drop = FALSE]
    if (nrow(df_grp) == 0) return(NULL)

    # Build named score vector
    scores <- stats::setNames(df_grp[[score.by]], df_grp[["gene"]])

    # Deduplicate: keep gene with max abs(score)
    if (anyDuplicated(names(scores))) {
      scores <- tapply(scores, names(scores), function(x) x[which.max(abs(x))])
    }

    # Remove NA
    scores <- scores[!is.na(scores)]
    if (length(scores) == 0) return(NULL)

    # Sort descending
    scores <- sort(scores, decreasing = TRUE)

    # Auto-detect scoreType
    local_scoreType <- scoreType
    if (all(scores > 0) && scoreType == "std") local_scoreType <- "pos"
    if (all(scores < 0) && scoreType == "std") local_scoreType <- "neg"

    # Run GSEA
    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList      = scores,
        TERM2GENE     = TERM2GENE,
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        pvalueCutoff  = pvalueCutoff,
        pAdjustMethod = p.adjust.method,
        scoreType     = local_scoreType,
        eps           = 0,
        nPermSimple   = 10000,
        by            = "fgsea",
        verbose       = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  [", grp, "] Error: ", conditionMessage(e))
      return(NULL)
    })

    if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
      if (verbose) message("  [", grp, "] No enriched pathways")
      return(NULL)
    }

    # Add group info
    gsea_res@result[["Group"]] <- grp

    if (verbose) {
      n_sig <- sum(gsea_res@result[["p.adjust"]] < 0.05)
      message("  [", grp, "] ", nrow(gsea_res@result),
              " pathways (", n_sig, " significant)")
    }

    gsea_res
  }

  # Execute: parallel or sequential
  cores <- max(1L, as.integer(cores))
  if (cores > 1 && .Platform$OS.type == "unix") {
    results_list <- parallel::mclapply(groups, .run_one_group,
                                       mc.cores = cores)
  } else {
    results_list <- lapply(groups, .run_one_group)
  }
  names(results_list) <- groups

  # Remove NULLs
  results_list <- results_list[!vapply(results_list, is.null, logical(1))]

  if (length(results_list) == 0) {
    warning("No GSEA results for any group.")
    return(list(enrichment = data.frame(), results = list(),
                input = input_df))
  }

  # --- clean pathway names ---
  if (clean.names) {
    for (nm in names(results_list)) {
      results_list[[nm]]@result[["Description"]] <- .clean_pathway_names(
        results_list[[nm]]@result[["Description"]], prefix = prefix
      )
    }
  }

  # --- merge enrichment table ---
  enrichment <- do.call(rbind, lapply(results_list, function(x) x@result))
  rownames(enrichment) <- NULL

  if (verbose) {
    message("Done. ", nrow(enrichment), " total pathway-group results, ",
            sum(enrichment[["p.adjust"]] < 0.05, na.rm = TRUE), " significant.")
  }

  list(
    enrichment = enrichment,
    results    = results_list,
    input      = input_df
  )
}


# ============================================================================
# PlotGsea — Visualise RunGsea results
# ============================================================================

#' Plot GSEA Results
#'
#' Generates various types of plots for Gene Set Enrichment Analysis (GSEA)
#' results produced by \code{\link{RunGsea}}.
#'
#' @param gsea_result A list returned by \code{RunGsea}, containing
#'   \code{enrichment} (data.frame), \code{results} (named list of
#'   \code{gseaResult} objects), and \code{input}.
#' @param type The type of plot to generate. One of:
#'   \describe{
#'     \item{\code{"line"}}{Classic 3-panel GSEA running score plot
#'       (enrichment score curve + hit marks + ranked metric).
#'       Requires \code{group_use} to specify one group.}
#'     \item{\code{"comparison"}}{NES bubble plot comparing pathways
#'       across all groups.}
#'     \item{\code{"bar"}}{Bidirectional NES bar plot per group.}
#'     \item{\code{"lollipop"}}{Lollipop chart showing NES with
#'       \code{-log10(p.adjust)} colour gradient.}
#'     \item{\code{"ridge"}}{Ridge plot of gene set enrichment scores
#'       via \code{enrichplot::ridgeplot()}.}
#'     \item{\code{"network"}}{Gene-pathway network graph (requires
#'       \pkg{igraph} and \pkg{ggrepel}).}
#'     \item{\code{"enrichmap"}}{Enrichment map: pathway nodes connected
#'       by shared genes, clustered (requires \pkg{igraph} and
#'       \pkg{ggforce}).}
#'     \item{\code{"wordcloud"}}{Word cloud of pathway terms or core
#'       enrichment genes (requires \pkg{ggwordcloud}).}
#'     \item{\code{"volcano_nes"}}{Volcano-style GSEA plot with
#'       \code{-log10(p.adjust)} on x-axis and NES on y-axis.
#'       Points coloured by significance (Activated / Repressed / ns).}
#'     \item{\code{"circle"}}{Circular GSEA plot showing running score
#'       curves, hit marks and ranked-list heatmap in circos tracks
#'       (requires \pkg{circlize}).}
#'     \item{\code{"sankey"}}{Sankey flow diagram linking core enrichment
#'       genes to their pathways (requires \pkg{ggsankey}).}
#'   }
#' @param group_use Character vector of group names to include.
#'   For \code{type = "line"}, specify \strong{one} group.
#'   For other types, defaults to all available groups.
#' @param id_use Character vector of pathway IDs to display. If \code{NULL}
#'   (default), pathways are selected automatically by \code{topTerm} and
#'   \code{direction}.
#' @param topTerm Number of top pathways to show per direction. Default 6.
#' @param direction Which enrichment direction to include:
#'   \code{"both"} (default), \code{"pos"}, or \code{"neg"}.
#' @param padjustCutoff Significance threshold for filtering pathways.
#'   Default 0.05.
#' @param line_width Line width for the running score curve. Default 1.5.
#' @param line_color Color(s) for running score lines. Default
#'   \code{"#6BB82D"}.
#' @param n_coregene Number of core enrichment genes to label in line plot.
#'   Default 0 (no labels).
#' @param features_label Character vector of specific genes to label in line
#'   plot. Overrides \code{n_coregene}.
#' @param label.size Size of gene labels. Default 3.5.
#' @param character_width Maximum character width before wrapping pathway
#'   names. Default 50.
#' @param word_type For \code{type = "wordcloud"}: source of words.
#'   \code{"term"} splits pathway names; \code{"feature"} uses core
#'   enrichment genes. Default \code{"term"}.
#' @param word_size For \code{type = "wordcloud"}: range of text sizes.
#'   Default \code{c(2, 8)}.
#' @param topWord For \code{type = "wordcloud"}: maximum number of words
#'   to display. Default 100.
#' @param network_layout For \code{type = "network"}: igraph layout
#'   algorithm name (e.g., \code{"fr"}, \code{"kk"}, \code{"circle"}).
#'   Default \code{"fr"}.
#' @param enrichmap_layout For \code{type = "enrichmap"}: igraph layout
#'   algorithm. Default \code{"fr"}.
#' @param enrichmap_cluster For \code{type = "enrichmap"}: igraph community
#'   detection algorithm (e.g., \code{"fast_greedy"}, \code{"louvain"}).
#'   Default \code{"fast_greedy"}.
#' @param enrichmap_nlabel For \code{type = "enrichmap"}: number of top
#'   terms per cluster to show as labels. Default 4.
#' @param nes_cutoff For \code{type = "volcano_nes"}: NES threshold for
#'   significance classification. Default 1.
#' @param circ_type For \code{type = "circle"}: inner track style.
#'   \code{"c"} = hit marks + heatmap, \code{"h"} = coloured segments +
#'   heatmap, \code{"m"} = tick marks only. Default \code{"c"}.
#' @param topGenes For \code{type = "sankey"}: number of top core genes
#'   to show per pathway. Default 5.
#' @param palette Color palette name for fills (passed to
#'   \code{palette_colors}). Default \code{"Spectral"}.
#' @param palcolor Custom color vector. Overrides \code{palette}.
#' @param combine Logical; whether to combine multiple panels into one plot
#'   using \pkg{patchwork}. Default \code{TRUE}.
#' @param ncol Number of columns when combining multiple plots. Default
#'   \code{NULL} (auto).
#' @param nrow Number of rows when combining multiple plots. Default
#'   \code{NULL} (auto).
#' @param seed Random seed. Default 11.
#'
#' @return A \code{ggplot} object or a list of \code{ggplot} objects.
#'
#' @details
#' The function dispatches to internal helpers depending on \code{type}:
#' \itemize{
#'   \item \code{.plot_gsea_line} — classic 3-panel GSEA plot using
#'     \pkg{patchwork}.
#'   \item \code{.plot_gsea_comparison} — bubble plot with NES colour, setSize
#'     size, and significance border.
#'   \item \code{.plot_gsea_bar} — bidirectional NES bar plot.
#'   \item \code{.plot_gsea_ridge} — ridge plot via
#'     \code{enrichplot::ridgeplot()}.
#' }
#'
#' @examples
#' \dontrun{
#' markers <- RunDE(seu, group.by = "celltype")
#' res <- RunGsea(markers)
#'
#' # Classic GSEA running score plot
#' PlotGsea(res, type = "line", group_use = "T_cell")
#'
#' # Comparison bubble across all groups
#' PlotGsea(res, type = "comparison", topTerm = 5)
#'
#' # NES bar plot for one group
#' PlotGsea(res, type = "bar", group_use = "T_cell", topTerm = 10)
#'
#' # Ridge plot
#' PlotGsea(res, type = "ridge", group_use = "T_cell")
#'
#' # Lollipop chart
#' PlotGsea(res, type = "lollipop", group_use = "T_cell", topTerm = 10)
#'
#' # Gene-pathway network
#' PlotGsea(res, type = "network", group_use = "T_cell", topTerm = 5)
#'
#' # Enrichment map
#' PlotGsea(res, type = "enrichmap", group_use = "T_cell", topTerm = 30)
#'
#' # Word cloud of pathway terms
#' PlotGsea(res, type = "wordcloud", group_use = "T_cell",
#'          word_type = "term")
#'
#' # Word cloud of core genes
#' PlotGsea(res, type = "wordcloud", group_use = "T_cell",
#'          word_type = "feature")
#'
#' # GSEA volcano (NES vs -log10 padj)
#' PlotGsea(res, type = "volcano_nes", group_use = "T_cell")
#'
#' # Circular GSEA plot
#' PlotGsea(res, type = "circle", group_use = "T_cell",
#'          id_use = c("HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1"))
#'
#' # Sankey gene-pathway flow
#' PlotGsea(res, type = "sankey", group_use = "T_cell", topGenes = 5)
#' }
#'
#' @export
PlotGsea <- function(gsea_result,
                     type = c("line", "comparison", "bar", "lollipop",
                              "ridge", "network", "enrichmap", "wordcloud",
                              "volcano_nes", "circle", "sankey"),
                     group_use = NULL,
                     id_use = NULL,
                     topTerm = 6,
                     direction = c("both", "pos", "neg"),
                     padjustCutoff = 0.05,
                     line_width = 1.5,
                     line_color = "#6BB82D",
                     n_coregene = 0,
                     features_label = NULL,
                     label.size = 3.5,
                     character_width = 50,
                     word_type = c("term", "feature"),
                     word_size = c(2, 8),
                     topWord = 100,
                     network_layout = "fr",
                     enrichmap_layout = "fr",
                     enrichmap_cluster = "fast_greedy",
                     enrichmap_nlabel = 4,
                     nes_cutoff = 1,
                     circ_type = c("c", "h", "m"),
                     topGenes = 5,
                     palette = "Spectral",
                     palcolor = NULL,
                     combine = TRUE,
                     ncol = NULL,
                     nrow = NULL,
                     seed = 11) {

  set.seed(seed)
  type <- match.arg(type)
  direction <- match.arg(direction)
  word_type <- match.arg(word_type)
  circ_type <- match.arg(circ_type)

  # --- validate input ---
  if (!is.list(gsea_result) ||
      !all(c("enrichment", "results") %in% names(gsea_result))) {
    stop("'gsea_result' must be a list from RunGsea() ",
         "containing 'enrichment' and 'results'.", call. = FALSE)
  }

  enrichment <- gsea_result[["enrichment"]]
  results    <- gsea_result[["results"]]

  if (nrow(enrichment) == 0 || length(results) == 0) {
    stop("No GSEA results to plot.", call. = FALSE)
  }

  # Default group_use
  if (is.null(group_use)) {
    group_use <- names(results)
  }

  # Validate group_use
  missing_groups <- setdiff(group_use, names(results))
  if (length(missing_groups) > 0 &&
      type %in% c("line", "bar", "lollipop", "ridge",
                   "network", "enrichmap", "wordcloud",
                   "volcano_nes", "circle", "sankey")) {
    warning("Groups not found in results: ",
            paste(missing_groups, collapse = ", "))
    group_use <- intersect(group_use, names(results))
  }

  switch(type,
    line       = .plot_gsea_line(results, enrichment, group_use, id_use,
                                topTerm, direction, padjustCutoff,
                                line_width, line_color, n_coregene,
                                features_label, label.size,
                                character_width, palette, palcolor,
                                combine, ncol, nrow),
    comparison = .plot_gsea_comparison(enrichment, group_use, id_use,
                                      topTerm, direction, padjustCutoff,
                                      character_width, palette, palcolor),
    bar        = .plot_gsea_bar(results, enrichment, group_use, id_use,
                                topTerm, direction, padjustCutoff,
                                character_width, palette, palcolor,
                                combine, ncol, nrow),
    lollipop   = .plot_gsea_lollipop(results, group_use, id_use,
                                      topTerm, direction, padjustCutoff,
                                      character_width, palette, palcolor,
                                      combine, ncol, nrow),
    ridge      = .plot_gsea_ridge(results, group_use, topTerm,
                                  padjustCutoff, palette, palcolor,
                                  combine, ncol, nrow),
    network    = .plot_gsea_network(results, group_use, id_use,
                                    topTerm, direction, padjustCutoff,
                                    character_width, network_layout,
                                    palette, palcolor,
                                    combine, ncol, nrow),
    enrichmap  = .plot_gsea_enrichmap(results, group_use, id_use,
                                      topTerm, direction, padjustCutoff,
                                      character_width, enrichmap_layout,
                                      enrichmap_cluster, enrichmap_nlabel,
                                      palette, palcolor,
                                      combine, ncol, nrow),
    wordcloud  = .plot_gsea_wordcloud(results, group_use, id_use,
                                      topTerm, direction, padjustCutoff,
                                      word_type, word_size, topWord,
                                      palette, palcolor,
                                      combine, ncol, nrow),
    volcano_nes = .plot_gsea_volcano_nes(results, group_use, id_use,
                                          topTerm, direction, padjustCutoff,
                                          nes_cutoff, character_width,
                                          palette, palcolor,
                                          combine, ncol, nrow),
    circle     = .plot_gsea_circle(results, group_use, id_use,
                                    topTerm, direction, padjustCutoff,
                                    circ_type, line_color,
                                    features_label, character_width,
                                    combine, ncol, nrow),
    sankey     = .plot_gsea_sankey(results, group_use, id_use,
                                   topTerm, direction, padjustCutoff,
                                   topGenes, character_width,
                                   palette, palcolor,
                                   combine, ncol, nrow)
  )
}


# ---------- internal: gsea running score extraction -------------------------

#' Extract GSEA running score data (internal)
#' @keywords internal
.gsea_scores <- function(geneList, geneSet, exponent = 1) {
  geneSet <- intersect(geneSet, names(geneList))
  N  <- length(geneList)
  Nh <- length(geneSet)
  Phit  <- Pmiss <- numeric(N)
  hits  <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  if (NR == 0) NR <- 1

  Phit  <- cumsum(Phit / NR)
  Pmiss[!hits] <- 1 / max(N - Nh, 1)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  data.frame(
    x = seq_along(runningES),
    runningScore = runningES,
    position = as.integer(hits),
    gene = names(geneList),
    stringsAsFactors = FALSE
  )
}

#' Extract gsInfo from a gseaResult object (internal)
#' @keywords internal
.gsInfo <- function(object, id_use) {
  geneList <- object@geneList
  if (is.numeric(id_use)) {
    id_use <- object@result[id_use, "ID"]
  }
  geneSet  <- object@geneSets[[id_use]]
  exponent <- object@params[["exponent"]]
  df <- .gsea_scores(geneList, geneSet, exponent)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList   <- geneList
  df$Description <- object@result[id_use, "Description"]
  if ("core_enrichment" %in% colnames(object@result)) {
    df$CoreGene <- object@result[id_use, "core_enrichment"]
  }
  if (length(object@gene2Symbol) == length(object@geneList)) {
    df$GeneName <- object@gene2Symbol
  } else {
    df$GeneName <- df$gene
  }
  df
}


# ---------- internal: select top pathway IDs --------------------------------

#' Select top pathway IDs by direction and significance (internal)
#' @keywords internal
.select_top_pathways <- function(res_df, topTerm, direction, padjustCutoff,
                                 id_use = NULL) {
  if (!is.null(id_use)) return(id_use)

  sig <- res_df[res_df[["p.adjust"]] < padjustCutoff, , drop = FALSE]
  sig <- sig[order(sig[["p.adjust"]]), , drop = FALSE]

  ids_up   <- sig[sig[["NES"]] > 0, "ID", drop = TRUE]
  ids_down <- sig[sig[["NES"]] < 0, "ID", drop = TRUE]

  switch(direction,
    "pos"  = utils::head(ids_up, topTerm),
    "neg"  = utils::head(ids_down, topTerm),
    "both" = unique(c(
      utils::head(ids_up,   ceiling(topTerm / 2)),
      utils::head(ids_down, ceiling(topTerm / 2))
    ))
  )
}


# ---------- internal: line plot (classic 3-panel GSEA) ----------------------

#' @keywords internal
.plot_gsea_line <- function(results, enrichment, group_use, id_use,
                            topTerm, direction, padjustCutoff,
                            line_width, line_color, n_coregene,
                            features_label, label.size,
                            character_width, palette, palcolor,
                            combine, ncol, nrow) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for line plots.",
         " Install with: install.packages('patchwork')")

  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    # Select pathways
    geneSetID_use <- .select_top_pathways(
      res_obj@result, topTerm, direction, padjustCutoff, id_use
    )
    if (length(geneSetID_use) == 0) next

    # Filter to IDs actually present
    geneSetID_use <- intersect(geneSetID_use, res_obj@result[["ID"]])
    if (length(geneSetID_use) == 0) next

    # Extract running score data
    gsdata_list <- lapply(geneSetID_use, function(id) .gsInfo(res_obj, id))
    gsdata <- do.call(rbind, gsdata_list)

    # Stats table
    stat <- res_obj@result[res_obj@result[["ID"]] %in% geneSetID_use, , drop = FALSE]
    stat <- stat[match(geneSetID_use, stat[["ID"]]), , drop = FALSE]
    stat[["p.sig"]] <- dplyr::case_when(
      stat[["p.adjust"]] > 0.05  ~ "ns",
      stat[["p.adjust"]] > 0.01  ~ "*",
      stat[["p.adjust"]] > 0.001 ~ "**",
      stat[["p.adjust"]] > 1e-4  ~ "***",
      TRUE                       ~ "****"
    )

    # Create description labels with stats
    stat[["DescLabel"]] <- paste0(
      stat[["Description"]],
      "\n(NES=", round(stat[["NES"]], 3),
      ", padj=", format(stat[["p.adjust"]], digits = 3, scientific = TRUE),
      ", ", stat[["p.sig"]], ")"
    )
    stat[["DescLabel"]] <- stringr::str_wrap(stat[["DescLabel"]],
                                              width = character_width)
    desc_map <- stats::setNames(stat[["DescLabel"]], stat[["Description"]])

    gsdata[["DescLabel"]] <- desc_map[gsdata[["Description"]]]
    gsdata[["DescLabel"]] <- factor(gsdata[["DescLabel"]],
                                     levels = unique(gsdata[["DescLabel"]]))

    # --- Panel 1: Enrichment Score Curve ---
    bg_dat <- data.frame(
      xmin = -Inf, xmax = Inf,
      ymin = c(0, -Inf), ymax = c(Inf, 0),
      fill = c(grDevices::adjustcolor("#C40003", 0.15),
               grDevices::adjustcolor("#1D008F", 0.15))
    )

    # Color mapping
    n_paths <- length(geneSetID_use)
    if (length(line_color) == n_paths) {
      color_use <- stats::setNames(line_color, levels(gsdata[["DescLabel"]]))
    } else {
      color_use <- palette_colors(levels(gsdata[["DescLabel"]]),
                                   palette = palette, palcolor = palcolor)
    }

    p1 <- ggplot2::ggplot(gsdata, ggplot2::aes(x = x)) +
      ggplot2::geom_rect(
        data = bg_dat,
        ggplot2::aes(xmin = xmin, xmax = xmax,
                     ymin = ymin, ymax = ymax, fill = I(fill)),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = 1, color = "grey40") +
      ggplot2::geom_line(
        ggplot2::aes(y = runningScore, color = DescLabel),
        linewidth = line_width
      ) +
      ggplot2::scale_color_manual(values = color_use) +
      ggplot2::scale_x_continuous(expand = c(0.01, 0)) +
      ggplot2::ylab("Enrichment Score") +
      ggplot2::xlab(NULL) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          color = "black", fill = "transparent", linewidth = 1
        ),
        plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, l = 0.2,
                                       unit = "cm"),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        panel.grid.major = ggplot2::element_line(colour = "grey90",
                                                  linetype = 2)
      )

    # Single pathway: add NES annotation + peak marker
    if (n_paths == 1) {
      peak_idx <- which.max(abs(gsdata$runningScore))
      nes_val  <- stat[["NES"]][1]
      p1 <- p1 +
        ggplot2::annotate(
          "segment",
          x = 0, xend = gsdata$x[peak_idx],
          y = gsdata$runningScore[peak_idx],
          yend = gsdata$runningScore[peak_idx],
          linetype = 2
        ) +
        ggplot2::annotate(
          "segment",
          x = gsdata$x[peak_idx], xend = gsdata$x[peak_idx],
          y = 0, yend = gsdata$runningScore[peak_idx],
          linetype = 2
        ) +
        ggplot2::annotate(
          "point",
          x = gsdata$x[peak_idx],
          y = gsdata$runningScore[peak_idx],
          fill = ifelse(nes_val < 0, "#5E34F5", "#F52323"),
          color = "black", size = 2.5,
          shape = ifelse(nes_val < 0, 25, 24)
        ) +
        ggplot2::ggtitle(
          stat[["Description"]][1],
          subtitle = paste0(
            "(NES=", round(nes_val, 3),
            ", padj=", format(stat[["p.adjust"]][1],
                              digits = 3, scientific = TRUE),
            ", ", stat[["p.sig"]][1], ")"
          )
        ) +
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(face = "italic"),
          legend.position = "none"
        )

      # Gene labels for single pathway
      if ((is.numeric(n_coregene) && n_coregene > 0) ||
          length(features_label) > 0) {
        if (length(features_label) == 0 && !is.null(gsdata$CoreGene)) {
          core_genes <- unlist(strsplit(gsdata$CoreGene[1], "/"))
          n_show <- min(n_coregene, length(core_genes))
          features_label_use <- gsdata$GeneName[
            gsdata$GeneName %in% core_genes
          ][seq_len(n_show)]
        } else {
          features_label_use <- features_label
        }
        df_gene <- gsdata[gsdata$position == 1 &
                           gsdata$GeneName %in% features_label_use, ,
                          drop = FALSE]
        if (nrow(df_gene) > 0) {
          if (!requireNamespace("ggrepel", quietly = TRUE))
            stop("Package 'ggrepel' is required for gene labels.")
          p1 <- p1 +
            ggplot2::geom_point(
              data = df_gene,
              ggplot2::aes(y = runningScore),
              color = "black", inherit.aes = TRUE
            ) +
            ggrepel::geom_text_repel(
              data = df_gene,
              ggplot2::aes(y = runningScore, label = GeneName),
              min.segment.length = 0, max.overlaps = 50,
              segment.colour = "grey40", size = label.size,
              color = "black", bg.color = "white", bg.r = 0.1
            )
        }
      }
    }

    # --- Panel 2: Hit marks ---
    # Remap ymin/ymax for stacked hit marks per pathway
    i <- 0
    for (term in rev(levels(gsdata$DescLabel))) {
      idx <- which(gsdata$ymin != 0 & gsdata$DescLabel == term)
      gsdata[idx, "ymin"] <- i
      gsdata[idx, "ymax"] <- i + 1
      i <- i + 1
    }

    p2 <- ggplot2::ggplot(gsdata, ggplot2::aes(x = x)) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = ymin, ymax = ymax, color = DescLabel)
      ) +
      ggplot2::scale_color_manual(values = color_use) +
      ggplot2::scale_x_continuous(expand = c(0.01, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(t = -0.1, b = 0, r = 0.2, l = 0.2,
                                       unit = "cm"),
        panel.border = ggplot2::element_rect(
          color = "black", fill = "transparent", linewidth = 1
        ),
        axis.line.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
      )

    # --- Panel 3: Ranked list metric ---
    df_rank <- data.frame(
      x = seq_along(res_obj@geneList),
      y = res_obj@geneList
    )
    # Clip extremes
    y_raw <- df_rank$y
    df_rank$y[df_rank$y > stats::quantile(y_raw, 0.98, na.rm = TRUE)] <-
      stats::quantile(y_raw, 0.98, na.rm = TRUE)
    df_rank$y[df_rank$y < stats::quantile(y_raw, 0.02, na.rm = TRUE)] <-
      stats::quantile(y_raw, 0.02, na.rm = TRUE)

    # Zero crossing
    min_y    <- df_rank$y[which.min(abs(df_rank$y))]
    cross_x  <- stats::median(df_rank$x[df_rank$y == min_y])

    p3 <- ggplot2::ggplot(df_rank, ggplot2::aes(x = x)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = x, xend = x, y = y, yend = 0),
        color = "grey30"
      ) +
      ggplot2::scale_x_continuous(expand = c(0.01, 0)) +
      ggplot2::ylab("Ranked List Metric") +
      ggplot2::xlab("Rank in Ordered Dataset") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2,
                                       unit = "cm"),
        axis.line = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          color = "black", fill = "transparent", linewidth = 1
        )
      )

    if (max(df_rank$y) > 0) {
      p3 <- p3 + ggplot2::annotate(
        "text", x = 0, y = Inf, vjust = 1.3, hjust = 0,
        color = "#C81A1F", size = 4, label = " Positively correlated"
      )
    }
    if (min(df_rank$y) < 0) {
      p3 <- p3 + ggplot2::annotate(
        "text", x = Inf, y = -Inf, vjust = -0.3, hjust = 1,
        color = "#3C298C", size = 4, label = "Negatively correlated "
      )
    }
    if (max(df_rank$y) > 0 && min(df_rank$y) < 0) {
      p3 <- p3 +
        ggplot2::geom_vline(xintercept = cross_x, linetype = 2) +
        ggplot2::annotate(
          "text", y = 0, x = cross_x,
          vjust = ifelse(diff(abs(range(df_rank$y))) > 0, -0.3, 1.3),
          size = 4, label = paste0("Zero cross at ", cross_x)
        )
    }

    # --- Assemble 3 panels with patchwork ---
    p_combined <- p1 / p2 / p3 +
      patchwork::plot_layout(heights = c(3, 1, 2))

    plist[[grp]] <- p_combined
  }

  if (length(plist) == 0) {
    warning("No significant pathways found for any group.")
    return(invisible(NULL))
  }

  if (combine && length(plist) > 1) {
    patchwork::wrap_plots(plist, ncol = ncol, nrow = nrow)
  } else if (length(plist) == 1) {
    plist[[1]]
  } else {
    plist
  }
}


# ---------- internal: comparison bubble plot --------------------------------

#' @keywords internal
.plot_gsea_comparison <- function(enrichment, group_use, id_use,
                                   topTerm, direction, padjustCutoff,
                                   character_width, palette, palcolor) {

  # Subset to requested groups
  enrichment <- enrichment[enrichment[["Group"]] %in% group_use, , drop = FALSE]

  # Collect pathway IDs to show
  if (!is.null(id_use)) {
    ids <- id_use
  } else {
    ids <- c()
    for (grp in group_use) {
      df_grp <- enrichment[enrichment[["Group"]] == grp, , drop = FALSE]
      ids <- unique(c(ids, .select_top_pathways(
        df_grp, topTerm, direction, padjustCutoff
      )))
    }
  }
  if (length(ids) == 0) {
    warning("No significant pathways for comparison plot.")
    return(invisible(NULL))
  }

  enrichment_sub <- enrichment[enrichment[["ID"]] %in% ids, , drop = FALSE]

  # Wrap descriptions
  enrichment_sub[["Description"]] <- stringr::str_wrap(
    enrichment_sub[["Description"]], width = character_width
  )
  # Order by NES
  term_order <- unique(enrichment_sub[["Description"]][
    order(enrichment_sub[["NES"]])
  ])
  enrichment_sub[["Description"]] <- factor(
    enrichment_sub[["Description"]], levels = term_order
  )
  enrichment_sub[["Significant"]] <- factor(
    enrichment_sub[["p.adjust"]] < padjustCutoff,
    levels = c("TRUE", "FALSE")
  )

  nes_max <- max(abs(enrichment_sub[["NES"]]), na.rm = TRUE)

  p <- ggplot2::ggplot(
    enrichment_sub,
    ggplot2::aes(x = Group, y = Description)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = setSize, fill = NES, color = Significant),
      shape = 21, stroke = 0.8
    ) +
    ggplot2::scale_size_area(name = "Gene Set Size", max_size = 6,
                              n.breaks = 4) +
    ggplot2::scale_fill_gradientn(
      name = "NES", n.breaks = 4,
      limits = c(-nes_max, nes_max),
      colors = palette_colors(palette = palette, palcolor = palcolor),
      guide = ggplot2::guide_colorbar(
        frame.colour = "black", ticks.colour = "black",
        title.hjust = 0, order = 1
      )
    ) +
    ggplot2::scale_color_manual(
      name = paste0("Significant\n(padj<", padjustCutoff, ")"),
      values = c("TRUE" = "black", "FALSE" = "grey85"),
      guide = ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        override.aes = list(fill = "grey30", shape = 21),
        order = 2
      )
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "grey80", linetype = 2),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(lineheight = 0.7),
      legend.position = "right"
    ) +
    ggplot2::labs(x = NULL, y = NULL)

  p
}


# ---------- internal: bar plot (bidirectional NES) --------------------------

#' @keywords internal
.plot_gsea_bar <- function(results, enrichment, group_use, id_use,
                            topTerm, direction, padjustCutoff,
                            character_width, palette, palcolor,
                            combine, ncol, nrow) {
  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    geneSetID_use <- .select_top_pathways(
      res_obj@result, topTerm, direction, padjustCutoff, id_use
    )
    if (length(geneSetID_use) == 0) next

    stat <- res_obj@result[res_obj@result[["ID"]] %in% geneSetID_use, ,
                            drop = FALSE]
    stat <- stat[order(stat[["NES"]]), , drop = FALSE]

    # Wrap description
    stat[["Description"]] <- stringr::str_wrap(
      stat[["Description"]], width = character_width
    )
    stat[["Description"]] <- factor(
      stat[["Description"]], levels = unique(stat[["Description"]])
    )
    stat[["Direction"]] <- factor(
      ifelse(stat[["NES"]] > 0, "Pos", "Neg"),
      levels = c("Pos", "Neg")
    )

    # Colors
    dir_colors <- palette_colors(
      x = c("Neg", "Pos"), palette = palette, palcolor = palcolor
    )

    p <- ggplot2::ggplot(
      stat,
      ggplot2::aes(x = .data[["NES"]], y = .data[["Description"]])
    ) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_col(
        ggplot2::aes(fill = .data[["Direction"]],
                     alpha = -log10(.data[["p.adjust"]])),
        color = "black"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          x = 0,
          y = .data[["Description"]],
          label = .data[["Description"]],
          hjust = ifelse(.data[["NES"]] > 0, 1, 0)
        ),
        nudge_x = ifelse(stat[["NES"]] > 0, -0.05, 0.05),
        lineheight = 0.7, size = 3.5
      ) +
      ggplot2::scale_fill_manual(values = dir_colors,
        guide = if (direction == "both") {
          ggplot2::guide_legend(order = 1)
        } else {
          "none"
        }
      ) +
      ggplot2::coord_cartesian(
        xlim = c(-max(abs(stat[["NES"]])), max(abs(stat[["NES"]])))
      ) +
      ggplot2::labs(x = "NES", y = NULL, title = grp) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "right"
      )

    plist[[grp]] <- p
  }

  if (length(plist) == 0) {
    warning("No significant pathways for bar plot.")
    return(invisible(NULL))
  }

  if (combine && length(plist) > 1) {
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
    patchwork::wrap_plots(plist, ncol = ncol, nrow = nrow)
  } else if (length(plist) == 1) {
    plist[[1]]
  } else {
    plist
  }
}


# ---------- internal: ridge plot --------------------------------------------

#' @keywords internal
.plot_gsea_ridge <- function(results, group_use, topTerm,
                              padjustCutoff, palette, palcolor,
                              combine, ncol, nrow) {
  if (!requireNamespace("enrichplot", quietly = TRUE))
    stop("Package 'enrichplot' is required for ridge plots.\n",
         "  BiocManager::install('enrichplot')")

  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    # Filter to significant
    sig_ids <- res_obj@result[["ID"]][
      res_obj@result[["p.adjust"]] < padjustCutoff
    ]
    if (length(sig_ids) == 0) next

    n_show <- min(topTerm, length(sig_ids))

    p <- tryCatch({
      enrichplot::ridgeplot(res_obj, showCategory = n_show) +
        ggplot2::labs(title = grp) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }, error = function(e) {
      warning("Ridge plot failed for group '", grp, "': ",
              conditionMessage(e))
      NULL
    })
    if (!is.null(p)) plist[[grp]] <- p
  }

  if (length(plist) == 0) {
    warning("No ridge plots generated.")
    return(invisible(NULL))
  }

  if (combine && length(plist) > 1) {
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Package 'patchwork' is required.")
    patchwork::wrap_plots(plist, ncol = ncol, nrow = nrow)
  } else if (length(plist) == 1) {
    plist[[1]]
  } else {
    plist
  }
}


# ---------- internal: small utilities for PlotGsea --------------------------

#' Capitalize first letter of each string (internal)
#' @keywords internal
.capitalize <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}

#' Unnest a list-column in a data.frame (internal)
#' @keywords internal
.unnest_cols <- function(df, col) {
  vals <- df[[col]]
  lens <- vapply(vals, length, integer(1))
  # drop rows with 0-length entries
  keep <- lens > 0
  df <- df[rep(seq_len(nrow(df)), lens[keep] * as.integer(keep)), , drop = FALSE]
  # actually we need to expand
  df2 <- df[rep(which(keep), lens[keep]), , drop = FALSE]
  df2[[col]] <- unlist(vals[keep])
  rownames(df2) <- NULL
  df2
}

#' Standard combine-or-return logic (internal)
#' @keywords internal
.combine_plist <- function(plist, combine, ncol, nrow) {
  if (length(plist) == 0) {
    warning("No plots generated.")
    return(invisible(NULL))
  }
  if (combine && length(plist) > 1) {
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
    patchwork::wrap_plots(plist, ncol = ncol, nrow = nrow)
  } else if (length(plist) == 1) {
    plist[[1]]
  } else {
    plist
  }
}


# ---------- internal: lollipop plot -----------------------------------------

#' @keywords internal
.plot_gsea_lollipop <- function(results, group_use, id_use,
                                 topTerm, direction, padjustCutoff,
                                 character_width, palette, palcolor,
                                 combine, ncol, nrow) {
  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    geneSetID_use <- .select_top_pathways(
      res_obj@result, topTerm, direction, padjustCutoff, id_use
    )
    if (length(geneSetID_use) == 0) next

    stat <- res_obj@result[res_obj@result[["ID"]] %in% geneSetID_use, ,
                            drop = FALSE]
    stat <- stat[order(stat[["NES"]]), , drop = FALSE]

    stat[["Description"]] <- .capitalize(stat[["Description"]])
    stat[["Description"]] <- stringr::str_wrap(
      stat[["Description"]], width = character_width
    )
    stat[["Description"]] <- factor(
      stat[["Description"]], levels = unique(stat[["Description"]])
    )
    stat[["neg_log10_padj"]] <- -log10(stat[["p.adjust"]])
    stat[["Direction"]] <- factor(
      ifelse(stat[["NES"]] > 0, "Pos", "Neg"), levels = c("Pos", "Neg")
    )

    # Count core enrichment genes
    if ("core_enrichment" %in% colnames(stat)) {
      stat[["Count"]] <- vapply(
        strsplit(stat[["core_enrichment"]], "/"),
        length, integer(1)
      )
    } else {
      stat[["Count"]] <- stat[["setSize"]]
    }

    p <- ggplot2::ggplot(
      stat,
      ggplot2::aes(x = .data[["Description"]],
                   y = .data[["NES"]])
    ) +
      ggplot2::geom_hline(yintercept = 0, color = "grey60") +
      # Outer thick segment (border effect)
      ggplot2::geom_segment(
        ggplot2::aes(y = 0,
                     xend = .data[["Description"]],
                     yend = .data[["NES"]]),
        color = "black", linewidth = 2
      ) +
      # Inner coloured segment
      ggplot2::geom_segment(
        ggplot2::aes(y = 0,
                     xend = .data[["Description"]],
                     yend = .data[["NES"]],
                     color = .data[["neg_log10_padj"]]),
        linewidth = 1
      ) +
      # Point at end
      ggplot2::geom_point(
        ggplot2::aes(size = .data[["Count"]],
                     fill = .data[["neg_log10_padj"]]),
        shape = 21, color = "black"
      ) +
      ggplot2::scale_size(name = "Gene Count", range = c(3, 6),
                           breaks = scales::breaks_extended(n = 4)) +
      ggplot2::guides(
        size = ggplot2::guide_legend(
          override.aes = list(fill = "grey30", shape = 21), order = 1
        )
      ) +
      ggplot2::scale_fill_gradientn(
        name = expression(-log[10](p.adjust)),
        n.breaks = 3,
        colours = palette_colors(palette = palette, palcolor = palcolor),
        na.value = "grey80",
        guide = ggplot2::guide_colorbar(
          frame.colour = "black", ticks.colour = "black",
          title.hjust = 0
        ),
        aesthetics = c("color", "fill")
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "NES", title = grp) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(colour = "grey80",
                                                  linetype = 2),
        axis.text.y = ggplot2::element_text(lineheight = 0.7),
        legend.position = "right"
      )

    plist[[grp]] <- p
  }

  .combine_plist(plist, combine, ncol, nrow)
}


# ---------- internal: network plot ------------------------------------------

#' @keywords internal
.plot_gsea_network <- function(results, group_use, id_use,
                                topTerm, direction, padjustCutoff,
                                character_width, network_layout,
                                palette, palcolor,
                                combine, ncol, nrow) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required for network plots.\n",
         "  install.packages('igraph')")
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Package 'ggrepel' is required for network plots.\n",
         "  install.packages('ggrepel')")

  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    geneSetID_use <- .select_top_pathways(
      res_obj@result, topTerm, direction, padjustCutoff, id_use
    )
    if (length(geneSetID_use) == 0) next

    df <- res_obj@result[res_obj@result[["ID"]] %in% geneSetID_use, ,
                          drop = FALSE]
    df[["neg_log10_padj"]] <- -log10(df[["p.adjust"]])
    df[["Description"]] <- .capitalize(df[["Description"]])
    df[["Description"]] <- stringr::str_wrap(
      df[["Description"]], width = character_width
    )
    # Add NES + p info to description label
    df[["p.sig"]] <- dplyr::case_when(
      df[["p.adjust"]] > 0.05  ~ "ns",
      df[["p.adjust"]] > 0.01  ~ "*",
      df[["p.adjust"]] > 0.001 ~ "**",
      df[["p.adjust"]] > 1e-4  ~ "***",
      TRUE                     ~ "****"
    )
    df[["DescLabel"]] <- paste0(
      df[["Description"]],
      "\n(NES=", round(df[["NES"]], 3), ", ", df[["p.sig"]], ")"
    )
    df[["DescLabel"]] <- factor(df[["DescLabel"]],
                                 levels = unique(df[["DescLabel"]]))

    # Build gene-pathway edges from core_enrichment
    if (!"core_enrichment" %in% colnames(df)) next
    df[["geneID"]] <- strsplit(df[["core_enrichment"]], "/")
    df_unnest <- .unnest_cols(df, "geneID")

    # Nodes: terms + genes
    nodes <- rbind(
      data.frame(ID = df[["DescLabel"]], class = "term",
                 metric = df[["neg_log10_padj"]],
                 stringsAsFactors = FALSE),
      data.frame(ID = unique(df_unnest[["geneID"]]), class = "gene",
                 metric = 0, stringsAsFactors = FALSE)
    )
    # Edges
    edges <- data.frame(
      from = df_unnest[["DescLabel"]],
      to   = df_unnest[["geneID"]],
      stringsAsFactors = FALSE
    )

    graph <- igraph::graph_from_data_frame(d = edges, vertices = nodes,
                                            directed = FALSE)
    # Layout
    if (network_layout %in% c("circle", "tree", "grid")) {
      layout_mat <- switch(network_layout,
        "circle" = igraph::layout_in_circle(graph),
        "tree"   = igraph::layout_as_tree(graph),
        "grid"   = igraph::layout_on_grid(graph)
      )
    } else {
      layout_fun <- utils::getFromNamespace(
        paste0("layout_with_", network_layout), "igraph"
      )
      layout_mat <- layout_fun(graph)
    }

    df_graph <- igraph::as_data_frame(graph, what = "both")
    df_nodes <- df_graph$vertices
    df_nodes[["dim1"]] <- layout_mat[, 1]
    df_nodes[["dim2"]] <- layout_mat[, 2]

    df_edges <- df_graph$edges
    df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
    df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
    df_edges[["to_dim1"]]   <- df_nodes[df_edges[["to"]], "dim1"]
    df_edges[["to_dim2"]]   <- df_nodes[df_edges[["to"]], "dim2"]

    # Colors per term
    term_colors <- palette_colors(levels(df[["DescLabel"]]),
                                   palette = palette, palcolor = palcolor)

    # Gene node color: blend from connected terms
    gene_nodes <- df_nodes[df_nodes$class == "gene", "name", drop = TRUE]
    gene_colors <- vapply(gene_nodes, function(g) {
      connected_terms <- df_unnest[["DescLabel"]][df_unnest[["geneID"]] == g]
      connected_terms <- unique(as.character(connected_terms))
      cols <- term_colors[connected_terms]
      cols <- cols[!is.na(cols)]
      if (length(cols) == 0) return("grey70")
      # Simple average blending
      rgb_mat <- grDevices::col2rgb(cols)
      avg_rgb <- rowMeans(rgb_mat)
      grDevices::rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue = 255)
    }, character(1))

    all_colors <- c(term_colors, stats::setNames(gene_colors, gene_nodes))
    df_nodes[["color"]] <- all_colors[df_nodes[["name"]]]

    # Edge colors from term
    df_edges[["color"]] <- all_colors[df_edges[["from"]]]

    # Label contrast
    rgb_sum <- colSums(grDevices::col2rgb(df_nodes[["color"]]))
    df_nodes[["label_color"]] <- ifelse(rgb_sum > 255 * 2, "black", "white")

    # Term number labels
    df_nodes[["label"]] <- NA
    term_names <- levels(df[["DescLabel"]])
    for (j in seq_along(term_names)) {
      df_nodes[df_nodes[["name"]] == term_names[j], "label"] <- j
    }

    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = df_edges,
        ggplot2::aes(x = from_dim1, y = from_dim2,
                     xend = to_dim1, yend = to_dim2,
                     color = color),
        alpha = 0.6, lineend = "round", show.legend = FALSE
      ) +
      ggplot2::geom_label(
        data = df_nodes[df_nodes$class == "gene", ],
        ggplot2::aes(x = dim1, y = dim2, label = name,
                     fill = color, color = label_color),
        size = 2.5, show.legend = FALSE
      ) +
      # Term nodes: outer ring
      ggplot2::geom_point(
        data = df_nodes[df_nodes$class == "term", ],
        ggplot2::aes(x = dim1, y = dim2),
        size = 8, color = "black", fill = "black",
        stroke = 1, shape = 21, show.legend = FALSE
      ) +
      # Term nodes: inner fill
      ggplot2::geom_point(
        data = df_nodes[df_nodes$class == "term", ],
        ggplot2::aes(x = dim1, y = dim2, fill = color),
        size = 7, color = "white", stroke = 1, shape = 21
      ) +
      ggrepel::geom_text_repel(
        data = df_nodes[df_nodes$class == "term", ],
        ggplot2::aes(x = dim1, y = dim2, label = label),
        fontface = "bold", min.segment.length = 0,
        segment.color = "black", point.size = NA,
        max.overlaps = 100, force = 0,
        color = "white", bg.color = "black", bg.r = 0.1, size = 5
      ) +
      ggplot2::scale_color_identity(guide = "none") +
      ggplot2::scale_fill_identity(
        name = "Term:",
        guide = "legend",
        labels = term_names,
        breaks = term_colors[term_names]
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(title = "Term:", byrow = TRUE)
      ) +
      ggplot2::labs(x = NULL, y = NULL, title = grp) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "right"
      )

    plist[[grp]] <- p
  }

  .combine_plist(plist, combine, ncol, nrow)
}


# ---------- internal: enrichment map ----------------------------------------

#' @keywords internal
.plot_gsea_enrichmap <- function(results, group_use, id_use,
                                  topTerm, direction, padjustCutoff,
                                  character_width, enrichmap_layout,
                                  enrichmap_cluster, enrichmap_nlabel,
                                  palette, palcolor,
                                  combine, ncol, nrow) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required for enrichment maps.\n",
         "  install.packages('igraph')")

  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    # Use larger default topTerm for enrichmap
    local_topTerm <- max(topTerm, 20)
    geneSetID_use <- .select_top_pathways(
      res_obj@result, local_topTerm, direction, padjustCutoff, id_use
    )
    if (length(geneSetID_use) < 2) next  # need >=2 for edges

    df <- res_obj@result[res_obj@result[["ID"]] %in% geneSetID_use, ,
                          drop = FALSE]
    df[["neg_log10_padj"]] <- -log10(df[["p.adjust"]])
    df[["Description"]] <- .capitalize(df[["Description"]])
    df[["Description"]] <- stringr::str_wrap(
      df[["Description"]], width = character_width
    )
    df[["Direction"]] <- factor(
      ifelse(df[["NES"]] > 0, "Pos", "Neg"), levels = c("Pos", "Neg")
    )

    # Parse core genes
    if (!"core_enrichment" %in% colnames(df)) next
    df[["geneID"]] <- strsplit(df[["core_enrichment"]], "/")
    df[["Count"]] <- vapply(df[["geneID"]], length, integer(1))
    rownames(df) <- df[["ID"]]

    # Build term-term edges by gene overlap
    nodes <- df
    combs <- utils::combn(nodes[["ID"]], 2)
    edge_weights <- apply(combs, 2, function(pair) {
      length(intersect(df[pair[1], "geneID"][[1]],
                       df[pair[2], "geneID"][[1]]))
    })
    edges <- data.frame(
      from   = combs[1, ],
      to     = combs[2, ],
      weight = edge_weights,
      stringsAsFactors = FALSE
    )
    edges <- edges[edges[["weight"]] > 0, , drop = FALSE]

    if (nrow(edges) == 0) {
      # No shared genes; fall back to isolated nodes
      edges <- data.frame(from = character(0), to = character(0),
                          weight = numeric(0))
    }

    graph <- igraph::graph_from_data_frame(d = edges, vertices = nodes,
                                            directed = FALSE)

    # Layout
    if (enrichmap_layout %in% c("circle", "tree", "grid")) {
      layout_mat <- switch(enrichmap_layout,
        "circle" = igraph::layout_in_circle(graph),
        "tree"   = igraph::layout_as_tree(graph),
        "grid"   = igraph::layout_on_grid(graph)
      )
    } else {
      layout_fun <- utils::getFromNamespace(
        paste0("layout_with_", enrichmap_layout), "igraph"
      )
      layout_mat <- layout_fun(graph)
    }

    # Community detection
    cluster_fun <- utils::getFromNamespace(
      paste0("cluster_", enrichmap_cluster), "igraph"
    )
    clusters <- cluster_fun(graph)

    df_graph <- igraph::as_data_frame(graph, what = "both")
    df_nodes <- df_graph$vertices
    df_nodes[["dim1"]] <- layout_mat[, 1]
    df_nodes[["dim2"]] <- layout_mat[, 2]
    df_nodes[["cluster"]] <- factor(
      paste0("C", clusters$membership),
      levels = paste0("C", sort(unique(clusters$membership)))
    )

    # Cluster labels: top terms per cluster
    cluster_labels <- tapply(seq_len(nrow(df_nodes)), df_nodes[["cluster"]],
      function(idx) {
        sub_df <- df_nodes[idx, , drop = FALSE]
        sub_df <- sub_df[order(sub_df[["neg_log10_padj"]], decreasing = TRUE), ]
        top_terms <- utils::head(sub_df[["Description"]], enrichmap_nlabel)
        paste(top_terms, collapse = "\n")
      }
    )
    df_nodes[["cluster_label"]] <- as.character(
      cluster_labels[as.character(df_nodes[["cluster"]])]
    )

    # Edges coordinates
    df_edges <- df_graph$edges
    if (nrow(df_edges) > 0) {
      df_edges[["from_dim1"]] <- df_nodes[df_edges[["from"]], "dim1"]
      df_edges[["from_dim2"]] <- df_nodes[df_edges[["from"]], "dim2"]
      df_edges[["to_dim1"]]   <- df_nodes[df_edges[["to"]], "dim1"]
      df_edges[["to_dim2"]]   <- df_nodes[df_edges[["to"]], "dim2"]
    }

    # Cluster colors
    cluster_colors <- palette_colors(levels(df_nodes[["cluster"]]),
                                      palette = palette, palcolor = palcolor)

    # Build plot
    p <- ggplot2::ggplot()

    # Try to add mark layer (ggforce ellipse)
    has_ggforce <- requireNamespace("ggforce", quietly = TRUE)
    if (has_ggforce && nlevels(df_nodes[["cluster"]]) > 1) {
      p <- p +
        ggforce::geom_mark_ellipse(
          data = df_nodes,
          ggplot2::aes(x = dim1, y = dim2,
                       color = cluster, fill = cluster,
                       label = cluster,
                       description = cluster_label),
          expand = grid::unit(3, "mm"),
          alpha = 0.1,
          label.margin = ggplot2::margin(1, 1, 1, 1, "mm"),
          label.fontsize = 8,
          label.fill = "grey95",
          con.size = 1, con.cap = 0
        )
    }

    # Edges
    if (nrow(df_edges) > 0) {
      p <- p +
        ggplot2::geom_segment(
          data = df_edges,
          ggplot2::aes(x = from_dim1, y = from_dim2,
                       xend = to_dim1, yend = to_dim2,
                       linewidth = weight),
          alpha = 0.15, lineend = "round"
        ) +
        ggplot2::scale_linewidth(
          name = "Shared Genes", range = c(0.3, 3),
          breaks = scales::breaks_extended(n = 4)
        ) +
        ggplot2::guides(
          linewidth = ggplot2::guide_legend(
            override.aes = list(alpha = 1, color = "grey"), order = 2
          )
        )
    }

    # Nodes
    p <- p +
      ggplot2::geom_point(
        data = df_nodes,
        ggplot2::aes(x = dim1, y = dim2,
                     size = Count, fill = cluster),
        color = "black", shape = 21
      ) +
      ggplot2::scale_size(
        name = "Gene Count", range = c(2, 6),
        breaks = scales::breaks_extended(n = 4)
      ) +
      ggplot2::guides(
        size = ggplot2::guide_legend(
          override.aes = list(fill = "grey30", shape = 21), order = 1
        )
      ) +
      ggplot2::scale_fill_manual(
        name = "Cluster",
        values = cluster_colors,
        aesthetics = c("colour", "fill")
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          override.aes = list(alpha = 1, color = "black", shape = 21,
                              size = 4),
          byrow = TRUE, order = 3
        ),
        color = "none"
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(0.15, 0)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(0.15, 0)) +
      ggplot2::labs(x = NULL, y = NULL, title = grp) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "right"
      )

    plist[[grp]] <- p
  }

  .combine_plist(plist, combine, ncol, nrow)
}


# ---------- internal: word cloud plot ---------------------------------------

#' @keywords internal
.plot_gsea_wordcloud <- function(results, group_use, id_use,
                                  topTerm, direction, padjustCutoff,
                                  word_type, word_size, topWord,
                                  palette, palcolor,
                                  combine, ncol, nrow) {
  if (!requireNamespace("ggwordcloud", quietly = TRUE))
    stop("Package 'ggwordcloud' is required for word cloud plots.\n",
         "  install.packages('ggwordcloud')")

  # Common stop-words to exclude
  words_excluded <- c(
    "of", "in", "to", "and", "the", "a", "an", "by", "for", "with",
    "from", "on", "at", "or", "as", "is", "via", "into", "its", "type",
    "i", "ii", "vs", "process", "regulation", "activity", "pathway",
    "positive", "negative", "response", "signaling", "signal", "cell",
    "cellular"
  )

  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    # Select all significant pathways (not just topTerm) for wordcloud
    geneSetID_use <- .select_top_pathways(
      res_obj@result, topTerm = Inf, direction, padjustCutoff, id_use
    )
    if (length(geneSetID_use) == 0) next

    df <- res_obj@result[res_obj@result[["ID"]] %in% geneSetID_use, ,
                          drop = FALSE]
    df[["neg_log10_padj"]] <- -log10(df[["p.adjust"]])

    if (word_type == "term") {
      # Split pathway descriptions into words
      word_list <- strsplit(
        tolower(as.character(df[["Description"]])), "\\s+"
      )
      df[["keyword"]] <- word_list
      df_words <- .unnest_cols(df, "keyword")

      # Aggregate: score = sum(-log10 padj) per word, count = appearances
      word_agg <- stats::aggregate(
        df_words[["neg_log10_padj"]],
        by = list(keyword = df_words[["keyword"]]),
        FUN = function(x) c(score = sum(x), count = length(x))
      )
      word_df <- data.frame(
        keyword = word_agg[[1]],
        score   = word_agg[[2]][, "score"],
        count   = as.integer(word_agg[[2]][, "count"]),
        stringsAsFactors = FALSE
      )

      # Filter stop words + short words
      word_df <- word_df[
        !tolower(word_df[["keyword"]]) %in% words_excluded &
        nchar(word_df[["keyword"]]) > 1 &
        !grepl("\\[.*\\]", word_df[["keyword"]]), ,
        drop = FALSE
      ]

    } else {
      # word_type == "feature": use core enrichment genes
      if (!"core_enrichment" %in% colnames(df)) next
      df[["keyword"]] <- strsplit(df[["core_enrichment"]], "/")
      df_words <- .unnest_cols(df, "keyword")

      word_agg <- stats::aggregate(
        df_words[["neg_log10_padj"]],
        by = list(keyword = df_words[["keyword"]]),
        FUN = function(x) c(score = sum(x), count = length(x))
      )
      word_df <- data.frame(
        keyword = word_agg[[1]],
        score   = word_agg[[2]][, "score"],
        count   = as.integer(word_agg[[2]][, "count"]),
        stringsAsFactors = FALSE
      )
    }

    if (nrow(word_df) == 0) next

    # Take top N words
    word_df <- word_df[
      utils::head(order(word_df[["score"]], decreasing = TRUE), topWord), ,
      drop = FALSE
    ]

    # Random angle for visual variety
    word_df[["angle"]] <- 90 * sample(
      c(0, 1), nrow(word_df), replace = TRUE, prob = c(0.65, 0.35)
    )

    p <- ggplot2::ggplot(
      word_df,
      ggplot2::aes(
        label = .data[["keyword"]],
        size  = .data[["count"]],
        color = .data[["score"]],
        angle = .data[["angle"]]
      )
    ) +
      ggwordcloud::geom_text_wordcloud(
        rm_outside = TRUE, eccentricity = 1,
        shape = "square", show.legend = TRUE, grid_margin = 3
      ) +
      ggplot2::scale_color_gradientn(
        name = "Score",
        colours = palette_colors(palette = palette, palcolor = palcolor),
        guide = ggplot2::guide_colorbar(
          frame.colour = "black", ticks.colour = "black",
          title.hjust = 0
        )
      ) +
      ggplot2::scale_size(
        name = "Count", range = word_size,
        breaks = ceiling(seq(
          min(word_df[["count"]], na.rm = TRUE),
          max(word_df[["count"]], na.rm = TRUE),
          length.out = 3
        ))
      ) +
      ggplot2::guides(
        size = ggplot2::guide_legend(
          override.aes = list(colour = "black", label = "A"), order = 1
        )
      ) +
      ggplot2::labs(title = grp) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

    plist[[grp]] <- p
  }

  .combine_plist(plist, combine, ncol, nrow)
}


# ---------- internal: GSEA volcano (NES vs -log10 padj) ---------------------

#' @keywords internal
.plot_gsea_volcano_nes <- function(results, group_use, id_use,
                                    topTerm, direction, padjustCutoff,
                                    nes_cutoff, character_width,
                                    palette, palcolor,
                                    combine, ncol, nrow) {
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Package 'ggrepel' is required. Install with: install.packages('ggrepel')")

  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    df <- as.data.frame(res_obj@result)
    if (nrow(df) == 0) next

    # Compute -log10(p.adjust) — cap p=0
    pvals <- df[["p.adjust"]]
    pvals[pvals == 0] <- min(pvals[pvals > 0], na.rm = TRUE) * 0.1
    df[["logp"]] <- -log10(pvals)

    # Classify significance
    df[["type"]] <- dplyr::case_when(
      df[["NES"]] >=  nes_cutoff & df[["p.adjust"]] < padjustCutoff ~ "Activated",
      df[["NES"]] <= -nes_cutoff & df[["p.adjust"]] < padjustCutoff ~ "Repressed",
      TRUE ~ "ns"
    )
    df[["type"]] <- factor(df[["type"]],
                            levels = c("Activated", "ns", "Repressed"))

    # Top terms to label
    topterm <- do.call(rbind, lapply(c("Activated", "Repressed"), function(tp) {
      sub <- df[df[["type"]] == tp, , drop = FALSE]
      sub <- sub[order(sub[["p.adjust"]]), , drop = FALSE]
      utils::head(sub, topTerm)
    }))

    # Use user-specified or default colours
    if (!is.null(palcolor) && length(palcolor) >= 3) {
      pt_colors <- palcolor[1:3]
    } else {
      pt_colors <- c("#CC3333", "#CCCCCC", "#0099CC")
    }
    names(pt_colors) <- c("Activated", "ns", "Repressed")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = logp, y = NES)) +
      ggplot2::geom_point(
        ggplot2::aes(color = type), alpha = 0.6, size = 2.5
      ) +
      ggplot2::geom_vline(
        xintercept = -log10(padjustCutoff),
        linewidth = 0.8, linetype = "solid", color = "grey75"
      ) +
      ggplot2::geom_hline(
        yintercept = c(-nes_cutoff, nes_cutoff),
        linewidth = 0.8, linetype = "dashed", color = "grey75"
      ) +
      ggplot2::scale_color_manual(
        name = "", values = pt_colors, drop = FALSE
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 4))
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(colour = "black"),
        legend.position = "top"
      ) +
      ggplot2::labs(
        x = expression(-log[10](p.adjust)),
        y = "Normalized Enrichment Score",
        title = grp
      )

    # Add labels for top terms
    if (nrow(topterm) > 0) {
      topterm[["label"]] <- stringr::str_wrap(
        topterm[["Description"]], width = character_width
      )
      p <- p +
        ggrepel::geom_text_repel(
          data = topterm,
          ggplot2::aes(x = logp, y = NES, label = label),
          fontface = "italic", max.overlaps = 50,
          force = 40, min.segment.length = ggplot2::unit(0.1, "cm"),
          size = 3.2, color = "black"
        )
    }

    plist[[grp]] <- p
  }

  .combine_plist(plist, combine, ncol, nrow)
}


# ---------- internal: circular GSEA plot (circlize) -------------------------
# Draws directly to the current graphics device using circlize (base graphics),
# matching the behaviour of GseaVis::circGsea.  Returns invisible(NULL).

#' @keywords internal
.plot_gsea_circle <- function(results, group_use, id_use,
                               topTerm, direction, padjustCutoff,
                               circ_type, curveCol,
                               markGene, character_width,
                               combine, ncol, nrow) {
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package 'circlize' is required for circle plots.\n",
         "  install.packages('circlize')")

  # Collect valid drawing tasks per group
 draw_tasks <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    geneSetID_use <- .select_top_pathways(
      res_obj@result, topTerm, direction, padjustCutoff, id_use
    )
    if (length(geneSetID_use) == 0) next
    geneSetID_use <- intersect(geneSetID_use, res_obj@result[["ID"]])
    if (length(geneSetID_use) == 0) next

    draw_tasks[[grp]] <- list(res_obj = res_obj,
                               geneSetID_use = geneSetID_use)
  }

  if (length(draw_tasks) == 0) {
    message("No significant pathways to plot for circle type.")
    return(invisible(NULL))
  }

  # Multi-panel layout when more than one group
  n_panels <- length(draw_tasks)
  if (n_panels > 1) {
    nr <- nrow %||% ceiling(sqrt(n_panels))
    nc <- ncol %||% ceiling(n_panels / nr)
    old_par <- graphics::par(mfrow = c(nr, nc), mar = c(1, 1, 2, 1))
    on.exit(graphics::par(old_par), add = TRUE)
  }

  for (grp in names(draw_tasks)) {
    task <- draw_tasks[[grp]]
    res_obj <- task$res_obj
    geneSetID_use <- task$geneSetID_use

    # --- prepare data ---
    gsdata <- do.call(rbind, lapply(geneSetID_use, function(sid) {
      df <- .gsInfo(res_obj, sid)
      df$id <- sid
      df
    }))

    gsdata1 <- gsdata[gsdata$position == 1, , drop = FALSE]

    htCol <- c("#08519C", "#A50F15")
    ht <- do.call(rbind, lapply(geneSetID_use, function(sid) {
      tmp <- gsdata[gsdata$id == sid, , drop = FALSE]
      v <- seq(1, max(sum(tmp$position), 1), length.out = 9)
      inv <- findInterval(rev(cumsum(tmp$position)), v)
      if (min(inv) == 0) inv <- inv + 1
      color <- grDevices::colorRampPalette(c(htCol[1], "white", htCol[2]))(10)
      xmin <- which(!duplicated(inv))
      xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
      data.frame(xmin = xmin, xmax = xmax, col = color[unique(inv)],
                 id = sid, stringsAsFactors = FALSE)
    }))

    df_lb <- res_obj@result[res_obj@result[["ID"]] %in% geneSetID_use, ,
                             drop = FALSE]

    if (!is.null(markGene) && length(markGene) > 0) {
      mgene_df <- gsdata[gsdata$GeneName %in% markGene &
                           gsdata$position == 1, , drop = FALSE]
    } else {
      mgene_df <- NULL
    }

    # --- draw directly using circlize ---
    circlize::circos.clear()
    circlize::circos.par(circle.margin = rep(0.1, 4))
    circlize::circos.initialize(
      sectors = gsdata1$id,
      xlim = c(1, max(gsdata$x))
    )

    # Track 1: Pathway ID + NES + pvalue label
    circlize::circos.track(
      sectors = geneSetID_use, ylim = c(0, 1),
      bg.col = "grey80",
      panel.fun = function(x, y) {
        tmp <- df_lb[df_lb[["ID"]] == circlize::CELL_META$sector.index, ]
        pval <- tmp[["pvalue"]]
        plabel <- if (pval < 0.001) "< 0.001"
                  else if (pval < 0.01) "< 0.01"
                  else if (pval < 0.05) "< 0.05"
                  else round(pval, 2)
        goName <- paste0(
          circlize::CELL_META$sector.index, "\n(",
          "NES:", round(tmp[["NES"]], 1), " | P:", plabel, ")"
        )
        circlize::circos.text(
          x = circlize::CELL_META$xcenter,
          y = circlize::CELL_META$ycenter,
          labels = goName, font = 2,
          niceFacing = TRUE, facing = "inside", cex = 0.7
        )
      }
    )

    # Track 2: Running score curve + segments
    segCol <- c("#0099CC", "#FF0033")
    circlize::circos.track(
      sectors = geneSetID_use, ylim = c(-1, 1),
      panel.fun = function(x, y) {
        sid <- circlize::CELL_META$sector.index
        tmp <- gsdata1[gsdata1$id == sid, , drop = FALSE]
        circlize::circos.segments(
          x0 = 0, x1 = max(gsdata$x), y0 = 0, y1 = 0,
          lty = "dashed", col = "black"
        )
        if (circ_type == "h") {
          colors <- grDevices::colorRampPalette(segCol)(max(nrow(tmp), 1))
          sorted_colors <- colors[match(
            tmp$runningScore, sort(tmp$runningScore)
          )]
          circlize::circos.lines(
            x = c(1, tmp$x, max(gsdata$x)),
            y = c(0, tmp$runningScore, 0),
            type = "h", baseline = 0, col = sorted_colors
          )
        } else if (circ_type == "m") {
          circlize::circos.segments(
            x0 = c(1, tmp$x, max(gsdata$x)),
            x1 = c(1, tmp$x, max(gsdata$x)),
            y0 = -0.25, y1 = 0.25, lty = "solid", col = "grey30"
          )
        }
        circlize::circos.lines(
          x = c(1, tmp$x, max(gsdata$x)),
          y = c(0, tmp$runningScore, 0),
          lwd = 2, col = curveCol
        )
      }
    )

    # Track 3: Heatmap + hit marks
    circlize::circos.track(
      sectors = geneSetID_use, ylim = c(0, 1),
      panel.fun = function(x, y) {
        sid <- circlize::CELL_META$sector.index
        if (circ_type == "c") {
          tmp <- gsdata1[gsdata1$id == sid, , drop = FALSE]
          circlize::circos.segments(
            x0 = tmp$x, x1 = tmp$x, y0 = 0, y1 = 1
          )
          ytop <- 0.5
        } else {
          ytop <- 1
        }
        tmpht <- ht[ht$id == sid, , drop = FALSE]
        circlize::circos.rect(
          xleft = tmpht$xmin, xright = tmpht$xmax,
          col = ggplot2::alpha(tmpht$col, 0.75),
          border = ggplot2::alpha(tmpht$col, 0.75),
          ybottom = rep(0, nrow(tmpht)),
          ytop = rep(ytop, nrow(tmpht))
        )
        circlize::circos.rect(
          xleft = 0, xright = max(tmpht$xmax),
          col = NA, border = "black", ybottom = 0, ytop = ytop
        )
      }
    )

    # Mark genes if provided
    if (!is.null(mgene_df) && nrow(mgene_df) > 0) {
      circlize::circos.labels(
        sectors = mgene_df$id, x = mgene_df$x,
        labels = mgene_df$GeneName, cex = 0.5,
        connection_height = circlize::mm_h(2.5)
      )
    }

    graphics::title(main = grp, line = -1)
  }

  circlize::circos.clear()
  invisible(NULL)
}


# ---------- internal: sankey gene-pathway flow plot -------------------------

#' @keywords internal
.plot_gsea_sankey <- function(results, group_use, id_use,
                               topTerm, direction, padjustCutoff,
                               topGenes, character_width,
                               palette, palcolor,
                               combine, ncol, nrow) {
  if (!requireNamespace("ggsankey", quietly = TRUE))
    stop("Package 'ggsankey' is required for sankey plots.\n",
         "  Install from GitHub: remotes::install_github('davidsjoberg/ggsankey')")

  plist <- list()

  for (grp in group_use) {
    res_obj <- results[[grp]]
    if (is.null(res_obj)) next

    geneSetID_use <- .select_top_pathways(
      res_obj@result, topTerm, direction, padjustCutoff, id_use
    )
    if (length(geneSetID_use) == 0) next

    df <- res_obj@result[res_obj@result[["ID"]] %in% geneSetID_use, ,
                          drop = FALSE]
    if (!"core_enrichment" %in% colnames(df) || nrow(df) == 0) next

    # Sort by pvalue (best first)
    df <- df[order(df[["pvalue"]]), , drop = FALSE]

    # Wrap description
    df[["Description"]] <- stringr::str_wrap(
      df[["Description"]], width = character_width
    )
    # Pathway order: reversed so first pathway at top in sankey
    term_levels <- rev(unique(df[["Description"]]))
    df[["Description"]] <- factor(df[["Description"]], levels = term_levels)

    # ---- Build gene-pathway pairs (ordered by pathway) ----
    # Iterate pathways in display order so genes are grouped by pathway
    gene_pairs <- do.call(rbind, lapply(term_levels, function(desc) {
      row_idx <- which(as.character(df[["Description"]]) == desc)
      if (length(row_idx) == 0) return(NULL)
      genes <- unlist(strsplit(df[["core_enrichment"]][row_idx[1]], "/"))
      genes <- utils::head(genes, topGenes)
      data.frame(
        geneID = genes,
        Description = desc,
        stringsAsFactors = FALSE
      )
    }))

    if (is.null(gene_pairs) || nrow(gene_pairs) == 0) next

    # Gene order: unique in pathway-traversal order (keeps genes grouped)
    gene_order <- unique(gene_pairs[["geneID"]])

    # Convert to sankey long format
    sankey_long <- ggsankey::make_long(gene_pairs, geneID, Description)

    # Set factor levels: pathways (right) + genes (left, in pathway order)
    all_nodes <- c(term_levels, gene_order)
    sankey_long[["node"]] <- factor(sankey_long[["node"]],
                                     levels = all_nodes)

    # Colors
    n_nodes <- length(all_nodes)
    node_colors <- palette_colors(
      all_nodes, palette = palette, palcolor = palcolor
    )

    # ---- Step 1: base sankey (no text yet) ----
    sankeyExpand <- c(0.5, 1)
    ps <- ggplot2::ggplot(
      sankey_long,
      ggplot2::aes(x = x, next_x = next_x,
                   node = node, next_node = next_node,
                   fill = factor(node), label = node)
    ) +
      ggsankey::geom_sankey(
        flow.alpha = 0.5, flow.fill = "grey",
        node.fill = node_colors
      ) +
      ggplot2::theme_void() +
      ggplot2::scale_x_discrete(
        expand = ggplot2::expansion(mult = sankeyExpand)
      ) +
      ggplot2::theme(legend.position = "none")

    # ---- Step 2: extract node positions for text placement ----
    ps_data <- ggplot2::ggplot_build(ps)
    # Node layer is typically the 2nd data element from geom_sankey
    ps_data_info <- tryCatch(
      ps_data$data[[2]],
      error = function(e) NULL
    )

    if (is.null(ps_data_info) || nrow(ps_data_info) == 0) {
      # Fallback: use geom_sankey_text directly
      ps2 <- ps + ggsankey::geom_sankey_text(size = 2.5, hjust = 1)
    } else {
      # Re-add text labels at exact node positions (left-aligned)
      ps2 <- ps +
        ggplot2::geom_text(
          data = ps_data_info,
          ggplot2::aes(
            next_x = 0, next_node = 0,
            x = xmin, y = (ymin + ymax) / 2,
            label = label, hjust = 1
          ),
          fontface = "bold", size = 2.5
        )
    }

    # ---- Step 3: build right-side dot plot ----
    dot_df <- df
    dot_df[["neg_log10_p"]] <- -log10(dot_df[["pvalue"]])
    dot_df[["neg_log10_padj"]] <- -log10(dot_df[["p.adjust"]])

    pp <- ggplot2::ggplot(dot_df) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = neg_log10_p,
          y = Description,
          size = neg_log10_padj,
          fill = neg_log10_padj
        ),
        color = "black", shape = 21
      ) +
      ggplot2::scale_size(
        name = expression(-log[10](p.adjust)),
        range = c(2, 6)
      ) +
      ggplot2::scale_fill_viridis_c(
        name = expression(-log[10](p.adjust)),
        option = "plasma", direction = -1,
        guide = ggplot2::guide_colorbar(
          frame.colour = "black", ticks.colour = "black"
        )
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(colour = "black"),
        plot.background = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = expression(-log[10](pvalue)), y = "")

    # ---- Step 4: overlay dot plot on right side ----
    if (!is.null(ps_data_info) && nrow(ps_data_info) > 0) {
      xmin_dot <- max(ps_data_info$xmax, na.rm = TRUE)
      # Get y-range from pathway nodes (right column, x == "2")
      y_range <- ps_data_info[ps_data_info$x == "2", , drop = FALSE]
      if (nrow(y_range) == 0) y_range <- ps_data_info
      ymin_dot <- min(c(y_range$ymin, y_range$ymax), na.rm = TRUE)
      ymax_dot <- max(c(y_range$ymin, y_range$ymax), na.rm = TRUE)

      p <- ps2 +
        ggplot2::annotation_custom(
          grob = ggplot2::ggplotGrob(pp),
          xmin = xmin_dot - 0.05,
          xmax = 2 + sankeyExpand[2],
          ymin = ymin_dot - 4.5,
          ymax = ymax_dot + 0.25
        ) +
        ggplot2::labs(title = grp) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
        )
    } else {
      # Fallback without overlay
      p <- ps2 +
        ggplot2::labs(title = grp) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
        )
    }

    plist[[grp]] <- p
  }

  .combine_plist(plist, combine, ncol, nrow)
}
