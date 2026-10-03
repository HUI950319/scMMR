# Run Monocle3 trajectory inference on a Seurat object

Thin wrapper around the public monocle3 API that converts a Seurat
object to a `cell_data_set`, runs the standard
`cluster_cells -> learn_graph -> order_cells` pipeline, and writes
clusters, partitions and pseudotime back into the Seurat object.

## Usage

``` r
RunMonocle3(
  srt,
  assay = NULL,
  layer = "counts",
  reduction = "umap",
  root_cells = NULL,
  root_pr_nodes = NULL,
  root_group = NULL,
  root_group_value = NULL,
  use_partition = TRUE,
  close_loop = TRUE,
  k = 50,
  cluster_method = c("louvain", "leiden"),
  resolution = NULL,
  partition_qval = 0.05,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Assay to pull counts from. Default `NULL` (active assay).

- layer:

  Layer/slot to pull counts from. Default `"counts"`.

- reduction:

  Name of the reduction in `srt` to use as the 2D embedding for
  Monocle3. Default `"umap"`.

- root_cells:

  Character vector of cell barcodes to use as the trajectory root.
  Either `root_cells`, `root_pr_nodes`, or `root_group` must be
  supplied.

- root_pr_nodes:

  Character vector of principal-graph node names (e.g. `"Y_1"`) to use
  as the root.

- root_group:

  A length-1 character: name of a metadata column. The most-stem-like
  value in that column will be used to pick root cells (cells whose
  values match `root_group_value`).

- root_group_value:

  The level of `root_group` to use as root (e.g. `"Stem"`).

- use_partition:

  Passed to `monocle3::learn_graph()`. If `TRUE` (default), each
  partition gets its own disjoint graph.

- close_loop:

  Passed to `monocle3::learn_graph()`. Default `TRUE`.

- k:

  Number of nearest neighbours for `cluster_cells`. Default 50.

- cluster_method:

  Clustering method for `cluster_cells` (`"louvain"` or `"leiden"`).
  Default `"louvain"`.

- resolution:

  Resolution parameter for `cluster_cells`. Default `NULL` (Monocle3
  auto-selects).

- partition_qval:

  q-value threshold for partitioning. Default 0.05.

- seed:

  Random seed. Default 11.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

The input Seurat object with:

- Metadata columns added: `Monocle3_clusters`, `Monocle3_partitions`,
  `Monocle3_Pseudotime`.

- `srt@tools$Monocle3` list containing: `cds` (the fitted
  `cell_data_set`), `edge_df` (principal-graph edges in 2D), `node_df`
  (branch / leaf / internal nodes in 2D), `trajectory` (a list of
  ggplot2 layers — `geom_segment` for the principal graph; just add to
  any ggplot), `milestones` (a list of ggplot2 layers — leaf nodes in
  black, branch nodes in red), `root_cells` or `root_pr_nodes`
  (whichever was used).

## Details

This wrapper deliberately uses only public Monocle3 functions (no
`monocle3:::` internals) and is non-interactive: roots must be supplied
as parameters.

**Choosing roots.** You must give the function a way to pick the
starting point of pseudotime. The three options, in order of preference:

1.  `root_cells`: explicit cell barcodes you trust as the start.

2.  `root_pr_nodes`: principal-graph node names. Run the function once
    without a root to see the available nodes in the returned `node_df`,
    then re-run with the chosen node(s).

3.  `root_group` + `root_group_value`: pick all cells whose metadata
    column `root_group` equals `root_group_value`.

If none is supplied, the function will run everything up to
`learn_graph` and return without ordering, so you can inspect `node_df`
and decide.

**Plotting the trajectory.** The returned `edge_df` and `node_df` are
plain data frames in the coordinates of `reduction`, so you can overlay
them on any `ggplot2` scatter:


      library(ggplot2)
      d <- cbind(srt@meta.data, srt[["umap"]]@cell.embeddings)
      ggplot(d, aes(umap_1, umap_2)) +
        geom_point(aes(colour = Monocle3_Pseudotime), size = .3) +
        geom_segment(data = srt@tools$Monocle3$edge_df,
                     aes(x = x, y = y, xend = xend, yend = yend),
                     inherit.aes = FALSE) +
        geom_point(data = srt@tools$Monocle3$node_df,
                   aes(x, y), inherit.aes = FALSE,
                   shape = 21, fill = "red", size = 2)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(scMMR)

seu <- qs::qread(system.file("extdata", "toy_test.qs", package = "scMMR"))

# First pass: no root, just to see available nodes
seu <- RunMonocle3(seu, reduction = "umap")
head(seu@tools$Monocle3$node_df)

# Second pass: pick a root by metadata
seu <- RunMonocle3(
  seu,
  reduction       = "umap",
  root_group      = "celltype",
  root_group_value = "B cells"
)
head(seu$Monocle3_Pseudotime)
} # }
```
