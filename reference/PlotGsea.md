# Plot GSEA Results

Generates various types of plots for Gene Set Enrichment Analysis (GSEA)
results produced by
[`RunGsea`](https://hui950319.github.io/scMMR/reference/RunGsea.md).

## Usage

``` r
PlotGsea(
  gsea_result,
  type = c("line", "comparison", "bar", "lollipop", "ridge", "network", "enrichmap",
    "wordcloud", "volcano_nes", "circle", "sankey"),
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
  seed = 11
)
```

## Arguments

- gsea_result:

  A list returned by `RunGsea`, containing `enrichment` (data.frame),
  `results` (named list of `gseaResult` objects), and `input`.

- type:

  The type of plot to generate. One of:

  `"line"`

  :   Classic 3-panel GSEA running score plot (enrichment score curve +
      hit marks + ranked metric). Requires `group_use` to specify one
      group.

  `"comparison"`

  :   NES bubble plot comparing pathways across all groups.

  `"bar"`

  :   Bidirectional NES bar plot per group.

  `"lollipop"`

  :   Lollipop chart showing NES with `-log10(p.adjust)` colour
      gradient.

  `"ridge"`

  :   Ridge plot of gene set enrichment scores via
      `enrichplot::ridgeplot()`.

  `"network"`

  :   Gene-pathway network graph (requires igraph and ggrepel).

  `"enrichmap"`

  :   Enrichment map: pathway nodes connected by shared genes, clustered
      (requires igraph and ggforce).

  `"wordcloud"`

  :   Word cloud of pathway terms or core enrichment genes (requires
      ggwordcloud).

  `"volcano_nes"`

  :   Volcano-style GSEA plot with `-log10(p.adjust)` on x-axis and NES
      on y-axis. Points coloured by significance (Activated / Repressed
      / ns).

  `"circle"`

  :   Circular GSEA plot showing running score curves, hit marks and
      ranked-list heatmap in circos tracks (requires circlize).

  `"sankey"`

  :   Sankey flow diagram linking core enrichment genes to their
      pathways (requires ggsankey).

- group_use:

  Character vector of group names to include. For `type = "line"`,
  specify **one** group. For other types, defaults to all available
  groups.

- id_use:

  Character vector of pathway IDs to display. If `NULL` (default),
  pathways are selected automatically by `topTerm` and `direction`.

- topTerm:

  Number of top pathways to show per direction. Default 6.

- direction:

  Which enrichment direction to include: `"both"` (default), `"pos"`, or
  `"neg"`.

- padjustCutoff:

  Significance threshold for filtering pathways. Default 0.05.

- line_width:

  Line width for the running score curve. Default 1.5.

- line_color:

  Color(s) for running score lines. Default `"#6BB82D"`.

- n_coregene:

  Number of core enrichment genes to label in line plot. Default 0 (no
  labels).

- features_label:

  Character vector of specific genes to label in line plot. Overrides
  `n_coregene`.

- label.size:

  Size of gene labels. Default 3.5.

- character_width:

  Maximum character width before wrapping pathway names. Default 50.

- word_type:

  For `type = "wordcloud"`: source of words. `"term"` splits pathway
  names; `"feature"` uses core enrichment genes. Default `"term"`.

- word_size:

  For `type = "wordcloud"`: range of text sizes. Default `c(2, 8)`.

- topWord:

  For `type = "wordcloud"`: maximum number of words to display. Default
  100.

- network_layout:

  For `type = "network"`: igraph layout algorithm name (e.g., `"fr"`,
  `"kk"`, `"circle"`). Default `"fr"`.

- enrichmap_layout:

  For `type = "enrichmap"`: igraph layout algorithm. Default `"fr"`.

- enrichmap_cluster:

  For `type = "enrichmap"`: igraph community detection algorithm (e.g.,
  `"fast_greedy"`, `"louvain"`). Default `"fast_greedy"`.

- enrichmap_nlabel:

  For `type = "enrichmap"`: number of top terms per cluster to show as
  labels. Default 4.

- nes_cutoff:

  For `type = "volcano_nes"`: NES threshold for significance
  classification. Default 1.

- circ_type:

  For `type = "circle"`: inner track style. `"c"` = hit marks + heatmap,
  `"h"` = coloured segments + heatmap, `"m"` = tick marks only. Default
  `"c"`.

- topGenes:

  For `type = "sankey"`: number of top core genes to show per pathway.
  Default 5.

- palette:

  Color palette name for fills (passed to `palette_colors`). Default
  `"Spectral"`.

- palcolor:

  Custom color vector. Overrides `palette`.

- combine:

  Logical; whether to combine multiple panels into one plot using
  patchwork. Default `TRUE`.

- ncol:

  Number of columns when combining multiple plots. Default `NULL`
  (auto).

- nrow:

  Number of rows when combining multiple plots. Default `NULL` (auto).

- seed:

  Random seed. Default 11.

## Value

A `ggplot` object or a list of `ggplot` objects.

## Details

The function dispatches to internal helpers depending on `type`:

- `.plot_gsea_line` — classic 3-panel GSEA plot using patchwork.

- `.plot_gsea_comparison` — bubble plot with NES colour, setSize size,
  and significance border.

- `.plot_gsea_bar` — bidirectional NES bar plot.

- `.plot_gsea_ridge` — ridge plot via `enrichplot::ridgeplot()`.

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- RunDE(seu, group.by = "celltype")
res <- RunGsea(markers)

# Classic GSEA running score plot
PlotGsea(res, type = "line", group_use = "T_cell")

# Comparison bubble across all groups
PlotGsea(res, type = "comparison", topTerm = 5)

# NES bar plot for one group
PlotGsea(res, type = "bar", group_use = "T_cell", topTerm = 10)

# Ridge plot
PlotGsea(res, type = "ridge", group_use = "T_cell")

# Lollipop chart
PlotGsea(res, type = "lollipop", group_use = "T_cell", topTerm = 10)

# Gene-pathway network
PlotGsea(res, type = "network", group_use = "T_cell", topTerm = 5)

# Enrichment map
PlotGsea(res, type = "enrichmap", group_use = "T_cell", topTerm = 30)

# Word cloud of pathway terms
PlotGsea(res, type = "wordcloud", group_use = "T_cell",
         word_type = "term")

# Word cloud of core genes
PlotGsea(res, type = "wordcloud", group_use = "T_cell",
         word_type = "feature")

# GSEA volcano (NES vs -log10 padj)
PlotGsea(res, type = "volcano_nes", group_use = "T_cell")

# Circular GSEA plot
PlotGsea(res, type = "circle", group_use = "T_cell",
         id_use = c("HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1"))

# Sankey gene-pathway flow
PlotGsea(res, type = "sankey", group_use = "T_cell", topGenes = 5)
} # }
```
