# Pathway Bubble Plot

Create a bubble plot visualizing pathway analysis results across cell
types. Automatically adapts visual encoding to GSEA or SCPA results
based on the `Method` column.

## Usage

``` r
PlotPathwayBubble(
  data,
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
  facet.by.sign = NULL
)
```

## Arguments

- data:

  A data.frame returned by
  [`RunPathwayAnalysis`](https://hui950319.github.io/scMMR/reference/RunPathwayAnalysis.md).
  Must contain columns: `Celltype`, `Pathway`, `Score`, `Sign`,
  `Method`, and at least one of `QValue` or `AdjPValue`.

- top.n:

  Number of top pathways to display. Default: 10.

- size.by:

  What to map to point size:

  `"pvalue"`

  :   Point size = `-log10(pvalue.col)`. Default for GSEA.

  `"score"`

  :   Point size = `|Score|` (i.e. \|NES\| or \|FC\|). Default for SCPA.

  `"qvalue"`

  :   Point size = raw `QValue`.

  Default: `NULL` (auto-selected based on method).

- pvalue.col:

  Column name for significance filtering and `size.by = "pvalue"`
  mapping. One of `"QValue"` or `"AdjPValue"`. Default: auto-detected
  (`"QValue"` for GSEA, `"AdjPValue"` for SCPA).

- pvalue.cutoff:

  Significance threshold for selecting top pathways. Default: 0.05.

- select.by:

  Strategy for choosing top pathways:

  `"frequency"`

  :   Pathways most frequently significant across cell types (default).

  `"score"`

  :   Pathways with highest average \|Score\|.

- pathways:

  Character vector of specific pathway names to display. Overrides
  `top.n` and `select.by` when provided.

- order.by:

  How to order pathways on the y-axis:

  `"frequency"`

  :   By count of significant appearances (default).

  `"score"`

  :   By mean Score value (low to high).

  `"sign"`

  :   By mean Score (high to low, Activated first).

- base.size:

  Base font size for the plot. Default: 15.

- colors.gsea:

  Color gradient vector for GSEA NES fill. Default:
  `c("blue", "white", "red")`.

- colors.scpa:

  Named vector for SCPA Sign fill. Default:
  `c(Activated = "red", Repressed = "blue")`.

- facet.by.sign:

  Logical; whether to facet by Activated/Repressed. Default: `TRUE` for
  GSEA, `FALSE` for SCPA. Set explicitly to override.

## Value

A `ggplot` object.

## Visual encoding

- GSEA:

  - **fill**: NES value (continuous blue-white-red gradient).

  - **size**: determined by `size.by` (default: `-log10(q-value)`).

  - **facet**: Activated / Repressed (by default).

- SCPA:

  - **fill**: Direction (Activated = red, Repressed = blue).

  - **size**: determined by `size.by` (default: `|FC|`).

  - **no facet** by default.

## See also

[`RunPathwayAnalysis`](https://hui950319.github.io/scMMR/reference/RunPathwayAnalysis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# GSEA: default (fill = NES, size = -log10(q))
PlotPathwayBubble(gsea_res)

# GSEA: size by |NES|
PlotPathwayBubble(gsea_res, size.by = "score")

# SCPA: default (fill = Sign, size = |FC|)
PlotPathwayBubble(scpa_res)

# SCPA: size by q-value
PlotPathwayBubble(scpa_res, size.by = "qvalue")

# Show specific pathways
PlotPathwayBubble(gsea_res, pathways = c("Oxidative Phosphorylation",
                                          "Interferon Alpha Response"))

# Top 15, selected by score, ordered by score
PlotPathwayBubble(gsea_res, top.n = 15, select.by = "score",
                  order.by = "score")
} # }
```
