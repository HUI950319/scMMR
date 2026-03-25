# Differential Abundance Testing via KNN Neighborhoods

Builds a KNN graph on the DNN's shared embedding, samples representative
neighborhoods, and tests for differential cell abundance between
conditions. Inspired by the Milo framework but operating on DNN-learned
representations.

## Usage

``` r
RankPercent(
  embedding,
  cell_meta,
  cell_type_col = "cell_type_pred",
  condition_col = "group",
  conditions = NULL,
  k = 30L,
  prop_sample = 0.1,
  test = c("fisher", "nb_glm"),
  min_cells = 20L,
  fdr_threshold = 0.1,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- embedding:

  Numeric matrix (n_cells x d), typically `pred$shared_embedding`.

- cell_meta:

  A data.frame. Rownames must be cell IDs.

- cell_type_col:

  Column name for cell type labels (default `"cell_type_pred"`).

- condition_col:

  Column name for condition/group labels (default `"group"`).

- conditions:

  Optional character vector of length 2. If provided, only these two
  conditions are tested. If `NULL` (default), all conditions are used;
  Fisher's or chi-squared test is applied.

- k:

  Number of nearest neighbors (default 30).

- prop_sample:

  Proportion of cells sampled as neighborhood centers (default 0.1).

- test:

  Statistical test: `"fisher"` (default) or `"nb_glm"` (negative
  binomial GLM via
  [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html)).

- min_cells:

  Minimum cells in a neighborhood to include (default 20).

- fdr_threshold:

  FDR threshold for summary (default 0.1).

- seed:

  Random seed (default 42).

- verbose:

  Logical. Print progress messages (default `TRUE`).

## Value

A named list:

- da_results:

  data.frame: nhood_idx, nhood_cell_id, n_cells, counts per condition,
  logFC, p_value, p_adj, cell_type_majority.

- cell_da_scores:

  Named numeric vector of per-cell DA scores. Ready for
  `Seurat::AddMetaData()`.

- cell_type_summary:

  data.frame: cell_type, n_nhoods, n_da_nhoods, mean_logFC, fraction_da.

- params:

  List of parameters used.

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- DNN_predict(query, model_path, return_embedding = TRUE)
q1 <- Seurat::AddMetaData(toy_test, pred$predictions)
da <- RankPercent(pred$shared_embedding, q1@meta.data,
                  conditions = c("PH", "SH"), k = 30)
q1 <- Seurat::AddMetaData(q1, da$cell_da_scores, col.name = "da_score")
PlotPercent(da)
} # }
```
