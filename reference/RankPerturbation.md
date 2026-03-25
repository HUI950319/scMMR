# Rank Cell Types by Transcriptional Perturbation

Measures how much each cell type's transcriptional state shifts between
two conditions in the DNN's learned embedding space. Cell types with
larger shifts are ranked higher (more "perturbed"). Inspired by Augur
but using multi-task DNN representations.

## Usage

``` r
RankPerturbation(
  embedding,
  cell_meta,
  cell_type_col = "cell_type_pred",
  condition_col = "group",
  conditions = NULL,
  method = c("wasserstein", "mmd", "energy", "classifier"),
  n_permutations = 500L,
  min_cells = 10L,
  n_pcs = 20L,
  seed = 42L,
  n_cores = 1L,
  balance_cells = FALSE,
  early_stop_alpha = 0.1,
  verbose = TRUE
)
```

## Arguments

- embedding:

  Numeric matrix (n_cells x d), typically `pred$shared_embedding` from
  `DNN_predict(return_embedding = TRUE)`. Rownames must be cell IDs.

- cell_meta:

  A data.frame (e.g. `seurat_obj@meta.data`). Rownames must be cell IDs.

- cell_type_col:

  Column name for cell type labels (default `"cell_type_pred"`).

- condition_col:

  Column name for condition/group labels (default `"group"`).

- conditions:

  Character vector of length 2 specifying which two conditions to
  compare (e.g. `c("PH", "SH")`). If `NULL` and exactly 2 unique
  conditions exist, they are used automatically.

- method:

  Distance metric: `"wasserstein"` (default, sliced Wasserstein distance
  with random projections), `"mmd"` (Maximum Mean Discrepancy with RBF
  kernel), `"energy"` (Energy distance, parameter-free), or
  `"classifier"` (LDA-based AUC, Augur-like).

- n_permutations:

  Number of permutations for p-value estimation (default 500). Set to 0
  to skip.

- min_cells:

  Minimum cells per cell-type-per-condition (default 10).

- n_pcs:

  Number of PCs for dimensionality reduction before distance computation
  (default 20). `NULL` to skip PCA.

- seed:

  Random seed (default 42).

- n_cores:

  Number of parallel cores for computing across cell types. Default 1
  (serial). On Unix, uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html); on
  Windows, uses
  [`parallel::parLapply`](https://rdrr.io/r/parallel/clusterApply.html).

- balance_cells:

  Logical. If `TRUE`, downsample the larger condition group to match the
  smaller one per cell type before distance computation (default
  `FALSE`).

- early_stop_alpha:

  Numeric. p-value threshold for adaptive early stopping of permutation
  tests. Permutations stop early when the running p-value is clearly
  above `2 * early_stop_alpha` (default 0.1). Set to 0 to disable early
  stopping.

- verbose:

  Logical. Print progress messages (default `TRUE`).

## Value

A named list:

- results:

  data.frame: cell_type, score, p_value, p_adj, n_cells_cond1,
  n_cells_cond2, rank.

- method:

  Method used.

- conditions:

  Two conditions compared.

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- DNN_predict(query, model_path, return_embedding = TRUE)
q1 <- Seurat::AddMetaData(toy_test, pred$predictions)

# Sliced Wasserstein (default)
res <- RankPerturbation(pred$shared_embedding, q1@meta.data,
                        conditions = c("PH", "SH"))

# Energy distance (fast, parameter-free)
res_e <- RankPerturbation(pred$shared_embedding, q1@meta.data,
                          conditions = c("PH", "SH"),
                          method = "energy")

# With parallel computation
res_p <- RankPerturbation(pred$shared_embedding, q1@meta.data,
                          conditions = c("PH", "SH"),
                          n_cores = 4)

print(res$results)
PlotPerturbation(res)
} # }
```
