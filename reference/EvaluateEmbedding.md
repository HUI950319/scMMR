# Evaluate DNN Embedding Quality

Assess the information content of the DNN shared embedding and compare
it with PCA space. Returns intrinsic dimensionality, KNN overlap,
distance correlation, RV coefficient, and silhouette comparison.

## Usage

``` r
EvaluateEmbedding(
  embedding,
  seurat_obj = NULL,
  pca_ref = NULL,
  cell_type_col = "cell_type",
  n_pcs = 50L,
  k_values = c(10L, 20L, 50L, 100L, 200L),
  n_pairs = 5000L,
  max_cells = 20000L,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- embedding:

  Numeric matrix (n_cells x d), typically `pred$shared_embedding` from
  `DNN_predict(return_embedding = TRUE)`. Rownames must be cell IDs.

- seurat_obj:

  A Seurat object (optional). If provided, PCA embeddings and cell type
  labels are extracted automatically. At least one of `seurat_obj` or
  `pca_ref` is needed for consistency metrics.

- pca_ref:

  Numeric matrix (n_cells x n_pcs), alternative to `seurat_obj`.
  Rownames must be cell IDs.

- cell_type_col:

  Character. Column name in Seurat metadata or a named character vector
  of cell type labels. Used for silhouette comparison. Default
  `"cell_type"`.

- n_pcs:

  Integer. Number of PCs to compute on the DNN embedding for intrinsic
  dimensionality assessment (default 50).

- k_values:

  Integer vector. Values of k for KNN overlap curve (default
  `c(10, 20, 50, 100, 200)`).

- n_pairs:

  Integer. Number of cell pairs for distance correlation (default 5000).

- max_cells:

  Integer. Downsample to this many cells for performance (default
  20000).

- seed:

  Integer. Random seed (default 42).

- verbose:

  Logical. Print progress messages (default TRUE).

## Value

A named list with components:

- pca_of_embedding:

  List: sdev, var_explained, cum_var, effective_dim

- consistency:

  List: knn_overlap (data.frame), dist_cor (list), rv_coef (numeric),
  silhouette (list). NULL if no PCA reference.

- summary:

  Character string summarising key findings.

- params:

  Parameters used.

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- DNN_predict(query, model_path, return_embedding = TRUE)
eval_res <- EvaluateEmbedding(pred$shared_embedding,
                               seurat_obj = seu,
                               cell_type_col = "cell_type")
print(eval_res$summary)
PlotEmbeddingEval(eval_res)
PlotEmbeddingEval(eval_res, which = "elbow")
} # }
```
