# Predict Cell Types Using a Trained Multi-Task DNN Model

Loads a trained multi-task model, runs predictions on query data, and
optionally computes gene importance via Integrated Gradients. When a GMT
file is provided, per-cell pathway importance scores are also computed
by aggregating gene-level IG attributions to pathways.

## Usage

``` r
DNN_predict(
  query,
  model_path,
  save_path = NULL,
  true_label_col = "cell_type",
  explain = FALSE,
  top_k_global = 15L,
  top_k_class = 10L,
  n_cells_explain = 50L,
  pathway_gmt = NULL,
  pathway_min_genes = 5L,
  pathway_n_steps = 20L,
  return_embedding = FALSE,
  device = "auto"
)
```

## Arguments

- query:

  Path to an h5ad file or a Seurat object containing the query
  single-cell data.

- model_path:

  Path to the saved model file (`.pt`).

- save_path:

  Path to save results as HDF5 (`.h5`). Set to `NULL` to skip saving.

- true_label_col:

  Column in `obs` with ground-truth labels for accuracy evaluation. Set
  to `NULL` if no labels available.

- explain:

  Logical. If `TRUE`, compute gene importance via Integrated Gradients
  (default `FALSE`).

- top_k_global:

  Number of top genes for global importance (default 15).

- top_k_class:

  Number of top genes per cell type (default 10).

- n_cells_explain:

  Number of cells for attribution baseline (speed tradeoff, default 50).

- pathway_gmt:

  Path to a GMT file for pathway-level scoring. Requires
  `explain = TRUE`. When provided, computes per-cell pathway importance
  via `|IG attributions| × GMT mask / pathway_size`. Built-in GMT files
  are available in the package
  (`system.file("extdata/gmt", package = "scMMR")`): `reactome.gmt`,
  `GO_bp.gmt`, `TF.gmt`, `immune.gmt` (human), plus `m_reactome.gmt`,
  `m_GO_bp.gmt`, `m_TF.gmt` (mouse). Set to `NULL` to skip pathway
  scoring.

- pathway_min_genes:

  Minimum number of overlapping genes between a pathway and the model's
  variable genes for the pathway to be retained (default 5).

- pathway_n_steps:

  Number of Integrated Gradients interpolation steps for pathway scoring
  (default 20). Lower values are faster but less precise. For gene
  importance `n_steps=50` is used; pathway scoring tolerates fewer steps
  since attributions are aggregated across many genes.

- return_embedding:

  Logical. If `TRUE`, include the 512-dim shared embedding in the output
  (default `FALSE`).

- device:

  Compute device: `"auto"` (default), `"cpu"`, or `"cuda"`.

## Value

A named list with:

- predictions:

  data.frame with columns: cell_id, cell_type_pred, confidence, is_ood,
  umap_1_pred, umap_2_pred. `cell_type_pred` always contains the model's
  best-guess label (never `"Unknown"`). Use `is_ood` or `confidence` for
  downstream filtering.

- shared_embedding:

  Matrix of shape (n_cells, 512) — shared representation from the model
  backbone. Only included when `return_embedding = TRUE`.

- imp_global:

  data.frame (gene, importance) if `explain=TRUE`, otherwise `NULL`.
  Computed via Integrated Gradients.

- imp_per_class:

  data.frame (cell_type, n_cells, rank, gene, importance) if
  `explain=TRUE`, otherwise `NULL`.

- pathway_scores:

  Matrix of shape (n_cells, n_pathways) with per-cell pathway importance
  scores. Rownames are cell IDs; column names are pathway names. Only
  included when `pathway_gmt` is provided and `explain = TRUE`. Scores
  represent each pathway's discriminative contribution (via Integrated
  Gradients), not expression enrichment.

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)
use_scMMR_python(condaenv = "scanpy-env")

# Basic prediction
pred <- DNN_predict(
  query      = "query.h5ad",
  model_path = "models/my_model.pt",
  explain    = TRUE,
  device     = "cpu"
)
head(pred$predictions)
pred$imp_global

# With pathway scoring (GMT file)
gmt <- system.file("extdata/gmt", "reactome.gmt", package = "scMMR")
pred2 <- DNN_predict(
  query       = "query.h5ad",
  model_path  = "models/my_model.pt",
  explain     = TRUE,
  pathway_gmt = gmt,
  device      = "cpu"
)
dim(pred2$pathway_scores)  # n_cells x n_pathways
} # }
```
