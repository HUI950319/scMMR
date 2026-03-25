# Train a Multi-Task DNN Model

One-step function that loads reference data, selects highly variable
genes, creates a multi-task neural network model, trains it, configures
the OOD confidence threshold, and saves the trained model. The model
jointly learns cell type classification and embedding (e.g. UMAP)
regression.

## Usage

``` r
DNN_train(
  input,
  label_column,
  embedding_key,
  save_path = "model.pt",
  n_top_genes = 6000L,
  batch_key = NULL,
  ood_threshold = 0.5,
  hidden_size = 512L,
  num_epochs = 50L,
  learning_rate = 0.001,
  batch_size = 256L,
  device = "auto",
  head_hidden_size = 256L,
  n_resnet_blocks = 4L,
  input_dropout_rate = 0.25,
  resnet_dropout_rate = 0.1,
  head_dropout_rate = 0.1,
  test_size = 0.2,
  early_stopping_patience = 10L,
  early_stopping_monitor = "cls",
  use_gradnorm = TRUE,
  gradnorm_alpha = 0.15,
  label_smoothing = 0.1,
  use_lr_scheduler = TRUE,
  use_uncertainty = FALSE,
  random_state = 42L
)
```

## Arguments

- input:

  Path to an h5ad file or a Seurat object containing the reference
  single-cell data. Must include cell type labels in `obs` and embedding
  coordinates in `obsm`.

- label_column:

  Column name in `obs` for cell type labels (e.g. `"cell_type"`).

- embedding_key:

  Key in `obsm` for embedding coordinates (e.g. `"umap"`, `"X_umap"`).

- save_path:

  Path to save the trained model (`.pt` file). Parent directories are
  created automatically.

- n_top_genes:

  Number of highly variable genes to select (default 6000).

- batch_key:

  Column in `obs` for batch-aware HVG selection. Set to `NULL` for no
  batch correction.

- ood_threshold:

  OOD (out-of-distribution) confidence threshold (default 0.5). Cells
  with confidence below this are flagged as `is_ood = TRUE` in
  predictions. Cell type labels always retain the model's best-guess;
  users can filter by `is_ood` or `confidence` downstream.

- hidden_size:

  Width of the shared backbone (default 512).

- num_epochs:

  Maximum training epochs (default 50).

- learning_rate:

  Learning rate (default 0.001).

- batch_size:

  Training batch size (default 256).

- device:

  Compute device: `"auto"` (default, uses GPU if available), `"cpu"`
  (force CPU), or `"cuda"` (require GPU).

- head_hidden_size:

  Width of task-specific heads (default 256).

- n_resnet_blocks:

  Number of ResNet blocks (default 4).

- input_dropout_rate:

  Dropout rate for input layer (default 0.25).

- resnet_dropout_rate:

  Dropout rate for ResNet blocks (default 0.1).

- head_dropout_rate:

  Dropout rate for task heads (default 0.1).

- test_size:

  Fraction of data for validation (default 0.2).

- early_stopping_patience:

  Epochs to wait before early stopping (default 10).

- early_stopping_monitor:

  Metric to monitor: `"cls"` or `"reg"` (default `"cls"`).

- use_gradnorm:

  Use GradNorm for multi-task balancing (default TRUE).

- gradnorm_alpha:

  GradNorm alpha parameter (default 0.15).

- label_smoothing:

  Label smoothing factor (default 0.1).

- use_lr_scheduler:

  Use cosine annealing LR scheduler (default TRUE).

- use_uncertainty:

  Use uncertainty-aware regression (default FALSE).

- random_state:

  Random seed (default 42).

## Value

A named list with:

- model:

  Trained Python MultiTaskModel object.

- var_genes:

  Character vector of selected variable genes.

- history:

  data.frame with per-epoch training metrics.

- best_val_acc:

  Best validation accuracy (%).

- best_val_rmse:

  Best validation RMSE.

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)
use_scMMR_python(condaenv = "scanpy-env")

result <- DNN_train(
  input         = "reference.h5ad",
  label_column  = "cell_type",
  embedding_key = "umap",
  save_path     = "models/my_model.pt",
  device        = "auto"
)
print(result$best_val_acc)
print(result$history)
} # }
```
