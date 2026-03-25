# Train a Deconvolution DNN Model

Trains a deep neural network to estimate cell type proportions from bulk
RNA-seq expression data. The model is trained on simulated pseudo-bulk
samples generated from a single-cell RNA-seq reference dataset.

## Usage

``` r
DNN_deconv_train(
  reference,
  label_column,
  save_path = "deconv_model.pt",
  n_pseudobulk = 5000L,
  n_cells_per_sample = 500L,
  n_top_genes = 6000L,
  batch_key = NULL,
  hidden_size = 512L,
  num_epochs = 50L,
  learning_rate = 0.001,
  batch_size = 256L,
  loss_function = "mse",
  device = "auto",
  head_hidden_size = 256L,
  n_resnet_blocks = 4L,
  input_dropout_rate = 0.25,
  resnet_dropout_rate = 0.1,
  head_dropout_rate = 0.1,
  test_size = 0.1,
  early_stopping_patience = 10L,
  use_lr_scheduler = TRUE,
  random_state = 42L,
  recon_weight = 0.1
)
```

## Arguments

- reference:

  Path to an h5ad file or a Seurat object containing the scRNA-seq
  reference. Must include raw counts and cell type labels.

- label_column:

  Column name in `obs`/`meta.data` for cell type labels (e.g.
  `"cell_type"`).

- save_path:

  Path to save the trained model (`.pt` file). Parent directories are
  created automatically.

- n_pseudobulk:

  Number of pseudo-bulk samples to generate for training (default 5000).
  More samples generally improve accuracy.

- n_cells_per_sample:

  Number of cells to mix per pseudo-bulk sample (default 500).

- n_top_genes:

  Number of highly variable genes to select (default 6000).

- batch_key:

  Column in `obs` for batch-aware HVG selection. Set to `NULL` for no
  batch correction.

- hidden_size:

  Width of the shared backbone (default 512).

- num_epochs:

  Maximum training epochs (default 50).

- learning_rate:

  Learning rate (default 0.001).

- batch_size:

  Training batch size (default 256).

- loss_function:

  Loss function: `"mse"` (default) or `"kl"` (KL divergence).

- device:

  Compute device: `"auto"` (default, uses GPU if available), `"cpu"`, or
  `"cuda"`.

- head_hidden_size:

  Width of the proportion head (default 256).

- n_resnet_blocks:

  Number of ResNet blocks (default 4).

- input_dropout_rate:

  Dropout rate for input layer (default 0.25).

- resnet_dropout_rate:

  Dropout rate for ResNet blocks (default 0.1).

- head_dropout_rate:

  Dropout rate for proportion head (default 0.1).

- test_size:

  Fraction of pseudo-bulk for validation (default 0.1).

- early_stopping_patience:

  Epochs to wait before early stopping (default 10).

- use_lr_scheduler:

  Use cosine annealing LR scheduler (default TRUE).

- random_state:

  Random seed (default 42).

- recon_weight:

  Weight for reconstruction loss relative to proportion loss (default
  0.1). Higher values push the decoder to reconstruct input expression
  more faithfully. Set to 0 to disable reconstruction loss.

## Value

A named list with:

- model:

  Trained Python DeconvModel object.

- var_genes:

  Character vector of selected variable genes.

- cell_types:

  Character vector of cell type names (output order).

- history:

  data.frame with per-epoch training metrics (train_loss, val_loss,
  train_corr, val_corr).

- best_val_loss:

  Best validation loss achieved.

- best_val_corr:

  Best validation Pearson correlation.

## Details

The model architecture reuses the same ResNet backbone as
[`DNN_train`](https://hui950319.github.io/scMMR/reference/DNN_train.md)
but with key differences:

- **Input**: continuous log-CPM expression (not binary)

- **Output**: cell type proportions (softmax, sum = 1)

- **Training data**: simulated pseudo-bulk from scRNA-seq

## See also

[`DNN_deconv_predict`](https://hui950319.github.io/scMMR/reference/DNN_deconv_predict.md)
for predicting proportions,
[`DNN_train`](https://hui950319.github.io/scMMR/reference/DNN_train.md)
for cell type classification.

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)
use_scMMR_python(condaenv = "scMMR")

# Train deconvolution model from scRNA-seq reference
result <- DNN_deconv_train(
  reference    = "reference.h5ad",
  label_column = "cell_type",
  save_path    = "models/deconv_model.pt",
  n_pseudobulk = 5000L,
  device       = "auto"
)
print(result$cell_types)
print(result$best_val_corr)
} # }
```
