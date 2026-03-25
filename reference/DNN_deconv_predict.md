# Predict Cell Type Proportions from Bulk RNA-seq

Uses a trained deconvolution DNN model (from
[`DNN_deconv_train`](https://hui950319.github.io/scMMR/reference/DNN_deconv_train.md))
to estimate cell type proportions from bulk RNA-seq expression data.

## Usage

``` r
DNN_deconv_predict(
  bulk_expr,
  model_path,
  device = "auto",
  adaptive = FALSE,
  adaptive_mode = "overall",
  adaptive_max_iter = 3L,
  adaptive_steps = 300L,
  adaptive_lr = 1e-04
)
```

## Arguments

- bulk_expr:

  Bulk expression data. Accepts multiple formats:

  - A **matrix** or **data.frame** (genes x samples). Rownames = gene
    names, colnames = sample names.

  - A **CSV file path** (genes x samples, first column = gene names or
    rownames).

  - An **h5ad file path**.

  - A **Seurat object** (uses RNA assay counts).

- model_path:

  Path to the trained deconvolution model (`.pt` file created by
  [`DNN_deconv_train`](https://hui950319.github.io/scMMR/reference/DNN_deconv_train.md)).

- device:

  Compute device: `"auto"` (default), `"cpu"`, or `"cuda"`.

- adaptive:

  Logical; if `TRUE`, run TAPE-style adaptive stage to refine
  proportions and extract cell-type-specific GEP (default `FALSE`).

- adaptive_mode:

  Adaptive mode: `"overall"` (default) adapts on all samples jointly
  (one shared GEP), or `"high-resolution"` adapts per-sample (individual
  GEPs, slower).

- adaptive_max_iter:

  Number of alternating decoder/encoder rounds (default 3).

- adaptive_steps:

  Gradient steps per phase per round (default 300).

- adaptive_lr:

  Learning rate for adaptive optimization (default 1e-4).

## Value

When `adaptive = FALSE`: a `data.frame` with `sample_id` and one column
per cell type (backward compatible).

When `adaptive = TRUE`: a named `list` with:

- proportions:

  data.frame of cell type proportions.

- sigmatrix:

  Cell-type-specific GEP. For `"overall"` mode: a data.frame (cell_types
  x genes). For `"high-resolution"` mode: a 3D array (samples x
  cell_types x genes).

## Details

The function handles gene alignment automatically: genes present in both
the bulk data and the model are used; missing genes are zero-filled.
Expression is normalized to log1p(CPM) to match the training procedure.

When `adaptive = TRUE`, the model runs a TAPE-style adaptive stage that
alternately fine-tunes the decoder and encoder on the target bulk data,
producing tissue-adapted cell-type-specific gene expression profiles
(GEP / sigmatrix).

## See also

[`DNN_deconv_train`](https://hui950319.github.io/scMMR/reference/DNN_deconv_train.md)
for training the model.

## Examples

``` r
if (FALSE) { # \dontrun{
library(scMMR)
use_scMMR_python(condaenv = "scMMR")

# Standard prediction (no adaptive)
props <- DNN_deconv_predict(
  bulk_expr  = bulk_matrix,
  model_path = "models/deconv_model.pt"
)

# Adaptive prediction with GEP extraction
result <- DNN_deconv_predict(
  bulk_expr      = bulk_matrix,
  model_path     = "models/deconv_model.pt",
  adaptive       = TRUE,
  adaptive_mode  = "overall"
)
result$proportions   # cell type proportions
result$sigmatrix     # cell-type-specific GEP (K x genes)
} # }
```
