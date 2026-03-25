# Plot Embedding Evaluation Results

Visualise results from
[`EvaluateEmbedding`](https://hui950319.github.io/scMMR/reference/EvaluateEmbedding.md).

## Usage

``` r
PlotEmbeddingEval(eval_result, which = "all", base_size = 12)
```

## Arguments

- eval_result:

  Output from
  [`EvaluateEmbedding()`](https://hui950319.github.io/scMMR/reference/EvaluateEmbedding.md).

- which:

  Character. Which plot(s) to draw: `"all"` (default, returns a list of
  all available plots), `"elbow"`, `"cumvar"`, `"knn"`, `"distcor"`, or
  `"silhouette"`.

- base_size:

  Numeric. Base font size (default 12).

## Value

A single ggplot object or a named list of ggplot objects (when
`which = "all"`).

## Examples

``` r
if (FALSE) { # \dontrun{
eval_res <- EvaluateEmbedding(pred$shared_embedding, seurat_obj = seu)

# All plots as a list
plots <- PlotEmbeddingEval(eval_res)
plots$elbow
plots$knn

# Single plot
PlotEmbeddingEval(eval_res, which = "silhouette")

# Combine with patchwork
library(patchwork)
plots$elbow + plots$cumvar + plots$knn + plots$distcor
} # }
```
