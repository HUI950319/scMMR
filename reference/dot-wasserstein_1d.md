# 1D Wasserstein distance (sort-based) (internal)

Uses sorted samples with evenly-spaced index selection, which is
significantly faster than the quantile-based approach for large vectors.

## Usage

``` r
.wasserstein_1d(x, y, n_quantiles = 200L)
```
