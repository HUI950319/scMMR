# MMD with RBF kernel (internal)

Uses multiple bandwidth values (median heuristic x scales) to avoid
sensitivity to a single sigma. Returns the raw (unbiased) MMD-squared
without clamping so that permutation tests remain valid.

## Usage

``` r
.mmd_rbf(mat1, mat2, sigma = NULL, max_cells = 2000L)
```

## Details

Optimisations vs. original:

- Downsamples groups larger than `max_cells` to control memory.

- Pre-computes squared distance matrices once (shared across bandwidths)
  instead of recomputing kernel matrices per sigma.
