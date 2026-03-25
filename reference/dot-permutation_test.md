# Permutation test with adaptive early stopping (internal)

Runs permutation test for a single cell type. After every
`check_interval` permutations, checks whether the running p-value is
clearly above `early_stop_alpha`. If so, terminates early to save
computation.

## Usage

``` r
.permutation_test(
  emb_ct,
  cond_ct,
  conditions,
  dist_fn,
  observed,
  n_permutations,
  early_stop_alpha = 0.1,
  check_interval = 100L
)
```

## Value

Numeric p-value.
