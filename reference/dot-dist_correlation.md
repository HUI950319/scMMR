# Distance correlation between two spaces (internal)

Randomly sample cell pairs, compute pairwise Euclidean distance in both
spaces, and report Spearman rank correlation.

## Usage

``` r
.dist_correlation(emb, ref, n_pairs = 5000L, seed = 42L)
```
