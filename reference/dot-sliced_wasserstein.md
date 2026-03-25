# Sliced Wasserstein distance with random projections (internal)

Projects high-dimensional distributions onto random 1D directions and
averages 1D Wasserstein distances. When `n_projections = NULL`, falls
back to coordinate-axis projections for backward compatibility.

## Usage

``` r
.sliced_wasserstein(mat1, mat2, n_projections = 50L)
```
