# Energy distance (E-distance) (internal)

A fast, parameter-free distributional distance. \$\$E(X,Y) = 2 E\\X -
Y\\ - E\\X - X'\\ - E\\Y - Y'\\\$\$

## Usage

``` r
.energy_distance(mat1, mat2, max_cells = 2000L)
```

## Details

For large populations, downsamples to `max_cells` per group. Cross-group
distances are computed in chunks to avoid memory issues.
