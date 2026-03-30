# Rasterised point geom

Convenience wrapper that creates a `geom_point` layer and rasterises it
via
[`rasterise_layer`](https://hui950319.github.io/scMMR/reference/rasterise_layer.md).

## Usage

``` r
geom_point_rast(..., raster.dpi = 300, raster.dev = NULL, raster.scale = 1)
```

## Arguments

- ...:

  Additional arguments passed to
  [`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html).

- raster.dpi:

  Numeric. Raster resolution. Default `300`.

- raster.dev:

  Character. Backend device. Default `NULL` (auto-detect).

- raster.scale:

  Numeric. Scale factor. Default `1`.

## Value

A rasterised ggplot2 layer.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
df <- data.frame(x = rnorm(1e5), y = rnorm(1e5))
ggplot(df, aes(x, y)) + geom_point_rast(alpha = 0.3)
} # }
```
