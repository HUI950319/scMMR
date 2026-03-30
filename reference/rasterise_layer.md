# Rasterise a ggplot2 layer

Wraps an existing ggplot2 layer so that its graphical output is rendered
to an off-screen device and converted to a bitmap (`rasterGrob`). This
dramatically reduces file size and rendering time for scatter plots with
many points, especially when saving to PDF/SVG.

## Usage

``` r
rasterise_layer(layer, dpi = 300, dev = NULL, scale = 1)
```

## Arguments

- layer:

  A ggplot2 layer object (e.g. the return value of `geom_point(...)`).

- dpi:

  Numeric scalar. Resolution in dots per inch. Default `300`.

- dev:

  Character. Rendering backend: `"cairo"`, `"ragg"`, or `"png"`. Default
  `NULL` (auto-detect best available).

- scale:

  Numeric scalar \> 0. Scaling factor applied to the rendered bitmap
  dimensions. Values \< 1 reduce resolution (faster), values \> 1
  increase it. Default `1`.

## Value

A modified ggplot2 layer whose `draw_geom` method produces rasterised
grobs.

## Details

The function is self-contained and does not depend on scattermore or
ggrastr. It supports three rendering backends with automatic fallback:

1.  **Cairo** (default if available) – in-memory raster device, fastest.

2.  **ragg** – high-quality AGG-based capture device.

3.  **png** – file-based fallback using
    [`grDevices::png()`](https://rdrr.io/r/grDevices/png.html) and
    [`png::readPNG()`](https://rdrr.io/pkg/png/man/readPNG.html).

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
df <- data.frame(x = rnorm(1e5), y = rnorm(1e5))
ggplot(df, aes(x, y)) + rasterise_layer(geom_point(alpha = 0.3))
} # }
```
