# Grid makeContext method for scMMR rasterised grobs

S3 method called by the grid rendering system. Renders the original grob
to an off-screen bitmap device and returns a `rasterGrob`.

## Usage

``` r
# S3 method for class 'scmmr_raster'
makeContext(x)
```

## Arguments

- x:

  A grob with class `"scmmr_raster"`.

## Value

A `rasterGrob`.
