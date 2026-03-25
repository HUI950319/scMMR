# Generate Named Color Vector from Palette

Create a named color vector for discrete or continuous data, similar to
`thisplot::palette_colors()`. Can be used directly with
`ggplot2::scale_color_manual(values = ...)`.

## Usage

``` r
palette_colors(
  x,
  n = 100,
  palette = "Paired",
  palcolor = NULL,
  type = c("auto", "discrete", "continuous"),
  matched = FALSE,
  reverse = FALSE,
  NA_keep = FALSE,
  NA_color = "grey80"
)
```

## Arguments

- x:

  A vector of character, factor, or numeric values. If missing, numeric
  values `1:n` are used.

- n:

  Integer. Number of colors to return for numeric values. Default 100.

- palette:

  Character. Name of a built-in palette. Use
  [`show_palettes()`](https://hui950319.github.io/scMMR/reference/show_palettes.md)
  to see all available names. Default `"Paired"`.

- palcolor:

  Optional character vector of custom colors. When provided, overrides
  `palette`.

- type:

  Character, one of `"auto"`, `"discrete"`, or `"continuous"`. Default
  `"auto"` detects from `x`.

- matched:

  Logical. If `TRUE`, return a color for each element of `x` (same
  length). If `FALSE` (default), return a named vector of unique colors.

- reverse:

  Logical. Reverse the color order. Default `FALSE`.

- NA_keep:

  Logical. Keep color for `NA` values. Default `FALSE`.

- NA_color:

  Character. Color for `NA`. Default `"grey80"`.

## Value

A named character vector of hex color codes.

## Examples

``` r
if (FALSE) { # \dontrun{
palette_colors(c("A", "B", "C"))
palette_colors(c("A", "B", "C"), palette = "Dark2")
palette_colors(c("A", "B", "C"), palcolor = c("red", "blue", "green"))
palette_colors(1:100, palette = "Spectral")
} # }
```
