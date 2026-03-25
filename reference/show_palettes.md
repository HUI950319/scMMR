# Display Available Color Palettes

Show all built-in color palettes as a stacked bar chart.

## Usage

``` r
show_palettes(
  palette_names = NULL,
  type = c("discrete", "continuous"),
  index = NULL,
  return_palettes = FALSE
)
```

## Arguments

- palette_names:

  Optional character vector of palette names to display. If `NULL`, all
  palettes are shown.

- type:

  Character vector. Filter by palette type: `"discrete"` and/or
  `"continuous"`. Default shows both.

- index:

  Optional integer vector. Show only palettes at these positions.

- return_palettes:

  Logical. If `TRUE`, return the palette list invisibly instead of just
  printing. Default `FALSE`.

## Value

If `return_palettes` is `TRUE`, returns a named list of color vectors.
Otherwise, prints the plot and returns palette names invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
show_palettes()
show_palettes(type = "discrete")
show_palettes(palette_names = c("Paired", "npg", "nejm", "Spectral"))
pals <- show_palettes(return_palettes = TRUE)
} # }
```
