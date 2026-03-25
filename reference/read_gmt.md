# Read a GMT (Gene Matrix Transposed) file

Parse a GMT (Gene Matrix Transposed) file and return gene sets as a
named list. Each element contains the gene symbols belonging to that
pathway.

## Usage

``` r
read_gmt(path)
```

## Arguments

- path:

  Path to the GMT file.

## Value

A named list of character vectors (pathway name → gene symbols).

## Examples

``` r
if (FALSE) { # \dontrun{
gmt <- read_gmt("h.all.v2022.1.Hs.symbols.gmt")
names(gmt)[1:5]
} # }
```
