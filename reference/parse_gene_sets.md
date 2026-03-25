# Parse gene sets from various input formats

Converts gene sets from multiple input formats into a standardized named
list. Supports named list (pass-through), data.frame (term-to-gene
mapping), or a file path to a GMT file.

## Usage

``` r
parse_gene_sets(gene.sets)
```

## Arguments

- gene.sets:

  Gene sets in one of three formats:

  - A named list of character vectors (returned as-is).

  - A data.frame with at least two columns: term (1st) and gene (2nd).

  - A file path to a GMT file (parsed via
    [`read_gmt`](https://hui950319.github.io/scMMR/reference/read_gmt.md)).

## Value

A named list of character vectors (pathway name → gene symbols).

## See also

[`read_gmt`](https://hui950319.github.io/scMMR/reference/read_gmt.md),
[`ComputeModuleScore`](https://hui950319.github.io/scMMR/reference/ComputeModuleScore.md),
[`RunPathwayAnalysis`](https://hui950319.github.io/scMMR/reference/RunPathwayAnalysis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# From GMT file
gs <- parse_gene_sets("hallmark.gmt")

# From named list
gs <- parse_gene_sets(list(PathA = c("TP53", "MDM2"), PathB = c("EGFR", "ERBB2")))

# From data.frame
df <- data.frame(term = c("PathA", "PathA"), gene = c("TP53", "MDM2"))
gs <- parse_gene_sets(df)
} # }
```
