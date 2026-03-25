# Clean pathway names by removing common prefixes and converting to Title Case

Removes database-specific prefixes such as `HALLMARK_`, `KEGG_`,
`REACTOME_`, etc., replaces underscores with spaces, and converts the
result to Title Case for cleaner visualization.

## Usage

``` r
.clean_pathway_names(names_vec, prefix = NULL)
```

## Arguments

- names_vec:

  Character vector of pathway names.

- prefix:

  Optional custom prefix to remove (regex pattern). If `NULL` (default),
  common prefixes are removed automatically.

## Value

Cleaned character vector.

## Examples

``` r
if (FALSE) { # \dontrun{
.clean_pathway_names(c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "KEGG_GLYCOLYSIS"))
# [1] "Tnfa Signaling Via Nfkb" "Glycolysis"
} # }
```
