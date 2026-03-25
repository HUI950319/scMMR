# Lazy-load Python helper (internal)

Called internally before any Python-dependent function. Only runs
`source_python()` once per session.

## Usage

``` r
.ensure_python()
```
