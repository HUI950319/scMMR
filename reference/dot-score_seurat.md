# Seurat-style module scoring (AddModuleScore-like)

Bins genes by mean expression, selects expression-matched control genes,
and computes score = mean(feature genes) - mean(control genes).

## Usage

``` r
.score_seurat(x, gene.sets, nbin, ctrl, cores)
```
