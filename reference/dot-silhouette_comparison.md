# Silhouette comparison between two spaces (internal)

Compute mean silhouette score per cell type in both embedding and
reference spaces. Falls back to a simple implementation if the 'cluster'
package is unavailable.

## Usage

``` r
.silhouette_comparison(emb, ref, labels, max_cells = 5000L)
```
