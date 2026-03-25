# Post-process marker results: add gene/group/p_val_adj/test_group columns

Post-process marker results: add gene/group/p_val_adj/test_group columns

## Usage

``` r
.post_process_markers(
  markers,
  group1_str,
  group2_str = "others",
  p.adjust.method = "bonferroni",
  group_levels = NULL
)
```
