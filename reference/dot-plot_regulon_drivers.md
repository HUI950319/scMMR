# Generate plots for one pathway in FindRegulonDrivers

Generate plots for one pathway in FindRegulonDrivers

## Usage

``` r
.plot_regulon_drivers(df, pw_name, label.n = 5, label.overlap = NULL)
```

## Arguments

- df:

  Data.frame for one pathway (from result_df).

- pw_name:

  Pathway name.

- label.n:

  Top N for rank barplot.

- label.overlap:

  Overlap threshold for scatter labeling.

## Value

List with rank_bar and scatter plots.
