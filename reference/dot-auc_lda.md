# AUC via stratified k-fold cross-validated LDA (internal)

Uses Linear Discriminant Analysis for faster and more stable
classification than logistic regression. Falls back to GLM if MASS is
unavailable.

## Usage

``` r
.auc_lda(mat1, mat2, n_folds = 5L)
```

## Details

When the number of features exceeds `min(n_per_class) / 2`, an
additional PCA step is applied to avoid the `p >> n` problem.
