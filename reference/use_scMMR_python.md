# Set the Python environment for scMMR

Configures which Python environment reticulate should use. Must be
called **before** any scMMR function that uses Python (i.e., before the
first call to `DNN_train` or `DNN_predict`).

## Usage

``` r
use_scMMR_python(condaenv = NULL, virtualenv = NULL, required = TRUE)
```

## Arguments

- condaenv:

  Name or path of a conda environment.

- virtualenv:

  Name or path of a virtualenv.

- required:

  If TRUE (default), error when the env is not found.
