# Install Python dependencies for scMMR

Creates a conda or virtualenv environment and installs the required
Python packages (torch, scanpy, anndata, etc.).

## Usage

``` r
install_scMMR_python(
  envname = "scMMR",
  method = c("conda", "virtualenv"),
  gpu = FALSE,
  python_version = "3.10"
)
```

## Arguments

- envname:

  Name of the conda/virtual environment (default "scMMR").

- method:

  Installation method: "conda" or "virtualenv" (default "conda").

- gpu:

  Logical. If `FALSE` (default), installs CPU-only PyTorch (~200 MB). If
  `TRUE`, installs CUDA-enabled PyTorch (~2 GB).

- python_version:

  Python version to install (default "3.10").
