# ══════════════════════════════════════════════════════════════════════════════
# python_env.R — Python environment management
# ══════════════════════════════════════════════════════════════════════════════

#' Install Python dependencies for scMMR
#'
#' Creates a conda or virtualenv environment and installs the required Python
#' packages (torch, scanpy, anndata, etc.).
#'
#' @param envname Name of the conda/virtual environment (default "scMMR").
#' @param method  Installation method: "conda" or "virtualenv" (default "conda").
#' @param gpu     Logical. If \code{FALSE} (default), installs CPU-only PyTorch
#'                (~200 MB). If \code{TRUE}, installs CUDA-enabled PyTorch (~2 GB).
#' @param python_version Python version to install (default "3.10").
#' @export
install_scMMR_python <- function(envname = "scMMR",
                                  method  = c("conda", "virtualenv"),
                                  gpu     = FALSE,
                                  python_version = "3.10") {
  method <- match.arg(method)

  # Base packages (always needed)
  packages <- c("anndata", "scanpy", "numpy", "pandas", "scipy", "h5py")

  if (method == "conda") {
    reticulate::conda_create(envname, python_version = python_version)

    # Install torch separately to control CPU vs GPU
    if (gpu) {
      reticulate::conda_install(envname, packages = c("pytorch", "torchvision",
                                                       "torchaudio"),
                                channel = "pytorch")
    } else {
      reticulate::conda_install(envname, packages = c("pytorch", "torchvision",
                                                       "torchaudio", "cpuonly"),
                                channel = "pytorch")
    }
    reticulate::conda_install(envname, packages = packages, pip = TRUE)

  } else {
    reticulate::virtualenv_create(envname, python = python_version)

    if (gpu) {
      reticulate::virtualenv_install(envname, packages = "torch")
    } else {
      reticulate::virtualenv_install(envname, packages = "torch",
                                      pip_options = c("--index-url",
                                                      "https://download.pytorch.org/whl/cpu"))
    }
    reticulate::virtualenv_install(envname, packages = packages)
  }

  message("\n\u2713 scMMR Python environment '", envname, "' installed successfully.")
  if (gpu) {
    message("  GPU (CUDA) support enabled.")
  } else {
    message("  CPU-only mode. Use gpu=TRUE to enable CUDA support.")
  }
  invisible(envname)
}


#' Set the Python environment for scMMR
#'
#' Configures which Python environment reticulate should use.
#' Must be called \strong{before} any scMMR function that uses Python
#' (i.e., before the first call to \code{DNN_train} or \code{DNN_predict}).
#'
#' @param condaenv   Name or path of a conda environment.
#' @param virtualenv Name or path of a virtualenv.
#' @param required   If TRUE (default), error when the env is not found.
#' @export
use_scMMR_python <- function(condaenv = NULL, virtualenv = NULL,
                              required = TRUE) {
  if (!is.null(condaenv)) {
    reticulate::use_condaenv(condaenv, required = required)
    message("\u2713 Using conda env: ", condaenv)
  } else if (!is.null(virtualenv)) {
    reticulate::use_virtualenv(virtualenv, required = required)
    message("\u2713 Using virtualenv: ", virtualenv)
  } else {
    stop("Provide either condaenv or virtualenv.")
  }
  invisible(NULL)
}


#' Lazy-load Python helper (internal)
#'
#' Called internally before any Python-dependent function.
#' Only runs \code{source_python()} once per session.
#'
#' @keywords internal
.ensure_python <- function() {
  if (.scMMR_env$python_loaded) return(invisible(NULL))

  py_path <- system.file("python", package = "scMMR")
  if (!nzchar(py_path)) {
    stop("Cannot find inst/python/ directory in scMMR package.")
  }

  # Add package python dir to sys.path so multi_task_model.py is importable
  # Use py_run_string because sys$path gets auto-converted to R vector
  reticulate::py_run_string(sprintf(
    "import sys; sys.path.insert(0, '%s') if '%s' not in sys.path else None",
    py_path, py_path
  ))

  helper_file <- file.path(py_path, "multitask_predict_helper.py")
  if (!file.exists(helper_file)) {
    stop("Python helper not found: ", helper_file)
  }

  reticulate::source_python(helper_file, envir = .scMMR_env)

  # Source deconv helper (deconvolution functions)
  deconv_file <- file.path(py_path, "deconv_helper.py")
  if (file.exists(deconv_file)) {
    reticulate::source_python(deconv_file, envir = .scMMR_env)
  }

  .scMMR_env$python_loaded <- TRUE
  invisible(NULL)
}
