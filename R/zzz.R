# ══════════════════════════════════════════════════════════════════════════════
# zzz.R — Package initialisation
# ══════════════════════════════════════════════════════════════════════════════

#' @import reticulate
#' @importFrom Matrix Matrix
#' @importFrom methods is

# Package-level environment (invisible to user)
.scMMR_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .scMMR_env$python_loaded <- FALSE
}
