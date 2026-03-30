# ============================================================================
# raster.R -- Lightweight layer rasterisation for ggplot2
# ============================================================================

#' Rasterise a ggplot2 layer
#'
#' Wraps an existing ggplot2 layer so that its graphical output is rendered to
#' an off-screen device and converted to a bitmap (\code{rasterGrob}).
#' This dramatically reduces file size and rendering time for scatter plots
#' with many points, especially when saving to PDF/SVG.
#'
#' The function is self-contained and does not depend on \pkg{scattermore} or
#' \pkg{ggrastr}. It supports three rendering backends with automatic
#' fallback:
#' \enumerate{
#'   \item \strong{Cairo} (default if available) -- in-memory raster device,
#'     fastest.
#'   \item \strong{ragg} -- high-quality AGG-based capture device.
#'   \item \strong{png} -- file-based fallback using \code{grDevices::png()}
#'     and \code{png::readPNG()}.
#' }
#'
#' @param layer A ggplot2 layer object (e.g. the return value of
#'   \code{geom_point(...)}).
#' @param dpi Numeric scalar. Resolution in dots per inch. Default \code{300}.
#' @param dev Character. Rendering backend: \code{"cairo"}, \code{"ragg"}, or
#'   \code{"png"}. Default \code{NULL} (auto-detect best available).
#' @param scale Numeric scalar > 0. Scaling factor applied to the rendered
#'   bitmap dimensions. Values < 1 reduce resolution (faster), values > 1
#'   increase it. Default \code{1}.
#'
#' @return A modified ggplot2 layer whose \code{draw_geom} method produces
#'   rasterised grobs.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' df <- data.frame(x = rnorm(1e5), y = rnorm(1e5))
#' ggplot(df, aes(x, y)) + rasterise_layer(geom_point(alpha = 0.3))
#' }
#'
#' @export
rasterise_layer <- function(layer, dpi = 300, dev = NULL, scale = 1) {
  if (!inherits(layer, "Layer")) {
    stop("'layer' must be a ggplot2 Layer object.", call. = FALSE)
  }
  if (!is.numeric(scale) || scale <= 0) {
    stop("'scale' must be a positive number.", call. = FALSE)
  }

  # Auto-detect best backend
  if (is.null(dev)) {
    dev <- if (requireNamespace("Cairo", quietly = TRUE)) {
      "cairo"
    } else if (requireNamespace("ragg", quietly = TRUE)) {
      "ragg"
    } else {
      "png"
    }
  }
  dev <- match.arg(dev, c("cairo", "ragg", "png"))

  # Validate backend availability
  if (dev == "cairo" && !requireNamespace("Cairo", quietly = TRUE)) {
    warning("Cairo package not available, falling back to 'png' backend.",
            call. = FALSE)
    dev <- "png"
  }
  if (dev == "ragg" && !requireNamespace("ragg", quietly = TRUE)) {
    warning("ragg package not available, falling back to 'png' backend.",
            call. = FALSE)
    dev <- "png"
  }
  if (dev == "png" && !requireNamespace("png", quietly = TRUE)) {
    stop("At least one of 'Cairo', 'ragg', or 'png' package is required ",
         "for rasterisation.", call. = FALSE)
  }

  force(layer)
  force(dpi)
  force(dev)
  force(scale)

  ggplot2::ggproto(NULL, layer,
    draw_geom = function(self, data, layout) {
      grobs <- ggplot2::ggproto_parent(layer, self)$draw_geom(data, layout)
      lapply(grobs, function(grob) {
        if (inherits(grob, "zeroGrob")) return(grob)
        class(grob) <- c("scmmr_raster", class(grob))
        grob$scmmr_dpi   <- dpi
        grob$scmmr_dev   <- dev
        grob$scmmr_scale <- scale
        grob
      })
    }
  )
}


#' Grid makeContext method for scMMR rasterised grobs
#'
#' S3 method called by the grid rendering system. Renders the original grob
#' to an off-screen bitmap device and returns a \code{rasterGrob}.
#'
#' @param x A grob with class \code{"scmmr_raster"}.
#' @return A \code{rasterGrob}.
#' @importFrom grid makeContext
#' @method makeContext scmmr_raster
#' @export
makeContext.scmmr_raster <- function(x) {

  # Extract and clean raster metadata
  dpi   <- x$scmmr_dpi   %||% 300
  dev   <- x$scmmr_dev   %||% "cairo"
  scale <- x$scmmr_scale %||% 1

  x$scmmr_dpi   <- NULL
  x$scmmr_dev   <- NULL
  x$scmmr_scale <- NULL
  class(x) <- setdiff(class(x), "scmmr_raster")

  # Viewport for coordinate context
  vp <- if (is.null(x$vp)) grid::viewport() else x$vp

  # Compute dimensions in inches
  width  <- grid::convertWidth(grid::unit(1, "npc"), "inch", valueOnly = TRUE)
  height <- grid::convertHeight(grid::unit(1, "npc"), "inch", valueOnly = TRUE)

  # Apply scale
  w_dev <- width  / scale
  h_dev <- height / scale

  # Save and restore current device
  dev_prev <- grDevices::dev.cur()
  on.exit(if (dev_prev > 1) grDevices::dev.set(dev_prev), add = TRUE)

  # --- Render to off-screen device ---
  if (dev == "cairo") {
    # Cairo in-memory raster: fastest path
    Cairo::Cairo(
      type = "raster", width = w_dev, height = h_dev,
      units = "in", dpi = dpi, bg = NA
    )
    grid::pushViewport(vp)
    grid::grid.draw(x)
    grid::popViewport()
    cap <- grid::grid.cap()
    grDevices::dev.off()

  } else if (dev == "ragg") {
    # ragg capture device: high-quality AGG backend
    ragg::agg_capture(
      width = w_dev, height = h_dev,
      units = "in", res = dpi, background = NA
    )
    grid::pushViewport(vp)
    grid::grid.draw(x)
    grid::popViewport()
    cap <- grid::grid.cap()
    grDevices::dev.off()

  } else {
    # File-based fallback: png device + png::readPNG
    tmpfile <- tempfile(fileext = ".png")
    on.exit(unlink(tmpfile), add = TRUE)

    grDevices::png(
      tmpfile, width = w_dev, height = h_dev,
      units = "in", res = dpi, bg = "transparent",
      type = if (capabilities("cairo")) "cairo" else "Xlib"
    )
    grid::pushViewport(vp)
    grid::grid.draw(x)
    grid::popViewport()
    grDevices::dev.off()

    # Read back as RGBA matrix → colour matrix
    rgba <- png::readPNG(tmpfile, native = FALSE)
    dm   <- dim(rgba)
    cap  <- matrix(
      grDevices::rgb(
        red   = as.vector(rgba[, , 1]),
        green = as.vector(rgba[, , 2]),
        blue  = as.vector(rgba[, , 3]),
        alpha = as.vector(rgba[, , 4])
      ),
      nrow = dm[1], ncol = dm[2]
    )
  }

  # Return rasterGrob at original (pre-scale) size
  grid::rasterGrob(
    cap,
    x = 0.5, y = 0.5,
    width  = grid::unit(width,  "inch"),
    height = grid::unit(height, "inch"),
    default.units = "npc",
    just = "center"
  )
}


#' Rasterised point geom
#'
#' Convenience wrapper that creates a \code{geom_point} layer and rasterises
#' it via \code{\link{rasterise_layer}}.
#'
#' @param raster.dpi Numeric. Raster resolution. Default \code{300}.
#' @param raster.dev Character. Backend device. Default \code{NULL}
#'   (auto-detect).
#' @param raster.scale Numeric. Scale factor. Default \code{1}.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{geom_point}}.
#'
#' @return A rasterised ggplot2 layer.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' df <- data.frame(x = rnorm(1e5), y = rnorm(1e5))
#' ggplot(df, aes(x, y)) + geom_point_rast(alpha = 0.3)
#' }
#'
#' @export
geom_point_rast <- function(..., raster.dpi = 300, raster.dev = NULL,
                            raster.scale = 1) {
  rasterise_layer(
    ggplot2::geom_point(...),
    dpi   = raster.dpi,
    dev   = raster.dev,
    scale = raster.scale
  )
}
