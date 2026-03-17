# ============================================================================
# palette.R — Color palette utilities for scMMR
# Simplified from thisplot::palette_colors(), independent of thisplot package
# ============================================================================

# --------------------------------------------------------------------------
# Internal: built-in palette list (RColorBrewer + ggsci hardcoded)
# --------------------------------------------------------------------------

#' @keywords internal
.get_palette_list <- function() {
  pal_list <- list()

  # --- RColorBrewer palettes ---
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    brewer_info <- RColorBrewer::brewer.pal.info
    for (pal_name in rownames(brewer_info)) {
      pal_n <- brewer_info[pal_name, "maxcolors"]
      pal_cat <- brewer_info[pal_name, "category"]
      cols <- RColorBrewer::brewer.pal(name = pal_name, n = pal_n)
      # reorder Paired for better default usage
      if (pal_name == "Paired") {
        cols <- cols[c(1:4, 7, 8, 5, 6, 9, 10, 11, 12)]
      }
      # reverse diverging palettes (so blue=low, red=high)
      if (pal_cat == "div") cols <- rev(cols)
      attr(cols, "type") <- ifelse(pal_cat == "qual", "discrete", "continuous")
      pal_list[[pal_name]] <- cols
    }
  }

  # --- ggsci palettes (hardcoded to avoid runtime dependency) ---
  ggsci_pals <- list(
    npg = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
            "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
            "#7E6148FF", "#B09C85FF"),
    nejm = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF",
             "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF"),
    lancet = c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF",
               "#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF",
               "#1B1919FF"),
    jama = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF",
             "#79AF97FF", "#6A6599FF", "#80796BFF"),
    aaas = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF",
             "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF",
             "#808180FF", "#1B1919FF"),
    d3 = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
            "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
            "#BCBD22FF", "#17BECFFF"),
    igv = c("#5050FFFF", "#CE3D32FF", "#749B58FF", "#F0E685FF",
            "#466983FF", "#BA6338FF", "#5DB1DDFF", "#802268FF",
            "#6BD76BFF", "#D595A7FF", "#924822FF", "#837B8DFF",
            "#C75127FF", "#D58F5CFF", "#7A65A5FF", "#E4AF69FF",
            "#3B1B53FF", "#CDDEB7FF", "#612A79FF", "#AE1F63FF"),
    simspec = c("#c22b86", "#f769a1", "#fcc5c1", "#253777",
                "#1d92c0", "#9ec9e1", "#015b33", "#42aa5e",
                "#d9f0a2", "#E66F00", "#f18c28", "#FFBB61")
  )
  for (nm in names(ggsci_pals)) {
    cols <- ggsci_pals[[nm]]
    attr(cols, "type") <- "discrete"
    pal_list[[nm]] <- cols
  }

  # --- continuous utility palettes ---
  jet_cols <- grDevices::colorRampPalette(
    c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
      "yellow", "#FF7F00", "red", "#7F0000")
  )(100)
  attr(jet_cols, "type") <- "continuous"
  pal_list[["jet"]] <- jet_cols

  GdRd <- c("gold", "red3")
  attr(GdRd, "type") <- "continuous"
  pal_list[["GdRd"]] <- GdRd

  pal_list
}


# ============================================================================
# palette_colors
# ============================================================================

#' Generate Named Color Vector from Palette
#'
#' Create a named color vector for discrete or continuous data, similar to
#' \code{thisplot::palette_colors()}. Can be used directly with
#' \code{ggplot2::scale_color_manual(values = ...)}.
#'
#' @param x A vector of character, factor, or numeric values.
#'   If missing, numeric values \code{1:n} are used.
#' @param n Integer. Number of colors to return for numeric values.
#'   Default 100.
#' @param palette Character. Name of a built-in palette.
#'   Use \code{show_palettes()} to see all available names.
#'   Default \code{"Paired"}.
#' @param palcolor Optional character vector of custom colors.
#'   When provided, overrides \code{palette}.
#' @param type Character, one of \code{"auto"}, \code{"discrete"}, or
#'   \code{"continuous"}. Default \code{"auto"} detects from \code{x}.
#' @param matched Logical. If \code{TRUE}, return a color for each element
#'   of \code{x} (same length). If \code{FALSE} (default), return a named
#'   vector of unique colors.
#' @param reverse Logical. Reverse the color order. Default \code{FALSE}.
#' @param NA_keep Logical. Keep color for \code{NA} values. Default \code{FALSE}.
#' @param NA_color Character. Color for \code{NA}. Default \code{"grey80"}.
#'
#' @return A named character vector of hex color codes.
#'
#' @examples
#' \dontrun{
#' palette_colors(c("A", "B", "C"))
#' palette_colors(c("A", "B", "C"), palette = "Dark2")
#' palette_colors(c("A", "B", "C"), palcolor = c("red", "blue", "green"))
#' palette_colors(1:100, palette = "Spectral")
#' }
#'
#' @export
palette_colors <- function(
    x,
    n         = 100,
    palette   = "Paired",
    palcolor  = NULL,
    type      = c("auto", "discrete", "continuous"),
    matched   = FALSE,
    reverse   = FALSE,
    NA_keep   = FALSE,
    NA_color  = "grey80"
) {
  palette_list <- .get_palette_list()

  if (missing(x)) {
    x <- 1:n
    type <- "continuous"
  }

  type <- match.arg(type)

  # resolve palcolor: custom colors override palette
  if (is.list(palcolor)) palcolor <- unlist(palcolor)
  if (is.null(palcolor) || length(palcolor) == 0 || all(palcolor == "")) {
    if (!palette %in% names(palette_list)) {
      stop("Palette '", palette, "' not found. ",
           "Use show_palettes() to see available palettes.")
    }
    palcolor <- palette_list[[palette]]
  }

  # named palcolor matching x values
  if (!is.null(names(palcolor))) {
    if (all(x %in% names(palcolor))) {
      palcolor <- palcolor[intersect(names(palcolor), x)]
    }
  }

  pal_n <- length(palcolor)

  if (type == "auto") {
    type <- if (is.numeric(x)) "continuous" else "discrete"
  }

  if (type == "discrete") {
    if (!is.factor(x)) x <- factor(x, levels = unique(x))
    n_x <- nlevels(x)

    if (isTRUE(attr(palcolor, "type") == "continuous")) {
      color <- grDevices::colorRampPalette(palcolor)(n_x)
    } else {
      color <- if (n_x <= pal_n) {
        palcolor[seq_len(n_x)]
      } else {
        grDevices::colorRampPalette(palcolor)(n_x)
      }
    }
    names(color) <- levels(x)

    if (any(is.na(x))) {
      color <- c(color, stats::setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      color <- color[as.character(x)]
      color[is.na(color)] <- NA_color
    }

  } else {
    # continuous
    if (!is.numeric(x) && all(!is.na(x))) {
      stop("x must be numeric for continuous color palettes.")
    }
    if (all(is.na(x))) {
      values <- as.factor(rep(0, n))
    } else if (length(unique(stats::na.omit(as.numeric(x)))) == 1) {
      values <- as.factor(rep(unique(stats::na.omit(as.numeric(x))), n))
    } else {
      breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                     length.out = n + 1)
      if (isTRUE(matched)) {
        values <- cut(x, breaks = breaks, include.lowest = TRUE)
      } else {
        values <- cut(1:100, breaks = breaks, include.lowest = TRUE)
      }
    }

    n_x <- nlevels(values)
    color <- if (n_x <= pal_n) {
      palcolor[seq_len(n_x)]
    } else {
      grDevices::colorRampPalette(palcolor)(n_x)
    }
    names(color) <- levels(values)

    if (any(is.na(x))) {
      color <- c(color, stats::setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      if (all(is.na(x))) {
        color <- NA_color
      } else {
        color <- color[as.character(values)]
        color[is.na(color)] <- NA_color
      }
    }
  }

  if (isTRUE(reverse)) color <- rev(color)
  if (isFALSE(NA_keep)) color <- color[names(color) != "NA"]

  color
}


# ============================================================================
# show_palettes
# ============================================================================

#' Display Available Color Palettes
#'
#' Show all built-in color palettes as a stacked bar chart.
#'
#' @param palette_names Optional character vector of palette names to display.
#'   If \code{NULL}, all palettes are shown.
#' @param type Character vector. Filter by palette type:
#'   \code{"discrete"} and/or \code{"continuous"}.
#'   Default shows both.
#' @param index Optional integer vector. Show only palettes at these
#'   positions.
#' @param return_palettes Logical. If \code{TRUE}, return the palette list
#'   invisibly instead of just printing. Default \code{FALSE}.
#'
#' @return If \code{return_palettes} is \code{TRUE}, returns a named list
#'   of color vectors. Otherwise, prints the plot and returns palette
#'   names invisibly.
#'
#' @examples
#' \dontrun{
#' show_palettes()
#' show_palettes(type = "discrete")
#' show_palettes(palette_names = c("Paired", "npg", "nejm", "Spectral"))
#' pals <- show_palettes(return_palettes = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual
#'   scale_x_continuous theme element_blank
#' @importFrom rlang .data
#' @export
show_palettes <- function(
    palette_names   = NULL,
    type            = c("discrete", "continuous"),
    index           = NULL,
    return_palettes = FALSE
) {
  palette_list <- .get_palette_list()

  # filter by type
  type <- match.arg(type, several.ok = TRUE)
  keep <- vapply(palette_list, function(x) {
    isTRUE(attr(x, "type") %in% type)
  }, logical(1))
  palette_list <- palette_list[keep]

  # filter by index
  if (!is.null(index)) {
    index <- index[index %in% seq_along(palette_list)]
    palette_list <- palette_list[index]
  }

  # filter by names
  if (!is.null(palette_names)) {
    missing_pals <- palette_names[!palette_names %in% names(palette_list)]
    if (length(missing_pals) > 0) {
      warning("Palettes not found: ", paste(missing_pals, collapse = ", "))
    }
    palette_list <- palette_list[palette_names[palette_names %in% names(palette_list)]]
  }

  if (length(palette_list) == 0) {
    message("No palettes to show.")
    return(invisible(NULL))
  }

  # build data.frame for plotting
  df <- data.frame(
    palette = rep(names(palette_list), vapply(palette_list, length, integer(1))),
    color   = unlist(palette_list, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  df$palette <- factor(df$palette, levels = rev(unique(df$palette)))
  df$color_order <- factor(seq_len(nrow(df)))
  df$proportion <- as.numeric(1 / table(df$palette)[as.character(df$palette)])

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      y    = .data$palette,
      x    = .data$proportion,
      fill = .data$color_order
    )
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = df$color) +
    ggplot2::scale_x_continuous(expand = c(0, 0), trans = "reverse") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )
  print(p)

  if (isTRUE(return_palettes)) {
    return(invisible(palette_list))
  }
  invisible(names(palette_list))
}
