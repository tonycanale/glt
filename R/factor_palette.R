#' Map numeric values to a diverging colour palette
#'
#' Assigns a hex colour to each element of a numeric vector using separate
#' colour ramps for negative and positive values, with a neutral midpoint colour
#' for zero.
#'
#' @param x Numeric vector of values to map.
#' @param neg_col Character. Colour for the most-negative values. Default
#'   `"blue"`.
#' @param mid_col Character. Colour for zero. Default `"white"`.
#' @param pos_col Character. Colour for the most-positive values. Default
#'   `"red"`.
#' @param n Integer. Number of steps in each colour ramp. Default `100`.
#'
#' @return A character vector of hex colour codes the same length as `x`.
#'   `NA` inputs produce `NA` outputs.
#'
#' @keywords internal
#' @noRd
signed_palette <- function(x,
                           neg_col = "red",
                           mid_col = "white",
                           pos_col = "blue",
                           n = 100,
                           log = FALSE,
                           max_abs = NULL) {
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }

  if (length(log) != 1L || is.na(log) || !is.logical(log)) {
    stop("`log` must be either TRUE or FALSE.")
  }

  if (length(n) != 1L || n < 2L) {
    stop("`n` must be an integer greater than one.")
  }

  cols <- rep(NA_character_, length(x))
  valid <- !is.na(x)

  if (!any(valid)) {
    return(cols)
  }

  if (any(!is.finite(x[valid]))) {
    stop("`x` must contain only finite values or NA.")
  }

  if (is.null(max_abs)) {
    max_abs <- max(abs(x[valid]))
  } else if (length(max_abs) != 1L ||
             !is.finite(max_abs) ||
             max_abs <= 0) {
    stop("`max_abs` must be a single positive finite number.")
  }

  # Return the midpoint colour when all entries are zero.
  if (max_abs == 0) {
    cols[valid] <- mid_col
    return(cols)
  }

  # Common magnitude scale for negative and positive values.
  magnitude <- pmin(abs(x[valid]), max_abs)

  if (log) {
    # log1p is defined at zero and preserves zero as the midpoint.
    scaled <- log1p(magnitude) / log1p(max_abs)
  } else {
    scaled <- magnitude / max_abs
  }

  position <- round(scaled * (n - 1L)) + 1L
  position <- pmax(1L, pmin(n, position))

  neg_pal <- grDevices::colorRampPalette(
    c(mid_col, neg_col)
  )(n)

  pos_pal <- grDevices::colorRampPalette(
    c(mid_col, pos_col)
  )(n)

  valid_values <- x[valid]
  valid_cols <- rep(mid_col, length(valid_values))

  negative <- valid_values < 0
  positive <- valid_values > 0

  valid_cols[negative] <- neg_pal[position[negative]]
  valid_cols[positive] <- pos_pal[position[positive]]

  cols[valid] <- valid_cols
  cols
}


#' Plot a real-valued matrix as a colour image
#'
#' Renders a numeric matrix as a raster image using a diverging colour palette
#' (via [signed_palette()] by default), with optional grid lines and axis labels.
#' Rows increase downward (matrix convention).
#'
#' @param mat A numeric matrix to plot.
#' @param palette_fun A function with signature `f(x, ...)` that maps a numeric
#'   vector to a character vector of hex colour codes. Defaults to
#'   [signed_palette()].
#' @param draw_grid Logical. Whether to draw cell grid lines. Default `TRUE`.
#' @param grid_col Character. Colour of the grid lines. Default `"grey80"`.
#' @param axes Logical. Whether to draw row and column axes. Default `TRUE`.
#' @param xlab Character. Label for the horizontal axis. Default `"column"`.
#' @param ylab Character. Label for the vertical axis. Default `"row"`.
#' @param asp Numeric. Aspect ratio passed to [plot()]. Default `1`.
#' @param main Character. Plot title. Default `""`.
#' @param legend Logical. Whether to draw a colour-scale legend to the right
#'   of the matrix image. Default `FALSE`.
#' @param legend_n Integer. Number of colour steps used to build the legend
#'   gradient. Default `100`.
#' @param legend_width Numeric. Width of the legend panel in inches.
#'   Default `1.4`.
#' @param ... Additional arguments forwarded to `palette_fun`.
#'
#' @return Invisibly `NULL`; called for its side effect of producing a plot.
#'
#' @seealso [signed_palette()]
#' @export
plot_real_matrix <- function(mat,
                             palette_fun = signed_palette,
                             draw_grid = TRUE,
                             grid_col = "grey80",
                             axes = TRUE,
                             xlab = "column",
                             ylab = "row",
                             asp = 1,
                             main ="",
                             legend = FALSE,
                             legend_n = 100,
                             legend_width = 1.4,
                             legend_digits = 2,
                             v_lines = NULL,
                             h_lines = NULL,
                             ...) {
  # ...existing code...
  if (!is.matrix(mat)) {
    stop("'mat' deve essere una matrice.")
  }
  if (!is.numeric(mat)) {
    stop("'mat' deve contenere valori numerici.")
  }
  
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  plot_asp <- if (nr == nc) asp else NA_real_
  
  cols <- palette_fun(as.vector(mat), ...)
  col_mat <- matrix(cols, nrow = nr, ncol = nc)
  
  col_mat[is.na(mat)] <- "#00000000"
  
  r <- as.raster(col_mat)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
if (legend) {
  device_width <- dev.size("in")[1]

  if (legend_width >= device_width) {
    stop("`legend_width` must be smaller than the device width.")
  }

  layout(
    matrix(1:2, nrow = 1),
    widths = c(device_width - legend_width, legend_width)
  )
}
  
  par(mar = c(4, 4, 2, 2) + 0.1)
  
  plot(NA,
       xlim = c(0.5, nc + 0.5),
       ylim = c(nr + 0.5, 0.5),
       xlab = xlab,
       ylab = ylab,
       xaxt = "n",
       yaxt = "n",
       bty = "n",
       main = main,
       asp = plot_asp)
  
  rasterImage(r,
              xleft = 0.5, ybottom = nr + 0.5,
              xright = nc + 0.5, ytop = 0.5,
              interpolate = FALSE)
  
  if (draw_grid) {
    abline(v = seq(0.5, nc + 0.5, by = 1), col = grid_col)
    abline(h = seq(0.5, nr + 0.5, by = 1), col = grid_col)
  }

  if (!is.null(h_lines)) {
    abline(h = h_lines, col = 1, lty = 2, lwd=1.5)
  }
  
  if (!is.null(v_lines)) {
    abline(v = v_lines, col = 1, lty = 2, lwd=1.5)
  }

  if (axes) {
    axis(1, at = 1:nc, labels = 1:nc)
    axis(2, at = 1:nr, labels = 1:nr, las = 1)
  }
  
  box()
  
  if (legend) {
    vals <- mat[!is.na(mat)]
    if (length(vals) == 0L) {
      warning("`mat` has no non-missing values; skipping legend.", call. = FALSE)
      return(invisible(NULL))
    }
    
    legend_vals <- seq(min(vals), max(vals), length.out = legend_n)
    legend_cols <- palette_fun(legend_vals, ...)
    legend_raster <- as.raster(matrix(rev(legend_cols), ncol = 1))
    
    par(mar = c(4, 1, 2, 3))
    
    plot(NA,
         xlim = c(0, 1),
         ylim = c(min(legend_vals), max(legend_vals)),
         xlab = "",
         ylab = "",
         xaxt = "n",
         yaxt = "n",
         bty = "n",
         main = "")
    
    rasterImage(legend_raster,
                xleft = 0, ybottom = min(legend_vals),
                xright = 1, ytop = max(legend_vals),
                interpolate = TRUE)
    
    legend_ticks <- pretty(range(legend_vals), n = 5)
legend_ticks <- legend_ticks[
  legend_ticks >= min(legend_vals) &
  legend_ticks <= max(legend_vals)
]

axis(
  4,
  at = legend_ticks,
  labels = formatC(
    legend_ticks,
    format = "f",
    digits = legend_digits
  ),
  las = 1
)
    box()
  }
  
  invisible(NULL)
}