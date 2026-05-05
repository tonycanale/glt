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
                           neg_col = "blue",
                           mid_col = "white",
                           pos_col = "red",
                           n = 100) {
  # ...existing code...
  cols <- character(length(x))
  
  neg_idx <- which(x < 0)
  pos_idx <- which(x > 0)
  zero_idx <- which(x == 0)
  
  if (length(neg_idx) > 0) {
    neg_vals <- x[neg_idx]
    neg_min <- min(neg_vals, na.rm = TRUE)
    neg_scaled <- (neg_vals - 0) / (neg_min - 0)
    neg_pal <- colorRampPalette(c(mid_col, neg_col))(n)
    neg_pos <- pmax(1, pmin(n, round(neg_scaled * (n - 1)) + 1))
    cols[neg_idx] <- neg_pal[neg_pos]
  }
  
  if (length(pos_idx) > 0) {
    pos_vals <- x[pos_idx]
    pos_max <- max(pos_vals, na.rm = TRUE)
    pos_scaled <- pos_vals / pos_max
    pos_pal <- colorRampPalette(c(mid_col, pos_col))(n)
    pos_pos <- pmax(1, pmin(n, round(pos_scaled * (n - 1)) + 1))
    cols[pos_idx] <- pos_pal[pos_pos]
  }
  
  cols[zero_idx] <- mid_col
  cols[is.na(x)] <- NA
  
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
  
  cols <- palette_fun(as.vector(mat), ...)
  col_mat <- matrix(cols, nrow = nr, ncol = nc)
  
  col_mat[is.na(mat)] <- "#00000000"
  
  r <- as.raster(col_mat)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
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
       asp = asp)
  
  rasterImage(r,
              xleft = 0.5, ybottom = nr + 0.5,
              xright = nc + 0.5, ytop = 0.5,
              interpolate = FALSE)
  
  if (draw_grid) {
    abline(v = seq(0.5, nc + 0.5, by = 1), col = grid_col)
    abline(h = seq(0.5, nr + 0.5, by = 1), col = grid_col)
  }
  
  if (axes) {
    axis(1, at = 1:nc, labels = 1:nc)
    axis(2, at = 1:nr, labels = 1:nr, las = 1)
  }
  
  box()
}