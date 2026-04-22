signed_palette <- function(x,
                           neg_col = "blue",
                           mid_col = "white",
                           pos_col = "red",
                           n = 100) {
  cols <- character(length(x))
  
  neg_idx <- which(x < 0)
  pos_idx <- which(x > 0)
  zero_idx <- which(x == 0)
  
  if (length(neg_idx) > 0) {
    neg_vals <- x[neg_idx]
    neg_min <- min(neg_vals, na.rm = TRUE)
    neg_scaled <- (neg_vals - 0) / (neg_min - 0)   # [neg_min, 0] -> [1, 0]
    neg_pal <- colorRampPalette(c(mid_col, neg_col))(n)
    neg_pos <- pmax(1, pmin(n, round(neg_scaled * (n - 1)) + 1))
    cols[neg_idx] <- neg_pal[neg_pos]
  }
  
  if (length(pos_idx) > 0) {
    pos_vals <- x[pos_idx]
    pos_max <- max(pos_vals, na.rm = TRUE)
    pos_scaled <- pos_vals / pos_max               # [0, pos_max] -> [0, 1]
    pos_pal <- colorRampPalette(c(mid_col, pos_col))(n)
    pos_pos <- pmax(1, pmin(n, round(pos_scaled * (n - 1)) + 1))
    cols[pos_idx] <- pos_pal[pos_pos]
  }
  
  cols[zero_idx] <- mid_col
  cols[is.na(x)] <- NA
  
  cols
}


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
  if (!is.matrix(mat)) {
    stop("'mat' deve essere una matrice.")
  }
  if (!is.numeric(mat)) {
    stop("'mat' deve contenere valori numerici.")
  }
  
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  # Mappa ogni elemento in un colore
  cols <- palette_fun(as.vector(mat), ...)
  col_mat <- matrix(cols, nrow = nr, ncol = nc)
  
  # Se ci sono NA, li rendo trasparenti
  col_mat[is.na(mat)] <- "#00000000"
  
  # as.raster interpreta la prima riga della matrice come la riga in alto:
  # esattamente ciò che vogliamo
  r <- as.raster(col_mat)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  par(mar = c(4, 4, 2, 2) + 0.1)
  
  plot(NA,
       xlim = c(0.5, nc + 0.5),
       ylim = c(nr + 0.5, 0.5),   # asse y invertito: righe crescono verso il basso
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