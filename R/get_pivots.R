#' Extract pivot row indices from a Delta matrix
#'
#' For each column of `Delta`, returns the row index of the first `1`, which
#' defines the pivot row for that factor.
#'
#' @param Delta An `m x H` integer matrix with entries in `{0, 1}`. Each column
#'   must contain at least one `1`.
#'
#' @return An integer vector of length `H`.
#'
#' @export
get_pivots <- function(Delta) {
  if (!is.matrix(Delta)) {
    stop("`Delta` must be a matrix.", call. = FALSE)
  }
  pivots <- apply(Delta, 2L, match, x = 1L)
  if (length(pivots) != length(unique(pivots))) {
    stop("Pivots need to be on different rows.", call. = FALSE)
  }
  if (anyNA(pivots)) {
    stop("Each column of `Delta` must contain at least one 1 to define a pivot.", call. = FALSE)
  }
  as.integer(pivots)
}
