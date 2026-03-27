#' @keywords internal
#' @noRd
update.pivots.new <- function(y, Delta, Eta, hyperpar) {
  Delta_int <- Delta
  storage.mode(Delta_int) <- "integer"

  update_pivots_cpp(
    y = y,
    Delta_in = Delta_int,
    Eta = Eta,
    hyperpar = hyperpar
  )
}
